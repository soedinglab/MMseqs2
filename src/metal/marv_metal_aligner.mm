// Metal gapless aligner behind MultiConfigScoreOnlyAlignerInterface plus the aligner factories the shared code calls.
// Short2 aligner: FP16 pass followed by a GPU-filtered int16 correction. Int aligner: int32 kernel.
#include "../aligners/make_aligner_headers.cuh"
#include "cuda_metal_launch.h"
#include "marv_metal_shared.h"

#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>

namespace libmarv {

namespace {

struct KernelChoice {
    std::string primary;
    std::string fallback;
    unsigned subjectsPerSimdgroup = 1;
    size_t threadgroupBytes = 0;
    unsigned threads = 1024;
};

struct TiledFamily { unsigned g, minNR, maxNR; };
const TiledFamily kTypedCache[] = {{4u, 4u, 48u}, {8u, 4u, 64u}, {16u, 4u, 96u}, {32u, 4u, 48u}};
const TiledFamily kInt8Cache[]  = {{16u, 48u, 96u}, {32u, 24u, 48u}};

template<size_t N>
bool pickTiled(const TiledFamily (&families)[N], uint32_t queryLength, uint32_t alphabetSize, size_t budget,
               size_t cacheElemBytes, unsigned forcedGroup, unsigned& groupSize, unsigned& numRegs, size_t& tgBytes) {
    for (const TiledFamily& f : families) {
        if (forcedGroup && forcedGroup != f.g) continue;
        unsigned nr = (((queryLength + f.g - 1u) / f.g) + 3u) & ~3u;
        nr = std::max(nr, f.minNR);
        if (nr > f.maxNR) continue;
        const size_t bytes = (size_t) alphabetSize * f.g * nr * cacheElemBytes;
        if (bytes > budget) continue;
        groupSize = f.g; numRegs = nr; tgBytes = bytes;
        return true;
    }
    return false;
}

unsigned untiledRegisters(uint32_t length) {
    const unsigned needed = (length + 127u) / 128u;
    for (unsigned n : {1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u, 10u, 11u, 12u, 13u, 14u, 15u, 16u, 20u, 24u, 32u, 40u, 48u, 64u})
        if (needed <= n) return n;
    throw std::runtime_error("libmarv/Metal: query exceeds the supported maximum of 8192");
}

std::string untiledKernel(uint32_t length, bool int32Mode) {
    return std::string(int32Mode ? "gl_untiled_i_r" : "gl_untiled_s_r") + std::to_string(untiledRegisters(length));
}

KernelChoice chooseKernel(uint32_t queryLength, uint32_t alphabetSize, bool int32Mode) {
    KernelChoice kc;
    if (queryLength > 32u * 4u * 64u) {
        throw std::runtime_error("libmarv/Metal: query length " + std::to_string(queryLength) + " exceeds the supported maximum of 8192 for the gapless scan");
    }
    if (int32Mode) {
        kc.primary = untiledKernel(queryLength, true);
        return kc;
    }
    const size_t budget = metalMaxThreadgroupMemory() - 1024;
    unsigned g = 0, nr = 0; size_t tgBytes = 0;
    // fewer registers per lane for medium queries (measured on M4)
    const unsigned preferredGroup = queryLength > 128 && queryLength <= 512 ? (queryLength <= 256 ? 8u : 16u) : 0u;
    if (queryLength > 256 && queryLength <= 512) kc.threads = 512;
    if (pickTiled(kTypedCache, queryLength, alphabetSize, budget, sizeof(uint16_t), preferredGroup, g, nr, tgBytes)) {
        kc.primary = "gl_tg_h_g" + std::to_string(g) + "_nr" + std::to_string(nr);
        kc.fallback = "gl_tg_s_g" + std::to_string(g) + "_nr" + std::to_string(nr);
    } else if (pickTiled(kInt8Cache, queryLength, alphabetSize, budget, 1, 0, g, nr, tgBytes)) {
        kc.primary = "gl_tg_h8_g" + std::to_string(g) + "_nr" + std::to_string(nr);
        kc.fallback = "gl_tg_s8_g" + std::to_string(g) + "_nr" + std::to_string(nr);
    } else {
        kc.primary = untiledKernel(queryLength, false);
        return kc;
    }
    kc.subjectsPerSimdgroup = 32u / g;
    kc.threadgroupBytes = tgBytes;
    return kc;
}

size_t resolveOrThrow(const void* ptr, const char* what) {
    size_t off = 0;
    if (!metalResolvePointer(ptr, &off)) {
        throw std::runtime_error(std::string("libmarv/Metal: ") + what + " is not device memory known to the Metal runtime");
    }
    return off;
}

class MetalGaplessAligner final : public MultiConfigScoreOnlyAlignerInterface {
public:
    explicit MetalGaplessAligner(bool int32Mode_) : int32Mode(int32Mode_) {}

    ~MetalGaplessAligner() override { marvBackendSync(); }

    bool isSingleTile(const OneToAllInputDataPSSM&) const override { return true; }
    int preferredKernelOutputArraySize() const override { return 1'000'000; }
    int getTileSize(const OneToAllInputDataPSSM&) const override { return 8192; }
    AlignmentAlgorithmE getAlignmentAlgorithm() const override { return AlignmentAlgorithmE::Gapless; }
    Datatype getScoreDatatype(int) const override { return int32Mode ? Datatype::Int : Datatype::Short2; }
    size_t getMinimumSuggestedTempBytes_multiTile(const OneToAllInputDataPSSM&) const override { return 0; }

    void makeGpuPssm(GpuConvertedPSSM& result, const cuda::std::mdspan<std::int8_t, cuda::std::dextents<int, 2>> hostPssmView,
                     cudaStream_t) override {
        const int rows = hostPssmView.extent(0);
        const int cols = hostPssmView.extent(1);
        result.kindOfPssm = GpuConvertedPSSM::Kind::Gapless;
        result.alphabetSize = rows;
        result.queryLength = cols;
        result.groupsize = 0;
        result.numRegs = 0;
        result.allocatedNumRegs = 0;
        result.columnstride = 128 * untiledRegisters(cols);
        result.gaplessScoreBound = 0;
        for (int c = 0; c < cols; ++c) {
            int best = 0;
            for (int r = 0; r < rows; ++r) best = std::max(best, (int) hostPssmView(r, c));
            result.gaplessScoreBound += best;
        }
        // padded once per query for the tiled cache and the fixed-register untiled kernels
        std::vector<int8_t> padded((size_t) rows * result.columnstride, 0);
        for (int r = 0; r < rows; ++r)
            memcpy(padded.data() + (size_t) r * result.columnstride, hostPssmView.data_handle() + (size_t) r * cols, cols);
        result.data.resize(padded.size());
        cudaMemcpy(result.data.data(), padded.data(), padded.size(), cudaMemcpyHostToDevice);
    }

    void scoreOnly_singleTile(int* const devAlignmentScores, const OneToAllInputDataPSSM& inputData,
                              const GpuConvertedPSSM& pssm, GapScoreArgs, cudaStream_t stream) override {
        run(devAlignmentScores, inputData, pssm, stream);
    }

    void scoreOnly_multiTile(char*, size_t, int* const devAlignmentScores, const OneToAllInputDataPSSM& inputData,
                             const GpuConvertedPSSM& pssm, GapScoreArgs, cudaStream_t stream) override {
        run(devAlignmentScores, inputData, pssm, stream);
    }

private:
    bool int32Mode;
    void run(int* const devAlignmentScores, const OneToAllInputDataPSSM& in, const GpuConvertedPSSM& pssm, cudaStream_t stream) {
        const uint32_t L = (uint32_t) in.queryLength;
        const int n = in.numAlignments;
        if (n <= 0 || L == 0) return;
        if (pssm.kindOfPssm != GpuConvertedPSSM::Kind::Gapless || pssm.queryLength != (SequenceLengthT) L) {
            throw std::runtime_error("libmarv/Metal: PSSM does not match the query");
        }
        const uint32_t alphabetSize = (uint32_t) pssm.alphabetSize;

        // d_subjectOffsets is transform(counting, {convertToZeroBased, parent})
        const auto& offsetOp = in.d_subjectOffsets.functor();
        const size_t* offsets = offsetOp.parent;

        MarvMetalGaplessParams params{};
        params.numItems = (uint32_t) n;
        params.queryLength = L;
        params.pssmStride = pssm.columnstride;
        params.alphabetSize = alphabetSize;
        params.charsBase = resolveOrThrow(in.d_subjects, "subject data");
        const size_t offBytes = resolveOrThrow(offsets, "subject offsets");
        const size_t lenBytes = resolveOrThrow(in.d_subjectLengths, "subject lengths");
        const size_t outBytes = resolveOrThrow(devAlignmentScores, "score output");
        params.pssmBase = resolveOrThrow(pssm.data.data(), "PSSM");
        if (offBytes % sizeof(size_t) || lenBytes % sizeof(int32_t) || outBytes % sizeof(int32_t)) {
            throw std::runtime_error("libmarv/Metal: misaligned offset, length or score arrays");
        }
        params.offsetsBase = offBytes / sizeof(size_t);
        if (offsetOp.convertToZeroBased) params.flags |= MARV_METAL_FLAG_ZERO_OFFSETS;
        params.lengthsBase = lenBytes / sizeof(int32_t);
        params.outBase = outBytes / sizeof(int32_t);
        const int32_t* idList = in.indexIndirection.ptr;
        if (idList != nullptr) {
            const size_t idBytes = resolveOrThrow(idList, "subject id list");
            if (idBytes % sizeof(int32_t)) throw std::runtime_error("libmarv/Metal: misaligned subject id list");
            params.idsBase = idBytes / sizeof(int32_t);
            params.flags |= MARV_METAL_FLAG_ID_LIST;
        } else {
            params.indexOffset = (uint32_t) in.indexIndirection.offset;
        }

        // setBytes copies the table into the command encoder
        unsigned char unusedTable[256]{};
        const unsigned char* table = in.letterTranslation ? in.letterTranslation : unusedTable;
        if (in.letterTranslation) params.flags |= MARV_METAL_FLAG_TRANSLATE;
        KernelChoice kc = chooseKernel(L, alphabetSize, int32Mode);
        if (L <= 128 && n < 4096) kc.threads = 128;
        const void* idArg = idList != nullptr ? (const void*) idList : (const void*) devAlignmentScores;
        MetalKernelArg args[8];
        args[0] = MetalKernelArg::buffer(in.d_subjects);
        args[1] = MetalKernelArg::buffer(offsets);
        args[2] = MetalKernelArg::buffer(in.d_subjectLengths);
        args[3] = MetalKernelArg::buffer(pssm.data.data());
        args[4] = MetalKernelArg::buffer(devAlignmentScores);
        args[5] = MetalKernelArg::bytes(params);
        args[6] = MetalKernelArg::buffer(idArg);
        args[7] = MetalKernelArg{table, 256, false};

        auto launch = [&](const std::string& name) {
            const size_t maxThreads = metalKernelMaxThreads(name.c_str());
            if (maxThreads == 0) throw std::runtime_error("libmarv/Metal: kernel not found: " + name);
            const unsigned sgPerTg = (unsigned) std::max<size_t>(1, std::min<size_t>(kc.threads / 32, maxThreads / 32));
            const unsigned itemsPerTg = sgPerTg * kc.subjectsPerSimdgroup;
            const unsigned grid = (unsigned)(((size_t) n + itemsPerTg - 1) / itemsPerTg);
            cudaError_t status = metalLaunchKernel(name.c_str(), args, 8, grid, sgPerTg * 32u, kc.threadgroupBytes, stream);
            if (status != cudaSuccess) throw std::runtime_error("libmarv/Metal: launch of " + name + " failed");
        };
        launch(kc.primary);
        if (!kc.fallback.empty() && pssm.gaplessScoreBound >= 2048) {
            // the correction reads the previous dispatch's scores and skips exact FP16 results
            params.flags |= MARV_METAL_FLAG_RETRY;
            launch(kc.fallback);
        }
    }
};

[[noreturn]] void notOnMetal(const char* what) {
    throw std::runtime_error(std::string("libmarv/Metal: ") + what + " is not implemented on the Metal backend");
}

}

std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_8bit() { return std::make_unique<MetalGaplessAligner>(false); }
std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_16bit() { return std::make_unique<MetalGaplessAligner>(false); }
std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_32bit() { return std::make_unique<MetalGaplessAligner>(true); }

std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_sw_aligner_32bit() { notOnMetal("Smith-Waterman"); }

std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_16bit() { notOnMetal("gapless end-position alignment"); }
std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_32bit() { notOnMetal("gapless end-position alignment"); }
std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_sw_endpos_aligner_32bit() { notOnMetal("Smith-Waterman end-position alignment"); }

}
