// SIMD (CPU) gapless aligner behind MultiConfigScoreOnlyAlignerInterface: query-striped kernels (Farrar layout) with
// 8-bit saturating (bias trick), 16-bit saturating and 32-bit scalar precision, parallel over subjects with OpenMP.
// GNU vector extensions with 16/32/64-byte vectors (64 needs AVX-512BW); MARV_SIMD_WIDTH=16|32|64 overrides.
#include "../aligners/make_aligner_headers.cuh"

#include <omp.h>
#include <functional>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#if defined(__SSE2__)
#include <immintrin.h>
#endif

namespace libmarv {

namespace {

template<int W> struct V;
template<> struct V<16> { typedef uint8_t u8 __attribute__((vector_size(16))); typedef int16_t i16 __attribute__((vector_size(16))); };
template<> struct V<32> { typedef uint8_t u8 __attribute__((vector_size(32))); typedef int16_t i16 __attribute__((vector_size(32))); };
#if defined(__AVX512BW__)
template<> struct V<64> { typedef uint8_t u8 __attribute__((vector_size(64))); typedef int16_t i16 __attribute__((vector_size(64))); };
#endif

#if defined(__has_builtin)
#if __has_builtin(__builtin_elementwise_add_sat) && __has_builtin(__builtin_elementwise_sub_sat) && __has_builtin(__builtin_elementwise_max)
#define MARV_SIMD_ELEMENTWISE_BUILTINS 1
#endif
#endif

template<int W> inline typename V<W>::u8 splat_u8(uint8_t x) { typename V<W>::u8 r; for (int i = 0; i < W; i++) r[i] = x; return r; }

#if defined(__SSE2__) && !defined(__AVX2__)
// hardware saturation for SSE2-only builds, also for emulated 32-byte vectors
template<class T, class Op> inline T sseBinary(T a, T b, Op op) {
    __m128i aa[sizeof(T) / 16], bb[sizeof(T) / 16], rr[sizeof(T) / 16];
    std::memcpy(aa, &a, sizeof(T));
    std::memcpy(bb, &b, sizeof(T));
    for (unsigned i = 0; i < sizeof(T) / 16; ++i) rr[i] = op(aa[i], bb[i]);
    T result;
    std::memcpy(&result, rr, sizeof(T));
    return result;
}
#endif

template<class T> inline T adds_u8(T a, T b) {
#if defined(__SSE2__)
#if defined(__AVX512BW__)
    if constexpr (sizeof(T) == 64) return (T) _mm512_adds_epu8((__m512i) a, (__m512i) b);
    else
#endif
    if constexpr (sizeof(T) == 16) return (T) _mm_adds_epu8((__m128i) a, (__m128i) b);
#if defined(__AVX2__)
    else return (T) _mm256_adds_epu8((__m256i) a, (__m256i) b);
#else
    else return sseBinary(a, b, [](__m128i x, __m128i y) { return _mm_adds_epu8(x, y); });
#endif
#endif
#ifdef MARV_SIMD_ELEMENTWISE_BUILTINS
    return __builtin_elementwise_add_sat(a, b);
#else
    T s = a + b; return s | (T)(s < a);
#endif
}
template<class T> inline T subs_u8(T a, T b) {
#if defined(__SSE2__)
#if defined(__AVX512BW__)
    if constexpr (sizeof(T) == 64) return (T) _mm512_subs_epu8((__m512i) a, (__m512i) b);
    else
#endif
    if constexpr (sizeof(T) == 16) return (T) _mm_subs_epu8((__m128i) a, (__m128i) b);
#if defined(__AVX2__)
    else return (T) _mm256_subs_epu8((__m256i) a, (__m256i) b);
#else
    else return sseBinary(a, b, [](__m128i x, __m128i y) { return _mm_subs_epu8(x, y); });
#endif
#endif
#ifdef MARV_SIMD_ELEMENTWISE_BUILTINS
    return __builtin_elementwise_sub_sat(a, b);
#else
    T d = a - b; return d & (T)(a >= b);
#endif
}
template<class T> inline T max_u8(T a, T b) {
#if defined(__SSE2__)
#if defined(__AVX512BW__)
    if constexpr (sizeof(T) == 64) return (T) _mm512_max_epu8((__m512i) a, (__m512i) b);
    else
#endif
    if constexpr (sizeof(T) == 16) return (T) _mm_max_epu8((__m128i) a, (__m128i) b);
#if defined(__AVX2__)
    else return (T) _mm256_max_epu8((__m256i) a, (__m256i) b);
#else
    else return sseBinary(a, b, [](__m128i x, __m128i y) { return _mm_max_epu8(x, y); });
#endif
#endif
#ifdef MARV_SIMD_ELEMENTWISE_BUILTINS
    return __builtin_elementwise_max(a, b);
#else
    T m = (T)(a < b); return (a & ~m) | (b & m);
#endif
}
template<class T> inline T adds_i16(T a, T b) {
#if defined(__SSE2__)
#if defined(__AVX512BW__)
    if constexpr (sizeof(T) == 64) return (T) _mm512_adds_epi16((__m512i) a, (__m512i) b);
    else
#endif
    if constexpr (sizeof(T) == 16) return (T) _mm_adds_epi16((__m128i) a, (__m128i) b);
#if defined(__AVX2__)
    else return (T) _mm256_adds_epi16((__m256i) a, (__m256i) b);
#else
    else return sseBinary(a, b, [](__m128i x, __m128i y) { return _mm_adds_epi16(x, y); });
#endif
#endif
#ifdef MARV_SIMD_ELEMENTWISE_BUILTINS
    return __builtin_elementwise_add_sat(a, b);
#else
    T s = a + b;
    T ovf = (a ^ s) & (b ^ s);          // negative where the sign flipped
    T sat = (a >> 15) ^ (int16_t) 0x7FFF;  // 0x8000 for negative a, 0x7FFF otherwise
    T m = ovf >> 15;
    return (s & ~m) | (sat & m);
#endif
}
template<class T> inline T max_i16(T a, T b) {
#if defined(__SSE2__)
#if defined(__AVX512BW__)
    if constexpr (sizeof(T) == 64) return (T) _mm512_max_epi16((__m512i) a, (__m512i) b);
    else
#endif
    if constexpr (sizeof(T) == 16) return (T) _mm_max_epi16((__m128i) a, (__m128i) b);
#if defined(__AVX2__)
    else return (T) _mm256_max_epi16((__m256i) a, (__m256i) b);
#else
    else return sseBinary(a, b, [](__m128i x, __m128i y) { return _mm_max_epi16(x, y); });
#endif
#endif
#ifdef MARV_SIMD_ELEMENTWISE_BUILTINS
    return __builtin_elementwise_max(a, b);
#else
    T m = a < b; return (a & ~m) | (b & m);
#endif
}

// lane shift by one; lane 0 becomes 0
inline V<16>::u8 shift1(V<16>::u8 v) { const V<16>::u8 z = {}; return __builtin_shufflevector(z, v, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30); }
inline V<16>::i16 shift1(V<16>::i16 v) { const V<16>::i16 z = {}; return __builtin_shufflevector(z, v, 7, 8, 9, 10, 11, 12, 13, 14); }
inline V<32>::u8 shift1(V<32>::u8 v) {
    const V<32>::u8 z = {};
    return __builtin_shufflevector(z, v, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46,
                                         47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62);
}
inline V<32>::i16 shift1(V<32>::i16 v) { const V<32>::i16 z = {}; return __builtin_shufflevector(z, v, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30); }

#if defined(__AVX512BW__)
inline V<64>::u8 shift1(V<64>::u8 v) { const V<64>::u8 z = {}; return __builtin_shufflevector(z, v, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126); }
inline V<64>::i16 shift1(V<64>::i16 v) { const V<64>::i16 z = {}; return __builtin_shufflevector(z, v, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62); }
#endif

template<class T> inline int hmax(T v) { int m = 0; for (unsigned i = 0; i < sizeof(T) / sizeof(v[0]); i++) m = std::max(m, (int) v[i]); return m; }

enum class Precision { U8, I16, I32 };

// TR: translate stored bytes through a 256-entry table before the PSSM lookup
template<int W, bool TR, int Segments = 0>
int score_u8(const int8_t* subject, int length, const typename V<W>::u8* profile, int segLen, uint8_t bias, int alphabetSize,
             const unsigned char* table, typename V<W>::u8* H) {
    typedef typename V<W>::u8 Vec;
    // Segments != 0: stripe held in registers
    Vec local[Segments ? Segments : 1];
    if constexpr (Segments != 0) { segLen = Segments; H = local; }
    for (int j = 0; j < segLen; j++) H[j] = Vec{};
    const Vec vBias = splat_u8<W>(bias);
    Vec vMax = {};
    for (int i = 0; i < length; i++) {
        int letter = TR ? (int) table[(uint8_t) subject[i]] : (int)(uint8_t) subject[i];
        if (letter >= alphabetSize) letter = alphabetSize - 1;
        const Vec* P = profile + (size_t) letter * segLen;
        Vec vH = shift1(H[segLen - 1]);
        for (int j = 0; j < segLen; j++) {
            vH = subs_u8(adds_u8(vH, P[j]), vBias);
            vMax = max_u8(vMax, vH);
            Vec prev = H[j]; H[j] = vH; vH = prev;
        }
    }
    return hmax(vMax);
}

template<int W, bool TR, int Segments = 0>
// keeps GCC from placing the generic loop across an instruction boundary
__attribute__((aligned(64))) int score_i16(const int8_t* subject, int length, const typename V<W>::i16* profile, int segLen, int alphabetSize,
              const unsigned char* table, typename V<W>::i16* H) {
    typedef typename V<W>::i16 Vec;
    // Segments != 0: stripe held in registers
    Vec local[Segments ? Segments : 1];
    if constexpr (Segments != 0) { segLen = Segments; H = local; }
    for (int j = 0; j < segLen; j++) H[j] = Vec{};
    const Vec zero = {};
    Vec vMax = {};
    for (int i = 0; i < length; i++) {
        int letter = TR ? (int) table[(uint8_t) subject[i]] : (int)(uint8_t) subject[i];
        if (letter >= alphabetSize) letter = alphabetSize - 1;
        const Vec* P = profile + (size_t) letter * segLen;
        Vec vH = shift1(H[segLen - 1]);
        for (int j = 0; j < segLen; j++) {
            vH = max_i16(adds_i16(vH, P[j]), zero);
            vMax = max_i16(vMax, vH);
            Vec prev = H[j]; H[j] = vH; vH = prev;
        }
    }
    return hmax(vMax);
}

template<bool TR>
int score_i32(const int8_t* subject, int length, const int8_t* pssm, int queryLength, int alphabetSize, const unsigned char* table, int* H) {
    std::fill(H, H + queryLength, 0);
    int best = 0;
    for (int i = 0; i < length; i++) {
        int letter = TR ? (int) table[(uint8_t) subject[i]] : (int)(uint8_t) subject[i];
        if (letter >= alphabetSize) letter = alphabetSize - 1;
        const int8_t* row = pssm + (size_t) letter * queryLength;
        int diag = 0;
        for (int c = 0; c < queryLength; c++) {
            int prev = H[c];
            int v = diag + row[c];
            H[c] = v > 0 ? v : 0;
            if (H[c] > best) best = H[c];
            diag = prev;
        }
    }
    return best;
}

int forcedWidth() {
    if (const char* e = getenv("MARV_SIMD_WIDTH")) {
        const int w = atoi(e);
        if (w == 16 || w == 32) return w;
#if defined(__AVX512BW__)
        if (w == 64) return w;
#endif
        throw std::runtime_error("libmarv/SIMD: MARV_SIMD_WIDTH must be 16 or 32 (64 requires an AVX-512BW build)");
    }
    return 0;
}

int widthForQuery(int queryLength, Precision precision, int forced) {
    if (forced) return forced;
#if defined(__AVX512BW__)
    if (queryLength > (precision == Precision::I16 ? 256 : 512)) return 64;
    return queryLength <= (precision == Precision::I16 ? 16 : 32) ? 16 : 32;
#elif defined(__AVX2__)
    return queryLength <= (precision == Precision::I16 ? 16 : 32) ? 16 : 32;
#elif defined(__SSE2__)
    (void) queryLength;
    (void) precision;
    return 16;
#elif defined(__aarch64__)
    // measured on M4
    const int crossover = precision == Precision::I16 ? 24 : 48;
    return queryLength > crossover ? 32 : 16;
#else
    (void) precision;
    return queryLength >= 160 ? 32 : 16;
#endif
}

template<class F> void dispatchStripe(int segments, F&& f) {
    switch (segments) {
        case 1: f.template operator()<1>(); break;
        case 2: f.template operator()<2>(); break;
        case 3: f.template operator()<3>(); break;
        case 4: f.template operator()<4>(); break;
        case 5: f.template operator()<5>(); break;
        case 6: f.template operator()<6>(); break;
        case 7: f.template operator()<7>(); break;
        case 8: f.template operator()<8>(); break;
        case 9: f.template operator()<9>(); break;
        case 10: f.template operator()<10>(); break;
        case 11: f.template operator()<11>(); break;
        case 12: f.template operator()<12>(); break;
        case 13: f.template operator()<13>(); break;
        case 14: f.template operator()<14>(); break;
        case 15: f.template operator()<15>(); break;
        case 16: f.template operator()<16>(); break;
        default: f.template operator()<0>(); break;
    }
}

class SimdGaplessAligner final : public MultiConfigScoreOnlyAlignerInterface {
public:
    explicit SimdGaplessAligner(Precision p) : precision(p), forced(forcedWidth()) {}

    bool isSingleTile(const OneToAllInputDataPSSM&) const override { return true; }
    int getTileSize(const OneToAllInputDataPSSM& in) const override { return in.queryLength; }
    int preferredKernelOutputArraySize() const override {
#if defined(__x86_64__)
        return 1'000'000;
#else
        return 10'000'000;
#endif
    }
    AlignmentAlgorithmE getAlignmentAlgorithm() const override { return AlignmentAlgorithmE::Gapless; }
    Datatype getScoreDatatype(int) const override {
        return precision == Precision::U8 ? Datatype::UInt8x4 : precision == Precision::I16 ? Datatype::Short2 : Datatype::Int;
    }
    size_t getMinimumSuggestedTempBytes_multiTile(const OneToAllInputDataPSSM&) const override { return 0; }

    // striped profile: lane l of segment j holds query column l*segLen + j; padding columns score 0
    void makeGpuPssm(GpuConvertedPSSM& result, const cuda::std::mdspan<std::int8_t, cuda::std::dextents<int, 2>> view, cudaStream_t) override {
        const int rows = view.extent(0), cols = view.extent(1);
        result.kindOfPssm = GpuConvertedPSSM::Kind::Gapless;
        result.alphabetSize = rows;
        result.queryLength = cols;
        result.allocatedNumRegs = 0;
        int minScore = 0;
        for (int r = 0; r < rows; r++) for (int c = 0; c < cols; c++) minScore = std::min(minScore, (int) view(r, c));
        const int width = widthForQuery(cols, precision, forced);
        if (precision == Precision::U8) {
            const int lanes = width, segLen = (cols + lanes - 1) / lanes;
            const uint8_t bias = (uint8_t)(-minScore);
            result.groupsize = lanes; result.numRegs = segLen; result.columnstride = lanes * segLen;
            result.substitutionScoreBias_optional = bias;
            result.data.resize((size_t) rows * segLen * lanes);
            uint8_t* out = (uint8_t*) result.data.data();
            for (int r = 0; r < rows; r++) for (int j = 0; j < segLen; j++) for (int l = 0; l < lanes; l++) {
                const int c = l * segLen + j;
                out[((size_t) r * segLen + j) * lanes + l] = (uint8_t)((c < cols ? (int) view(r, c) : 0) + bias);
            }
        } else if (precision == Precision::I16) {
            const int lanes = width / 2, segLen = (cols + lanes - 1) / lanes;
            result.groupsize = lanes; result.numRegs = segLen; result.columnstride = lanes * segLen;
            result.data.resize((size_t) rows * segLen * lanes * 2);
            int16_t* out = (int16_t*) result.data.data();
            for (int r = 0; r < rows; r++) for (int j = 0; j < segLen; j++) for (int l = 0; l < lanes; l++) {
                const int c = l * segLen + j;
                out[((size_t) r * segLen + j) * lanes + l] = (int16_t)(c < cols ? (int) view(r, c) : 0);
            }
        } else {
            result.groupsize = 1; result.numRegs = cols; result.columnstride = cols;
            result.data.resize((size_t) rows * cols);
            memcpy(result.data.data(), view.data_handle(), (size_t) rows * cols);
        }
    }

    void scoreOnly_singleTile(int* const scores, const OneToAllInputDataPSSM& in, const GpuConvertedPSSM& pssm, GapScoreArgs, cudaStream_t) override {
        run(scores, in, pssm);
    }
    void scoreOnly_multiTile(char*, size_t, int* const scores, const OneToAllInputDataPSSM& in, const GpuConvertedPSSM& pssm, GapScoreArgs, cudaStream_t) override {
        run(scores, in, pssm);
    }

    // every worker applies all queries to one block of subjects while it is in cache
    bool scoreOnly_batch(cuda::std::span<const ScoreOnlyBatchQuery> queries, GapScoreArgs, cudaStream_t) override {
        if (queries.size() < 2) return false;
        const int blockSubjects = 256;
        const int n = queries.front().input.numAlignments;
        for (const auto& q : queries)
            if (q.input.numAlignments != n || q.input.queryLength == 0) return false;
        std::vector<RangeWorker> workers;
        workers.reserve(queries.size());
        for (const auto& q : queries) workers.push_back(batchWorker(q));
        const int blocks = n / blockSubjects + (n % blockSubjects != 0);
        #pragma omp parallel for schedule(dynamic, 1)
        for (int block = 0; block < blocks; ++block) {
            const int begin = block * blockSubjects;
            const int end = begin + std::min(blockSubjects, n - begin);
            for (const auto& worker : workers) worker(begin, end);
        }
        return true;
    }

    using RangeWorker = std::function<void(int, int)>;
    RangeWorker batchWorker(const ScoreOnlyBatchQuery& query) {
        const auto& p = *query.pssm;
        const auto& in = query.input;
        if (p.kindOfPssm != GpuConvertedPSSM::Kind::Gapless || p.queryLength != in.queryLength)
            throw std::runtime_error("libmarv/SIMD: PSSM does not match the query");
        if (precision == Precision::I32) {
            return [query](int begin, int end) {
                const auto& in = query.input;
                std::vector<int> H(in.queryLength);
                const auto* raw = reinterpret_cast<const int8_t*>(query.pssm->data.data());
                for (int i = begin; i < end; ++i) {
                    if (in.letterTranslation)
                        query.scores[i] = score_i32<true>(in.getSubject(i), in.getSubjectLength(i), raw,
                            in.queryLength, query.pssm->alphabetSize, in.letterTranslation, H.data());
                    else
                        query.scores[i] = score_i32<false>(in.getSubject(i), in.getSubjectLength(i), raw,
                            in.queryLength, query.pssm->alphabetSize, nullptr, H.data());
                }
            };
        }
        const int width = precision == Precision::U8 ? p.groupsize : 2 * p.groupsize;
        RangeWorker result;
        auto bindWidth = [&]<int W>() {
            if (in.letterTranslation) result = typedBatchWorker<W, true>(query);
            else result = typedBatchWorker<W, false>(query);
        };
#if defined(__AVX512BW__)
        if (width == 64) bindWidth.template operator()<64>();
        else
#endif
        if (width == 32) bindWidth.template operator()<32>();
        else bindWidth.template operator()<16>();
        return result;
    }

private:
    template<int W, bool TR>
    RangeWorker typedBatchWorker(const ScoreOnlyBatchQuery& query) {
        RangeWorker result;
        dispatchStripe(query.pssm->numRegs, [&]<int S>() {
            if (precision == Precision::U8) {
                result = [query](int begin, int end) {
                    using Vec = typename V<W>::u8;
                    const auto& p = *query.pssm;
                    const auto& in = query.input;
                    std::vector<Vec> H(S ? 0 : p.numRegs);
                    const auto* profile = reinterpret_cast<const Vec*>(p.data.data());
                    for (int i = begin; i < end; ++i)
                        query.scores[i] = score_u8<W, TR, S>(in.getSubject(i), in.getSubjectLength(i), profile,
                            p.numRegs, uint8_t(p.substitutionScoreBias_optional.value_or(0)), p.alphabetSize,
                            in.letterTranslation, H.data());
                };
            } else {
                result = [query](int begin, int end) {
                    using Vec = typename V<W>::i16;
                    const auto& p = *query.pssm;
                    const auto& in = query.input;
                    std::vector<Vec> H(S ? 0 : p.numRegs);
                    const auto* profile = reinterpret_cast<const Vec*>(p.data.data());
                    for (int i = begin; i < end; ++i)
                        query.scores[i] = score_i16<W, TR, S>(in.getSubject(i), in.getSubjectLength(i), profile,
                            p.numRegs, p.alphabetSize, in.letterTranslation, H.data());
                };
            }
        });
        return result;
    }
    Precision precision;
    int forced;

    void run(int* const scores, const OneToAllInputDataPSSM& in, const GpuConvertedPSSM& pssm) {
        if (in.numAlignments <= 0 || in.queryLength == 0) return;
        if (pssm.kindOfPssm != GpuConvertedPSSM::Kind::Gapless || pssm.queryLength != in.queryLength) throw std::runtime_error("libmarv/SIMD: PSSM does not match the query");
        const unsigned char* table = in.letterTranslation;
        const int width = precision == Precision::U8 ? pssm.groupsize : precision == Precision::I16 ? 2 * pssm.groupsize : 16;
#if defined(__AVX512BW__)
        if (width == 64) { if (table) runAll<64, true>(scores, in, pssm, table); else runAll<64, false>(scores, in, pssm, table); }
        else
#endif
        if (width == 32) { if (table) runAll<32, true>(scores, in, pssm, table); else runAll<32, false>(scores, in, pssm, table); }
        else             { if (table) runAll<16, true>(scores, in, pssm, table); else runAll<16, false>(scores, in, pssm, table); }
    }

    template<int W, bool TR>
    void runAll(int* const scores, const OneToAllInputDataPSSM& in, const GpuConvertedPSSM& pssm, const unsigned char* table) {
        const int n = in.numAlignments;
        const int L = in.queryLength;
        const int alphabetSize = pssm.alphabetSize;
        const int segLen = pssm.numRegs;
        if (precision == Precision::U8) {
            typedef typename V<W>::u8 Vec;
            const Vec* profile = (const Vec*) pssm.data.data();
            const uint8_t bias = (uint8_t) pssm.substitutionScoreBias_optional.value_or(0);
            auto run = [&]<int S>() {
                #pragma omp parallel
                {
                    std::vector<Vec> H(S ? 0 : segLen);
                    #pragma omp for schedule(dynamic, 64)
                    for (int i = 0; i < n; i++) scores[i] = score_u8<W, TR, S>(in.getSubject(i), in.getSubjectLength(i), profile, segLen, bias, alphabetSize, table, H.data());
                }
            };
            dispatchStripe(segLen, run);
        } else if (precision == Precision::I16) {
            typedef typename V<W>::i16 Vec;
            const Vec* profile = (const Vec*) pssm.data.data();
            auto run = [&]<int S>() {
                #pragma omp parallel
                {
                    std::vector<Vec> H(S ? 0 : segLen);
                    #pragma omp for schedule(dynamic, 64)
                    for (int i = 0; i < n; i++) scores[i] = score_i16<W, TR, S>(in.getSubject(i), in.getSubjectLength(i), profile, segLen, alphabetSize, table, H.data());
                }
            };
            dispatchStripe(segLen, run);
        } else {
            const int8_t* raw = (const int8_t*) pssm.data.data();
            #pragma omp parallel
            {
                std::vector<int> H(L);
                #pragma omp for schedule(dynamic, 64)
                for (int i = 0; i < n; i++) scores[i] = score_i32<TR>(in.getSubject(i), in.getSubjectLength(i), raw, L, alphabetSize, table, H.data());
            }
        }
    }
};

[[noreturn]] void notOnSimd(const char* what) {
    throw std::runtime_error(std::string("libmarv/SIMD: ") + what + " is not implemented on the SIMD backend");
}

} // namespace

std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_8bit() { return std::make_unique<SimdGaplessAligner>(Precision::U8); }
std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_16bit() { return std::make_unique<SimdGaplessAligner>(Precision::I16); }
std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_32bit() { return std::make_unique<SimdGaplessAligner>(Precision::I32); }
std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_sw_aligner_32bit() { notOnSimd("Smith-Waterman"); }
std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_16bit() { notOnSimd("gapless end-position alignment"); }
std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_32bit() { notOnSimd("gapless end-position alignment"); }
std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_sw_endpos_aligner_32bit() { notOnSimd("Smith-Waterman end-position alignment"); }

} // namespace libmarv
