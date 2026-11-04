// Check every score while limiting scratch to one block: this forces repeated
// grid-stride iterations and partial final blocks, independent of GPU size.
#include "aligners/make_aligner_headers.cuh"
#include <algorithm>
#include <cstdio>
#include <cstdint>
#include <stdexcept>
#include <vector>
using namespace libmarv;

static void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

static void run(int bits, int queryLength) {
    constexpr int count = 521;
    std::vector<int> lengths(count);
    std::vector<size_t> offsets(count + 1);
    std::vector<int8_t> subjects;
    uint32_t rng = 1234567;
    auto random = [&]() { rng = rng * 1664525u + 1013904223u; return rng >> 16; };
    for (int s = 0; s < count; ++s) {
        offsets[s] = subjects.size();
        lengths[s] = 17 * s / (count - 1);
        for (int j = 0; j < lengths[s]; ++j) subjects.push_back(random() % 21);
        while (subjects.size() % 4) subjects.push_back(20);
    }
    offsets[count] = subjects.size();
    subjects.resize(subjects.size() + 16, 20);
    std::vector<int8_t> pssm(21 * queryLength);
    for (auto& score : pssm) score = int(random() % 11) - 5;
    // Short subjects keep scores below every direct kernel's saturation limit.
    std::vector<int> expected(count);
    std::vector<int> previous(queryLength + 1), current(queryLength + 1);
    for (int s = 0; s < count; ++s) {
        std::fill(previous.begin(), previous.end(), 0);
        for (int j = 0; j < lengths[s]; ++j) {
            current[0] = 0;
            const int letter = subjects[offsets[s] + j];
            for (int k = 1; k <= queryLength; ++k) {
                current[k] = std::max(0, previous[k - 1] + pssm[letter * queryLength + k - 1]);
                expected[s] = std::max(expected[s], current[k]);
            }
            previous.swap(current);
        }
    }
    helpers::SimpleAllocationDevice<int8_t, 0> d_subjects(subjects.size());
    helpers::SimpleAllocationDevice<size_t, 0> d_offsets(count);
    helpers::SimpleAllocationDevice<int, 0> d_lengths(count), d_scores(count + 1);
    CUDACHECK(cudaMemcpy(d_subjects.data(), subjects.data(), subjects.size(), cudaMemcpyHostToDevice));
    CUDACHECK(cudaMemcpy(d_offsets.data(), offsets.data(), count * sizeof(size_t), cudaMemcpyHostToDevice));
    CUDACHECK(cudaMemcpy(d_lengths.data(), lengths.data(), count * sizeof(int), cudaMemcpyHostToDevice));
    auto aligner = bits == 8 ? make_multiconfig_gapless_aligner_8bit()
        : bits == 16 ? make_multiconfig_gapless_aligner_16bit() : make_multiconfig_gapless_aligner_32bit();
    GpuConvertedPSSM profile;
    aligner->makeGpuPssm(profile, cuda::std::mdspan<int8_t, cuda::std::dextents<int, 2>>(pssm.data(), 21, queryLength), nullptr);
    OneToAllInputDataPSSM in;
    in.d_subjects = d_subjects.data();
    in.d_subjectOffsets = makeCustomOffsetIterator(d_offsets.data(), false);
    in.d_subjectLengths = d_lengths.data();
    in.queryLength = queryLength;
    in.maximumSubjectLength = lengths.back();
    in.numAlignments = 1;
    check(!aligner->isSingleTile(in), "test query must require multiple tiles");
    const size_t scratchBytes = aligner->getMinimumSuggestedTempBytes_multiTile(in);
    helpers::SimpleAllocationDevice<char, 0> scratch(scratchBytes);
    for (int n : {count, 259, count}) {
        in.numAlignments = n;
        CUDACHECK(cudaMemset(d_scores.data(), 0xff, (count + 1) * sizeof(int)));
        aligner->scoreOnly_multiTile(scratch.data(), scratchBytes, d_scores.data(), in, profile, {0, 0}, nullptr);
        std::vector<int> actual(count + 1);
        CUDACHECK(cudaMemcpy(actual.data(), d_scores.data(), actual.size() * sizeof(int), cudaMemcpyDeviceToHost));
        check(std::equal(expected.begin(), expected.begin() + n, actual.begin()), "multi-tile score differs from scalar reference");
        check(std::all_of(actual.begin() + n, actual.end(), [](int score) { return score == -1; }), "kernel wrote beyond alignment count");
    }
    std::printf("multi-tile bits=%d query=%d one_block_scratch=%zu: all scores match\n", bits, queryLength, scratchBytes);
}

int main() try {
    for (int bits : {8, 16, 32})
        for (int length : {8192, 8193}) run(bits, length);
    return 0;
} catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
}
