// The threaded staging copy must cover remainders and OpenMP teams smaller than requested.
#include "dbbatching.cuh"
#include <cstdio>
#include <numeric>

int main() {
    using namespace libmarv;
    for (int bytes : {0, 1, 7, 32, 103}) {
        std::vector<char> chars(bytes + 8);
        std::iota(chars.begin(), chars.end(), char(1));
        const std::vector<int> lengths{3, bytes, 5};
        const std::vector<size_t> offsets{0, 3, size_t(bytes + 3), size_t(bytes + 8)};
        std::vector<DBdataView> partitions;
        partitions.emplace_back(0, 3, 0, chars.data(), lengths.data(), offsets.data(), nullptr, nullptr);
        DeviceBatchCopyToPinnedPlan plan;
        plan.copyRanges = {{0, 0, 1, 1}, {0, 0, 0, 1}};
        std::vector<char> expected(bytes + 5, char(-1));
        std::vector<int> expectedLengths(4, -1);
        std::vector<size_t> expectedOffsets(5, size_t(-1));
        executePinnedCopyPlanSerial(plan, expected.data() + 1, expectedLengths.data() + 1,
                                   expectedOffsets.data() + 1, partitions);
        for (int threads : {1, 3, 4, 8}) {
            std::vector<char> actual(bytes + 5, char(-1));
            std::vector<int> actualLengths(4, -1);
            std::vector<size_t> actualOffsets(5, size_t(-1));
            executePinnedCopyPlanParallel(plan, actual.data() + 1, actualLengths.data() + 1,
                                         actualOffsets.data() + 1, partitions, threads);
            if (actual != expected || actualLengths != expectedLengths || actualOffsets != expectedOffsets) {
                fprintf(stderr, "staging copy differs: bytes=%d threads=%d\n", bytes, threads);
                return 1;
            }
        }
    }
    puts("database staging copy checks passed");
}
