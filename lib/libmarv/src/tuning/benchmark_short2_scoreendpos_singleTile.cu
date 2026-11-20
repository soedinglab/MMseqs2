
#include "benchmarkdata.hpp"
#include "common.hpp"
#include "configtuples.hpp"
#include "benchmark_scoreendpos.cuh"

#include <vector>

#include "../util.cuh"

std::vector<BenchmarkData> benchmark_short2_scoreEndpos_singleTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults
){
#if defined(__CUDACC__)
    using ApproachAndTypeTuple = std::tuple<
        libmarv::ApproachAndType<Approach::hardcodedzero, Datatype::Short2>,
        libmarv::ApproachAndType<Approach::kernelparamzero, Datatype::Short2>
    >;

    constexpr bool singleTile = true;

    std::vector<BenchmarkData> vec = score_endpos::ForEachCombination_call_singleconfigFunction<
        BlockSizeTuple,
        GroupSizeTuple,
        NumItemsTuple,
        ApproachAndTypeTuple,
        singleTile
    >{}(timingIterations, sleepAmount, checkResults, 1);

    return vec;
#else
    return {};
#endif
}