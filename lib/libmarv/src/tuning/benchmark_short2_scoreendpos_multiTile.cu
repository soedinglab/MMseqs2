
#include "benchmarkdata.hpp"
#include "common.hpp"
#include "configtuples.hpp"
#include "benchmark_scoreendpos.cuh"

#include <vector>

#include "../util.cuh"

std::vector<BenchmarkData> benchmark_short2_scoreEndpos_multiTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults,
    int tiles
){
#if defined(__CUDACC__)
    using ApproachAndTypeTuple = std::tuple<
        libmarv::ApproachAndType<Approach::hardcodedzero, Datatype::Short2>,
        libmarv::ApproachAndType<Approach::kernelparamzero, Datatype::Short2>
    >;

    assert(tiles > 1);
    if(tiles <= 1){
        throw std::runtime_error("tiles <= 1");
    }
    constexpr bool singleTile = false;

    std::vector<BenchmarkData> vec = score_endpos::ForEachCombination_call_singleconfigFunction<
        BlockSizeTuple,
        GroupSizeTuple,
        NumItemsTuple,
        ApproachAndTypeTuple,
        singleTile
    >{}(timingIterations, sleepAmount, checkResults, tiles);

    return vec;
#else
    return {};
#endif
}