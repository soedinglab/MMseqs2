
#include "benchmarkdata.hpp"
#include "common.hpp"
#include "configtuples.hpp"
#include "benchmark_scoreendpos.cuh"

#include <vector>

#include "../util.cuh"

std::vector<BenchmarkData> benchmark_half2_scoreEndpos_singleTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults
){
    using ApproachAndTypeTuple = std::tuple<
        libmarv::ApproachAndType<Approach::hardcodedzero, Datatype::Half2>,
        libmarv::ApproachAndType<Approach::kernelparamzero, Datatype::Half2>
        // libmarv::ApproachAndType<Approach::half2fmarelu, Datatype::Half2>        
    >;

    constexpr bool singleTile = true;

    using BlockSizeTupleFoo = std::tuple<IC<128>>;
    //using GroupSizeTupleFoo = std::tuple<IC<4>, IC<8>, IC<16>, IC<32>>;
    using GroupSizeTupleFoo = std::tuple<IC<4>, IC<8>>;
    using NumItemsTupleFoo = std::tuple<IC<1>, IC<2>, IC<3>, IC<4>, IC<5>, IC<6>, IC<7>, IC<8>, IC<10>, IC<12>, IC<14>, IC<16>, IC<18>>;
    // using NumItemsTupleFoo = std::tuple<IC<2>>;

    std::vector<BenchmarkData> vec = score_endpos::ForEachCombination_call_singleconfigFunction<
        BlockSizeTuple,
        GroupSizeTuple,
        NumItemsTuple,
        // BlockSizeTupleFoo,
        // GroupSizeTupleFoo,
        // NumItemsTupleFoo,
        ApproachAndTypeTuple,
        singleTile
    >{}(timingIterations, sleepAmount, checkResults, 1);

    return vec;
}