
#include "benchmarkdata.hpp"
#include "common.hpp"
#include "configtuples.hpp"
#include "benchmark_scoreonly.cuh"

#include <vector>

#include "../util.cuh"

std::vector<BenchmarkData> benchmark_float_scoreOnly_singleTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults
){
    using ApproachAndTypeTuple = std::tuple<
        libmarv::ApproachAndType<Approach::unused, Datatype::Float>        
    >;

    constexpr bool singleTile = true;

    using BlockSizeTupleFoo = std::tuple<IC<128>, IC<256>, IC<512>>;
    //using GroupSizeTupleFoo = std::tuple<IC<4>, IC<8>, IC<16>, IC<32>>;
    using GroupSizeTupleFoo = std::tuple<IC<4>>;
    using NumItemsTupleFoo = std::tuple<IC<16>, IC<32>, IC<35>, IC<64>>;
    // using NumItemsTupleFoo = std::tuple<IC<2>>;

    std::vector<BenchmarkData> vec = scoreonly::ForEachCombination_call_singleconfigFunction<  
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