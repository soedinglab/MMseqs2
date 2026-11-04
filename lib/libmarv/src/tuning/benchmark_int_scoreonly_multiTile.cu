
#include "benchmarkdata.hpp"
#include "common.hpp"
#include "configtuples.hpp"
#include "benchmark_scoreonly.cuh"

#include <vector>

#include "../util.cuh"

std::vector<BenchmarkData> benchmark_int_scoreOnly_multiTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults,
    int tiles
){
    using ApproachAndTypeTuple = std::tuple<
        libmarv::ApproachAndType<Approach::unused, Datatype::Int>        
    >;

    assert(tiles > 1);
    if(tiles <= 1){
        throw std::runtime_error("tiles <= 1");
    }
    constexpr bool singleTile = false;

    using BlockSizeTupleFoo = std::tuple<IC<128>>;
    //using GroupSizeTupleFoo = std::tuple<IC<4>, IC<8>, IC<16>, IC<32>>;
    using GroupSizeTupleFoo = std::tuple<IC<32>>;
    using NumItemsTupleFoo = std::tuple<IC<24>>;
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
    >{}(timingIterations, sleepAmount, checkResults, tiles);

    return vec;
}