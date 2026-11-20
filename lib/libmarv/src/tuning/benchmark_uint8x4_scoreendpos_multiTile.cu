
#include "benchmarkdata.hpp"
#include "common.hpp"
#include "configtuples.hpp"
#include "benchmark_scoreendpos.cuh"

#include <vector>

#include "../util.cuh"

std::vector<BenchmarkData> benchmark_uint8x4_scoreEndpos_multiTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults,
    int tiles
){
#if 1

    return {};

#else

#if defined(__CUDACC__)
    using ApproachAndTypeTuple = std::tuple<
        libmarv::ApproachAndType<Approach::unused, Datatype::UInt8x4>
    >;

    assert(tiles > 1);
    if(tiles <= 1){
        throw std::runtime_error("tiles <= 1");
    }
    constexpr bool singleTile = false;

    using BlockSizeTupleFoo = std::tuple<IC<128>>;
    //using GroupSizeTupleFoo = std::tuple<IC<4>, IC<8>, IC<16>, IC<32>>;
    using GroupSizeTupleFoo = std::tuple<IC<4>, IC<8>>;
    using NumItemsTupleFoo = std::tuple<IC<1>, IC<2>, IC<3>, IC<4>, IC<5>, IC<6>, IC<7>, IC<8>, IC<10>, IC<12>, IC<14>, IC<16>, IC<18>>;
    // using NumItemsTupleFoo = std::tuple<IC<1>>;

    // using NumItemsTupleFoo = std::tuple<IC<4>, IC<8>, IC<12>, IC<16>>;

    std::vector<BenchmarkData> vec = score_endpos::ForEachCombination_call_singleconfigFunction<
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
#else
    return {};
#endif

#endif
}