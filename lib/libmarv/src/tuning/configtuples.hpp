
#ifndef CONFIGTUPLES_HPP
#define CONFIGTUPLES_HPP

#include <utility>

template<int N>
using IC = std::integral_constant<int, N>;

using BlockSizeTuple = std::tuple<IC<128>, IC<256>, IC<512>>;
using GroupSizeTuple = std::tuple<IC<2>, IC<4>, IC<8>, IC<16>, IC<32>>;
using NumItemsTuple = std::tuple<IC<4>, IC<8>, IC<12>, IC<16>, IC<20>, IC<24>, IC<28>, IC<32>,
                                IC<36>, IC<40>, IC<44>, IC<48>, IC<52>, IC<56>, IC<60>, IC<64>>;


// using BlockSizeTuple = std::tuple<IC<128>, IC<256>>;
// using GroupSizeTuple = std::tuple<IC<4>, IC<8>>;
// using NumItemsTuple = std::tuple<IC<32>>;


#endif