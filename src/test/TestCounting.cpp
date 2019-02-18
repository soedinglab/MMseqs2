#include "Util.h"

#include <iostream>
#include <chrono>

const char* binary_name = "test_counting";

//u_int64_t revcomp64_v2 (const u_int64_t& x, size_t sizeKmer)
//{
//    u_int64_t res = x;
//
//    res = ((res>> 2 & 0x3333333333333333) | (res & 0x3333333333333333) <<  2);
//    res = ((res>> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) <<  4);
//    res = ((res>> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) <<  8);
//    res = ((res>>16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
//    res = ((res>>32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
//    res = res ^ 0xAAAAAAAAAAAAAAAA;
//
//    return (res >> (2*( 32 - sizeKmer))) ;
//}

int main(int, char **){
    auto start = std::chrono::system_clock::now();
    size_t revComp = 0;
    for(u_int64_t i = 0; i < 10001010101011; i++){
        revComp += Util::revComplement(i, 21);
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = end - start;
    std::cout << revComp << "\t" << elapsed.count() << '\n';
    return 0;
}
