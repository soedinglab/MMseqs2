#include <iostream>
#include <list>
#include <algorithm>
#include <math.h>
#include <CacheFriendlyOperations.h>
#include <map>

const char* binary_name = "test_counting";

#include <iostream>
#include <chrono>
#include <ctime>
#include <immintrin.h> // AVX

void print128_num(__m128i var)
{
    uint8_t *val = (uint8_t*) &var;
    printf("Numerical: %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i \n",
           val[0], val[1], val[2], val[3], val[4], val[5],
           val[6], val[7], val[8], val[9], val[10], val[11],
           val[12], val[13], val[14], val[15]);
}
// Compute reverse complement of k-mer in 2-bit-per-nucleotide encoding (A: 00, C: 01, T: 10, G: 11)
uint64_t revComplement(const u_int64_t kmer, const int k) {
    // broadcast 64bit to 128 bit
    __m128i x = _mm_cvtsi64_si128(kmer);

    // create lookup (set 16 bytes in 128 bit)
    // a lookup entry at the index of two nucleotides (4 bit) describes the reverse
    // complement of these two nucleotides in the higher 4 bits (lookup1) or in the
    // lower 4 bits (lookup2)
    __m128i lookup1 = _mm_set_epi8(0x50,0x10,0xD0,0x90,0x40,0x00,0xC0,0x80,0x70,
                                   0x30,0xF0,0xB0,0x60,0x20,0xE0,0xA0);
    __m128i lookup2 = _mm_set_epi8(0x05,0x01,0x0D,0x09,0x04,0x00,0x0C,0x08,0x07,
                                   0x03,0x0F,0x0B,0x06,0x02,0x0E,0x0A);
    // set upper 8 bytes to 0 and revert order of lower 8 bytes
    __m128i upper = _mm_set_epi8(0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0,1,2,3,4,5,6,7);

    __m128i kmer1 = _mm_and_si128(x, _mm_set1_epi8(0x0F)); // get lower 4 bits
    __m128i kmer2 = _mm_and_si128(x, _mm_set1_epi8(0xF0)); // get higher 4 bits

//    // shift right by 2 nucleotides
//    print128_num(kmer2);
    kmer2 >>= 4;
//    kmer2 = _mm_srli_epi64(kmer2, 4);
//
//    print128_num(kmer2);
//    std::cout << kmer2 << std::endl;
    // use _mm_shuffle_epi8 to look up reverse complement
    kmer1 =_mm_shuffle_epi8(lookup1, kmer1);
    kmer2 = _mm_shuffle_epi8(lookup2, kmer2);

    // _mm_or_si128: bitwise OR
    x = _mm_or_si128(kmer1, kmer2);

    // set upper 8 bytes to 0 and revert order of lower 8 bytes
    x = _mm_shuffle_epi8(x, upper);

    // shift out the unused nucleotide positions (1 <= k <=32 )
    // broadcast 128 bit to 64 bit
    return (((uint64_t)_mm_cvtsi128_si64(x)) >> (uint64_t)(64-2*k));
}

u_int64_t revcomp64_v2 (const u_int64_t& x, size_t sizeKmer)
{
    u_int64_t res = x;

    res = ((res>> 2 & 0x3333333333333333) | (res & 0x3333333333333333) <<  2);
    res = ((res>> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) <<  4);
    res = ((res>> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) <<  8);
    res = ((res>>16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
    res = ((res>>32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
    res = res ^ 0xAAAAAAAAAAAAAAAA;

    return (res >> (2*( 32 - sizeKmer))) ;
}

int main(int argc, char ** avgv){
    auto start = std::chrono::system_clock::now();
    size_t revComp = 0;
    for(u_int64_t i = 0; i < 10001010101011; i++){
        revComp += revComplement(i, 21);
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = end - start;
    std::cout << revComp << "\t" << elapsed.count() << '\n';
    return 0;
}
