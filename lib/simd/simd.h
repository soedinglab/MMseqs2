// SIMD helper
// optimze based on technolegy double, float and integer (32) SIMD instructions
// writen by Martin Steinegger

#ifndef SIMD_H
#define SIMD_H
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <iostream>

#define AVX512_ALIGN_DOUBLE		64
#define AVX512_VECSIZE_DOUBLE	8
#define AVX512_ALIGN_FLOAT		64
#define AVX512_VECSIZE_FLOAT	16
#define AVX512_ALIGN_INT		64
#define AVX512_VECSIZE_INT		16

#define AVX_ALIGN_DOUBLE		32
#define AVX_VECSIZE_DOUBLE		4
#define AVX_ALIGN_FLOAT			32
#define AVX_VECSIZE_FLOAT		8
#define AVX2_ALIGN_INT			32
#define AVX2_VECSIZE_INT		8

#define SSE_ALIGN_DOUBLE		16
#define SSE_VECSIZE_DOUBLE		2
#define SSE_ALIGN_FLOAT			16
#define SSE_VECSIZE_FLOAT		4
#define SSE_ALIGN_INT			16
#define SSE_VECSIZE_INT			4

#define MAX_ALIGN_DOUBLE	AVX512_ALIGN_DOUBLE
#define MAX_VECSIZE_DOUBLE	AVX512_VECSIZE_DOUBLE
#define MAX_ALIGN_FLOAT		AVX512_ALIGN_FLOAT
#define MAX_VECSIZE_FLOAT	AVX512_VECSIZE_FLOAT
#define MAX_ALIGN_INT		AVX512_ALIGN_INT
#define MAX_VECSIZE_INT		AVX512_VECSIZE_INT

#define SIMDE_ENABLE_NATIVE_ALIASES
#include <simde/simde-features.h>

// FIXME: Finish AVX512 implementation
//#if defined(SIMDE_X86_AVX512F_NATIVE) && defined(SIMDE_X86_AVX512BW_NATIVE)
//#define AVX512
//#endif

#if defined(AVX512) || defined(SIMDE_X86_AVX2_NATIVE)
#define AVX2
#endif

#ifdef AVX512
#include <simde/x86/avx512f.h>
#include <simde/x86/avx512bw.h>

// double support
#ifndef SIMD_DOUBLE
#define SIMD_DOUBLE
#define ALIGN_DOUBLE        AVX512_ALIGN_DOUBLE
#define VECSIZE_DOUBLE      AVX512_VECSIZE_DOUBLE
typedef __m512d simd_double;
#define simdf64_add(x,y)    _mm512_add_pd(x,y)
#define simdf64_sub(x,y)    _mm512_sub_pd(x,y)
#define simdf64_mul(x,y)    _mm512_mul_pd(x,y)
#define simdf64_div(x,y)    _mm512_div_pd(x,y)
#define simdf64_max(x,y)    _mm512_max_pd(x,y)
#define simdf64_load(x)     _mm512_load_pd(x)
#define simdf64_store(x,y)  _mm512_store_pd(x,y)
#define simdf64_set(x)      _mm512_set1_pd(x)
#define simdf64_setzero(x)  _mm512_setzero_pd()
#define simdf64_gt(x,y)     _mm512_cmpnle_pd_mask(x,y)
#define simdf64_lt(x,y)     _mm512_cmplt_pd_mask(x,y)
#define simdf64_or(x,y)     _mm512_or_si512(x,y)
#define simdf64_and(x,y)    _mm512_and_si512 (x,y)
#define simdf64_andnot(x,y) _mm512_andnot_si512(x,y)
#define simdf64_xor(x,y)    _mm512_xor_si512(x,y)
#endif //SIMD_DOUBLE
// float support
#ifndef SIMD_FLOAT
#define SIMD_FLOAT
#define ALIGN_FLOAT         AVX512_ALIGN_FLOAT
#define VECSIZE_FLOAT       AVX512_VECSIZE_FLOAT
typedef __m512  simd_float;
#define simdf32_add(x,y)    _mm512_add_ps(x,y)
#define simdf32_sub(x,y)    _mm512_sub_ps(x,y)
#define simdf32_mul(x,y)    _mm512_mul_ps(x,y)
#define simdf32_div(x,y)    _mm512_div_ps(x,y)
#define simdf32_rcp(x)      _mm512_rcp_ps(x)
#define simdf32_max(x,y)    _mm512_max_ps(x,y)
#define simdf32_min(x,y)    _mm512_min_ps(x,y)
#define simdf32_load(x)     _mm512_load_ps(x)
#define simdf32_store(x,y)  _mm512_store_ps(x,y)
#define simdf32_set(x)      _mm512_set1_ps(x)
#define simdf32_setzero(x)  _mm512_setzero_ps()
#define simdf32_gt(x,y)     _mm512_cmpnle_ps_mask(x,y)
#define simdf32_eq(x,y)     _mm512_cmpeq_ps_mask(x,y)
#define simdf32_lt(x,y)     _mm512_cmplt_ps_mask(x,y)
#define simdf32_or(x,y)     _mm512_or_si512(x,y)
#define simdf32_and(x,y)    _mm512_and_si512(x,y)
#define simdf32_andnot(x,y) _mm512_andnot_si512(x,y)
#define simdf32_xor(x,y)    _mm512_xor_si512(x,y)
#define simdf32_f2i(x) 	    _mm512_cvtps_epi32(x)  // convert s.p. float to integer
#define simdf_f2icast(x)    _mm512_castps_si512(x)
#endif //SIMD_FLOAT
// integer support 
#ifndef SIMD_INT
#define SIMD_INT
#define ALIGN_INT           AVX512_ALIGN_INT
#define VECSIZE_INT         AVX512_VECSIZE_INT
typedef __m512i simd_int;
#define simdi32_add(x,y)    _mm512_add_epi32(x,y)
#define simdi16_add(x,y)    _mm512_add_epi16(x,y)
#define simdi16_adds(x,y)   _mm512_adds_epi16(x,y)
#define simdui8_adds(x,y)   _mm512_adds_epu8()
#define simdi32_sub(x,y)    _mm512_sub_epi32(x,y)
#define simdui8_subs(x,y)   _mm512_subs_epu8()
#define simdi32_mul(x,y)    _mm512_mullo_epi32(x,y)
#define simdui8_max(x,y)    _mm512_max_epu8()
#define simdi16_max(x,y)    _mm512_max_epi16(x,y)
#define simdi32_max(x,y)    _mm512_max_epi32(x,y)
#define simdi_load(x)       _mm512_load_si512(x)
#define simdi_streamload(x) _mm512_stream_load_si512(x)
#define simdi_store(x,y)    _mm512_store_si512(x,y)
#define simdi_storeu(x,y)   _mm512_storeu_si512(x,y)
#define simdi32_set(x)      _mm512_set1_epi32(x)
#define simdi16_set(x)      _mm512_set1_epi16(x)
#define simdi8_set(x)       _mm512_set1_epi8(x)
#define simdi32_shuffle(x,y) _mm512_shuffle_epi32(x,y)
#define simdi16_shuffle(x,y) NOT_YET_IMP(x,y)
#define simdi8_shuffle(x,y)  _mm512_shuffle_epi8(x,y)
#define simdi_setzero()     _mm512_setzero_si512()
#define simdi32_gt(x,y)     _mm512_cmpgt_epi32(x,y)
#define simdi8_gt(x,y)      NOT_YET_IMP()
#define simdi16_gt(x,y)     NOT_YET_IMP()
#define simdi8_eq(x,y)      NOT_YET_IMP()
#define simdi32_lt(x,y)     NOT_YET_IMP()
#define simdi16_lt(x,y)     NOT_YET_IMP()
#define simdi8_lt(x,y)      NOT_YET_IMP()

#define simdi_or(x,y)       _mm512_or_si512(x,y)
#define simdi_and(x,y)      _mm512_and_si512(x,y)
#define simdi_andnot(x,y)   _mm512_andnot_si512(x,y)
#define simdi_xor(x,y)      _mm512_xor_si512(x,y)
#define simdi8_shiftl(x,y)  NOT_YET_IMP()
#define simdi8_shiftr(x,y)  NOT_YET_IMP()
#define simdi8_movemask(x)  NOT_YET_IMP()
#define simdi16_extract(x,y) NOT_YET_IMP()
#define simdi16_slli(x,y)	_mm512_slli_epi16(x,y) // shift integers in a left by y
#define simdi16_srli(x,y)	_mm512_srli_epi16(x,y) // shift integers in a right by y
#define simdi32_slli(x,y)	_mm512_slli_epi32(x,y) // shift integers in a left by y
#define simdi32_srli(x,y)	_mm512_srli_epi32(x,y) // shift integers in a right by y
#define simdi32_i2f(x) 	    _mm512_cvtepi32_ps(x)  // convert integer to s.p. float
#define simdi_i2fcast(x)    _mm512_castsi512_ps(x)
#endif //SIMD_INT
#endif //AVX512_SUPPORT


#ifdef AVX2
#include <simde/x86/avx2.h>
// integer support  (usable with AVX2)
#ifndef SIMD_INT
#define SIMD_INT
#define ALIGN_INT           AVX2_ALIGN_INT
#define VECSIZE_INT         AVX2_VECSIZE_INT
uint16_t simd_hmax16_sse(const __m128i buffer);
uint8_t simd_hmax8_sse(const __m128i buffer);
inline uint16_t simd_hmax16_avx(const __m256i buffer) {
    const __m128i abcd = _mm256_castsi256_si128(buffer);
    const uint16_t first = simd_hmax16_sse(abcd);
    const __m128i efgh = _mm256_extracti128_si256(buffer, 1);
    const uint16_t second = simd_hmax16_sse(efgh);
    return std::max(first, second);
}

inline uint8_t simd_hmax8_avx(const __m256i buffer) {
    const __m128i abcd = _mm256_castsi256_si128(buffer);
    const uint8_t first = simd_hmax8_sse(abcd);
    const __m128i efgh = _mm256_extracti128_si256(buffer, 1);
    const uint8_t second = simd_hmax8_sse(efgh);
    return std::max(first, second);
}

template  <unsigned int N>
inline __m256i _mm256_shift_left(__m256i a) {
    __m256i mask = _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0) );
    return _mm256_alignr_epi8(a,mask,16-N);
}

inline unsigned short extract_epi16(__m256i v, int pos) {
    switch(pos){
        case 0: return _mm256_extract_epi16(v, 0);
        case 1: return _mm256_extract_epi16(v, 1);
        case 2: return _mm256_extract_epi16(v, 2);
        case 3: return _mm256_extract_epi16(v, 3);
        case 4: return _mm256_extract_epi16(v, 4);
        case 5: return _mm256_extract_epi16(v, 5);
        case 6: return _mm256_extract_epi16(v, 6);
        case 7: return _mm256_extract_epi16(v, 7);
        case 8: return _mm256_extract_epi16(v, 8);
        case 9: return _mm256_extract_epi16(v, 9);
        case 10: return _mm256_extract_epi16(v, 10);
        case 11: return _mm256_extract_epi16(v, 11);
        case 12: return _mm256_extract_epi16(v, 12);
        case 13: return _mm256_extract_epi16(v, 13);
        case 14: return _mm256_extract_epi16(v, 14);
        case 15: return _mm256_extract_epi16(v, 15);
    }
    return 0;
}

typedef __m256i simd_int;
#define simdi32_add(x,y)    _mm256_add_epi32(x,y)
#define simdi16_add(x,y)    _mm256_add_epi16(x,y)
#define simdi16_adds(x,y)   _mm256_adds_epi16(x,y)
#define simdui8_adds(x,y)   _mm256_adds_epu8(x,y)
#define simdi32_sub(x,y)    _mm256_sub_epi32(x,y)
#define simdui16_subs(x,y)  _mm256_subs_epu16(x,y)
#define simdui8_subs(x,y)   _mm256_subs_epu8(x,y)
#define simdi32_mul(x,y)    _mm256_mullo_epi32(x,y)
#define simdi32_max(x,y)    _mm256_max_epi32(x,y)
#define simdi16_max(x,y)    _mm256_max_epi16(x,y)
#define simdi16_hmax(x)     simd_hmax16_avx(x)
#define simdui8_max(x,y)    _mm256_max_epu8(x,y)
#define simdi8_hmax(x)      simd_hmax8_avx(x)
#define simdi_load(x)       _mm256_load_si256(x)
#define simdi_loadu(x)       _mm256_loadu_si256(x)
#define simdi_streamload(x) _mm256_stream_load_si256(x)
#define simdi_store(x,y)    _mm256_store_si256(x,y)
#define simdi_storeu(x,y)   _mm256_storeu_si256(x,y)
#define simdi32_set(x)      _mm256_set1_epi32(x)
#define simdi16_set(x)      _mm256_set1_epi16(x)
#define simdi8_set(x)       _mm256_set1_epi8(x)
#define simdi32_shuffle(x,y) _mm256_shuffle_epi32(x,y)
#define simdi16_shuffle(x,y) _mm256_shuffle_epi16(x,y)
#define simdi8_shuffle(x,y)  _mm256_shuffle_epi8(x,y)
#define simdi_setzero()     _mm256_setzero_si256()
#define simdi32_gt(x,y)     _mm256_cmpgt_epi32(x,y)
#define simdi8_gt(x,y)      _mm256_cmpgt_epi8(x,y)
#define simdi16_gt(x,y)     _mm256_cmpgt_epi16(x,y)
#define simdi8_eq(x,y)      _mm256_cmpeq_epi8(x,y)
#define simdi16_eq(x,y)     _mm256_cmpeq_epi16(x,y)
#define simdi32_eq(x,y)     _mm256_cmpeq_epi32(x,y)
#define simdi32_lt(x,y)     _mm256_cmpgt_epi32(y,x) // inverse
#define simdi16_lt(x,y)     _mm256_cmpgt_epi16(y,x) // inverse
#define simdi8_lt(x,y)      _mm256_cmpgt_epi8(y,x)
#define simdi_or(x,y)       _mm256_or_si256(x,y)
#define simdi_and(x,y)      _mm256_and_si256(x,y)
#define simdi_andnot(x,y)   _mm256_andnot_si256(x,y)
#define simdi_xor(x,y)      _mm256_xor_si256(x,y)
#define simdi8_shiftl(x,y)  _mm256_shift_left<y>(x)
//TODO fix like shift_left
#define simdi8_shiftr(x,y)  _mm256_srli_si256(x,y)
#define SIMD_MOVEMASK_MAX   0xffffffff
#define simdi8_movemask(x)  _mm256_movemask_epi8(x)
#define simdi16_extract(x,y) extract_epi16(x,y)
#define simdi16_slli(x,y)	_mm256_slli_epi16(x,y) // shift integers in a left by y
#define simdi16_srli(x,y)	_mm256_srli_epi16(x,y) // shift integers in a right by y
#define simdi32_slli(x,y)   _mm256_slli_epi32(x,y) // shift integers in a left by y
#define simdi32_srli(x,y)   _mm256_srli_epi32(x,y) // shift integers in a right by y
#define simdi32_i2f(x) 	    _mm256_cvtepi32_ps(x)  // convert integer to s.p. float
#define simdi_i2fcast(x)    _mm256_castsi256_ps(x)
#endif

#include <simde/x86/avx.h>
// double support (usable with AVX1)
#ifndef SIMD_DOUBLE
#define SIMD_DOUBLE
#define ALIGN_DOUBLE        AVX_ALIGN_DOUBLE
#define VECSIZE_DOUBLE      AVX_VECSIZE_DOUBLE
typedef __m256d simd_double;
#define simdf64_add(x,y)    _mm256_add_pd(x,y)
#define simdf64_sub(x,y)    _mm256_sub_pd(x,y)
#define simdf64_mul(x,y)    _mm256_mul_pd(x,y)
#define simdf64_div(x,y)    _mm256_div_pd(x,y)
#define simdf64_max(x,y)    _mm256_max_pd(x,y)
#define simdf64_load(x)     _mm256_load_pd(x)
#define simdf64_store(x,y)  _mm256_store_pd(x,y)
#define simdf64_set(x)      _mm256_set1_pd(x)
#define simdf64_setzero(x)  _mm256_setzero_pd()
#define simdf64_gt(x,y)     _mm256_cmp_pd(x,y,_CMP_GT_OS)
#define simdf64_lt(x,y)     _mm256_cmp_pd(x,y,_CMP_LT_OS)
#define simdf64_or(x,y)     _mm256_or_pd(x,y)
#define simdf64_and(x,y)    _mm256_and_pd(x,y)
#define simdf64_andnot(x,y) _mm256_andnot_pd(x,y)
#define simdf64_xor(x,y)    _mm256_xor_pd(x,y)
#endif //SIMD_DOUBLE
// float support (usable with AVX1)
#ifndef SIMD_FLOAT
#define SIMD_FLOAT
#define ALIGN_FLOAT         AVX_ALIGN_FLOAT
#define VECSIZE_FLOAT       AVX_VECSIZE_FLOAT
typedef __m256 simd_float;
#define simdf32_add(x,y)    _mm256_add_ps(x,y)
#define simdf32_sub(x,y)    _mm256_sub_ps(x,y)
#define simdf32_mul(x,y)    _mm256_mul_ps(x,y)
#define simdf32_div(x,y)    _mm256_div_ps(x,y)
#define simdf32_rcp(x)      _mm256_rcp_ps(x)
#define simdf32_max(x,y)    _mm256_max_ps(x,y)
#define simdf32_min(x,y)    _mm256_min_ps(x,y)
#define simdf32_load(x)     _mm256_load_ps(x)
#define simdf32_store(x,y)  _mm256_store_ps(x,y)
#define simdf32_set(x)      _mm256_set1_ps(x)
#define simdf32_setzero(x)  _mm256_setzero_ps()
#define simdf32_gt(x,y)     _mm256_cmp_ps(x,y,_CMP_GT_OS)
#define simdf32_eq(x,y)     _mm256_cmp_ps(x,y,_CMP_EQ_OS)
#define simdf32_lt(x,y)     _mm256_cmp_ps(x,y,_CMP_LT_OS)
#define simdf32_or(x,y)     _mm256_or_ps(x,y)
#define simdf32_and(x,y)    _mm256_and_ps(x,y)
#define simdf32_andnot(x,y) _mm256_andnot_ps(x,y)
#define simdf32_xor(x,y)    _mm256_xor_ps(x,y)
#define simdf32_f2i(x) 	    _mm256_cvtps_epi32(x)  // convert s.p. float to integer
#define simdf_f2icast(x)    _mm256_castps_si256(x) // compile time cast
#endif
#endif

#include <simde/x86/sse4.1.h>
inline uint16_t simd_hmax16_sse(const __m128i buffer) {
    __m128i tmp1 = _mm_subs_epu16(_mm_set1_epi16((short)65535), buffer);
    __m128i tmp3 = _mm_minpos_epu16(tmp1);
    return (65535 - _mm_cvtsi128_si32(tmp3));
}

inline uint8_t simd_hmax8_sse(const __m128i buffer) {
    __m128i tmp1 = _mm_subs_epu8(_mm_set1_epi8((char)255), buffer);
    __m128i tmp2 = _mm_min_epu8(tmp1, _mm_srli_epi16(tmp1, 8));
    __m128i tmp3 = _mm_minpos_epu16(tmp2);
    return (int8_t)(255 -(int8_t) _mm_cvtsi128_si32(tmp3));
}

// double support
#ifndef SIMD_DOUBLE
#define SIMD_DOUBLE
#define ALIGN_DOUBLE        SSE_ALIGN_DOUBLE
#define VECSIZE_DOUBLE      SSE_VECSIZE_DOUBLE
typedef __m128d simd_double;
#define simdf64_add(x,y)    _mm_add_pd(x,y)
#define simdf64_sub(x,y)    _mm_sub_pd(x,y)
#define simdf64_mul(x,y)    _mm_mul_pd(x,y)
#define simdf64_div(x,y)    _mm_div_pd(x,y)
#define simdf64_max(x,y)    _mm_max_pd(x,y)
#define simdf64_load(x)     _mm_load_pd(x)
#define simdf64_store(x,y)  _mm_store_pd(x,y)
#define simdf64_set(x)      _mm_set1_pd(x)
#define simdf64_setzero(x)  _mm_setzero_pd()
#define simdf64_gt(x,y)     _mm_cmpgt_pd(x,y)
#define simdf64_lt(x,y)     _mm_cmplt_pd(x,y)
#define simdf64_or(x,y)     _mm_or_pd(x,y)
#define simdf64_and(x,y)    _mm_and_pd(x,y)
#define simdf64_andnot(x,y) _mm_andnot_pd(x,y)
#define simdf64_xor(x,y)    _mm_xor_pd(x,y)
#endif //SIMD_DOUBLE

// float support
#ifndef SIMD_FLOAT
#define SIMD_FLOAT
#define ALIGN_FLOAT         SSE_ALIGN_FLOAT
#define VECSIZE_FLOAT       SSE_VECSIZE_FLOAT
typedef __m128  simd_float;
#define simdf32_add(x,y)    _mm_add_ps(x,y)
#define simdf32_sub(x,y)    _mm_sub_ps(x,y)
#define simdf32_mul(x,y)    _mm_mul_ps(x,y)
#define simdf32_div(x,y)    _mm_div_ps(x,y)
#define simdf32_rcp(x)      _mm_rcp_ps(x)
#define simdf32_max(x,y)    _mm_max_ps(x,y)
#define simdf32_min(x,y)    _mm_min_ps(x,y)
#define simdf32_load(x)     _mm_load_ps(x)
#define simdf32_store(x,y)  _mm_store_ps(x,y)
#define simdf32_set(x)      _mm_set1_ps(x)
#define simdf32_setzero(x)  _mm_setzero_ps()
#define simdf32_gt(x,y)     _mm_cmpgt_ps(x,y)
#define simdf32_eq(x,y)     _mm_cmpeq_ps(x,y)
#define simdf32_lt(x,y)     _mm_cmplt_ps(x,y)
#define simdf32_or(x,y)     _mm_or_ps(x,y)
#define simdf32_and(x,y)    _mm_and_ps(x,y)
#define simdf32_andnot(x,y) _mm_andnot_ps(x,y)
#define simdf32_xor(x,y)    _mm_xor_ps(x,y)
#define simdf32_f2i(x) 	    _mm_cvtps_epi32(x)  // convert s.p. float to integer
#define simdf_f2icast(x)    _mm_castps_si128(x) // compile time cast
#endif //SIMD_FLOAT

// integer support
#ifndef SIMD_INT
#define SIMD_INT
inline unsigned short extract_epi16(__m128i v, int pos) {
    switch(pos){
        case 0: return _mm_extract_epi16(v, 0);
        case 1: return _mm_extract_epi16(v, 1);
        case 2: return _mm_extract_epi16(v, 2);
        case 3: return _mm_extract_epi16(v, 3);
        case 4: return _mm_extract_epi16(v, 4);
        case 5: return _mm_extract_epi16(v, 5);
        case 6: return _mm_extract_epi16(v, 6);
        case 7: return _mm_extract_epi16(v, 7);
    }
    return 0;
}
#define ALIGN_INT           SSE_ALIGN_INT
#define VECSIZE_INT         SSE_VECSIZE_INT
typedef __m128i simd_int;
#define simdi32_add(x,y)    _mm_add_epi32(x,y)
#define simdi16_add(x,y)    _mm_add_epi16(x,y)
#define simdi16_adds(x,y)   _mm_adds_epi16(x,y)
#define simdui8_adds(x,y)   _mm_adds_epu8(x,y)
#define simdi32_sub(x,y)    _mm_sub_epi32(x,y)
#define simdui16_subs(x,y)  _mm_subs_epu16(x,y)
#define simdui8_subs(x,y)   _mm_subs_epu8(x,y)
#define simdi32_mul(x,y)    _mm_mullo_epi32(x,y) // SSE4.1
#define simdi32_max(x,y)    _mm_max_epi32(x,y) // SSE4.1
#define simdi16_max(x,y)    _mm_max_epi16(x,y)
#define simdi16_hmax(x)     simd_hmax16_sse(x)
#define simdui8_max(x,y)    _mm_max_epu8(x,y)
#define simdi8_hmax(x)      simd_hmax8_sse(x)
#define simdi_load(x)       _mm_load_si128(x)
#define simdi_loadu(x)      _mm_loadu_si128(x)
#define simdi_streamload(x) _mm_stream_load_si128(x)
#define simdi_storeu(x,y)   _mm_storeu_si128(x,y)
#define simdi_store(x,y)    _mm_store_si128(x,y)
#define simdi32_set(x)      _mm_set1_epi32(x)
#define simdi16_set(x)      _mm_set1_epi16(x)
#define simdi8_set(x)       _mm_set1_epi8(x)
#define simdi32_shuffle(x,y) _mm_shuffle_epi32(x,y)
#define simdi16_shuffle(x,y) _mm_shuffle_epi16(x,y)
#define simdi8_shuffle(x,y)  _mm_shuffle_epi8(x,y)
#define simdi_setzero()     _mm_setzero_si128()
#define simdi32_gt(x,y)     _mm_cmpgt_epi32(x,y)
#define simdi8_gt(x,y)      _mm_cmpgt_epi8(x,y)
#define simdi32_eq(x,y)     _mm_cmpeq_epi32(x,y)
#define simdi16_eq(x,y)     _mm_cmpeq_epi16(x,y)
#define simdi8_eq(x,y)      _mm_cmpeq_epi8(x,y)
#define simdi32_lt(x,y)     _mm_cmplt_epi32(x,y)
#define simdi16_lt(x,y)     _mm_cmplt_epi16(x,y)
#define simdi8_lt(x,y)      _mm_cmplt_epi8(x,y)
#define simdi16_gt(x,y)     _mm_cmpgt_epi16(x,y)
#define simdi_or(x,y)       _mm_or_si128(x,y)
#define simdi_and(x,y)      _mm_and_si128(x,y)
#define simdi_andnot(x,y)   _mm_andnot_si128(x,y)
#define simdi_xor(x,y)      _mm_xor_si128(x,y)
#define simdi8_shiftl(x,y)  _mm_slli_si128(x,y)
#define simdi8_shiftr(x,y)  _mm_srli_si128(x,y)
#define SIMD_MOVEMASK_MAX   0xffff
#define simdi8_movemask(x)  _mm_movemask_epi8(x)
#define simdi16_extract(x,y) extract_epi16(x,y)
#define simdi16_slli(x,y)	_mm_slli_epi16(x,y) // shift integers in a left by y
#define simdi16_srli(x,y)	_mm_srli_epi16(x,y) // shift integers in a right by y
#define simdi32_slli(x,y)	_mm_slli_epi32(x,y) // shift integers in a left by y
#define simdi32_srli(x,y)	_mm_srli_epi32(x,y) // shift integers in a right by y
#define simdi32_i2f(x) 	    _mm_cvtepi32_ps(x)  // convert integer to s.p. float
#define simdi_i2fcast(x)    _mm_castsi128_ps(x)
#endif //SIMD_INT

inline void *mem_align(size_t boundary, size_t size) {
    void *pointer;
    if (posix_memalign(&pointer, boundary, size) != 0) {
#define MEM_ALIGN_ERROR "mem_align could not allocate memory.\n"
        fwrite(MEM_ALIGN_ERROR, sizeof(MEM_ALIGN_ERROR), 1, stderr);
#undef MEM_ALIGN_ERROR
        exit(3);
    }
    return pointer;
}
#ifdef SIMD_FLOAT
inline simd_float * malloc_simd_float(const size_t size) {
    return (simd_float *) mem_align(ALIGN_FLOAT, size);
}
#endif
#ifdef SIMD_DOUBLE
inline simd_double * malloc_simd_double(const size_t size) {
    return (simd_double *) mem_align(ALIGN_DOUBLE, size);
}
#endif
#ifdef SIMD_INT
inline simd_int * malloc_simd_int(const size_t size) {
    return (simd_int *) mem_align(ALIGN_INT, size);
}
#endif

template <typename T>
T** malloc_matrix(int dim1, int dim2) {
#define ICEIL(x_int, fac_int) ((x_int + fac_int - 1) / fac_int) * fac_int
    // Compute mem sizes rounded up to nearest multiple of ALIGN_FLOAT
    size_t size_pointer_array = ICEIL(dim1*sizeof(T*), ALIGN_FLOAT);
    size_t dim2_padded = ICEIL(dim2*sizeof(T), ALIGN_FLOAT)/sizeof(T);

    T** matrix = (T**) mem_align( ALIGN_FLOAT, size_pointer_array + dim1*dim2_padded*sizeof(T) );
    if (matrix == NULL)
        return matrix;

    T* ptr = (T*) (matrix + (size_pointer_array/sizeof(T*)) );
    for (int i=0; i<dim1; ++i) {
        matrix[i] = ptr;
        ptr += dim2_padded;
    }
#undef ICEIL
    return matrix;
}


inline float ScalarProd20(const float* qi, const float* tj) {
//#ifdef AVX
//  float __attribute__((aligned(ALIGN_FLOAT))) res;
//  __m256 P; // query 128bit SSE2 register holding 4 floats
//  __m256 S; // aux register
//  __m256 R; // result
//  __m256* Qi = (__m256*) qi;
//  __m256* Tj = (__m256*) tj;
//
//  R = _mm256_mul_ps(*(Qi++),*(Tj++));
//  P = _mm256_mul_ps(*(Qi++),*(Tj++));
//  S = _mm256_mul_ps(*Qi,*Tj); // floats A, B, C, D, ?, ?, ? ,?
//  R = _mm256_add_ps(R,P);     // floats 0, 1, 2, 3, 4, 5, 6, 7
//  P = _mm256_permute2f128_ps(R, R, 0x01); // swap hi and lo 128 bits: 4, 5, 6, 7, 0, 1, 2, 3
//  R = _mm256_add_ps(R,P);     // 0+4, 1+5, 2+6, 3+7, 0+4, 1+5, 2+6, 3+7
//  R = _mm256_add_ps(R,S);     // 0+4+A, 1+5+B, 2+6+C, 3+7+D, ?, ?, ? ,?
//  R = _mm256_hadd_ps(R,R);    // 04A15B, 26C37D, ?, ?, 04A15B, 26C37D, ?, ?
//  R = _mm256_hadd_ps(R,R);    // 01234567ABCD, ?, 01234567ABCD, ?, 01234567ABCD, ?, 01234567ABCD, ?
//  _mm256_store_ps(&res, R);
//  return res;
//#else
//
//
//TODO fix this
    float __attribute__((aligned(16))) res;
    __m128 P; // query 128bit SSE2 register holding 4 floats
    __m128 R;// result
    __m128* Qi = (__m128*) qi;
    __m128* Tj = (__m128*) tj;

    __m128 P1 = _mm_mul_ps(*(Qi),*(Tj));
    __m128 P2 = _mm_mul_ps(*(Qi+1),*(Tj+1));
    __m128 R1 = _mm_add_ps(P1, P2);

    __m128 P3 = _mm_mul_ps(*(Qi + 2), *(Tj + 2));
    __m128 P4 = _mm_mul_ps(*(Qi + 3), *(Tj + 3));
    __m128 R2 = _mm_add_ps(P3, P4);
    __m128 P5 = _mm_mul_ps(*(Qi+4), *(Tj+4));

    R = _mm_add_ps(R1, R2);
    R = _mm_add_ps(R,P5);

//    R = _mm_hadd_ps(R,R);
//    R = _mm_hadd_ps(R,R);
    P = _mm_shuffle_ps(R, R, _MM_SHUFFLE(2,0,2,0));
    R = _mm_shuffle_ps(R, R, _MM_SHUFFLE(3,1,3,1));
    R = _mm_add_ps(R,P);
    P = _mm_shuffle_ps(R, R, _MM_SHUFFLE(2,0,2,0));
    R = _mm_shuffle_ps(R, R, _MM_SHUFFLE(3,1,3,1));
    R = _mm_add_ps(R,P);
    _mm_store_ss(&res, R);
    return res;
//#endif
//    return tj[0] * qi[0] + tj[1] * qi[1] + tj[2] * qi[2] + tj[3] * qi[3]
//            + tj[4] * qi[4] + tj[5] * qi[5] + tj[6] * qi[6] + tj[7] * qi[7]
//            + tj[8] * qi[8] + tj[9] * qi[9] + tj[10] * qi[10] + tj[11] * qi[11]
//            + tj[12] * qi[12] + tj[13] * qi[13] + tj[14] * qi[14]
//            + tj[15] * qi[15] + tj[16] * qi[16] + tj[17] * qi[17]
//            + tj[18] * qi[18] + tj[19] * qi[19];
}

#endif //SIMD_H
