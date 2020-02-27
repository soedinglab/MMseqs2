// sse2wasm is still very incomplete
// licensed under GPLv3 see LICENCE file
#ifndef SSE2WASM
#define SSE2WASM

#define __wasm_unimplemented_simd128__
#include <wasm_simd128.h>
#define SSE 1
typedef v128_t __m128d;
typedef v128_t __m128i;
typedef v128_t __m128;

#define _mm_add_ps        wasm_f32x4_add
#define _mm_sub_ps        wasm_f32x4_sub
#define _mm_mul_ps        wasm_f32x4_mul
#define _mm_div_ps        wasm_f32x4_div
#define _mm_rcp_ps(x)     wasm_f32x4_div(wasm_f32x4_splat(1), (x))
#define _mm_max_ps        wasm_f32x4_max
#define _mm_min_ps        wasm_f32x4_min
#define _mm_load_ps       wasm_v128_load
#define _mm_load_ss       wasm_v128_load
#define _mm_store_ps      wasm_v128_store
#define _mm_store_ss      wasm_v128_store
#define _mm_set1_ps       wasm_f32x4_splat
#define _mm_setzero_ps(x) wasm_f32x4_splat(0)
#define _mm_cmpgt_ps      wasm_f32x4_gt
#define _mm_cmpeq_ps      wasm_f32x4_eq
#define _mm_cmplt_ps      wasm_f32x4_lt
#define _mm_or_ps         wasm_v128_or
#define _mm_and_ps        wasm_v128_and
#define _mm_andnot_ps     wasm_v128_andnot
#define _mm_xor_ps        wasm_v128_xor
#define _mm_cvtps_epi32   NOT_YET_IMP
#define _mm_castps_si128  NOT_YET_IMP

#define _mm_add_epi32     wasm_i32x4_add
#define _mm_add_epi16     wasm_i16x8_add
#define _mm_add_epi8      wasm_i8x16_add
#define _mm_adds_epi16    wasm_i16x8_add_saturate
#define _mm_adds_epu8     wasm_u8x16_add_saturate
#define _mm_sub_epi32     wasm_i32x4_sub
#define _mm_sub_epi16     wasm_i16x8_sub
#define _mm_sub_epi8      wasm_i8x16_sub
#define _mm_subs_epu16    wasm_u16x8_sub_saturate
#define _mm_subs_epu8     wasm_u8x16_sub_saturate
#define _mm_mullo_epi32   NOT_YET_IMP
#define _mm_max_epi32     wasm_i32x4_max_s
#define _mm_max_epi16     wasm_i16x8_max_s
#define _mm_max_epu8      wasm_i16x8_max_u
#define _mm_min_epu8      wasm_i16x8_min_u
#define _mm_load_si128    wasm_v128_load
#define _mm_loadu_si128   wasm_v128_load
#define _mm_stream_load_si128 NOT_YET_IMP
#define _mm_storeu_si128  wasm_v128_store
#define _mm_store_si128   wasm_v128_store
#define _mm_set1_epi32    wasm_i32x4_splat
#define _mm_set1_epi16    wasm_i16x8_splat
#define _mm_set1_epi8     wasm_i8x16_splat
#define _mm_set_epi32     wasm_i32x4_make
#define _mm_set_epi16     wasm_i16x8_make
#define _mm_set_epi8      wasm_i8x16_make
#define _mm_setzero_si128(x) wasm_i64x2_splat(0)
#define _mm_cmpgt_epi32   wasm_i32x4_gt
#define _mm_cmpgt_epi16   wasm_i16x8_gt
#define _mm_cmpgt_epi8    wasm_i8x16_gt
#define _mm_cmpeq_epi32   wasm_i32x4_eq
#define _mm_cmpeq_epi16   wasm_i16x8_eq
#define _mm_cmpeq_epi8    wasm_i8x16_eq
#define _mm_cmplt_epi32   wasm_i32x4_lt
#define _mm_cmplt_epi16   wasm_i16x8_lt
#define _mm_cmplt_epi8    wasm_i8x16_lt
#define _mm_or_si128      wasm_v128_or
#define _mm_and_si128     wasm_v128_and
#define _mm_andnot_si128  wasm_v128_andnot
#define _mm_xor_si128     wasm_v128_xor
#define _mm_extract_epi16 wasm_i16x8_extract_lane
#define _mm_extract_epi8  wasm_i8x16_extract_lane
#define _mm_slli_epi16    NOT_YET_IMP
#define _mm_srli_epi16    NOT_YET_IMP
#define _mm_slli_epi32    NOT_YET_IMP
#define _mm_srli_epi32    NOT_YET_IMP
#define _mm_cvtepi32_ps(x) (__m128)((__m128i)(x))
#define _mm_castsi128_ps  NOT_YET_IMP

static inline void _mm_storel_epi64(__m128i* mem_addr, __m128i a) {
    return;
}

static inline __m128i _mm_setr_epi32(int e3, int e2, int e1, int e0) {
    union {
        int32_t as_i32[4];
        __m128i  as_vec;
    } t;
    t.as_i32[0] = e3;
    t.as_i32[1] = e2;
    t.as_i32[2] = e1;
    t.as_i32[3] = e0;
    return t.as_vec;
}

//__builtin_convertvector

static inline __m128i _mm_slli_si128(__m128i a, int imm8) {
    return _mm_setzero_si128();
}

static inline __m128i _mm_srli_si128(__m128i a, int imm8) {
    return _mm_setzero_si128();
}

static inline unsigned short _mm_movemask_epi8(__m128i a) {
    unsigned int result=0;

    union {
        __m128i si;
        char as_char[16];
    } t;
    t.si = a;
    result |= (t.as_char[15] & 0x80) << (15-7);
    result |= (t.as_char[14] & 0x80) << (14-7);
    result |= (t.as_char[13] & 0x80) << (13-7);
    result |= (t.as_char[12] & 0x80) << (12-7);
    result |= (t.as_char[11] & 0x80) << (11-7);
    result |= (t.as_char[10] & 0x80) << (10-7);
    result |= (t.as_char[9]  & 0x80) <<  (9-7);
    result |= (t.as_char[8]  & 0x80) <<  (8-7);
    result |= (t.as_char[7]  & 0x80);
    result |= (t.as_char[6]  & 0x80) >>  (7-6);
    result |= (t.as_char[5]  & 0x80) >>  (7-5);
    result |= (t.as_char[4]  & 0x80) >>  (7-4);
    result |= (t.as_char[3]  & 0x80) >>  (7-3);
    result |= (t.as_char[2]  & 0x80) >>  (7-2);
    result |= (t.as_char[1]  & 0x80) >>  (7-1);
    result |= (t.as_char[0]  & 0x80) >>   7;

    return result;
}

#define _MM_SHUFFLE(a, b, c, d) _mm_setzero_si128()


static inline __m128i _mm_shuffle_epi32(__m128i a, __m128i b) {
    return _mm_setzero_si128();
}

static inline __m128i _mm_shuffle_epi16(__m128i a, __m128i b) {
    return _mm_setzero_si128();
}
// wasm_v8x16_shuffle
static inline __m128i _mm_shuffle_epi8(__m128i a, __m128i b) {
    return _mm_setzero_si128();
}

static inline int64_t _mm_cvtsi128_si64(__m128 a) {
    return wasm_i64x2_extract_lane(a, 0);
}

static inline __m128i _mm_cvtsi64_si128(int64_t a) {
    return wasm_i64x2_make(a, 0);
}

static inline __m128i _mm_cvtsi32_si128(int a) {
    return wasm_i64x2_make((int64_t)a, 0);
}

#endif
