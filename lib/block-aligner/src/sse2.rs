#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

pub type Simd = __m128i; // use for storing DP scores
pub type HalfSimd = __m128i; // used for storing bytes (sequence or scoring matrix)
pub type LutSimd = __m128i; // used for storing a row in a scoring matrix (always 128 bits)
pub type TraceType = i16;
/// Number of 16-bit lanes in a SIMD vector.
pub const L: usize = 8;
pub const L_BYTES: usize = L * 2;
pub const HALFSIMD_MUL: usize = 2;
// using min = 0 is faster, but restricts range of scores (and restricts the max block size)
pub const ZERO: i16 = 1 << 14;
pub const MIN: i16 = 0;

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn store_trace(ptr: *mut TraceType, trace: TraceType) { *ptr = trace; }

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn simd_adds_i16(a: Simd, b: Simd) -> Simd { _mm_adds_epi16(a, b) }

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn simd_subs_i16(a: Simd, b: Simd) -> Simd { _mm_subs_epi16(a, b) }

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn simd_max_i16(a: Simd, b: Simd) -> Simd { _mm_max_epi16(a, b) }

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn simd_cmpeq_i16(a: Simd, b: Simd) -> Simd { _mm_cmpeq_epi16(a, b) }

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn simd_cmpgt_i16(a: Simd, b: Simd) -> Simd { _mm_cmpgt_epi16(a, b) }

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn simd_blend_i8(a: Simd, b: Simd, mask: Simd) -> Simd { _mm_or_si128(_mm_andnot_si128(mask, a), _mm_and_si128(mask, b)) }

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn simd_load(ptr: *const Simd) -> Simd { _mm_load_si128(ptr) }

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn simd_loadu(ptr: *const Simd) -> Simd { _mm_loadu_si128(ptr) }

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn simd_store(ptr: *mut Simd, a: Simd) { _mm_store_si128(ptr, a) }

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn simd_set1_i16(v: i16) -> Simd { _mm_set1_epi16(v) }

#[macro_export]
#[doc(hidden)]
macro_rules! simd_extract_i16 {
    ($a:expr, $num:expr) => {
        {
            debug_assert!($num < L);
            #[cfg(target_arch = "x86")]
            use std::arch::x86::*;
            #[cfg(target_arch = "x86_64")]
            use std::arch::x86_64::*;
            _mm_extract_epi16($a, $num as i32) as i16
        }
    };
}

#[macro_export]
#[doc(hidden)]
macro_rules! simd_insert_i16 {
    ($a:expr, $v:expr, $num:expr) => {
        {
            debug_assert!($num < L);
            #[cfg(target_arch = "x86")]
            use std::arch::x86::*;
            #[cfg(target_arch = "x86_64")]
            use std::arch::x86_64::*;
            _mm_insert_epi16($a, $v as i32, $num as i32)
        }
    };
}

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn simd_movemask_i8(a: Simd) -> u16 { _mm_movemask_epi8(a) as u16 }

#[macro_export]
#[doc(hidden)]
macro_rules! simd_sl_i16 {
    ($a:expr, $b:expr, $num:expr) => {
        {
            debug_assert!($num <= L);
            #[cfg(target_arch = "x86")]
            use std::arch::x86::*;
            #[cfg(target_arch = "x86_64")]
            use std::arch::x86_64::*;
            _mm_or_si128(_mm_slli_si128($a, (2 * $num) as i32), _mm_srli_si128($b, ((L - $num) * 2) as i32))
        }
    };
}

#[macro_export]
#[doc(hidden)]
macro_rules! simd_sr_i16 {
    ($a:expr, $b:expr, $num:expr) => {
        {
            debug_assert!($num <= L);
            #[cfg(target_arch = "x86")]
            use std::arch::x86::*;
            #[cfg(target_arch = "x86_64")]
            use std::arch::x86_64::*;
            _mm_or_si128(_mm_slli_si128($a, ((L - $num) * 2) as i32), _mm_srli_si128($b, (2 * $num) as i32))
        }
    };
}

// hardcoded to STEP = 8
#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn simd_step(a: Simd, b: Simd) -> Simd {
    a
}

// shift in zeros
macro_rules! simd_sllz_i16 {
    ($a:expr, $num:expr) => {
        {
            debug_assert!($num < L);
            #[cfg(target_arch = "x86")]
            use std::arch::x86::*;
            #[cfg(target_arch = "x86_64")]
            use std::arch::x86_64::*;
            _mm_slli_si128($a, ($num * 2) as i32)
        }
    };
}

// broadcast last 16-bit element to the whole vector
#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn simd_broadcasthi_i16(v: Simd) -> Simd {
    let v = _mm_shufflehi_epi16(v, 0b11111111);
    _mm_shuffle_epi32(v, 0b11111111)
}

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn simd_slow_extract_i16(v: Simd, i: usize) -> i16 {
    debug_assert!(i < L);

    #[repr(align(16))]
    struct A([i16; L]);

    let mut a = A([0i16; L]);
    simd_store(a.0.as_mut_ptr() as *mut Simd, v);
    *a.0.as_ptr().add(i)
}

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn simd_hmax_i16(v: Simd) -> i16 {
    let mut v2 = _mm_max_epi16(v, _mm_srli_si128(v, 2));
    v2 = _mm_max_epi16(v2, _mm_srli_si128(v2, 4));
    v2 = _mm_max_epi16(v2, _mm_srli_si128(v2, 8));
    simd_extract_i16!(v2, 0)
}

#[macro_export]
#[doc(hidden)]
macro_rules! simd_prefix_hadd_i16 {
    ($a:expr, $num:expr) => {
        {
            debug_assert!($num <= L);
            #[cfg(target_arch = "x86")]
            use std::arch::x86::*;
            #[cfg(target_arch = "x86_64")]
            use std::arch::x86_64::*;
            let mut v = _mm_subs_epi16($a, _mm_set1_epi16(ZERO));
            if $num > 4 {
                v = _mm_adds_epi16(v, _mm_srli_si128(v, 8));
            }
            if $num > 2 {
                v = _mm_adds_epi16(v, _mm_srli_si128(v, 4));
            }
            if $num > 1 {
                v = _mm_adds_epi16(v, _mm_srli_si128(v, 2));
            }
            simd_extract_i16!(v, 0)
        }
    };
}

#[macro_export]
#[doc(hidden)]
macro_rules! simd_prefix_hmax_i16 {
    ($a:expr, $num:expr) => {
        {
            debug_assert!($num <= L);
            #[cfg(target_arch = "x86")]
            use std::arch::x86::*;
            #[cfg(target_arch = "x86_64")]
            use std::arch::x86_64::*;
            let mut v = $a;
            if $num > 4 {
                v = _mm_max_epi16(v, _mm_srli_si128(v, 8));
            }
            if $num > 2 {
                v = _mm_max_epi16(v, _mm_srli_si128(v, 4));
            }
            if $num > 1 {
                v = _mm_max_epi16(v, _mm_srli_si128(v, 2));
            }
            simd_extract_i16!(v, 0)
        }
    };
}

#[macro_export]
#[doc(hidden)]
macro_rules! simd_suffix_hmax_i16 {
    ($a:expr, $num:expr) => {
        {
            debug_assert!($num <= L);
            #[cfg(target_arch = "x86")]
            use std::arch::x86::*;
            #[cfg(target_arch = "x86_64")]
            use std::arch::x86_64::*;
            let mut v = $a;
            if $num > 4 {
                v = _mm_max_epi16(v, _mm_slli_si128(v, 8));
            }
            if $num > 2 {
                v = _mm_max_epi16(v, _mm_slli_si128(v, 4));
            }
            if $num > 1 {
                v = _mm_max_epi16(v, _mm_slli_si128(v, 2));
            }
            simd_extract_i16!(v, 7)
        }
    };
}

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn simd_hargmax_i16(v: Simd, max: i16) -> usize {
    let v2 = _mm_cmpeq_epi16(v, _mm_set1_epi16(max));
    (simd_movemask_i8(v2).trailing_zeros() as usize) / 2
}

#[target_feature(enable = "sse2")]
#[inline]
#[allow(non_snake_case)]
#[allow(dead_code)]
pub unsafe fn simd_naive_prefix_scan_i16(R_max: Simd, gap_cost: Simd, _gap_cost_lane: PrefixScanConsts) -> Simd {
    let mut curr = R_max;

    for _i in 0..(L - 1) {
        let prev = curr;
        curr = simd_sl_i16!(curr, _mm_setzero_si128(), 1);
        curr = _mm_adds_epi16(curr, gap_cost);
        curr = _mm_max_epi16(curr, prev);
    }

    curr
}

pub type PrefixScanConsts = ();

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn get_prefix_scan_consts(gap: Simd) -> (Simd, PrefixScanConsts) {
    let mut shift1 = simd_sllz_i16!(gap, 1);
    shift1 = _mm_adds_epi16(shift1, gap);
    let mut shift2 = simd_sllz_i16!(shift1, 2);
    shift2 = _mm_adds_epi16(shift2, shift1);
    let mut shift4 = simd_sllz_i16!(shift2, 4);
    shift4 = _mm_adds_epi16(shift4, shift2);

    (shift4, ())
}

#[target_feature(enable = "sse2")]
#[inline]
#[allow(non_snake_case)]
pub unsafe fn simd_prefix_scan_i16(R_max: Simd, gap_cost: Simd, _gap_cost_lane: PrefixScanConsts) -> Simd {
    // Optimized prefix add and max for every eight elements
    // Note: be very careful to avoid lane-crossing which has a large penalty.
    // Also, make sure to use as little registers as possible to avoid
    // memory loads (latencies really matter since this is critical path).
    // Keep the CPU busy with instructions!
    // Note: relies on min score = 0 for speed!
    let mut shift1 = simd_sllz_i16!(R_max, 1);
    shift1 = _mm_adds_epi16(shift1, gap_cost);
    shift1 = _mm_max_epi16(R_max, shift1);
    let mut shift2 = simd_sllz_i16!(shift1, 2);
    shift2 = _mm_adds_epi16(shift2, _mm_slli_epi16(gap_cost, 1));
    shift2 = _mm_max_epi16(shift1, shift2);
    let mut shift4 = simd_sllz_i16!(shift2, 4);
    shift4 = _mm_adds_epi16(shift4, _mm_slli_epi16(gap_cost, 2));
    shift4 = _mm_max_epi16(shift2, shift4);

    shift4
}

// lookup two 128-bit tables
#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn halfsimd_lookup2_i16(lut1: LutSimd, lut2: LutSimd, v: HalfSimd) -> Simd {
    #[repr(align(16))]
    struct A([i8; L * HALFSIMD_MUL]);
    #[repr(align(16))]
    struct T([i8; L * HALFSIMD_MUL * 2]);
    #[repr(align(16))]
    struct A2([i16; L]);

    let mut idx = A([0i8; L * HALFSIMD_MUL]);
    simd_store(idx.0.as_mut_ptr() as *mut HalfSimd, v);
    let idx_ptr = idx.0.as_ptr();

    let mut table = T([0i8; L * HALFSIMD_MUL * 2]);
    simd_store(table.0.as_mut_ptr() as *mut LutSimd, lut1);
    simd_store(table.0.as_mut_ptr().add(L * HALFSIMD_MUL) as *mut LutSimd, lut2);
    let table_ptr = table.0.as_ptr();

    let mut res = A2([0i16; L]);

    let mut i = 0;
    while i < L {
        *res.0.as_mut_ptr().add(i) = *table_ptr.add(*idx_ptr.add(i) as usize) as i16;
        i += 1;
    }

    simd_load(res.0.as_ptr() as *const Simd)
}

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn halfsimd_lookup1_i16(lut: LutSimd, v: HalfSimd) -> Simd {
    #[repr(align(16))]
    struct A([i8; L * HALFSIMD_MUL]);
    #[repr(align(16))]
    struct A2([i16; L]);

    let mut idx = A([0i8; L * HALFSIMD_MUL]);
    simd_store(idx.0.as_mut_ptr() as *mut HalfSimd, v);
    let idx_ptr = idx.0.as_ptr();

    let mut table = A([0i8; L * HALFSIMD_MUL]);
    simd_store(table.0.as_mut_ptr() as *mut LutSimd, lut);
    let table_ptr = table.0.as_ptr();

    let mut res = A2([0i16; L]);

    let mut i = 0;
    while i < L {
        *res.0.as_mut_ptr().add(i) = *table_ptr.add((*idx_ptr.add(i) as usize) & 0b1111) as i16;
        i += 1;
    }

    simd_load(res.0.as_ptr() as *const Simd)
}

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn halfsimd_lookup_bytes_i16(match_scores: HalfSimd, mismatch_scores: HalfSimd, a: HalfSimd, b: HalfSimd) -> Simd {
    let mask = _mm_cmpeq_epi8(a, b);
    let c = simd_blend_i8(mismatch_scores, match_scores, mask);

    #[repr(align(16))]
    struct A([i8; L * HALFSIMD_MUL]);
    #[repr(align(16))]
    struct A2([i16; L]);

    let mut a = A([0i8; L * HALFSIMD_MUL]);
    simd_store(a.0.as_mut_ptr() as *mut HalfSimd, c);
    let a_ptr = a.0.as_ptr();

    let mut res = A2([0i16; L]);

    let mut i = 0;
    while i < L {
        *res.0.as_mut_ptr().add(i) = *a_ptr.add(i) as i16;
        i += 1;
    }

    simd_load(res.0.as_ptr() as *const Simd)
}

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn halfsimd_load(ptr: *const HalfSimd) -> HalfSimd { _mm_load_si128(ptr) }

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn halfsimd_loadu(ptr: *const HalfSimd) -> HalfSimd { _mm_loadu_si128(ptr) }

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn lutsimd_load(ptr: *const LutSimd) -> LutSimd { _mm_load_si128(ptr) }

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn lutsimd_loadu(ptr: *const LutSimd) -> LutSimd { _mm_loadu_si128(ptr) }

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn halfsimd_store(ptr: *mut HalfSimd, a: HalfSimd) { _mm_store_si128(ptr, a) }

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn halfsimd_sub_i8(a: HalfSimd, b: HalfSimd) -> HalfSimd { _mm_sub_epi8(a, b) }

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn halfsimd_set1_i8(v: i8) -> HalfSimd { _mm_set1_epi8(v) }

#[target_feature(enable = "sse2")]
#[inline]
pub unsafe fn halfsimd_get_idx(i: usize) -> usize { i + i / L * L }

#[macro_export]
#[doc(hidden)]
macro_rules! halfsimd_sr_i8 {
    ($a:expr, $b:expr, $num:expr) => {
        {
            debug_assert!($num <= L);
            #[cfg(target_arch = "x86")]
            use std::arch::x86::*;
            #[cfg(target_arch = "x86_64")]
            use std::arch::x86_64::*;
            let mask = _mm_srli_si128(_mm_set1_epi32(-1i32), (L + $num) as i32);
            _mm_or_si128(_mm_slli_si128($a, (L - $num) as i32), _mm_and_si128(_mm_srli_si128($b, $num as i32), mask))
        }
    };
}

#[target_feature(enable = "sse2")]
#[allow(dead_code)]
pub unsafe fn simd_dbg_i16(v: Simd) {
    #[repr(align(16))]
    struct A([i16; L]);

    let mut a = A([0i16; L]);
    simd_store(a.0.as_mut_ptr() as *mut Simd, v);

    for i in (0..a.0.len()).rev() {
        print!("{:6} ", a.0[i]);
    }
    println!();
}

#[target_feature(enable = "sse2")]
#[allow(dead_code)]
pub unsafe fn halfsimd_dbg_i8(v: HalfSimd) {
    #[repr(align(16))]
    struct A([i8; L * HALFSIMD_MUL]);

    let mut a = A([0i8; L * HALFSIMD_MUL]);
    halfsimd_store(a.0.as_mut_ptr() as *mut HalfSimd, v);

    for i in (0..a.0.len()).rev() {
        print!("{:3} ", a.0[i]);
    }
    println!();
}

#[target_feature(enable = "sse2")]
#[allow(dead_code)]
pub unsafe fn simd_assert_vec_eq(a: Simd, b: [i16; L]) {
    #[repr(align(16))]
    struct A([i16; L]);

    let mut arr = A([0i16; L]);
    simd_store(arr.0.as_mut_ptr() as *mut Simd, a);
    assert_eq!(arr.0, b);
}

#[target_feature(enable = "sse2")]
#[allow(dead_code)]
pub unsafe fn halfsimd_assert_vec_eq(a: HalfSimd, b: [i8; L]) {
    #[repr(align(16))]
    struct A([i8; L * HALFSIMD_MUL]);

    let mut arr = A([0i8; L * HALFSIMD_MUL]);
    halfsimd_store(arr.0.as_mut_ptr() as *mut HalfSimd, a);
    assert_eq!(&arr.0[..L], b);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_smoke() {
        #[target_feature(enable = "sse2")]
        unsafe fn inner() {
            #[repr(align(16))]
            struct A([i16; L]);

            let test = A([1, 2, 3, 4, 5, 6, 7, 8]);
            let test_rev = A([8, 7, 6, 5, 4, 3, 2, 1]);
            let test_mask = A([0, -1, 0, -1, 0, -1, 0, -1]);
            let vec0 = simd_load(test.0.as_ptr() as *const Simd);
            let vec0_rev = simd_load(test_rev.0.as_ptr() as *const Simd);
            let vec0_mask = simd_load(test_mask.0.as_ptr() as *const Simd);

            let mut vec1 = simd_sl_i16!(vec0, vec0, 1);
            simd_assert_vec_eq(vec1, [8, 1, 2, 3, 4, 5, 6, 7]);

            vec1 = simd_sr_i16!(vec0, vec0, 1);
            simd_assert_vec_eq(vec1, [2, 3, 4, 5, 6, 7, 8, 1]);

            vec1 = simd_adds_i16(vec0, vec0);
            simd_assert_vec_eq(vec1, [2, 4, 6, 8, 10, 12, 14, 16]);

            vec1 = simd_subs_i16(vec0, vec0);
            simd_assert_vec_eq(vec1, [0, 0, 0, 0, 0, 0, 0, 0]);

            vec1 = simd_max_i16(vec0, vec0_rev);
            simd_assert_vec_eq(vec1, [8, 7, 6, 5, 5, 6, 7, 8]);

            vec1 = simd_cmpeq_i16(vec0, vec0_rev);
            simd_assert_vec_eq(vec1, [0, 0, 0, 0, 0, 0, 0, 0]);

            vec1 = simd_cmpeq_i16(vec0, vec0);
            simd_assert_vec_eq(vec1, [-1, -1, -1, -1, -1, -1, -1, -1]);

            vec1 = simd_cmpgt_i16(vec0, vec0_rev);
            simd_assert_vec_eq(vec1, [0, 0, 0, 0, -1, -1, -1, -1]);

            vec1 = simd_blend_i8(vec0, vec0_rev, vec0_mask);
            simd_assert_vec_eq(vec1, [1, 7, 3, 5, 5, 3, 7, 1]);

            let mut val = simd_extract_i16!(vec0, 0);
            assert_eq!(val, 1);

            val = simd_slow_extract_i16(vec0, 0);
            assert_eq!(val, 1);

            vec1 = simd_insert_i16!(vec0, 0, 2);
            simd_assert_vec_eq(vec1, [1, 2, 0, 4, 5, 6, 7, 8]);

            let val1 = simd_movemask_i8(vec0_mask);
            assert_eq!(val1, 0b1100110011001100);

            vec1 = simd_sllz_i16!(vec0, 1);
            simd_assert_vec_eq(vec1, [0, 1, 2, 3, 4, 5, 6, 7]);

            vec1 = simd_broadcasthi_i16(vec0);
            simd_assert_vec_eq(vec1, [8, 8, 8, 8, 8, 8, 8, 8]);

            val = simd_hmax_i16(vec0);
            assert_eq!(val, 8);

            let zeros = simd_set1_i16(ZERO);
            val = simd_prefix_hadd_i16!(simd_adds_i16(vec0, zeros), 4);
            assert_eq!(val, 10);

            val = simd_prefix_hmax_i16!(vec0, 4);
            assert_eq!(val, 4);

            val = simd_suffix_hmax_i16!(vec0, 4);
            assert_eq!(val, 8);

            let val2 = simd_hargmax_i16(vec0, 4);
            assert_eq!(val2, 3);
        }
        unsafe { inner(); }
    }

    #[test]
    fn test_prefix_scan() {
        #[target_feature(enable = "sse2")]
        unsafe fn inner() {
            #[repr(align(16))]
            struct A([i16; L]);

            let vec = A([8, 9, 10, 15, 12, 13, 14, 11]);
            let gap = simd_set1_i16(0);
            let (_, consts) = get_prefix_scan_consts(gap);
            let res = simd_prefix_scan_i16(simd_load(vec.0.as_ptr() as *const Simd), gap, consts);
            simd_assert_vec_eq(res, [8, 9, 10, 15, 15, 15, 15, 15]);

            let vec = A([8, 9, 10, 15, 12, 13, 14, 11]);
            let gap = simd_set1_i16(-1);
            let (_, consts) = get_prefix_scan_consts(gap);
            let res = simd_prefix_scan_i16(simd_load(vec.0.as_ptr() as *const Simd), gap, consts);
            simd_assert_vec_eq(res, [8, 9, 10, 15, 14, 13, 14, 13]);
        }
        unsafe { inner(); }
    }
}
