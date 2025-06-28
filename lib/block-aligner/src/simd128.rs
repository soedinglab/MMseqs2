use std::arch::wasm32::*;

pub type Simd = v128;
// no v64 type, so HalfSimd is just v128 with upper half ignored
pub type HalfSimd = v128;
pub type LutSimd = v128;
pub type TraceType = i16;
/// Number of 16-bit lanes in a SIMD vector.
pub const L: usize = 8;
pub const L_BYTES: usize = L * 2;
pub const HALFSIMD_MUL: usize = 2;
// using min = 0 is faster, but restricts range of scores (and restricts the max block size)
pub const ZERO: i16 = 1 << 14;
pub const MIN: i16 = 0;

// Note: SIMD vectors treated as little-endian

// No non-temporal store in WASM
#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn store_trace(ptr: *mut TraceType, trace: TraceType) { *ptr = trace; }

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn simd_adds_i16(a: Simd, b: Simd) -> Simd { i16x8_add_sat(a, b) }

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn simd_subs_i16(a: Simd, b: Simd) -> Simd { i16x8_sub_sat(a, b) }

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn simd_max_i16(a: Simd, b: Simd) -> Simd { i16x8_max(a, b) }

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn simd_cmpeq_i16(a: Simd, b: Simd) -> Simd { i16x8_eq(a, b) }

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn simd_cmpgt_i16(a: Simd, b: Simd) -> Simd { i16x8_gt(a, b) }

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn simd_blend_i8(a: Simd, b: Simd, mask: Simd) -> Simd { v128_bitselect(b, a, mask) }

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn simd_load(ptr: *const Simd) -> Simd { v128_load(ptr) }

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn simd_loadu(ptr: *const Simd) -> Simd { v128_load(ptr) }

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn simd_store(ptr: *mut Simd, a: Simd) { v128_store(ptr, a) }

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn simd_set1_i16(v: i16) -> Simd { i16x8_splat(v) }

#[macro_export]
#[doc(hidden)]
macro_rules! simd_extract_i16 {
    ($a:expr, $num:expr) => {
        {
            debug_assert!($num < L);
            use std::arch::wasm32::*;
            i16x8_extract_lane::<{ $num }>($a)
        }
    };
}

#[macro_export]
#[doc(hidden)]
macro_rules! simd_insert_i16 {
    ($a:expr, $v:expr, $num:expr) => {
        {
            debug_assert!($num < L);
            use std::arch::wasm32::*;
            i16x8_replace_lane::<{ $num }>($a, $v)
        }
    };
}

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn simd_movemask_i8(a: Simd) -> u16 {
    i8x16_bitmask(a) as u16
    /*const MUL: i64 = {
        let mut m = 0u64;
        m |= 1u64 << (0 - 0);
        m |= 1u64 << (8 - 1);
        m |= 1u64 << (16 - 2);
        m |= 1u64 << (24 - 3);
        m |= 1u64 << (32 - 4);
        m |= 1u64 << (40 - 5);
        m |= 1u64 << (48 - 6);
        m |= 1u64 << (56 - 7);
        m as i64
    };
    let b = i64x2_mul(v128_and(a, i8x16_splat(0b10000000u8 as i8)), i64x2_splat(MUL));
    let res1 = i8x16_extract_lane::<{ L * 2 - 1 }>(b) as u32;
    let res2 = i8x16_extract_lane::<{ L - 1 }>(b) as u32;
    (res1 << 8) | res2*/
}

#[macro_export]
#[doc(hidden)]
macro_rules! simd_sl_i16 {
    ($a:expr, $b:expr, $num:expr) => {
        {
            debug_assert!($num <= L);
            use std::arch::wasm32::*;
            i16x8_shuffle::<{ 8 - $num }, { 9 - $num }, { 10 - $num }, { 11 - $num }, { 12 - $num }, { 13 - $num }, { 14 - $num }, { 15 - $num }>($b, $a)
        }
    };
}

#[macro_export]
#[doc(hidden)]
macro_rules! simd_sr_i16 {
    ($a:expr, $b:expr, $num:expr) => {
        {
            debug_assert!($num <= L);
            use std::arch::wasm32::*;
            i16x8_shuffle::<{ 0 + $num }, { 1 + $num }, { 2 + $num }, { 3 + $num }, { 4 + $num }, { 5 + $num }, { 6 + $num }, { 7 + $num }>($b, $a)
        }
    };
}

// hardcoded to STEP = 8
#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn simd_step(a: Simd, _b: Simd) -> Simd {
    a
}

macro_rules! simd_sllz_i16 {
    ($a:expr, $num:expr) => {
        {
            simd_sl_i16!($a, simd_set1_i16(0), $num)
        }
    };
}

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn simd_broadcasthi_i16(v: Simd) -> Simd {
    i16x8_shuffle::<7, 7, 7, 7, 7, 7, 7, 7>(v, v)
}

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn simd_slow_extract_i16(v: Simd, i: usize) -> i16 {
    debug_assert!(i < L);

    #[repr(align(16))]
    struct A([i16; L]);

    let mut a = A([0i16; L]);
    simd_store(a.0.as_mut_ptr() as *mut Simd, v);
    *a.0.as_ptr().add(i)
}

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn simd_hmax_i16(v: Simd) -> i16 {
    let mut v2 = i16x8_max(v, simd_sr_i16!(v, v, 1));
    v2 = i16x8_max(v2, simd_sr_i16!(v2, v2, 2));
    v2 = i16x8_max(v2, simd_sr_i16!(v2, v2, 4));
    simd_extract_i16!(v2, 0)
}

#[macro_export]
#[doc(hidden)]
macro_rules! simd_prefix_hadd_i16 {
    ($a:expr, $num:expr) => {
        {
            debug_assert!($num <= L);
            use std::arch::wasm32::*;
            let mut v = i16x8_sub_sat($a, i16x8_splat(ZERO));
            if $num > 4 {
                v = i16x8_add_sat(v, simd_sr_i16!(v, v, 4));
            }
            if $num > 2 {
                v = i16x8_add_sat(v, simd_sr_i16!(v, v, 2));
            }
            if $num > 1 {
                v = i16x8_add_sat(v, simd_sr_i16!(v, v, 1));
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
            use std::arch::wasm32::*;
            let mut v = $a;
            if $num > 4 {
                v = i16x8_max(v, simd_sr_i16!(v, v, 4));
            }
            if $num > 2 {
                v = i16x8_max(v, simd_sr_i16!(v, v, 2));
            }
            if $num > 1 {
                v = i16x8_max(v, simd_sr_i16!(v, v, 1));
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
            use std::arch::wasm32::*;
            let mut v = $a;
            if $num > 4 {
                v = i16x8_max(v, simd_sl_i16!(v, v, 4));
            }
            if $num > 2 {
                v = i16x8_max(v, simd_sl_i16!(v, v, 2));
            }
            if $num > 1 {
                v = i16x8_max(v, simd_sl_i16!(v, v, 1));
            }
            simd_extract_i16!(v, 7)
        }
    };
}

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn simd_hargmax_i16(v: Simd, max: i16) -> usize {
    let v2 = i16x8_eq(v, i16x8_splat(max));
    (simd_movemask_i8(v2).trailing_zeros() as usize) / 2
}

#[target_feature(enable = "simd128")]
#[inline]
#[allow(non_snake_case)]
#[allow(dead_code)]
pub unsafe fn simd_naive_prefix_scan_i16(R_max: Simd, gap_cost: Simd, _gap_cost_lane: PrefixScanConsts) -> Simd {
    let mut curr = R_max;

    for _i in 0..(L - 1) {
        let prev = curr;
        curr = simd_sllz_i16!(curr, 1);
        curr = i16x8_add_sat(curr, gap_cost);
        curr = i16x8_max(curr, prev);
    }

    curr
}

pub type PrefixScanConsts = ();

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn get_prefix_scan_consts(gap: Simd) -> (Simd, PrefixScanConsts) {
    let mut shift1 = simd_sllz_i16!(gap, 1);
    shift1 = i16x8_add_sat(shift1, gap);
    let mut shift2 = simd_sllz_i16!(shift1, 2);
    shift2 = i16x8_add_sat(shift2, shift1);
    let mut shift4 = simd_sllz_i16!(shift2, 4);
    shift4 = i16x8_add_sat(shift4, shift2);

    (shift4, ())
}

#[target_feature(enable = "simd128")]
#[inline]
#[allow(non_snake_case)]
pub unsafe fn simd_prefix_scan_i16(R_max: Simd, gap_cost: Simd, _gap_cost_lane: PrefixScanConsts) -> Simd {
    let mut shift1 = simd_sllz_i16!(R_max, 1);
    shift1 = i16x8_add_sat(shift1, gap_cost);
    shift1 = i16x8_max(shift1, R_max);
    let mut shift2 = simd_sllz_i16!(shift1, 2);
    shift2 = i16x8_add_sat(shift2, i16x8_shl(gap_cost, 1));
    shift2 = i16x8_max(shift1, shift2);
    let mut shift4 = simd_sllz_i16!(shift2, 4);
    shift4 = i16x8_add_sat(shift4, i16x8_shl(gap_cost, 2));
    shift4 = i16x8_max(shift2, shift4);

    shift4
}

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn halfsimd_lookup2_i16(lut1: LutSimd, lut2: LutSimd, v: HalfSimd) -> Simd {
    // must use a mask to avoid zeroing lanes that are too large
    let mask = i8x16_splat(0b1111);
    let v_mask = v128_and(v, mask);
    let a = i8x16_swizzle(lut1, v_mask);
    let b = i8x16_swizzle(lut2, v_mask);
    let lut_mask = i8x16_gt(v, mask);
    let c = v128_bitselect(b, a, lut_mask);
    i16x8_extend_low_i8x16(c)
}

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn halfsimd_lookup1_i16(lut: LutSimd, v: HalfSimd) -> Simd {
    i16x8_extend_low_i8x16(i8x16_swizzle(lut, v128_and(v, i8x16_splat(0b1111))))
}

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn halfsimd_lookup_bytes_i16(match_scores: HalfSimd, mismatch_scores: HalfSimd, a: HalfSimd, b: HalfSimd) -> Simd {
    let mask = i8x16_eq(a, b);
    let c = v128_bitselect(match_scores, mismatch_scores, mask);
    i16x8_extend_low_i8x16(c)
}

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn halfsimd_load(ptr: *const HalfSimd) -> HalfSimd { v128_load(ptr) }

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn halfsimd_loadu(ptr: *const HalfSimd) -> HalfSimd { v128_load(ptr) }

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn lutsimd_load(ptr: *const LutSimd) -> LutSimd { v128_load(ptr) }

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn lutsimd_loadu(ptr: *const LutSimd) -> LutSimd { v128_load(ptr) }

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn halfsimd_store(ptr: *mut HalfSimd, a: HalfSimd) { v128_store(ptr, a) }

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn halfsimd_sub_i8(a: HalfSimd, b: HalfSimd) -> HalfSimd { i8x16_sub(a, b) }

#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn halfsimd_set1_i8(v: i8) -> HalfSimd { i8x16_splat(v) }

// only the low 8 bytes are out of each v128 for halfsimd
#[target_feature(enable = "simd128")]
#[inline]
pub unsafe fn halfsimd_get_idx(i: usize) -> usize { i + i / L * L }

#[macro_export]
#[doc(hidden)]
macro_rules! halfsimd_sr_i8 {
    ($a:expr, $b:expr, $num:expr) => {
        {
            debug_assert!($num <= L);
            use std::arch::wasm32::*;
            // special indexing to skip over the high 8 bytes that are unused
            const fn get_idx(i: usize) -> usize { if i >= L { i + L } else { i } }
            i8x16_shuffle::<
                { get_idx(0 + $num) }, { get_idx(1 + $num) }, { get_idx(2 + $num) }, { get_idx(3 + $num) },
                { get_idx(4 + $num) }, { get_idx(5 + $num) }, { get_idx(6 + $num) }, { get_idx(7 + $num) },
                8, 9, 10, 11,
                12, 13, 14, 15
            >($b, $a)
        }
    };
}

#[target_feature(enable = "simd128")]
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

#[target_feature(enable = "simd128")]
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

#[target_feature(enable = "simd128")]
#[allow(dead_code)]
pub unsafe fn simd_assert_vec_eq(a: Simd, b: [i16; L]) {
    #[repr(align(16))]
    struct A([i16; L]);

    let mut arr = A([0i16; L]);
    simd_store(arr.0.as_mut_ptr() as *mut Simd, a);
    assert_eq!(arr.0, b);
}

#[target_feature(enable = "simd128")]
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
        #[target_feature(enable = "simd128")]
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
    fn test_endianness() {
        #[target_feature(enable = "simd128")]
        unsafe fn inner() {
            #[repr(align(16))]
            struct A([i16; L]);

            let vec = A([1, 2, 3, 4, 5, 6, 7, 8]);
            let vec = simd_load(vec.0.as_ptr() as *const Simd);
            let res = simd_sl_i16!(vec, vec, 1);
            simd_assert_vec_eq(res, [8, 1, 2, 3, 4, 5, 6, 7]);

            let vec = A([1, 2, 3, 4, 5, 6, 7, 8]);
            let vec = simd_load(vec.0.as_ptr() as *const Simd);
            let res = simd_sr_i16!(vec, vec, 1);
            simd_assert_vec_eq(res, [2, 3, 4, 5, 6, 7, 8, 1]);

            #[repr(align(16))]
            struct B([i8; L * HALFSIMD_MUL]);

            let vec = B([1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0]);
            let vec = halfsimd_load(vec.0.as_ptr() as *const HalfSimd);
            let res = halfsimd_sr_i8!(vec, vec, 1);
            halfsimd_assert_vec_eq(res, [2, 3, 4, 5, 6, 7, 8, 1]);

            simd_assert_vec_eq(simd_adds_i16(simd_set1_i16(i16::MIN), simd_set1_i16(i16::MIN)), [i16::MIN; 8]);
            simd_assert_vec_eq(simd_adds_i16(simd_set1_i16(i16::MAX), simd_set1_i16(i16::MIN)), [-1; 8]);
            simd_assert_vec_eq(simd_subs_i16(simd_set1_i16(i16::MAX), simd_set1_i16(i16::MIN)), [i16::MAX; 8]);
        }
        unsafe { inner(); }
    }

    #[test]
    fn test_prefix_scan() {
        #[target_feature(enable = "simd128")]
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
