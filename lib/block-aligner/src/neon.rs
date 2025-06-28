use std::arch::aarch64::*;

pub type Simd = int16x8_t;
pub type HalfSimd = int8x8_t;
pub type LutSimd = int8x16_t;
pub type TraceType = i16;
/// Number of 16-bit lanes in a SIMD vector.
pub const L: usize = 8;
pub const L_BYTES: usize = L * 2;
pub const HALFSIMD_MUL: usize = 1;
// using min = 0 is faster, but restricts range of scores (and restricts the max block size)
pub const ZERO: i16 = 1 << 14;
pub const MIN: i16 = 0;

// No non-temporal store in Neon
#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn store_trace(ptr: *mut TraceType, trace: TraceType) { *ptr = trace; }

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn simd_adds_i16(a: Simd, b: Simd) -> Simd { vqaddq_s16(a, b) }

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn simd_subs_i16(a: Simd, b: Simd) -> Simd { vqsubq_s16(a, b) }

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn simd_max_i16(a: Simd, b: Simd) -> Simd { vmaxq_s16(a, b) }

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn simd_cmpeq_i16(a: Simd, b: Simd) -> Simd { vreinterpretq_s16_u16(vceqq_s16(a, b)) }

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn simd_cmpgt_i16(a: Simd, b: Simd) -> Simd { vreinterpretq_s16_u16(vcgtq_s16(a, b)) }

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn simd_blend_i8(a: Simd, b: Simd, mask: Simd) -> Simd {
    // assume that each element in mask is either 0 or -1
    vbslq_s16(vreinterpretq_u16_s16(mask), b, a)
}

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn simd_load(ptr: *const Simd) -> Simd { vld1q_s16(ptr as *const i16) }

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn simd_loadu(ptr: *const Simd) -> Simd { vld1q_s16(ptr as *const i16) }

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn simd_store(ptr: *mut Simd, a: Simd) { vst1q_s16(ptr as *mut i16, a) }

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn simd_set1_i16(v: i16) -> Simd { vdupq_n_s16(v) }

#[macro_export]
#[doc(hidden)]
macro_rules! simd_extract_i16 {
    ($a:expr, $num:expr) => {
        {
            debug_assert!($num < L);
            #[cfg(target_arch = "aarch64")]
            use std::arch::aarch64::*;
            vgetq_lane_s16($a, $num as i32)
        }
    };
}

#[macro_export]
#[doc(hidden)]
macro_rules! simd_insert_i16 {
    ($a:expr, $v:expr, $num:expr) => {
        {
            debug_assert!($num < L);
            #[cfg(target_arch = "aarch64")]
            use std::arch::aarch64::*;
            vsetq_lane_s16($v, $a, $num as i32)
        }
    };
}

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn simd_movemask_i8(a: Simd) -> u16 {
    // assume that each byte is either 0 or -1
    static POW2: [u8; 16] = [
        1 << 0, 1 << 1, 1 << 2, 1 << 3, 1 << 4, 1 << 5, 1 << 6, 1 << 7,
        1 << 0, 1 << 1, 1 << 2, 1 << 3, 1 << 4, 1 << 5, 1 << 6, 1 << 7
    ];
    let mask = vld1q_u8(POW2.as_ptr());
    let masked = vandq_u8(vreinterpretq_u8_s16(a), mask);
    let lo = vaddv_u8(vget_low_u8(masked)) as u16;
    let hi = vaddv_u8(vget_high_u8(masked)) as u16;
    (hi << 8) | lo
}

#[macro_export]
#[doc(hidden)]
macro_rules! simd_sl_i16 {
    ($a:expr, $b:expr, $num:expr) => {
        {
            debug_assert!($num <= L);
            #[cfg(target_arch = "aarch64")]
            use std::arch::aarch64::*;
            vextq_s16($b, $a, (L - $num) as i32)
        }
    };
}

#[macro_export]
#[doc(hidden)]
macro_rules! simd_sr_i16 {
    ($a:expr, $b:expr, $num:expr) => {
        {
            debug_assert!($num <= L);
            #[cfg(target_arch = "aarch64")]
            use std::arch::aarch64::*;
            if $num == L {
                $a
            } else {
                vextq_s16($b, $a, $num as i32)
            }
        }
    };
}

// hardcoded to STEP = 8
#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn simd_step(a: Simd, b: Simd) -> Simd {
    a
}

macro_rules! simd_sllz_i16 {
    ($a:expr, $num:expr) => {
        {
            simd_sl_i16!($a, simd_set1_i16(0), $num)
        }
    };
}

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn simd_broadcasthi_i16(v: Simd) -> Simd {
    vdupq_laneq_s16(v, 7)
}

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn simd_slow_extract_i16(v: Simd, i: usize) -> i16 {
    debug_assert!(i < L);

    #[repr(align(16))]
    struct A([i16; L]);

    let mut a = A([0i16; L]);
    simd_store(a.0.as_mut_ptr() as *mut Simd, v);
    *a.0.as_ptr().add(i)
}

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn simd_hmax_i16(v: Simd) -> i16 { vmaxvq_s16(v) }

#[macro_export]
#[doc(hidden)]
macro_rules! simd_prefix_hadd_i16 {
    ($a:expr, $num:expr) => {
        {
            debug_assert!($num <= L);
            let mut v = simd_subs_i16($a, simd_set1_i16(ZERO));
            if $num > 4 {
                v = simd_adds_i16(v, simd_sr_i16!(v, v, 4));
            }
            if $num > 2 {
                v = simd_adds_i16(v, simd_sr_i16!(v, v, 2));
            }
            if $num > 1 {
                v = simd_adds_i16(v, simd_sr_i16!(v, v, 1));
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
            let mut v = $a;
            if $num > 4 {
                v = simd_max_i16(v, simd_sr_i16!(v, v, 4));
            }
            if $num > 2 {
                v = simd_max_i16(v, simd_sr_i16!(v, v, 2));
            }
            if $num > 1 {
                v = simd_max_i16(v, simd_sr_i16!(v, v, 1));
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
            let mut v = $a;
            if $num > 4 {
                v = simd_max_i16(v, simd_sl_i16!(v, v, 4));
            }
            if $num > 2 {
                v = simd_max_i16(v, simd_sl_i16!(v, v, 2));
            }
            if $num > 1 {
                v = simd_max_i16(v, simd_sl_i16!(v, v, 1));
            }
            simd_extract_i16!(v, 7)
        }
    };
}

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn simd_hargmax_i16(v: Simd, max: i16) -> usize {
    let v2 = simd_cmpeq_i16(v, simd_set1_i16(max));
    (simd_movemask_i8(v2).trailing_zeros() as usize) / 2
}

#[target_feature(enable = "neon")]
#[inline]
#[allow(non_snake_case)]
#[allow(dead_code)]
pub unsafe fn simd_naive_prefix_scan_i16(R_max: Simd, gap_cost: Simd, _gap_cost_lane: PrefixScanConsts) -> Simd {
    let mut curr = R_max;

    for _i in 0..(L - 1) {
        let prev = curr;
        curr = simd_sllz_i16!(curr, 1);
        curr = simd_adds_i16(curr, gap_cost);
        curr = simd_max_i16(curr, prev);
    }

    curr
}

pub type PrefixScanConsts = ();

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn get_prefix_scan_consts(gap: Simd) -> (Simd, PrefixScanConsts) {
    let mut shift1 = simd_sllz_i16!(gap, 1);
    shift1 = simd_adds_i16(shift1, gap);
    let mut shift2 = simd_sllz_i16!(shift1, 2);
    shift2 = simd_adds_i16(shift2, shift1);
    let mut shift4 = simd_sllz_i16!(shift2, 4);
    shift4 = simd_adds_i16(shift4, shift2);

    (shift4, ())
}

#[target_feature(enable = "neon")]
#[inline]
#[allow(non_snake_case)]
pub unsafe fn simd_prefix_scan_i16(R_max: Simd, gap_cost: Simd, _gap_cost_lane: PrefixScanConsts) -> Simd {
    let mut shift1 = simd_sllz_i16!(R_max, 1);
    shift1 = simd_adds_i16(shift1, gap_cost);
    shift1 = simd_max_i16(shift1, R_max);
    let mut shift2 = simd_sllz_i16!(shift1, 2);
    shift2 = simd_adds_i16(shift2, vshlq_n_s16(gap_cost, 1));
    shift2 = simd_max_i16(shift1, shift2);
    let mut shift4 = simd_sllz_i16!(shift2, 4);
    shift4 = simd_adds_i16(shift4, vshlq_n_s16(gap_cost, 2));
    shift4 = simd_max_i16(shift2, shift4);

    shift4
}

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn halfsimd_lookup2_i16(lut1: LutSimd, lut2: LutSimd, v: HalfSimd) -> Simd {
    let v2 = vcombine_u8(vreinterpret_u8_s8(v), vdup_n_u8(0));
    let c = vget_low_s8(vqtbl2q_s8(int8x16x2_t(lut1, lut2), v2));
    vmovl_s8(c)
}

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn halfsimd_lookup1_i16(lut: LutSimd, v: HalfSimd) -> Simd {
    let v2 = vcombine_u8(vand_u8(vreinterpret_u8_s8(v), vdup_n_u8(0b1111)), vdup_n_u8(0));
    let c = vget_low_s8(vqtbl1q_s8(lut, v2));
    vmovl_s8(c)
}

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn halfsimd_lookup_bytes_i16(match_scores: HalfSimd, mismatch_scores: HalfSimd, a: HalfSimd, b: HalfSimd) -> Simd {
    let mask = vceq_s8(a, b);
    let c = vbsl_s8(mask, match_scores, mismatch_scores);
    vmovl_s8(c)
}

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn halfsimd_load(ptr: *const HalfSimd) -> HalfSimd { vld1_s8(ptr as *const i8) }

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn halfsimd_loadu(ptr: *const HalfSimd) -> HalfSimd { vld1_s8(ptr as *const i8) }

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn lutsimd_load(ptr: *const LutSimd) -> LutSimd { vld1q_s8(ptr as *const i8) }

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn lutsimd_loadu(ptr: *const LutSimd) -> LutSimd { vld1q_s8(ptr as *const i8) }

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn halfsimd_store(ptr: *mut HalfSimd, a: HalfSimd) { vst1_s8(ptr as *mut i8, a) }

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn halfsimd_sub_i8(a: HalfSimd, b: HalfSimd) -> HalfSimd { vsub_s8(a, b) }

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn halfsimd_set1_i8(v: i8) -> HalfSimd { vdup_n_s8(v) }

#[target_feature(enable = "neon")]
#[inline]
pub unsafe fn halfsimd_get_idx(i: usize) -> usize { i }

#[macro_export]
#[doc(hidden)]
macro_rules! halfsimd_sr_i8 {
    ($a:expr, $b:expr, $num:expr) => {
        {
            debug_assert!($num <= L);
            #[cfg(target_arch = "aarch64")]
            use std::arch::aarch64::*;
            vext_s8($b, $a, $num as i32)
        }
    };
}

#[target_feature(enable = "neon")]
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

#[target_feature(enable = "neon")]
#[allow(dead_code)]
pub unsafe fn halfsimd_dbg_i8(v: HalfSimd) {
    #[repr(align(8))]
    struct A([i8; L]);

    let mut a = A([0i8; L]);
    halfsimd_store(a.0.as_mut_ptr() as *mut HalfSimd, v);

    for i in (0..a.0.len()).rev() {
        print!("{:3} ", a.0[i]);
    }
    println!();
}

#[target_feature(enable = "neon")]
#[allow(dead_code)]
pub unsafe fn simd_assert_vec_eq(a: Simd, b: [i16; L]) {
    #[repr(align(16))]
    struct A([i16; L]);

    let mut arr = A([0i16; L]);
    simd_store(arr.0.as_mut_ptr() as *mut Simd, a);
    assert_eq!(arr.0, b);
}

#[target_feature(enable = "neon")]
#[allow(dead_code)]
pub unsafe fn halfsimd_assert_vec_eq(a: HalfSimd, b: [i8; L]) {
    #[repr(align(8))]
    struct A([i8; L]);

    let mut arr = A([0i8; L]);
    halfsimd_store(arr.0.as_mut_ptr() as *mut HalfSimd, a);
    assert_eq!(arr.0, b);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_smoke() {
        #[target_feature(enable = "neon")]
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
        #[target_feature(enable = "neon")]
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
