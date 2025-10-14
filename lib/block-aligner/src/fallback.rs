#![allow(dead_code)]

use core::{cmp, ptr};

pub type Simd = u128;     // 8 × i16 packed little-endian
pub type HalfSimd = u128; // 16 × i8 packed little-endian
pub type LutSimd = u128;  // 16 × i8 packed little-endian
pub type TraceType = i16;

/// Number of 16-bit lanes in a SIMD vector.
pub const L: usize = 8;
pub const L_BYTES: usize = L * 2;
pub const HALFSIMD_MUL: usize = 2;
// using min = 0 is faster, but restricts range of scores (and restricts the max block size)
pub const ZERO: i16 = 1 << 14;
pub const MIN: i16 = 0;

#[inline]
pub fn to_i16x8(x: Simd) -> [i16; L] {
    let b = x.to_le_bytes();
    let mut a = [0i16; L];
    let mut i = 0;
    while i < L {
        a[i] = i16::from_le_bytes([b[2 * i], b[2 * i + 1]]);
        i += 1;
    }
    a
}

#[inline]
pub fn from_i16x8(a: [i16; L]) -> Simd {
    let mut b = [0u8; 16];
    let mut i = 0;
    while i < L {
        let bytes = a[i].to_le_bytes();
        b[2 * i] = bytes[0];
        b[2 * i + 1] = bytes[1];
        i += 1;
    }
    Simd::from_le_bytes(b)
}

#[inline]
pub fn to_u8x16(x: u128) -> [u8; 16] {
    x.to_le_bytes()
}

#[inline]
pub fn from_u8x16(a: [u8; 16]) -> u128 {
    u128::from_le_bytes(a)
}

#[inline]
pub fn bitselect(a: u128, b: u128, mask: u128) -> u128 {
    (a & !mask) | (b & mask)
}

#[inline]
pub fn i16_shl_wrap(v: i16, n: u32) -> i16 {
    (((v as u16) << n) as u16) as i16
}

#[inline]
pub unsafe fn store_trace(ptr_: *mut TraceType, trace: TraceType) { *ptr_ = trace; }

#[inline]
pub unsafe fn simd_load(ptr_: *const Simd) -> Simd { ptr::read(ptr_) }

#[inline]
pub unsafe fn simd_loadu(ptr_: *const Simd) -> Simd { ptr::read_unaligned(ptr_) }

#[inline]
pub unsafe fn simd_store(ptr_: *mut Simd, a: Simd) { ptr::write(ptr_, a) }

#[inline]
pub unsafe fn halfsimd_load(ptr_: *const HalfSimd) -> HalfSimd { ptr::read(ptr_) }

#[inline]
pub unsafe fn halfsimd_loadu(ptr_: *const HalfSimd) -> HalfSimd { ptr::read_unaligned(ptr_) }

#[inline]
pub unsafe fn lutsimd_load(ptr_: *const LutSimd) -> LutSimd { ptr::read(ptr_) }

#[inline]
pub unsafe fn lutsimd_loadu(ptr_: *const LutSimd) -> LutSimd { ptr::read_unaligned(ptr_) }

#[inline]
pub unsafe fn halfsimd_store(ptr_: *mut HalfSimd, a: HalfSimd) { ptr::write(ptr_, a) }

#[inline]
pub unsafe fn simd_set1_i16(v: i16) -> Simd {
    from_i16x8([v; L])
}

#[inline]
pub unsafe fn simd_adds_i16(a: Simd, b: Simd) -> Simd {
    let aa = to_i16x8(a);
    let bb = to_i16x8(b);
    let mut r = [0i16; L];
    let mut i = 0;
    while i < L {
        r[i] = aa[i].saturating_add(bb[i]);
        i += 1;
    }
    from_i16x8(r)
}

#[inline]
pub unsafe fn simd_subs_i16(a: Simd, b: Simd) -> Simd {
    let aa = to_i16x8(a);
    let bb = to_i16x8(b);
    let mut r = [0i16; L];
    let mut i = 0;
    while i < L {
        r[i] = aa[i].saturating_sub(bb[i]);
        i += 1;
    }
    from_i16x8(r)
}

#[inline]
pub unsafe fn simd_max_i16(a: Simd, b: Simd) -> Simd {
    let aa = to_i16x8(a);
    let bb = to_i16x8(b);
    let mut r = [0i16; L];
    let mut i = 0;
    while i < L {
        r[i] = cmp::max(aa[i], bb[i]);
        i += 1;
    }
    from_i16x8(r)
}

#[inline]
pub unsafe fn simd_cmpeq_i16(a: Simd, b: Simd) -> Simd {
    let aa = to_i16x8(a);
    let bb = to_i16x8(b);
    let mut r = [0i16; L];
    let mut i = 0;
    while i < L {
        r[i] = if aa[i] == bb[i] { -1 } else { 0 };
        i += 1;
    }
    from_i16x8(r)
}

#[inline]
pub unsafe fn simd_cmpgt_i16(a: Simd, b: Simd) -> Simd {
    let aa = to_i16x8(a);
    let bb = to_i16x8(b);
    let mut r = [0i16; L];
    let mut i = 0;
    while i < L {
        r[i] = if aa[i] > bb[i] { -1 } else { 0 };
        i += 1;
    }
    from_i16x8(r)
}

// Per-byte blend: (~mask & a) | (mask & b)
#[inline]
pub unsafe fn simd_blend_i8(a: Simd, b: Simd, mask: Simd) -> Simd {
    bitselect(a, b, mask)
}

#[inline]
pub unsafe fn simd_movemask_i8(a: Simd) -> u16 {
    let bytes = to_u8x16(a);
    let mut m = 0u16;
    let mut i = 0;
    while i < 16 {
        // Take the sign bit of each byte (like SSE2)
        let bit = ((bytes[i] as i8) < 0) as u16;
        m |= bit << i;
        i += 1;
    }
    m
}

// hardcoded to STEP = 8
#[inline]
pub unsafe fn simd_step(a: Simd, _b: Simd) -> Simd { a }

// Shift left by lanes; bring in from b's high end
#[macro_export]
#[doc(hidden)]
macro_rules! simd_sl_i16 {
    ($a:expr, $b:expr, $num:expr) => {{
        debug_assert!($num <= L);
        let aa = crate::fallback::to_i16x8($a);
        let bb = crate::fallback::to_i16x8($b);
        let mut r = [0i16; L];
        let mut i = 0usize;
        while i < L {
            r[i] = if i < $num { bb[L - $num + i] } else { aa[i - $num] };
            i += 1;
        }
        crate::fallback::from_i16x8(r)
    }};
}

// Shift right by lanes; bring in from b's low end
#[macro_export]
#[doc(hidden)]
macro_rules! simd_sr_i16 {
    ($a:expr, $b:expr, $num:expr) => {{
        debug_assert!($num <= L);
        let aa = crate::fallback::to_i16x8($a);
        let bb = crate::fallback::to_i16x8($b);
        let mut r = [0i16; L];
        let mut i = 0usize;
        while i < L {
            r[i] = if i + $num < L { aa[i + $num] } else { bb[i + $num - L] };
            i += 1;
        }
        crate::fallback::from_i16x8(r)
    }};
}

// shift in zeros (by lanes)
macro_rules! simd_sllz_i16 {
    ($a:expr, $num:expr) => {{
        debug_assert!($num < L);
        simd_sl_i16!($a, unsafe { simd_set1_i16(0) }, $num)
    }};
}

// broadcast highest 16-bit element to the whole vector
#[inline]
pub unsafe fn simd_broadcasthi_i16(v: Simd) -> Simd {
    let a = to_i16x8(v);
    from_i16x8([a[7]; L])
}

#[inline]
pub unsafe fn simd_slow_extract_i16(v: Simd, i: usize) -> i16 {
    debug_assert!(i < L);
    to_i16x8(v)[i]
}

#[inline]
pub unsafe fn simd_hmax_i16(v: Simd) -> i16 {
    let a = to_i16x8(v);
    let mut m = a[0];
    let mut i = 1;
    while i < L {
        if a[i] > m { m = a[i]; }
        i += 1;
    }
    m
}

#[macro_export]
#[doc(hidden)]
macro_rules! simd_extract_i16 {
    ($a:expr, $num:expr) => {{
        debug_assert!($num < L);
        $crate::fallback::to_i16x8($a)[$num]
    }};
}

#[macro_export]
#[doc(hidden)]
macro_rules! simd_insert_i16 {
    ($a:expr, $v:expr, $num:expr) => {{
        debug_assert!($num < L);
        let mut tmp = $crate::fallback::to_i16x8($a);
        tmp[$num] = $v as i16;
        $crate::fallback::from_i16x8(tmp)
    }};
}

#[macro_export]
#[doc(hidden)]
macro_rules! simd_prefix_hadd_i16 {
    ($a:expr, $num:expr) => {{
        debug_assert!($num <= L);
        let mut s: i32 = 0;
        let aa = $crate::fallback::to_i16x8($a);
        let mut i = 0usize;
        while i < $num {
            s += (aa[i] as i32) - (ZERO as i32);
            i += 1;
        }
        s.clamp(i16::MIN as i32, i16::MAX as i32) as i16
    }};
}

#[macro_export]
#[doc(hidden)]
macro_rules! simd_prefix_hmax_i16 {
    ($a:expr, $num:expr) => {{
        debug_assert!($num <= L);
        let aa = $crate::fallback::to_i16x8($a);
        let mut m = aa[0];
        let mut i = 1usize;
        while i < $num {
            if aa[i] > m { m = aa[i]; }
            i += 1;
        }
        m
    }};
}

#[macro_export]
#[doc(hidden)]
macro_rules! simd_suffix_hmax_i16 {
    ($a:expr, $num:expr) => {{
        debug_assert!($num <= L);
        let aa = $crate::fallback::to_i16x8($a);
        let start = L - $num;
        let mut m = aa[start];
        let mut i = start + 1;
        while i < L {
            if aa[i] > m { m = aa[i]; }
            i += 1;
        }
        m
    }};
}

#[inline]
pub unsafe fn simd_hargmax_i16(v: Simd, max: i16) -> usize {
    let a = to_i16x8(v);
    let mut i = 0usize;
    while i < L {
        if a[i] == max { return i; }
        i += 1;
    }
    L 
}

pub type PrefixScanConsts = ();

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

#[inline]
pub unsafe fn simd_naive_prefix_scan_i16(R_max: Simd, gap_cost: Simd, _gap_cost_lane: PrefixScanConsts) -> Simd {
    let mut curr = R_max;
    let mut i = 0usize;
    while i < (L - 1) {
        let prev = curr;
        curr = simd_sllz_i16!(curr, 1);
        curr = simd_adds_i16(curr, gap_cost);
        curr = simd_max_i16(curr, prev);
        i += 1;
    }
    curr
}

#[inline]
pub unsafe fn simd_prefix_scan_i16(R_max: Simd, gap_cost: Simd, _c: PrefixScanConsts) -> Simd {
    let mut shift1 = simd_sllz_i16!(R_max, 1);
    shift1 = simd_adds_i16(shift1, gap_cost);
    shift1 = simd_max_i16(shift1, R_max);

    let gc_arr = to_i16x8(gap_cost);
    let mut gc_shl1 = [0i16; L];
    let mut gc_shl2 = [0i16; L];
    let mut i = 0usize;
    while i < L {
        gc_shl1[i] = i16_shl_wrap(gc_arr[i], 1);
        gc_shl2[i] = i16_shl_wrap(gc_arr[i], 2);
        i += 1;
    }
    let gc_shl1 = from_i16x8(gc_shl1);
    let gc_shl2 = from_i16x8(gc_shl2);

    let mut shift2 = simd_sllz_i16!(shift1, 2);
    shift2 = simd_adds_i16(shift2, gc_shl1);
    shift2 = simd_max_i16(shift2, shift1);

    let mut shift4 = simd_sllz_i16!(shift2, 4);
    shift4 = simd_adds_i16(shift4, gc_shl2);
    shift4 = simd_max_i16(shift4, shift2);

    shift4
}

#[inline]
pub unsafe fn halfsimd_sub_i8(a: HalfSimd, b: HalfSimd) -> HalfSimd {
    let aa = to_u8x16(a);
    let bb = to_u8x16(b);
    let mut r = [0u8; 16];
    let mut i = 0usize;
    while i < 16 {
        r[i] = aa[i].wrapping_sub(bb[i]);
        i += 1;
    }
    from_u8x16(r)
}

#[inline]
pub unsafe fn halfsimd_set1_i8(v: i8) -> HalfSimd {
    let b = v as u8;
    from_u8x16([b; 16])
}

#[inline]
pub unsafe fn halfsimd_get_idx(i: usize) -> usize { i + i / L * L }

#[macro_export]
#[doc(hidden)]
macro_rules! halfsimd_sr_i8 {
    ($a:expr, $b:expr, $num:expr) => {{
        debug_assert!($num <= L);
        let aa = $crate::to_u8x16($a);
        let bb = $crate::to_u8x16($b);
        let mut r = [0u8; 16];
        let mut i = 0usize;
        while i < L {
            r[i] = if i + $num < L { aa[i + $num] } else { bb[i + $num - L] };
            i += 1;
        }
        $crate::from_u8x16(r)
    }};
}

#[inline]
pub unsafe fn halfsimd_lookup2_i16(lut1: LutSimd, lut2: LutSimd, v: HalfSimd) -> Simd {
    let idx = to_u8x16(v);
    let t1 = to_u8x16(lut1);
    let t2 = to_u8x16(lut2);
    let mut table = [0u8; 32];
    table[..16].copy_from_slice(&t1);
    table[16..].copy_from_slice(&t2);

    let mut out = [0i16; L];
    let mut i = 0usize;
    while i < L {
        let j = idx[i] as usize;
        out[i] = table[j] as i8 as i16;
        i += 1;
    }
    from_i16x8(out)
}

#[inline]
pub unsafe fn halfsimd_lookup1_i16(lut: LutSimd, v: HalfSimd) -> Simd {
    let idx = to_u8x16(v);
    let t = to_u8x16(lut);
    let mut out = [0i16; L];
    let mut i = 0usize;
    while i < L {
        out[i] = t[(idx[i] & 0x0F) as usize] as i8 as i16;
        i += 1;
    }
    from_i16x8(out)
}

#[inline]
pub unsafe fn halfsimd_lookup_bytes_i16(
    match_scores: HalfSimd,
    mismatch_scores: HalfSimd,
    a: HalfSimd,
    b: HalfSimd
) -> Simd {
    let aa = to_u8x16(a);
    let bb = to_u8x16(b);
    let ms = to_u8x16(match_scores);
    let mms = to_u8x16(mismatch_scores);
    let mut c = [0u8; 16];
    let mut i = 0usize;
    while i < 16 {
        c[i] = if aa[i] == bb[i] { ms[i] } else { mms[i] };
        i += 1;
    }
    let mut out = [0i16; L];
    i = 0;
    while i < L {
        out[i] = c[i] as i8 as i16;
        i += 1;
    }
    from_i16x8(out)
}

#[inline]
pub unsafe fn simd_dbg_i16(v: Simd) {
    let a = to_i16x8(v);
    for i in (0..a.len()).rev() {
        print!("{:6} ", a[i]);
    }
    println!();
}

#[inline]
pub unsafe fn halfsimd_dbg_i8(v: HalfSimd) {
    let a = to_u8x16(v);
    for i in (0..a.len()).rev() {
        print!("{:3} ", a[i] as i8);
    }
    println!();
}

#[inline]
pub unsafe fn simd_assert_vec_eq(a: Simd, b: [i16; L]) {
    assert_eq!(to_i16x8(a), b);
}

#[inline]
pub unsafe fn halfsimd_assert_vec_eq(a: HalfSimd, b: [i8; L]) {
    let arr = to_u8x16(a);
    let mut low = [0i8; L];
    let mut i = 0usize;
    while i < L {
        low[i] = arr[i] as i8;
        i += 1;
    }
    assert_eq!(low, b);
}
