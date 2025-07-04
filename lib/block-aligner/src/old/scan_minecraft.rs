#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use crate::avx2::*;

#[cfg(target_arch = "wasm32")]
use crate::simd128::*;

use crate::scores::*;

use std::{cmp, ptr, i16};
use std::marker::PhantomData;

const NULL: u8 = b'A' + 26u8; // this null byte value works for both amino acids and nucleotides

// Notes:
//
// BLOSUM62 matrix max = 11, min = -4; gap open = -11 (includes extension), gap extend = -1
//
// R[i][j] = max(R[i - 1][j] + gap_extend, D[i - 1][j] + gap_open)
// C[i][j] = max(C[i][j - 1] + gap_extend, D[i][j - 1] + gap_open)
// D[i][j] = max(D[i - 1][j - 1] + matrix[query[i]][reference[j]], R[i][j], C[i][j])
//
// indexing (we want to calculate D11):
//      x0   x1
//    +--------
// 0x | 00   01
// 1x | 10   11
//
// note that 'x' represents any bit
//
// Each block is made up of vertical SIMD vectors of length 8 or 16 16-bit integers.

// TODO: create matrices with const fn

pub struct Block<'a, P: ScoreParams, M: 'a + Matrix, const TRACE: bool, const X_DROP: bool> {
    res: AlignResult,
    trace: Trace,
    query: &'a PaddedBytes,
    i: usize,
    reference: &'a PaddedBytes,
    j: usize,
    matrix: &'a M,
    x_drop: i32,
    _phantom: PhantomData<P>
}

impl<'a, P: ScoreParams, M: 'a + Matrix, const TRACE: bool, const X_DROP: bool> Block<'a, P, M, { TRACE }, { X_DROP }> {
    const EVEN_BITS: u32 = 0x55555555u32;

    /// Adaptive banded alignment.
    ///
    /// The x drop option indicates whether to terminate the alignment process early when
    /// the max score in the current band drops below the max score encountered so far. If
    /// x drop is not enabled, then the band will keep shifting until the end of the reference
    /// string is reached.
    ///
    /// Limitations:
    /// 1. Requires x86 AVX2 or WASM SIMD support.
    /// 2. The reference and the query can only contain uppercase alphabetical characters.
    /// 3. The actual size of the band is K + 1 rounded up to the next multiple of the
    ///    vector length of 16 (for x86 AVX2) or 8 (for WASM SIMD).
    pub fn align(query: &'a PaddedBytes, reference: &'a PaddedBytes, matrix: &'a M, x_drop: i32) -> Self {
        assert!(P::GAP_OPEN <= P::GAP_EXTEND);

        if X_DROP {
            assert!(x_drop >= 0);
        }

        let mut a = Self {
            res: AlignResult { score: 0, query_idx: 0, reference_idx: 0 },
            trace: if TRACE { Trace::new(query.len(), reference.len()) } else { Trace::new(0, 0) },
            query,
            i: 0,
            reference,
            j: 0,
            matrix,
            x_drop,
            _phantom: PhantomData
        };

        unsafe { a.align_core(); }
        a
    }

    #[cfg_attr(any(target_arch = "x86", target_arch = "x86_64"), target_feature(enable = "avx2"))]
    #[cfg_attr(target_arch = "wasm32", target_feature(enable = "simd128"))]
    #[allow(non_snake_case)]
    unsafe fn align_core(&mut self) {
        let neg_inf = simd_set1_i16(i16::MIN);

        let mut best_max = i32::MIN;
        let mut best_argmax_i = 0usize;
        let mut best_argmax_j = 0usize;

        let mut dir = Direction::Diagonal;
        let mut prev_dir = Direction::Diagonal;

        let mut off = 0i32;
        let mut prev_off = 0i32;

        let mut corner1 = i16::MIN as i32;
        let mut corner2 = 0i32;

        let mut D = simd_insert_i16::<{ L - 1 }>(neg_inf, 0i16);
        let mut C = neg_inf;

        let mut D_buf = Aligned([i16::MIN; L]);
        D_buf.0[L - 1] = 0;
        let mut R_buf = Aligned([i16::MIN; L]);

        loop {
            let off_add = simd_set1_i16(clamp(prev_off - off));

            #[cfg(feature = "debug")]
            {
                println!("i: {}", self.i);
                println!("j: {}", self.j);
                println!("{:?}", dir);
                println!("off: {}", off);
            }

            let (new_D, new_C, D_max, D_argmax) = match dir {
                Direction::Diagonal => {
                    let off_add = prev_off - off;

                    self.place_block_diag(
                        clamp(corner2 - off),
                        clamp((D_buf.0[L - 2] as i32) + off_add),
                        clamp((simd_extract_i16::<{ L - 2 }>(D) as i32) + off_add),
                        clamp((R_buf.0[L - 1] as i32) + off_add),
                        clamp((simd_extract_i16::<{ L - 1 }>(C) as i32) + off_add),
                        D_buf.0.as_mut_ptr(),
                        R_buf.0.as_mut_ptr()
                    )
                },
                Direction::Right => {
                    let corner = if prev_dir == Direction::Down { clamp(corner1 - off) } else { i16::MIN };

                    self.place_block_rd::<true>(
                        simd_adds_i16(D, off_add),
                        simd_adds_i16(C, off_add),
                        corner,
                        D_buf.0.as_mut_ptr(),
                        R_buf.0.as_mut_ptr()
                    )
                },
                Direction::Down => {
                    let corner = if prev_dir == Direction::Right { clamp(corner1 - off) } else { i16::MIN };
                    let D_buf_ptr = D_buf.0.as_mut_ptr();
                    let C_buf_ptr = R_buf.0.as_mut_ptr();
                    simd_store(D_buf_ptr as _, simd_adds_i16(simd_load(D_buf_ptr as _), off_add));
                    simd_store(C_buf_ptr as _, simd_adds_i16(simd_load(C_buf_ptr as _), off_add));

                    self.place_block_rd::<false>(
                        neg_inf,
                        neg_inf,
                        corner,
                        D_buf.0.as_mut_ptr(),
                        R_buf.0.as_mut_ptr()
                    )
                }
            };
            D = new_D;
            C = new_C;

            let right_max = simd_hmax_i16(D);
            let down_max = simd_hmax_i16(simd_load(D_buf.0.as_ptr() as _));
            prev_dir = dir;

            if X_DROP {
                let max = simd_hmax_i16(D_max);

                if off + (max as i32) > best_max {
                    let lane_idx = (simd_movemask_i8(
                            simd_cmpeq_i16(D_max, simd_set1_i16(max))).trailing_zeros() / 2) as usize;
                    best_argmax_i = self.i + lane_idx;
                    best_argmax_j = self.j + simd_slow_extract_i16(D_argmax, lane_idx) as usize;
                    best_max = off + max as i32;
                }

                if off + (cmp::max(right_max, down_max) as i32) < best_max - self.x_drop {
                    // x drop termination
                    break;
                }
            }

            // first check if the shift direction is "forced"
            if self.i + L > self.query.len() && self.j + L > self.reference.len() {
                // reached the end of the strings
                break;
            } else if self.j + L > self.reference.len() {
                self.i += L;
                dir = Direction::Down;
            } else if self.i + L > self.query.len() {
                self.j += L;
                dir = Direction::Right;
            } else {
                // move according to max
                if down_max > right_max {
                    self.i += L;
                    dir = Direction::Down;
                } else if right_max > down_max {
                    self.j += L;
                    dir = Direction::Right;
                } else if right_max == down_max && down_max == D_buf.0[L - 1] {
                    self.i += L - 1;
                    self.j += L - 1;
                    dir = Direction::Diagonal;
                } else {
                    // arbitrary
                    self.j += L;
                    dir = Direction::Right;
                }
            }

            corner1 = corner2;
            corner2 = off + D_buf.0[L - 1] as i32;
            prev_off = off;
            off += simd_extract_i16::<0>(D) as i32;
        }

        self.res = if X_DROP {
            AlignResult {
                score: best_max,
                query_idx: best_argmax_i,
                reference_idx: best_argmax_j
            }
        } else {
            debug_assert!(self.i <= self.query.len());
            AlignResult {
                score: off + simd_slow_extract_i16(D, self.query.len() - self.i) as i32,
                query_idx: self.query.len(),
                reference_idx: self.reference.len()
            }
        };
    }

    // Place block diagonally, overlapping the previous block's lower right corner element.
    //
    // Assumes all inputs are already relative to the current offset.
    #[cfg_attr(any(target_arch = "x86", target_arch = "x86_64"), target_feature(enable = "avx2"))]
    #[cfg_attr(target_arch = "wasm32", target_feature(enable = "simd128"))]
    #[allow(non_snake_case)]
    #[cold]
    unsafe fn place_block_diag(&mut self,
                               corner11: i16,
                               corner10: i16,
                               corner01: i16,
                               R_corner: i16,
                               C_corner: i16,
                               D_buf: *mut i16,
                               R_buf: *mut i16) -> (Simd, Simd, Simd, Simd) {
        let (neg_inf, gap_open, gap_extend) = self.get_const_simd();
        let query = halfsimd_convert_char(halfsimd_loadu(self.query.as_ptr(self.i) as _), M::NUC);
        let mut D00 = simd_sl_i16!(neg_inf, simd_set1_i16(corner10), 2);
        let mut D10 = neg_inf;
        let mut C10 = neg_inf;
        let mut R_insert = simd_set1_i16(R_corner);
        let mut D_insert = simd_set1_i16(corner01);
        let mut D_max = neg_inf;
        let mut D_argmax = simd_set1_i16(0);
        let mut curr_i = simd_set1_i16(0);

        if TRACE {
            self.trace.dir(0b00);
        }

        for i in 0..L {
            let matrix_ptr = self.matrix.as_ptr(convert_char(self.reference.get(self.j + i), M::NUC) as usize);
            let scores1 = halfsimd_load(matrix_ptr as *const HalfSimd);
            let scores2 = if M::NUC {
                halfsimd_set1_i8(0) // unused, should be optimized out
            } else {
                halfsimd_load((matrix_ptr as *const HalfSimd).add(1))
            };

            // efficiently lookup scores for each query character
            let scores = if M::NUC {
                halfsimd_lookup1_i16(scores1, query)
            } else {
                halfsimd_lookup2_i16(scores1, scores2, query)
            };

            let mut D11 = simd_adds_i16(D00, scores);
            let mut C11 = simd_max_i16(simd_adds_i16(C10, gap_extend), simd_adds_i16(D10, gap_open));
            D11 = simd_max_i16(D11, C11);

            if i == 0 {
                D11 = simd_insert_i16::<0>(D11, corner11);
                C11 = simd_insert_i16::<0>(C11, C_corner);
            }

            let trace_D_C = if TRACE {
                simd_movemask_i8(simd_cmpeq_i16(D11, C11))
            } else {
                0 // should be optimized out
            };

            let D11_open = simd_adds_i16(D11, gap_open);
            let mut R11 = simd_sl_i16!(D11_open, R_insert, 1);
            R_insert = neg_inf;
            // avoid doing prefix scan if possible!
            if simd_movemask_i8(simd_cmpgt_i16(R11, D11_open)) != 0 {
                R11 = simd_prefix_scan_i16(R11, P::GAP_EXTEND as i16);
                D11 = simd_max_i16(D11, R11);
            }

            if TRACE {
                let trace_D_R = simd_movemask_i8(simd_cmpeq_i16(D11, R11));
                self.trace.add(((trace_D_R & Self::EVEN_BITS) << 1) | (trace_D_C & Self::EVEN_BITS));
            }

            if X_DROP {
                D_max = simd_max_i16(D_max, D11);
                let mask = simd_cmpeq_i16(D_max, D11);
                D_argmax = simd_blend_i8(D_argmax, curr_i, mask);
                curr_i = simd_adds_i16(curr_i, simd_set1_i16(1));
            }

            #[cfg(feature = "debug")]
            {
                print!("s:   ");
                simd_dbg_i16(scores);
                print!("C11: ");
                simd_dbg_i16(C11);
                print!("R11: ");
                simd_dbg_i16(R11);
                print!("D11: ");
                simd_dbg_i16(D11);
            }

            D00 = simd_sl_i16!(D11, D_insert, 1);
            D_insert = neg_inf;

            ptr::write(D_buf.add(i), simd_extract_i16::<{ L - 1 }>(D11));
            let R_buf_val = {
                let R_last = simd_max_i16(D11_open, simd_adds_i16(R11, gap_extend));
                simd_extract_i16::<{ L - 1 }>(R_last)
            };
            ptr::write(R_buf.add(i), R_buf_val);

            D10 = D11;
            C10 = C11;

            if !X_DROP && self.i + L > self.query.len()
                && self.j + i >= self.reference.len() {
                break;
            }
        }

        (D10, C10, D_max, D_argmax)
    }

    // Place block right or down.
    //
    // Assumes all inputs are already relative to the current offset.
    #[cfg_attr(any(target_arch = "x86", target_arch = "x86_64"), target_feature(enable = "avx2"))]
    #[cfg_attr(target_arch = "wasm32", target_feature(enable = "simd128"))]
    #[allow(non_snake_case)]
    #[inline]
    unsafe fn place_block_rd<const RIGHT: bool>(&mut self,
                                                mut D10: Simd,
                                                mut C10: Simd,
                                                corner: i16,
                                                D_buf: *mut i16,
                                                R_buf: *mut i16) -> (Simd, Simd, Simd, Simd) {
        let (neg_inf, gap_open, gap_extend) = self.get_const_simd();
        let query = halfsimd_convert_char(halfsimd_loadu(self.query.as_ptr(self.i) as _), M::NUC);
        let mut D00 = simd_sl_i16!(D10, simd_set1_i16(corner), 1);
        let mut D_max = neg_inf;
        let mut D_argmax = simd_set1_i16(0);
        let mut curr_i = simd_set1_i16(0);

        if TRACE {
            self.trace.dir(if RIGHT { 0b01 } else { 0b10 });
        }

        for i in 0..L {
            let matrix_ptr = self.matrix.as_ptr(convert_char(self.reference.get(self.j + i), M::NUC) as usize);
            let scores1 = halfsimd_load(matrix_ptr as *const HalfSimd);
            let scores2 = if M::NUC {
                halfsimd_set1_i8(0) // unused, should be optimized out
            } else {
                halfsimd_load((matrix_ptr as *const HalfSimd).add(1))
            };

            // efficiently lookup scores for each query character
            let scores = if M::NUC {
                halfsimd_lookup1_i16(scores1, query)
            } else {
                halfsimd_lookup2_i16(scores1, scores2, query)
            };

            let mut D11 = simd_adds_i16(D00, scores);
            let C11 = simd_max_i16(simd_adds_i16(C10, gap_extend), simd_adds_i16(D10, gap_open));
            D11 = simd_max_i16(D11, C11);

            let trace_D_C = if TRACE {
                simd_movemask_i8(simd_cmpeq_i16(D11, C11))
            } else {
                0 // should be optimized out
            };

            let D11_open = simd_adds_i16(D11, gap_open);
            let R_insert = if RIGHT { neg_inf } else { simd_set1_i16(*R_buf.add(i)) };
            let mut R11 = simd_sl_i16!(D11_open, R_insert, 1);
            // avoid doing prefix scan if possible!
            if simd_movemask_i8(simd_cmpgt_i16(R11, D11_open)) != 0 {
                R11 = simd_prefix_scan_i16(R11, P::GAP_EXTEND as i16);
                D11 = simd_max_i16(D11, R11);
            }

            if TRACE {
                let trace_D_R = simd_movemask_i8(simd_cmpeq_i16(D11, R11));
                self.trace.add(((trace_D_R & Self::EVEN_BITS) << 1) | (trace_D_C & Self::EVEN_BITS));
            }

            if X_DROP {
                D_max = simd_max_i16(D_max, D11);
                let mask = simd_cmpeq_i16(D_max, D11);
                D_argmax = simd_blend_i8(D_argmax, curr_i, mask);
                curr_i = simd_adds_i16(curr_i, simd_set1_i16(1));
            }

            #[cfg(feature = "debug")]
            {
                print!("s:   ");
                simd_dbg_i16(scores);
                print!("C11: ");
                simd_dbg_i16(C11);
                print!("R11: ");
                simd_dbg_i16(R11);
                print!("D11: ");
                simd_dbg_i16(D11);
            }

            let D_insert = if RIGHT { neg_inf } else { simd_set1_i16(*D_buf.add(i)) };
            D00 = simd_sl_i16!(D11, D_insert, 1);

            ptr::write(D_buf.add(i), simd_extract_i16::<{ L - 1 }>(D11));
            let R_buf_val = {
                let R_last = simd_max_i16(D11_open, simd_adds_i16(R11, gap_extend));
                simd_extract_i16::<{ L - 1 }>(R_last)
            };
            ptr::write(R_buf.add(i), R_buf_val);

            D10 = D11;
            C10 = C11;

            if !X_DROP && self.i + L > self.query.len()
                && self.j + i >= self.reference.len() {
                break;
            }
        }

        (D10, C10, D_max, D_argmax)
    }

    #[inline(always)]
    pub fn res(&self) -> AlignResult {
        self.res
    }

    #[inline(always)]
    pub fn trace(&self) -> &Trace {
        assert!(TRACE);
        &self.trace
    }

    #[cfg_attr(any(target_arch = "x86", target_arch = "x86_64"), target_feature(enable = "avx2"))]
    #[cfg_attr(target_arch = "wasm32", target_feature(enable = "simd128"))]
    #[inline]
    unsafe fn get_const_simd(&self) -> (Simd, Simd, Simd) {
        // some useful constant simd vectors
        let neg_inf = simd_set1_i16(i16::MIN);
        let gap_open = simd_set1_i16(P::GAP_OPEN as i16);
        let gap_extend = simd_set1_i16(P::GAP_EXTEND as i16);
        (neg_inf, gap_open, gap_extend)
    }
}

#[inline(always)]
fn convert_char(c: u8, nuc: bool) -> u8 {
    debug_assert!(c >= b'A' && c <= NULL);
    if nuc { c } else { c - b'A' }
}

#[cfg_attr(any(target_arch = "x86", target_arch = "x86_64"), target_feature(enable = "avx2"))]
#[cfg_attr(target_arch = "wasm32", target_feature(enable = "simd128"))]
#[inline]
unsafe fn halfsimd_convert_char(v: HalfSimd, nuc: bool) -> HalfSimd {
    if nuc { v } else { halfsimd_sub_i8(v, halfsimd_set1_i8(b'A' as i8)) }
}

#[inline(always)]
fn clamp(x: i32) -> i16 {
    cmp::min(cmp::max(x, i16::MIN as i32), i16::MAX as i32) as i16
}

#[inline(always)]
fn div_ceil(n: usize, d: usize) -> usize {
    (n + d - 1) / d
}

#[derive(Clone, PartialEq, Debug)]
pub struct PaddedBytes {
    s: Vec<u8>
}

impl PaddedBytes {
    #[inline(always)]
    pub fn from_bytes(b: &[u8]) -> Self {
        let mut v = b.to_owned();
        v.insert(0, NULL);
        v.resize(v.len() + L, NULL);
        Self { s: v }
    }

    #[inline(always)]
    pub fn from_str(s: &str) -> Self {
        Self::from_bytes(s.as_bytes())
    }

    #[inline(always)]
    pub fn from_string(s: String) -> Self {
        let mut v = s.into_bytes();
        v.insert(0, NULL);
        v.resize(v.len() + L, NULL);
        Self { s: v }
    }

    #[inline(always)]
    pub fn get(&self, i: usize) -> u8 {
        unsafe { *self.s.get_unchecked(i) }
    }

    #[inline(always)]
    pub fn set(&mut self, i: usize, c: u8) {
        unsafe { *self.s.get_unchecked_mut(i) = c; }
    }

    #[inline(always)]
    pub fn as_ptr(&self, i: usize) -> *const u8 {
        unsafe { self.s.as_ptr().add(i) }
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.s.len() - L - 1
    }
}

#[derive(Copy, Clone, PartialEq, Debug)]
pub struct AlignResult {
    pub score: i32,
    pub query_idx: usize,
    pub reference_idx: usize
}

#[derive(Copy, Clone, PartialEq, Debug)]
enum Direction {
    Right,
    Down,
    Diagonal
}

#[derive(Clone)]
pub struct Trace {
    trace: Vec<u32>,
    shift_dir: Vec<u32>,
    idx: usize
}

impl Trace {
    #[inline(always)]
    pub fn new(query_len: usize, reference_len: usize) -> Self {
        let len = query_len + reference_len;
        Self {
            trace: vec![0; div_ceil(len, 16)],
            shift_dir: vec![0; div_ceil(div_ceil(len, L), 16)],
            idx: 0
        }
    }

    #[inline(always)]
    pub fn add(&mut self, t: u32) {
        unsafe { *self.trace.get_unchecked_mut(self.idx) = t; }
        self.idx += 1;
    }

    #[inline(always)]
    pub fn dir(&mut self, d: u32) {
        let i = self.idx / L;
        unsafe {
            *self.shift_dir.get_unchecked_mut(i / 16) |= d << (i % 16);
        }
    }

    #[inline(always)]
    pub fn clear(&mut self) {
        self.trace.fill(0);
        self.shift_dir.fill(0);
        self.idx = 0;
    }
}

#[cfg(test)]
mod tests {
    use crate::scores::*;

    use super::*;

    #[test]
    fn test_no_x_drop() {
        type TestParams = GapParams<-11, -1>;

        let r = PaddedBytes::from_bytes(b"AAAA");
        let q = PaddedBytes::from_bytes(b"AARA");
        let a = Block::<TestParams, _, false, false>::align(&q, &r, &BLOSUM62, 0);
        assert_eq!(a.res().score, 11);

        let r = PaddedBytes::from_bytes(b"AAAA");
        let q = PaddedBytes::from_bytes(b"AAAA");
        let a = Block::<TestParams, _, false, false>::align(&q, &r, &BLOSUM62, 0);
        assert_eq!(a.res().score, 16);

        let r = PaddedBytes::from_bytes(b"AAAA");
        let q = PaddedBytes::from_bytes(b"AARA");
        let a = Block::<TestParams, _, false, false>::align(&q, &r, &BLOSUM62, 0);
        assert_eq!(a.res().score, 11);

        let r = PaddedBytes::from_bytes(b"AAAA");
        let q = PaddedBytes::from_bytes(b"RRRR");
        let a = Block::<TestParams, _, false, false>::align(&q, &r, &BLOSUM62, 0);
        assert_eq!(a.res().score, -4);

        let r = PaddedBytes::from_bytes(b"AAAA");
        let q = PaddedBytes::from_bytes(b"AAA");
        let a = Block::<TestParams, _, false, false>::align(&q, &r, &BLOSUM62, 0);
        assert_eq!(a.res().score, 1);

        type TestParams2 = GapParams<-1, -1>;

        let r = PaddedBytes::from_bytes(b"AAAN");
        let q = PaddedBytes::from_bytes(b"ATAA");
        let a = Block::<TestParams2, _, false, false>::align(&q, &r, &NW1, 0);
        assert_eq!(a.res().score, 1);

        let r = PaddedBytes::from_bytes(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        let q = PaddedBytes::from_bytes(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        let a = Block::<TestParams2, _, false, false>::align(&q, &r, &NW1, 0);
        assert_eq!(a.res().score, 32);

        let r = PaddedBytes::from_bytes(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        let q = PaddedBytes::from_bytes(b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
        let a = Block::<TestParams2, _, false, false>::align(&q, &r, &NW1, 0);
        assert_eq!(a.res().score, -32);

        let r = PaddedBytes::from_bytes(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        let q = PaddedBytes::from_bytes(b"TATATATATATATATATATATATATATATATA");
        let a = Block::<TestParams2, _, false, false>::align(&q, &r, &NW1, 0);
        assert_eq!(a.res().score, 0);

        let r = PaddedBytes::from_bytes(b"TTAAAAAAATTTTTTTTTTTT");
        let q = PaddedBytes::from_bytes(b"TTTTTTTTAAAAAAATTTTTTTTT");
        let a = Block::<TestParams2, _, false, false>::align(&q, &r, &NW1, 0);
        assert_eq!(a.res().score, 9);

        let r = PaddedBytes::from_bytes(b"AAAA");
        let q = PaddedBytes::from_bytes(b"C");
        let a = Block::<TestParams2, _, false, false>::align(&q, &r, &NW1, 0);
        assert_eq!(a.res().score, -4);
        let a = Block::<TestParams2, _, false, false>::align(&r, &q, &NW1, 0);
        assert_eq!(a.res().score, -4);
    }

    #[test]
    fn test_x_drop() {
        type TestParams = GapParams<-11, -1>;

        let r = PaddedBytes::from_bytes(b"AAARRA");
        let q = PaddedBytes::from_bytes(b"AAAAAA");
        let a = Block::<TestParams, _, false, true>::align(&q, &r, &BLOSUM62, 1);
        assert_eq!(a.res(), AlignResult { score: 14, query_idx: 6, reference_idx: 6 });

        let r = PaddedBytes::from_bytes(b"AAAAAAAAAAAAAAARRRRRRRRRRRRRRRRAAAAAAAAAAAAA");
        let q = PaddedBytes::from_bytes(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        let a = Block::<TestParams, _, false, true>::align(&q, &r, &BLOSUM62, 1);
        assert_eq!(a.res(), AlignResult { score: 60, query_idx: 15, reference_idx: 15 });
    }
}
