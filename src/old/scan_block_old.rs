#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use crate::avx2::*;

#[cfg(target_arch = "wasm32")]
use crate::simd128::*;

use crate::scores::*;

use std::{cmp, ptr, i16, alloc};
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

pub struct Block<'a, P: ScoreParams, M: 'a + Matrix, const B: usize, const TRACE: bool, const X_DROP: bool> {
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

impl<'a, P: ScoreParams, M: 'a + Matrix, const B: usize, const TRACE: bool, const X_DROP: bool> Block<'a, P, M, { B }, { TRACE }, { X_DROP }> {
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

        let mut dir = Direction::Down;

        let mut off = 0i32;
        let mut prev_off;

        let mut D10_buf = Aligned::new(B);
        let mut C10_buf = Aligned::new(B);

        let mut D_buf = Aligned::new(B);
        let mut R_buf = Aligned::new(B);
        for j in 0..B * L {
            let D_insert = if j == 0 { 0 } else { (P::GAP_OPEN as i16) + ((j - 1) as i16) * (P::GAP_EXTEND as i16) };
            D_buf.set(j, D_insert);
            R_buf.set(j, D_insert + (P::GAP_OPEN as i16));
        }
        self.i += 1;

        let mut start = 0usize;

        let mut temp_buf1 = Aligned::new(1);
        let mut temp_buf2 = Aligned::new(1);

        loop {
            prev_off = off;
            let mut D_max = neg_inf;
            let mut D_argmax = simd_set1_i16(0);
            let mut curr_right_max = neg_inf;
            let mut curr_down_max = neg_inf;

            #[cfg(feature = "debug")]
            {
                println!("i: {}", self.i);
                println!("j: {}", self.j);
                println!("start: {}", start);
                println!("{:?}", dir);
                println!("off: {}", off);
            }

            match dir {
                Direction::Right => {
                    off += D10_buf.get(0) as i32;
                    let off_add = simd_set1_i16(clamp(prev_off - off));
                    // don't care about previous corner cell, even though it might have been computed
                    // "cutting corners" in the implementation
                    let mut corner = i16::MIN;
                    let mut curr_i = simd_set1_i16(0);

                    for i in (0..B * L).step_by(L) {
                        let next_corner = clamp((D10_buf.get(i + L - 1) as i32) + prev_off - off);

                        let D10_buf_ptr = D10_buf.as_mut_ptr().add(i);
                        simd_store(D10_buf_ptr as _, simd_adds_i16(simd_load(D10_buf_ptr as _), off_add));
                        let C10_buf_ptr = C10_buf.as_mut_ptr().add(i);
                        simd_store(C10_buf_ptr as _, simd_adds_i16(simd_load(C10_buf_ptr as _), off_add));

                        let (early_exit, curr_D_max, curr_D_argmax) = self.place_block(
                            self.i + i,
                            self.j + (B - 1) * L,
                            D10_buf_ptr,
                            C10_buf_ptr,
                            corner,
                            temp_buf1.as_mut_ptr(),
                            temp_buf2.as_mut_ptr()
                        );

                        corner = next_corner;
                        curr_right_max = simd_max_i16(curr_right_max, simd_load(D10_buf_ptr as _));

                        if X_DROP {
                            D_max = simd_max_i16(D_max, curr_D_max);
                            let mask = simd_cmpeq_i16(D_max, curr_D_max);
                            D_argmax = simd_blend_i8(D_argmax, simd_adds_i16(curr_D_argmax, curr_i), mask);
                            curr_i = simd_adds_i16(curr_i, simd_set1_i16(L as i16));
                        }

                        if early_exit {
                            break;
                        }
                    }

                    // shift and offset bottom row
                    curr_down_max = self.shift_and_offset(
                        D_buf.as_mut_ptr(),
                        R_buf.as_mut_ptr(),
                        temp_buf1.as_mut_ptr(),
                        temp_buf2.as_mut_ptr(),
                        off_add
                    );
                },
                Direction::Down => {
                    off += D_buf.get(0) as i32;
                    let off_add = simd_set1_i16(clamp(prev_off - off));
                    // don't care about previous corner cell, even though it might have been computed
                    // "cutting corners" in the implementation
                    let mut corner = i16::MIN;
                    let mut curr_j = simd_set1_i16(0);

                    for j in (0..B * L).step_by(L) {
                        let next_corner = clamp((D_buf.get(j + L - 1) as i32) + prev_off - off);

                        let D_buf_ptr = D_buf.as_mut_ptr().add(j);
                        simd_store(D_buf_ptr as _, simd_adds_i16(simd_load(D_buf_ptr as _), off_add));
                        let R_buf_ptr = R_buf.as_mut_ptr().add(j);
                        simd_store(R_buf_ptr as _, simd_adds_i16(simd_load(R_buf_ptr as _), off_add));

                        let (early_exit, curr_D_max, curr_D_argmax) = self.place_block(
                            self.i + start * L,
                            self.j + j,
                            temp_buf1.as_mut_ptr(),
                            temp_buf2.as_mut_ptr(),
                            corner,
                            D_buf_ptr,
                            R_buf_ptr
                        );

                        corner = next_corner;
                        curr_down_max = simd_max_i16(curr_down_max, simd_load(D_buf_ptr as _));

                        if X_DROP {
                            D_max = simd_max_i16(D_max, curr_D_max);
                            let mask = simd_cmpeq_i16(D_max, curr_D_max);
                            D_argmax = simd_blend_i8(D_argmax, simd_adds_i16(curr_D_argmax, curr_j), mask);
                            curr_j = simd_adds_i16(curr_j, simd_set1_i16(L as i16));
                        }

                        if early_exit {
                            break;
                        }
                    }

                    // shift and offset right column
                    curr_right_max = self.shift_and_offset(
                        D10_buf.as_mut_ptr(),
                        C10_buf.as_mut_ptr(),
                        temp_buf1.as_mut_ptr(),
                        temp_buf2.as_mut_ptr(),
                        off_add
                    );
                }
            }

            let right_max = simd_hmax_i16(curr_right_max);
            let down_max = simd_hmax_i16(curr_down_max);

            if X_DROP {
                let max = simd_hmax_i16(D_max);

                if off + (max as i32) > best_max {
                    let lane_idx = (simd_movemask_i8(
                            simd_cmpeq_i16(D_max, simd_set1_i16(max))).trailing_zeros() / 2) as usize;
                    let idx = simd_slow_extract_i16(D_argmax, lane_idx) as usize;
                    match dir {
                        Direction::Right => {
                            best_argmax_i = self.i + (idx / L) * L + lane_idx;
                            best_argmax_j = self.j + (B - 1) * L + idx % L;
                        },
                        Direction::Down => {
                            best_argmax_i = self.i + start * L + lane_idx;
                            best_argmax_j = self.j + idx;
                        }
                    }
                    best_max = off + max as i32;
                }

                if off + (cmp::max(right_max, down_max) as i32) < best_max - self.x_drop {
                    // x drop termination
                    break;
                }
            }

            // first check if the shift direction is "forced"
            if self.i + (start + 1) * L > self.query.len() && self.j + B * L > self.reference.len() {
                // reached the end of the strings
                break;
            } else if start < B - 1 {
                start += 1;
            } else if self.j + B * L > self.reference.len() {
                self.i += L;
                dir = Direction::Down;
            } else if self.i + B * L > self.query.len() {
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
                } else {
                    // arbitrary
                    self.j += L;
                    dir = Direction::Right;
                }
            }
        }

        self.res = if X_DROP {
            AlignResult {
                score: best_max,
                query_idx: best_argmax_i,
                reference_idx: best_argmax_j
            }
        } else {
            debug_assert!(self.i <= self.query.len());
            debug_assert!(self.query.len() - self.i < B * L);
            AlignResult {
                score: off + D10_buf.get(self.query.len() - self.i + (B - 1 - start) * L) as i32,
                query_idx: self.query.len(),
                reference_idx: self.reference.len()
            }
        };
    }

    #[cfg_attr(any(target_arch = "x86", target_arch = "x86_64"), target_feature(enable = "avx2"))]
    #[cfg_attr(target_arch = "wasm32", target_feature(enable = "simd128"))]
    #[allow(non_snake_case)]
    #[inline]
    unsafe fn shift_and_offset(&self, buf1: *mut i16, buf2: *mut i16, temp_buf1: *mut i16, temp_buf2: *mut i16, off_add: Simd) -> Simd {
        let neg_inf = simd_set1_i16(i16::MIN);
        let mut curr_max = neg_inf;

        for i in (0..(B - 1) * L).step_by(L) {
            let next1 = simd_adds_i16(simd_load(buf1.add(i + L) as _), off_add);
            let next2 = simd_adds_i16(simd_load(buf2.add(i + L) as _), off_add);
            simd_store(buf1.add(i) as _, next1);
            simd_store(buf2.add(i) as _, next2);
            curr_max = simd_max_i16(curr_max, next1);
        }

        let next1 = simd_load(temp_buf1 as _);
        let next2 = simd_load(temp_buf2 as _);
        simd_store(buf1.add((B - 1) * L) as _, next1);
        simd_store(buf2.add((B - 1) * L) as _, next2);
        simd_store(temp_buf1 as _, neg_inf);
        simd_store(temp_buf2 as _, neg_inf);
        simd_max_i16(curr_max, next1)
    }

    // Place block right or down.
    //
    // Assumes all inputs are already relative to the current offset.
    #[cfg_attr(any(target_arch = "x86", target_arch = "x86_64"), target_feature(enable = "avx2"))]
    #[cfg_attr(target_arch = "wasm32", target_feature(enable = "simd128"))]
    #[allow(non_snake_case)]
    #[inline]
    unsafe fn place_block(&mut self,
                          start_i: usize,
                          start_j: usize,
                          D10_buf: *mut i16,
                          C10_buf: *mut i16,
                          corner: i16,
                          D_buf: *mut i16,
                          R_buf: *mut i16) -> (bool, Simd, Simd) {
        let (neg_inf, gap_open, gap_extend) = self.get_const_simd();
        let query = halfsimd_loadu(self.query.as_ptr(start_i) as _);
        let mut D10 = simd_load(D10_buf as _);
        let mut C10 = simd_load(C10_buf as _);
        let mut D00 = simd_sl_i16!(D10, simd_set1_i16(corner), 1);
        let mut D_max = neg_inf;
        let mut D_argmax = simd_set1_i16(0);
        let mut curr_i = simd_set1_i16(0);
        let mut early_exit = false;

        // TODO: trace direction

        // hottest loop in the whole program
        for j in 0..L {
            // efficiently lookup scores for each query character
            let matrix_ptr = self.matrix.as_ptr(self.reference.get(start_j + j) as usize);
            let scores1 = halfsimd_load(matrix_ptr as *const HalfSimd);
            let scores2 = if M::NUC {
                halfsimd_set1_i8(0) // unused, should be optimized out
            } else {
                halfsimd_load((matrix_ptr as *const HalfSimd).add(1))
            };
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
            let R_insert = simd_set1_i16(*R_buf.add(j));
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

            let D_insert = simd_set1_i16(*D_buf.add(j));
            D00 = simd_sl_i16!(D11, D_insert, 1);

            ptr::write(D_buf.add(j), simd_extract_i16::<{ L - 1 }>(D11));
            let R_buf_val = {
                let R_last = simd_max_i16(D11_open, simd_adds_i16(R11, gap_extend));
                simd_extract_i16::<{ L - 1 }>(R_last)
            };
            ptr::write(R_buf.add(j), R_buf_val);

            D10 = D11;
            C10 = C11;

            if !X_DROP && start_i + L > self.query.len()
                && start_j + j >= self.reference.len() {
                early_exit = true;
                break;
            }
        }

        simd_store(D10_buf as _, D10);
        simd_store(C10_buf as _, C10);
        (early_exit, D_max, D_argmax)
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
    let c = c.to_ascii_uppercase();
    debug_assert!(c >= b'A' && c <= NULL);
    if nuc { c } else { c - b'A' }
}

#[inline(always)]
fn clamp(x: i32) -> i16 {
    cmp::min(cmp::max(x, i16::MIN as i32), i16::MAX as i32) as i16
}

#[inline(always)]
fn div_ceil(n: usize, d: usize) -> usize {
    (n + d - 1) / d
}

pub struct Aligned {
    layout: alloc::Layout,
    ptr: *const i16
}

impl Aligned {
    #[cfg_attr(any(target_arch = "x86", target_arch = "x86_64"), target_feature(enable = "avx2"))]
    #[cfg_attr(target_arch = "wasm32", target_feature(enable = "simd128"))]
    #[inline]
    pub unsafe fn new(blocks: usize) -> Self {
        let layout = alloc::Layout::from_size_align_unchecked(blocks * L * 2, L_BYTES);
        let ptr = alloc::alloc(layout) as *const i16;
        let neg_inf = simd_set1_i16(i16::MIN);
        for i in (0..blocks * L).step_by(L) {
            simd_store(ptr.add(i) as _, neg_inf);
        }
        Self { layout, ptr }
    }

    #[inline(always)]
    pub fn get(&self, i: usize) -> i16 {
        unsafe { *self.ptr.add(i) }
    }

    #[inline(always)]
    pub fn set(&mut self, i: usize, v: i16) {
        unsafe { ptr::write(self.ptr.add(i) as _, v); }
    }

    #[inline(always)]
    pub fn as_mut_ptr(&mut self) -> *mut i16 {
        self.ptr as _
    }

    #[inline(always)]
    pub fn as_ptr(&self) -> *const i16 {
        self.ptr
    }
}

impl Drop for Aligned {
    fn drop(&mut self) {
        unsafe { alloc::dealloc(self.ptr as _, self.layout); }
    }
}

#[derive(Clone, PartialEq, Debug)]
pub struct PaddedBytes {
    s: Vec<u8>,
    len: usize
}

impl PaddedBytes {
    #[inline(always)]
    pub fn from_bytes(b: &[u8], blocks: usize, nuc: bool) -> Self {
        let mut v = b.to_owned();
        let len = v.len();
        v.insert(0, NULL);
        v.resize(v.len() + blocks * L, NULL);
        v.iter_mut().for_each(|c| *c = convert_char(*c, nuc));
        Self { s: v, len }
    }

    #[inline(always)]
    pub fn from_str(s: &str, blocks: usize, nuc: bool) -> Self {
        Self::from_bytes(s.as_bytes(), blocks, nuc)
    }

    #[inline(always)]
    pub fn from_string(s: String, blocks: usize, nuc: bool) -> Self {
        let mut v = s.into_bytes();
        let len = v.len();
        v.insert(0, NULL);
        v.resize(v.len() + blocks * L, NULL);
        v.iter_mut().for_each(|c| *c = convert_char(*c, nuc));
        Self { s: v, len }
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
        self.len
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
    Down
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

        let r = PaddedBytes::from_bytes(b"AAAA", 2, false);
        let q = PaddedBytes::from_bytes(b"AARA", 2, false);
        let a = Block::<TestParams, _, 2, false, false>::align(&q, &r, &BLOSUM62, 0);
        assert_eq!(a.res().score, 11);

        let r = PaddedBytes::from_bytes(b"AAAA", 2, false);
        let q = PaddedBytes::from_bytes(b"AAAA", 2, false);
        let a = Block::<TestParams, _, 2, false, false>::align(&q, &r, &BLOSUM62, 0);
        assert_eq!(a.res().score, 16);

        let r = PaddedBytes::from_bytes(b"AAAA", 2, false);
        let q = PaddedBytes::from_bytes(b"AARA", 2, false);
        let a = Block::<TestParams, _, 2, false, false>::align(&q, &r, &BLOSUM62, 0);
        assert_eq!(a.res().score, 11);

        let r = PaddedBytes::from_bytes(b"AAAA", 2, false);
        let q = PaddedBytes::from_bytes(b"RRRR", 2, false);
        let a = Block::<TestParams, _, 2, false, false>::align(&q, &r, &BLOSUM62, 0);
        assert_eq!(a.res().score, -4);

        let r = PaddedBytes::from_bytes(b"AAAA", 2, false);
        let q = PaddedBytes::from_bytes(b"AAA", 2, false);
        let a = Block::<TestParams, _, 2, false, false>::align(&q, &r, &BLOSUM62, 0);
        assert_eq!(a.res().score, 1);

        type TestParams2 = GapParams<-1, -1>;

        let r = PaddedBytes::from_bytes(b"AAAN", 2, true);
        let q = PaddedBytes::from_bytes(b"ATAA", 2, true);
        let a = Block::<TestParams2, _, 2, false, false>::align(&q, &r, &NW1, 0);
        assert_eq!(a.res().score, 1);

        let r = PaddedBytes::from_bytes(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 2, true);
        let q = PaddedBytes::from_bytes(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 2, true);
        let a = Block::<TestParams2, _, 2, false, false>::align(&q, &r, &NW1, 0);
        assert_eq!(a.res().score, 32);

        let r = PaddedBytes::from_bytes(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 2, true);
        let q = PaddedBytes::from_bytes(b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", 2, true);
        let a = Block::<TestParams2, _, 2, false, false>::align(&q, &r, &NW1, 0);
        assert_eq!(a.res().score, -32);

        let r = PaddedBytes::from_bytes(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 2, true);
        let q = PaddedBytes::from_bytes(b"TATATATATATATATATATATATATATATATA", 2, true);
        let a = Block::<TestParams2, _, 2, false, false>::align(&q, &r, &NW1, 0);
        assert_eq!(a.res().score, 0);

        let r = PaddedBytes::from_bytes(b"TTAAAAAAATTTTTTTTTTTT", 2, true);
        let q = PaddedBytes::from_bytes(b"TTTTTTTTAAAAAAATTTTTTTTT", 2, true);
        let a = Block::<TestParams2, _, 2, false, false>::align(&q, &r, &NW1, 0);
        assert_eq!(a.res().score, 9);

        let r = PaddedBytes::from_bytes(b"AAAA", 2, true);
        let q = PaddedBytes::from_bytes(b"C", 2, true);
        let a = Block::<TestParams2, _, 2, false, false>::align(&q, &r, &NW1, 0);
        assert_eq!(a.res().score, -4);
        let a = Block::<TestParams2, _, 2, false, false>::align(&r, &q, &NW1, 0);
        assert_eq!(a.res().score, -4);
    }

    #[test]
    fn test_x_drop() {
        type TestParams = GapParams<-11, -1>;

        let r = PaddedBytes::from_bytes(b"AAARRA", 2, false);
        let q = PaddedBytes::from_bytes(b"AAAAAA", 2, false);
        let a = Block::<TestParams, _, 2, false, true>::align(&q, &r, &BLOSUM62, 1);
        assert_eq!(a.res(), AlignResult { score: 14, query_idx: 6, reference_idx: 6 });

        let r = PaddedBytes::from_bytes(b"AAAAAAAAAAAAAAARRRRRRRRRRRRRRRRAAAAAAAAAAAAA", 2, false);
        let q = PaddedBytes::from_bytes(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 2, false);
        let a = Block::<TestParams, _, 2, false, true>::align(&q, &r, &BLOSUM62, 1);
        assert_eq!(a.res(), AlignResult { score: 60, query_idx: 15, reference_idx: 15 });
    }
}
