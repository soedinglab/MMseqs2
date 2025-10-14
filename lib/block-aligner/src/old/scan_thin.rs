#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use crate::avx2::*;

#[cfg(target_arch = "wasm32")]
use crate::simd128::*;

use crate::scores::*;

use std::{alloc, cmp, ptr, i16};
use std::marker::PhantomData;

const NULL: u8 = b'A' + 26u8; // this null byte value works for both amino acids and nucleotides

#[inline(always)]
fn convert_char(c: u8, nuc: bool) -> u8 {
    debug_assert!(c >= b'A' && c <= NULL);
    if nuc { c } else { c - b'A' }
}

#[inline(always)]
fn clamp(x: i32) -> i16 {
    cmp::min(cmp::max(x, i16::MIN as i32), i16::MAX as i32) as i16
}

#[derive(Copy, Clone, PartialEq, Debug)]
pub struct EndIndex {
    pub query_idx: usize,
    pub ref_idx: usize
}

#[derive(Copy, Clone, PartialEq, Debug)]
enum Direction {
    Right,
    Down(usize)
}

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
// Each band is made up of strided SIMD vectors of length 8 or 16 16-bit integers.

#[allow(non_snake_case)]
pub struct ScanAligner<'a, P: ScoreParams, M: 'a + Matrix, const K: usize, const TRACE: bool, const X_DROP: bool> {
    query_buf_layout: alloc::Layout,
    query_buf_ptr: *mut HalfSimd,
    delta_Dx0_layout: alloc::Layout,
    delta_Dx0_ptr: *mut Simd,
    delta_Cx0_layout: alloc::Layout,
    delta_Cx0_ptr: *mut Simd,
    abs_A00: i32,

    trace: Vec<u32>,

    ring_buf_idx: usize,
    ref_idx: usize,

    best_max: i32,
    best_argmax_i: usize,
    best_argmax_j: usize,

    shift_dir: Direction,

    query: &'a [u8],
    matrix: &'a M,

    _phantom: PhantomData<P>
}

impl<'a, P: ScoreParams, M: 'a + Matrix, const K: usize, const TRACE: bool, const X_DROP: bool> ScanAligner<'a, P, M, { K }, { TRACE }, { X_DROP }> {
    // round K up to multiple of L
    // add one to K to make shifting down easier
    const CEIL_K: usize = ((K + 1 + L - 1) / L) * L;
    const STRIDE: usize = Self::CEIL_K / L;

    const EVEN_BITS: u32 = 0x55555555u32;

    #[cfg_attr(any(target_arch = "x86", target_arch = "x86_64"), target_feature(enable = "avx2"))]
    #[cfg_attr(target_arch = "wasm32", target_feature(enable = "simd128"))]
    #[allow(non_snake_case)]
    pub unsafe fn new(query: &'a [u8], matrix: &'a M) -> Self {
        assert!(Self::CEIL_K <= P::I);
        assert!(P::GAP_OPEN <= P::GAP_EXTEND);
        assert!(P::I % L == 0);

        // These chunks of memory are contiguous ring buffers that represent the current band
        let query_buf_layout = alloc::Layout::from_size_align_unchecked(Self::CEIL_K * HALFSIMD_MUL, L_BYTES);
        let query_buf_ptr = alloc::alloc(query_buf_layout) as *mut u8;

        let delta_Dx0_layout = alloc::Layout::from_size_align_unchecked(Self::CEIL_K * 2, L_BYTES);
        let delta_Dx0_ptr = alloc::alloc(delta_Dx0_layout) as *mut i16;

        let delta_Cx0_layout = alloc::Layout::from_size_align_unchecked(Self::CEIL_K * 2, L_BYTES);
        let delta_Cx0_ptr = alloc::alloc(delta_Cx0_layout) as *mut i16;

        // Initialize DP columns (the first band)
        // Not extremely optimized, since it only runs once
        {
            for i in 0..Self::CEIL_K {
                let buf_idx = (i % Self::STRIDE) * L + i / Self::STRIDE;
                debug_assert!(buf_idx < Self::CEIL_K);

                if i <= query.len() {
                    ptr::write(query_buf_ptr.add(halfsimd_get_idx(buf_idx)), convert_char(if i > 0 {
                        *query.get_unchecked(i - 1) } else { NULL }, M::NUC));

                    let val = if i > 0 {
                        (P::GAP_OPEN as i32) + ((i as i32) - 1) * (P::GAP_EXTEND as i32)
                    } else {
                        0
                    };

                    ptr::write(delta_Dx0_ptr.add(buf_idx), val as i16);
                } else {
                    ptr::write(query_buf_ptr.add(halfsimd_get_idx(buf_idx)), convert_char(NULL, M::NUC));
                    ptr::write(delta_Dx0_ptr.add(buf_idx), i16::MIN);
                }

                ptr::write(delta_Cx0_ptr.add(buf_idx), i16::MIN);
            }
        }

        Self {
            query_buf_layout,
            query_buf_ptr: query_buf_ptr as *mut HalfSimd,
            delta_Dx0_layout,
            delta_Dx0_ptr: delta_Dx0_ptr as *mut Simd,
            delta_Cx0_layout,
            delta_Cx0_ptr: delta_Cx0_ptr as *mut Simd,
            abs_A00: 0i32,

            trace: vec![],

            ring_buf_idx: 0,
            ref_idx: 0,

            best_max: 0, // max of first column
            best_argmax_i: 0,
            best_argmax_j: 0,

            shift_dir: Direction::Right,

            query,
            matrix,

            _phantom: PhantomData
        }
    }

    // TODO: deal with trace when shifting down
    // TODO: count number of down/right shifts for profiling

    #[inline(always)]
    fn shift_idx(&self) -> usize {
        self.ring_buf_idx
    }

    #[inline(always)]
    fn query_idx(&self) -> usize {
        self.ring_buf_idx + Self::CEIL_K - 1
    }

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
    #[cfg_attr(any(target_arch = "x86", target_arch = "x86_64"), target_feature(enable = "avx2"))]
    #[cfg_attr(target_arch = "wasm32", target_feature(enable = "simd128"))]
    #[allow(non_snake_case)]
    pub unsafe fn align(&mut self, reference: &[u8], x_drop: i32) {
        if X_DROP {
            assert!(x_drop >= 0);
        }

        // optional 32-bit traceback
        // 0b00 = up and left, 0b10 or 0b11 = up, 0b01 = left
        if TRACE {
            self.trace.resize(self.trace.len() + (reference.len() + 1) * Self::CEIL_K / L, Self::EVEN_BITS << 1);
        }

        let gap_open = simd_set1_i16(P::GAP_OPEN as i16);
        let gap_extend = simd_set1_i16(P::GAP_EXTEND as i16);
        let neg_inf = simd_set1_i16(i16::MIN);

        let stride_gap_scalar = (Self::STRIDE as i16) * (P::GAP_EXTEND as i16);
        let stride_gap = simd_set1_i16(stride_gap_scalar);
        let stride_gap1234 = simd_set4_i16(stride_gap_scalar * 4,
                                           stride_gap_scalar * 3,
                                           stride_gap_scalar * 2,
                                           stride_gap_scalar * 1);

        // values that are "shared" between the code for shifting down and shifting right
        let mut delta_D00 = simd_sl_i16!(simd_load({
            // get last stride vector and shift it
            let idx = (self.ring_buf_idx + Self::STRIDE - 1) % Self::STRIDE;
            self.delta_Dx0_ptr.add(idx)
        }), neg_inf, 1);
        let mut abs_R_band = i32::MIN;
        let mut abs_D_band = i16::MIN as i32;
        let mut j = 0usize;

        'outer: while j < reference.len() {
            match self.shift_dir {
                Direction::Down(shift_iter) => {
                    // fixed number of shift iterations because newly calculated D values are
                    // decreasing due to gap penalties
                    for _i in 0..shift_iter {
                        // Don't go past the end of the query
                        if self.shift_idx() >= self.query.len() {
                            self.shift_dir = Direction::Right;
                            continue 'outer;
                        }

                        let shift_vec_idx = self.ring_buf_idx % Self::STRIDE;
                        debug_assert!(shift_vec_idx < Self::CEIL_K / L);

                        // Update ring buffers to slide current band down
                        // the benefit of using ring buffers is apparent here: shifting down
                        // only requires shifting one simd vector and incrementing an index
                        let shift_D_ptr = self.delta_Dx0_ptr.add(shift_vec_idx);
                        let shift_query_ptr = self.query_buf_ptr.add(shift_vec_idx);
                        let shift_C_ptr = self.delta_Cx0_ptr.add(shift_vec_idx);

                        let c = if self.query_idx() < self.query.len() {
                            *self.query.get_unchecked(self.query_idx())
                        } else {
                            NULL
                        };
                        let query_insert = halfsimd_set1_i8(convert_char(c, M::NUC) as i8);

                        // abs_R_band is only used for the first iteration
                        // it already has the gap extend cost included
                        abs_D_band = cmp::max(abs_D_band + P::GAP_OPEN as i32, abs_R_band);
                        abs_R_band = i32::MIN;

                        let delta_Dx0_insert = simd_set1_i16(clamp(abs_D_band - self.abs_A00));

                        // Now shift in new values for each band
                        halfsimd_store(shift_query_ptr, halfsimd_sr_i8!(query_insert, halfsimd_load(shift_query_ptr), 1));
                        delta_D00 = simd_load(shift_D_ptr);
                        simd_store(shift_D_ptr, simd_sr_i16!(delta_Dx0_insert, delta_D00, 1));
                        simd_store(shift_C_ptr, simd_sr_i16!(neg_inf, simd_load(shift_C_ptr), 1));

                        self.ring_buf_idx += 1;
                    }

                    self.shift_dir = Direction::Right;
                },
                Direction::Right => {
                    // Load scores for the current reference character
                    let matrix_ptr = self.matrix.as_ptr(convert_char(*reference.get_unchecked(j), M::NUC) as usize);
                    let scores1 = halfsimd_load(matrix_ptr as *const HalfSimd);
                    let scores2 = if M::NUC {
                        halfsimd_set1_i8(0) // unused, should be optimized out
                    } else {
                        halfsimd_load((matrix_ptr as *const HalfSimd).add(1))
                    };

                    // Vector for prefix scan calculations
                    let mut delta_R_max = neg_inf;
                    // add the first D value of the previous band to the absolute A value of the
                    // previous band to get the absolute A value of the current band
                    let abs_band = self.abs_A00.saturating_add({
                        let ptr = self.delta_Dx0_ptr.add(self.ring_buf_idx % Self::STRIDE);
                        simd_extract_i16::<0>(simd_load(ptr)) as i32
                    });
                    // need to offset the values from the previous band
                    let abs_offset = simd_set1_i16(clamp(self.abs_A00 - abs_band));

                    delta_D00 = simd_adds_i16(delta_D00, abs_offset);

                    // Begin initial pass
                    {
                        let mut extend_to_end = stride_gap;

                        for i in 0..Self::STRIDE {
                            let idx = (self.ring_buf_idx + i) % Self::STRIDE;
                            debug_assert!(idx < Self::CEIL_K / L);

                            // efficiently lookup scores for each query character
                            let scores = if M::NUC {
                                halfsimd_lookup1_i16(scores1, halfsimd_load(self.query_buf_ptr.add(idx)))
                            } else {
                                halfsimd_lookup2_i16(scores1, scores2, halfsimd_load(self.query_buf_ptr.add(idx)))
                            };

                            let mut delta_D11 = simd_adds_i16(delta_D00, scores);

                            let delta_D10 = simd_adds_i16(simd_load(self.delta_Dx0_ptr.add(idx)), abs_offset);
                            let delta_C10 = simd_adds_i16(simd_load(self.delta_Cx0_ptr.add(idx)), abs_offset);
                            let delta_C11 = simd_max_i16(
                                simd_adds_i16(delta_C10, gap_extend), simd_adds_i16(delta_D10, gap_open));

                            delta_D11 = simd_max_i16(delta_D11, delta_C11);

                            if TRACE {
                                let trace_idx = (Self::CEIL_K / L) * (j + 1) + i;
                                debug_assert!(trace_idx < self.trace.len());
                                *self.trace.get_unchecked_mut(trace_idx) =
                                    simd_movemask_i8(simd_cmpeq_i16(delta_C11, delta_D11));
                            }

                            extend_to_end = simd_subs_i16(extend_to_end, gap_extend);
                            delta_R_max = simd_max_i16(delta_R_max, simd_adds_i16(delta_D11, extend_to_end));

                            // Slide band right by directly overwriting the previous band
                            simd_store(self.delta_Dx0_ptr.add(idx), delta_D11);
                            simd_store(self.delta_Cx0_ptr.add(idx), delta_C11);

                            delta_D00 = delta_D10;
                        }
                    }
                    // End initial pass

                    // Begin prefix scan
                    {
                        let prev_delta_R_max_last = simd_extract_i16::<{ L - 1 }>(delta_R_max) as i32;

                        delta_R_max = simd_sl_i16!(delta_R_max, neg_inf, 1);
                        delta_R_max = simd_prefix_scan_i16(delta_R_max, stride_gap, stride_gap1234, neg_inf);

                        let curr_delta_R_max_last = simd_extract_i16::<{ L - 1 }>(simd_adds_i16(delta_R_max, stride_gap)) as i32;
                        // this is the absolute R value for the last cell of the band, plus
                        // the gap open cost
                        abs_R_band = abs_band.saturating_add(
                            cmp::max(prev_delta_R_max_last, curr_delta_R_max_last) + (P::GAP_OPEN as i32));
                    }
                    // End prefix scan

                    let mut delta_D_max = neg_inf;
                    let mut delta_D_argmax = simd_set1_i16(0);

                    // Begin final pass
                    {
                        let mut delta_R01 = simd_adds_i16(simd_subs_i16(delta_R_max, gap_extend), gap_open);
                        let mut delta_D01 = neg_inf;
                        let mut curr_i = simd_set1_i16(0);

                        for i in 0..Self::STRIDE {
                            let idx = (self.ring_buf_idx + i) % Self::STRIDE;
                            debug_assert!(idx < Self::CEIL_K / L);

                            let delta_R11 = simd_max_i16(
                                simd_adds_i16(delta_R01, gap_extend), simd_adds_i16(delta_D01, gap_open));
                            let mut delta_D11 = simd_load(self.delta_Dx0_ptr.add(idx));
                            delta_D11 = simd_max_i16(delta_D11, delta_R11);

                            if TRACE {
                                let trace_idx = (Self::CEIL_K / L) * (j + 1) + i;
                                debug_assert!(trace_idx < self.trace.len());
                                let prev_trace = *self.trace.get_unchecked(trace_idx);
                                let curr_trace = simd_movemask_i8(simd_cmpeq_i16(delta_R11, delta_D11));
                                *self.trace.get_unchecked_mut(trace_idx) =
                                    (prev_trace & Self::EVEN_BITS) | ((curr_trace & Self::EVEN_BITS) << 1);
                            }

                            // consistently update the max D value for each stride vector
                            delta_D_max = simd_max_i16(delta_D_max, delta_D11);
                            let mask = simd_cmpeq_i16(delta_D_max, delta_D11);
                            delta_D_argmax = simd_blend_i8(delta_D_argmax, curr_i, mask);
                            curr_i = simd_adds_i16(curr_i, simd_set1_i16(1));

                            simd_store(self.delta_Dx0_ptr.add(idx), delta_D11);

                            delta_D01 = delta_D11;
                            delta_R01 = delta_R11;
                        }

                        // this is the absolute D value for the last cell of the band
                        abs_D_band = abs_band.saturating_add(simd_extract_i16::<{ L - 1 }>(delta_D01) as i32);

                        // updating delta_D00 is important if the band shifts right
                        delta_D00 = simd_sl_i16!(delta_D01, neg_inf, 1);
                    }
                    // End final pass

                    let (max, lane_idx) = simd_hargmax_i16(delta_D_max);
                    let max = (max as i32).saturating_add(abs_band);
                    // "slow" because it allows an index only known at run time
                    let stride_idx = simd_slow_extract_i16(delta_D_argmax, lane_idx) as u16 as usize;
                    let argmax = stride_idx + lane_idx * Self::STRIDE;

                    self.abs_A00 = abs_band;

                    if X_DROP && max < self.best_max - x_drop {
                        break;
                    }

                    // if not x drop, then keep track of values only for the current band
                    let cond = !X_DROP || max > self.best_max;
                    self.best_argmax_i = if cond { argmax + self.shift_idx() } else { self.best_argmax_i };
                    self.best_argmax_j = if cond { j + self.ref_idx + 1 } else { self.best_argmax_j };
                    self.best_max = if cond { max } else { self.best_max };

                    // high threshold for starting to shift down, to prevent switching back and
                    // forth between down and right all time
                    self.shift_dir = if argmax > Self::CEIL_K * 5 / 8 {
                        Direction::Down(argmax - Self::CEIL_K / 2)
                    } else {
                        Direction::Right
                    };

                    j += 1;
                }
            }
        }

        self.ref_idx += reference.len();
    }

    pub fn score(&self) -> i32 {
        self.best_max
    }

    pub fn end_idx(&self) -> EndIndex {
        EndIndex {
            query_idx: self.best_argmax_i,
            ref_idx: self.best_argmax_j
        }
    }

    pub fn raw_trace(&self) -> &[u32] {
        assert!(TRACE);
        &self.trace
    }
}

impl<'a, P: ScoreParams, M: 'a + Matrix, const K: usize, const TRACE: bool, const X_DROP: bool> Drop for ScanAligner<'a, P, M, { K }, { TRACE }, { X_DROP }> {
    fn drop(&mut self) {
        unsafe {
            alloc::dealloc(self.query_buf_ptr as *mut u8, self.query_buf_layout);
            alloc::dealloc(self.delta_Dx0_ptr as *mut u8, self.delta_Dx0_layout);
            alloc::dealloc(self.delta_Cx0_ptr as *mut u8, self.delta_Cx0_layout);
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::scores::*;

    use super::*;

    #[test]
    fn test_scan_align() {
        type TestParams = Params<-11, -1, 1024>;

        unsafe {
            let r = b"AAAA";
            let q = b"AARA";
            let mut a = ScanAligner::<TestParams, _, 2, false, false>::new(q, &BLOSUM62);
            a.align(r, 0);
            assert_eq!(a.score(), 11);

            let r = b"AAAA";
            let q = b"AARA";
            let mut a = ScanAligner::<TestParams, _, 6, false, false>::new(q, &BLOSUM62);
            a.align(r, 0);
            assert_eq!(a.score(), 11);

            let r = b"AAAA";
            let q = b"AAAA";
            let mut a = ScanAligner::<TestParams, _, 2, false, false>::new(q, &BLOSUM62);
            a.align(r, 0);
            assert_eq!(a.score(), 16);

            let r = b"AAAA";
            let q = b"AARA";
            let mut a = ScanAligner::<TestParams, _, 1, false, false>::new(q, &BLOSUM62);
            a.align(r, 0);
            assert_eq!(a.score(), 11);

            let r = b"AAAA";
            let q = b"RRRR";
            let mut a = ScanAligner::<TestParams, _, 8, false, false>::new(q, &BLOSUM62);
            a.align(r, 0);
            assert_eq!(a.score(), -4);

            let r = b"AAAA";
            let q = b"AAA";
            let mut a = ScanAligner::<TestParams, _, 2, false, false>::new(q, &BLOSUM62);
            a.align(r, 0);
            assert_eq!(a.score(), 1);

            type TestParams2 = Params<-1, -1, 2048>;

            let r = b"AAAN";
            let q = b"ATAA";
            let mut a = ScanAligner::<TestParams2, _, 4, false, false>::new(q, &NW1);
            a.align(r, 0);
            assert_eq!(a.score(), 1);

            let r = b"AAAA";
            let q = b"C";
            let mut a = ScanAligner::<TestParams2, _, 8, false, false>::new(q, &NW1);
            a.align(r, 0);
            assert_eq!(a.score(), -4);
            let mut a = ScanAligner::<TestParams2, _, 8, false, false>::new(r, &NW1);
            a.align(q, 0);
            assert_eq!(a.score(), -1);
        }
    }

    #[test]
    fn test_x_drop() {
        type TestParams = Params<-11, -1, 1024>;

        unsafe {
            let r = b"AAARRA";
            let q = b"AAAAAA";
            let mut a = ScanAligner::<TestParams, _, 3, false, true>::new(q, &BLOSUM62);
            a.align(r, 1);
            assert_eq!(a.score(), 12);
            assert_eq!(a.end_idx(), EndIndex { query_idx: 3, ref_idx: 3 });

            let r = b"AAARRA";
            let q = b"AAAAAA";
            let mut a = ScanAligner::<TestParams, _, 20, false, true>::new(q, &BLOSUM62);
            a.align(r, 1);
            assert_eq!(a.score(), 12);
            assert_eq!(a.end_idx(), EndIndex { query_idx: 3, ref_idx: 3 });
        }
    }
}
