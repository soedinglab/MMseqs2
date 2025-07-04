#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use crate::avx2::*;

#[cfg(target_arch = "wasm32")]
use crate::simd128::*;

use crate::scores::*;

use std::{alloc, cmp, ptr, i16};
use std::marker::PhantomData;

const NULL: u8 = b'A' + 26u8; // this null byte value works for both amino acids and nucleotides

#[inline]
unsafe fn convert_char(c: u8, nuc: bool) -> u8 {
    debug_assert!(c >= b'A' && c <= NULL);
    if nuc { c } else { c - b'A' }
}

#[inline]
unsafe fn clamp(x: i32) -> i16 {
    cmp::min(cmp::max(x, i16::MIN as i32), i16::MAX as i32) as i16
}

#[derive(Copy, Clone, PartialEq, Debug)]
pub struct EndIndex {
    pub query_idx: usize,
    pub ref_idx: usize
}

#[derive(Copy, Clone, PartialEq, Debug)]
pub enum Direction {
    Right,
    Down
}

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
// A band consists of multiple intervals. Each interval is made up of strided vectors.
//
// TODO: update stuff to match later adaptive banding code

#[allow(non_snake_case)]
pub struct ScanAligner<'a, P: ScoreParams, M: 'a + Matrix, const K_HALF: usize, const TRACE: bool, const X_DROP: bool> {
    query_buf_layout: alloc::Layout,
    query_buf_ptr: *mut HalfSimd,
    delta_Dx0_layout: alloc::Layout,
    delta_Dx0_ptr: *mut Simd,
    delta_Cx0_layout: alloc::Layout,
    delta_Cx0_ptr: *mut Simd,
    abs_Ax0_layout: alloc::Layout,
    abs_Ax0_ptr: *mut i32,

    trace: Vec<u32>,

    query_idx: usize,
    shift_idx: isize,
    ring_buf_idx: usize,
    ref_idx: usize,

    best_max: i32,
    best_argmax_i: isize,
    best_argmax_j: usize,

    shift_dir: Direction,

    query: &'a [u8],
    matrix: &'a M,

    _phantom: PhantomData<P>
}

impl<'a, P: ScoreParams, M: 'a + Matrix, const K_HALF: usize, const TRACE: bool, const X_DROP: bool> ScanAligner<'a, P, M, { K_HALF }, { TRACE }, { X_DROP }> {
    const K: usize = K_HALF * 2 + 1;
    const CEIL_K: usize = ((Self::K + L - 1) / L) * L; // round up to multiple of L
    const NUM_INTERVALS: usize = (Self::CEIL_K + P::I - 1) / P::I;

    // Use precomputed strides so compiler can avoid division/modulo instructions
    const STRIDE_I: usize = P::I / L;
    const STRIDE_LAST: usize = (Self::CEIL_K - ((Self::CEIL_K - 1) / P::I) * P::I) / L;

    const EVEN_BITS: u32 = 0x55555555u32;

    #[cfg_attr(any(target_arch = "x86", target_arch = "x86_64"), target_feature(enable = "avx2"))]
    #[cfg_attr(target_arch = "wasm32", target_feature(enable = "simd128"))]
    #[allow(non_snake_case)]
    pub unsafe fn new(query: &'a [u8], matrix: &'a M) -> Self {
        assert!(P::GAP_OPEN <= P::GAP_EXTEND);
        assert!(P::I % L == 0);

        // These chunks of memory are contiguous ring buffers that represent every interval in the current band
        let query_buf_layout = alloc::Layout::from_size_align_unchecked(Self::CEIL_K * HALFSIMD_MUL, L_BYTES);
        let query_buf_ptr = alloc::alloc(query_buf_layout) as *mut u8;

        let delta_Dx0_layout = alloc::Layout::from_size_align_unchecked(Self::CEIL_K * 2, L_BYTES);
        let delta_Dx0_ptr = alloc::alloc(delta_Dx0_layout) as *mut i16;

        let delta_Cx0_layout = alloc::Layout::from_size_align_unchecked(Self::CEIL_K * 2, L_BYTES);
        let delta_Cx0_ptr = alloc::alloc(delta_Cx0_layout) as *mut i16;

        // 32-bit absolute values
        let abs_Ax0_layout = alloc::Layout::array::<i32>(Self::NUM_INTERVALS).unwrap();
        let abs_Ax0_ptr = alloc::alloc(abs_Ax0_layout) as *mut i32;

        // Initialize DP columns
        // Not extremely optimized, since it only runs once
        {
            let mut abs_prev = 0;

            for idx in 0..Self::CEIL_K {
                let i = (idx as isize) - (K_HALF as isize);
                let interval_idx = idx / P::I;
                let stride = cmp::min(P::I, Self::CEIL_K - interval_idx * P::I) / L;
                let buf_idx = interval_idx * P::I + (((idx % P::I) % stride) * L + (idx % P::I) / stride);
                debug_assert!(buf_idx < Self::CEIL_K);

                if i >= 0 && i <= query.len() as isize {
                    ptr::write(query_buf_ptr.add(halfsimd_get_idx(buf_idx)), convert_char(if i > 0 {
                        *query.get_unchecked(i as usize - 1) } else { NULL }, M::NUC));

                    let val = if i > 0 {
                        (P::GAP_OPEN as i32) + ((i as i32) - 1) * (P::GAP_EXTEND as i32)
                    } else {
                        0
                    };

                    if idx % P::I == 0 {
                        ptr::write(abs_Ax0_ptr.add(interval_idx), val);
                        abs_prev = val;
                    }

                    ptr::write(delta_Dx0_ptr.add(buf_idx), (val - abs_prev) as i16);
                } else {
                    if idx % P::I == 0 {
                        ptr::write(abs_Ax0_ptr.add(interval_idx), 0);
                    }

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
            abs_Ax0_layout,
            abs_Ax0_ptr,

            trace: vec![],

            query_idx: Self::CEIL_K - K_HALF - 1,
            shift_idx: -(K_HALF as isize),
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

    /// Banded alignment.
    ///
    /// Limitations:
    /// 1. Requires x86 AVX2 or WASM SIMD support.
    /// 2. The reference and the query can only contain uppercase alphabetical characters.
    /// 3. The actual size of the band is K_HALF * 2 + 1 rounded up to the next multiple of the
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

        let stride_gap_I_scalar = (Self::STRIDE_I as i16) * (P::GAP_EXTEND as i16);
        let stride_gap_I = simd_set1_i16(stride_gap_I_scalar);
        let stride_gap_last_scalar = (Self::STRIDE_LAST as i16) * (P::GAP_EXTEND as i16);
        let stride_gap_last = simd_set1_i16(stride_gap_last_scalar);
        let stride_gap1234_I = simd_set4_i16(stride_gap_I_scalar * 4,
                                             stride_gap_I_scalar * 3,
                                             stride_gap_I_scalar * 2,
                                             stride_gap_I_scalar * 1);
        let stride_gap1234_last = simd_set4_i16(stride_gap_last_scalar * 4,
                                                stride_gap_last_scalar * 3,
                                                stride_gap_last_scalar * 2,
                                                stride_gap_last_scalar * 1);
        for j in 0..reference.len() {
            // Load scores for the current reference character
            let matrix_ptr = self.matrix.as_ptr(convert_char(*reference.get_unchecked(j), M::NUC) as usize);
            let scores1 = halfsimd_load(matrix_ptr as *const HalfSimd);
            let scores2 = if M::NUC {
                halfsimd_set1_i8(0) // unused, should be optimized out
            } else {
                halfsimd_load((matrix_ptr as *const HalfSimd).add(1))
            };

            let mut band_idx = 0usize;
            let mut abs_R_interval = i16::MIN as i32;
            let mut abs_D_interval = i16::MIN as i32;
            let mut abs_D_max = i32::MIN;
            let mut abs_D_argmax = 0isize;

            while band_idx < Self::CEIL_K {
                let last_interval = (band_idx + P::I) >= Self::CEIL_K;
                let stride = if last_interval { Self::STRIDE_LAST } else { Self::STRIDE_I };
                let stride_gap = if last_interval { stride_gap_last } else { stride_gap_I };
                let mut delta_D00;
                let mut abs_interval = *self.abs_Ax0_ptr.add(band_idx / P::I);

                // Update ring buffers to slide current band down
                {
                    let idx = band_idx / L + if last_interval {
                        self.ring_buf_idx % Self::STRIDE_LAST
                    } else {
                        self.ring_buf_idx % Self::STRIDE_I
                    };
                    let delta_Dx0_idx = self.delta_Dx0_ptr.add(idx);
                    // Save first vector of the previous interval before it is replaced
                    delta_D00 = simd_load(delta_Dx0_idx);

                    if self.shift_idx + (band_idx as isize) >= 0 {
                        abs_interval = abs_interval.saturating_add(simd_extract_i16::<0>(delta_D00) as i32);
                    }

                    let query_buf_idx = self.query_buf_ptr.add(idx);
                    let delta_Cx0_idx = self.delta_Cx0_ptr.add(idx);

                    if last_interval {
                        // This must be the last interval
                        let c = if self.query_idx < self.query.len() {
                            *self.query.get_unchecked(self.query_idx)
                        } else {
                            NULL
                        };
                        let query_insert = halfsimd_set1_i8(convert_char(c, M::NUC) as i8);

                        // Now shift in new values for each interval
                        halfsimd_store(query_buf_idx, halfsimd_sr_i8!(query_insert, halfsimd_load(query_buf_idx), 1));
                        simd_store(delta_Dx0_idx, simd_sr_i16!(neg_inf, delta_D00, 1));
                        simd_store(delta_Cx0_idx, simd_sr_i16!(neg_inf, simd_load(delta_Cx0_idx), 1));
                    } else {
                        // Not the last interval; need to shift in a value from the next interval
                        let next_band_idx = band_idx + P::I;
                        let next_last_interval = (next_band_idx + P::I) >= Self::CEIL_K;
                        let next_idx = next_band_idx / L + if next_last_interval {
                            self.ring_buf_idx % Self::STRIDE_LAST
                        } else {
                            self.ring_buf_idx % Self::STRIDE_I
                        };
                        let next_abs_interval = *self.abs_Ax0_ptr.add(next_band_idx / P::I);
                        let abs_offset = simd_set1_i16(clamp(next_abs_interval - abs_interval));
                        debug_assert!(next_idx < Self::CEIL_K / L);

                        let query_insert = halfsimd_load(self.query_buf_ptr.add(next_idx));
                        let delta_Dx0_insert = simd_adds_i16(simd_load(self.delta_Dx0_ptr.add(next_idx)), abs_offset);
                        let delta_Cx0_insert = simd_adds_i16(simd_load(self.delta_Cx0_ptr.add(next_idx)), abs_offset);

                        // Now shift in new values for each interval
                        halfsimd_store(query_buf_idx, halfsimd_sr_i8!(query_insert, halfsimd_load(query_buf_idx), 1));
                        simd_store(delta_Dx0_idx, simd_sr_i16!(delta_Dx0_insert, delta_D00, 1));
                        simd_store(delta_Cx0_idx, simd_sr_i16!(delta_Cx0_insert, simd_load(delta_Cx0_idx), 1));
                    }
                }

                // Vector for prefix scan calculations
                let mut delta_R_max = neg_inf;
                let abs_offset = simd_set1_i16(clamp(*self.abs_Ax0_ptr.add(band_idx / P::I) - abs_interval));
                delta_D00 = simd_adds_i16(delta_D00, abs_offset);

                // Begin initial pass
                {
                    let mut extend_to_end = stride_gap;

                    for i in 0..stride {
                        let idx = {
                            let mut idx = self.ring_buf_idx + 1 + i;
                            idx = if last_interval { idx % Self::STRIDE_LAST } else { idx % Self::STRIDE_I };
                            band_idx / L + idx
                        };
                        debug_assert!(idx < Self::CEIL_K / L);

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
                            let trace_idx = (Self::CEIL_K / L) * (j + 1) + band_idx / L + i;
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
                    delta_R_max = simd_insert_i16::<0>(delta_R_max, clamp(abs_R_interval - abs_interval));

                    let stride_gap1234 = if last_interval { stride_gap1234_last } else { stride_gap1234_I };
                    delta_R_max = simd_prefix_scan_i16(delta_R_max, stride_gap, stride_gap1234, neg_inf);

                    let curr_delta_R_max_last = simd_extract_i16::<{ L - 1 }>(simd_adds_i16(delta_R_max, stride_gap)) as i32;
                    abs_R_interval = abs_interval.saturating_add(cmp::max(prev_delta_R_max_last, curr_delta_R_max_last));
                }
                // End prefix scan

                let mut delta_D_max = neg_inf;
                let mut delta_D_argmax = simd_set1_i16(0);

                // Begin final pass
                {
                    let mut delta_R01 = simd_adds_i16(simd_subs_i16(delta_R_max, gap_extend), gap_open);
                    let mut delta_D01 = simd_insert_i16::<0>(neg_inf, clamp(abs_D_interval - abs_interval));
                    let mut curr_i = simd_set1_i16(0);

                    for i in 0..stride {
                        let idx = {
                            let mut idx = self.ring_buf_idx + 1 + i;
                            idx = if last_interval { idx % Self::STRIDE_LAST } else { idx % Self::STRIDE_I };
                            band_idx / L + idx
                        };
                        debug_assert!(idx < Self::CEIL_K / L);

                        let delta_R11 = simd_max_i16(
                            simd_adds_i16(delta_R01, gap_extend), simd_adds_i16(delta_D01, gap_open));
                        let mut delta_D11 = simd_load(self.delta_Dx0_ptr.add(idx));
                        delta_D11 = simd_max_i16(delta_D11, delta_R11);

                        if TRACE {
                            let trace_idx = (Self::CEIL_K / L) * (j + 1) + band_idx / L + i;
                            debug_assert!(trace_idx < self.trace.len());
                            let prev_trace = *self.trace.get_unchecked(trace_idx);
                            let curr_trace = simd_movemask_i8(simd_cmpeq_i16(delta_R11, delta_D11));
                            *self.trace.get_unchecked_mut(trace_idx) =
                                (prev_trace & Self::EVEN_BITS) | ((curr_trace & Self::EVEN_BITS) << 1);
                        }

                        if X_DROP {
                            delta_D_max = simd_max_i16(delta_D_max, delta_D11);
                            let mask = simd_cmpeq_i16(delta_D_max, delta_D11);
                            delta_D_argmax = simd_blend_i8(delta_D_argmax, curr_i, mask);
                            curr_i = simd_adds_i16(curr_i, simd_set1_i16(1));
                        }

                        simd_store(self.delta_Dx0_ptr.add(idx), delta_D11);

                        delta_D01 = delta_D11;
                        delta_R01 = delta_R11;
                    }

                    abs_D_interval = abs_interval.saturating_add(simd_extract_i16::<{ L - 1 }>(delta_D01) as i32);
                }
                // End final pass

                if X_DROP {
                    let (max, lane_idx) = simd_hmax_i16(delta_D_max);
                    let max = (max as i32).saturating_add(abs_interval);
                    let stride_idx = simd_slow_extract_i16(delta_D_argmax, lane_idx) as u16 as usize;
                    let max_idx = stride_idx + lane_idx * stride + band_idx;

                    if max > abs_D_max {
                        abs_D_max = max;
                        abs_D_argmax = max_idx as isize;
                    }
                }

                debug_assert!(band_idx / P::I < Self::NUM_INTERVALS);
                *self.abs_Ax0_ptr.add(band_idx / P::I) = abs_interval;
                band_idx += P::I;
            }

            self.ring_buf_idx += 1;
            self.query_idx += 1;
            self.shift_idx += 1;

            if X_DROP {
                if abs_D_max < self.best_max - x_drop {
                    break;
                } else if abs_D_max > self.best_max {
                    self.best_max = abs_D_max;
                    self.best_argmax_i = abs_D_argmax + self.shift_idx;
                    self.best_argmax_j = j + self.ref_idx + 1;
                }

                self.shift_dir = if abs_D_argmax > K_HALF { Direction::Down } else { Direction::Right };
            }
        }

        self.ref_idx += reference.len();
    }

    #[cfg_attr(any(target_arch = "x86", target_arch = "x86_64"), target_feature(enable = "avx2"))]
    #[cfg_attr(target_arch = "wasm32", target_feature(enable = "simd128"))]
    pub unsafe fn score(&self) -> i32 {
        if X_DROP {
            self.best_max
        } else {
            // Extract the score from the last band
            assert!((self.query.len() as isize) - self.shift_idx >= 0);
            let res_i = ((self.query.len() as isize) - self.shift_idx) as usize;
            let band_idx = (res_i / P::I) * P::I;
            let stride = cmp::min(P::I, Self::CEIL_K - band_idx) / L;
            let idx = band_idx / L + (self.ring_buf_idx + (res_i % P::I)) % stride;
            debug_assert!(idx < Self::CEIL_K / L);

            let delta = simd_slow_extract_i16(simd_load(self.delta_Dx0_ptr.add(idx)), (res_i % P::I) / stride) as i32;
            let abs = *self.abs_Ax0_ptr.add(res_i / P::I);

            delta + abs
        }
    }

    pub unsafe fn end_idx(&self) -> EndIndex {
        if X_DROP {
            assert!(self.best_argmax_i >= 0);
            EndIndex {
                query_idx: self.best_argmax_i as usize,
                ref_idx: self.best_argmax_j
            }
        } else {
            EndIndex {
                query_idx: self.query.len(),
                ref_idx: self.ref_idx
            }
        }
    }

    pub fn raw_trace(&self) -> &[u32] {
        assert!(TRACE);
        &self.trace
    }
}

impl<'a, P: ScoreParams, M: 'a + Matrix, const K_HALF: usize, const TRACE: bool, const X_DROP: bool> Drop for ScanAligner<'a, P, M, { K_HALF }, { TRACE }, { X_DROP }> {
    fn drop(&mut self) {
        unsafe {
            alloc::dealloc(self.query_buf_ptr as *mut u8, self.query_buf_layout);
            alloc::dealloc(self.delta_Dx0_ptr as *mut u8, self.delta_Dx0_layout);
            alloc::dealloc(self.delta_Cx0_ptr as *mut u8, self.delta_Cx0_layout);
            alloc::dealloc(self.abs_Ax0_ptr as *mut u8, self.abs_Ax0_layout);
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
            let mut a = ScanAligner::<TestParams, _, 1, false, false>::new(q, &BLOSUM62);
            a.align(r, 0);
            assert_eq!(a.score(), 11);

            let r = b"AAAA";
            let q = b"AARA";
            let mut a = ScanAligner::<TestParams, _, 3, false, false>::new(q, &BLOSUM62);
            a.align(r, 0);
            assert_eq!(a.score(), 11);

            let r = b"AAAA";
            let q = b"AAAA";
            let mut a = ScanAligner::<TestParams, _, 1, false, false>::new(q, &BLOSUM62);
            a.align(r, 0);
            assert_eq!(a.score(), 16);

            let r = b"AAAA";
            let q = b"AARA";
            let mut a = ScanAligner::<TestParams, _, 0, false, false>::new(q, &BLOSUM62);
            a.align(r, 0);
            assert_eq!(a.score(), 11);

            let r = b"AAAA";
            let q = b"RRRR";
            let mut a = ScanAligner::<TestParams, _, 4, false, false>::new(q, &BLOSUM62);
            a.align(r, 0);
            assert_eq!(a.score(), -4);

            let r = b"AAAA";
            let q = b"AAA";
            let mut a = ScanAligner::<TestParams, _, 1, false, false>::new(q, &BLOSUM62);
            a.align(r, 0);
            assert_eq!(a.score(), 1);

            type TestParams2 = Params<-1, -1, 2048>;

            let r = b"AAAN";
            let q = b"ATAA";
            let mut a = ScanAligner::<TestParams2, _, 2, false, false>::new(q, &NW1);
            a.align(r, 0);
            assert_eq!(a.score(), 1);

            let r = b"AAAA";
            let q = b"C";
            let mut a = ScanAligner::<TestParams2, _, 4, false, false>::new(q, &NW1);
            a.align(r, 0);
            assert_eq!(a.score(), -4);
            let mut a = ScanAligner::<TestParams2, _, 4, false, false>::new(r, &NW1);
            a.align(q, 0);
            assert_eq!(a.score(), -4);
        }
    }

    #[test]
    fn test_x_drop() {
        type TestParams = Params<-11, -1, 1024>;

        unsafe {
            let r = b"AAARRA";
            let q = b"AAAAAA";
            let mut a = ScanAligner::<TestParams, _, 1, false, true>::new(q, &BLOSUM62);
            a.align(r, 1);
            assert_eq!(a.score(), 12);
            assert_eq!(a.end_idx(), EndIndex { query_idx: 3, ref_idx: 3 });

            let r = b"AAARRA";
            let q = b"AAAAAA";
            let mut a = ScanAligner::<TestParams, _, 10, false, true>::new(q, &BLOSUM62);
            a.align(r, 1);
            assert_eq!(a.score(), 12);
            assert_eq!(a.end_idx(), EndIndex { query_idx: 3, ref_idx: 3 });
        }
    }
}
