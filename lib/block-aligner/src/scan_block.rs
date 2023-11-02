//! Main block aligner algorithm and supporting data structures.

#[cfg(feature = "simd_sse2")]
use crate::sse2::*;

#[cfg(feature = "simd_avx2")]
use crate::avx2::*;

#[cfg(feature = "simd_wasm")]
use crate::simd128::*;

#[cfg(feature = "simd_neon")]
use crate::neon::*;

use crate::scores::*;
use crate::cigar::*;

use std::{cmp, ptr, i16, alloc};
use std::ops::RangeInclusive;

#[cfg(feature = "mca")]
use std::arch::asm;

// Notes:
//
// R means row, C means column (typically stands for the DP tables)
//
// BLOSUM62 matrix max = 11, min = -4; gap open = -11 (includes extension), gap extend = -1
//
// Dynamic programming formula:
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
// The term "block" gets used for two different things (unfortunately):
//
// 1. A square region of the DP matrix that shifts, grows, and shrinks.
// This is helpful for conceptually visualizing the algorithm.
//
// 2. A rectangular region representing only cells in the DP matrix that are calculated
// due to shifting or growing. Since the step size is smaller than the block size, the
// square blocks overlap. Only the non-overlapping new cells (a rectangular region) are
// computed in each step. This usage applies for the "block" in the "place_block" function
// (this is also sometimes known as the "compute rect" function).

/// Keeps track of internal state and some parameters for block aligner.
///
/// This does not describe the whole state. The allocated scratch spaces
/// and other local variables are also needed.
struct State<'a, M: Matrix> {
    query: &'a PaddedBytes,
    i: usize,
    reference: &'a PaddedBytes,
    j: usize,
    min_size: usize,
    max_size: usize,
    matrix: &'a M,
    gaps: Gaps,
    x_drop: i32
}

/// Keeps track of internal state and some parameters for block aligner for
/// sequence to profile alignment.
///
/// This does not describe the whole state. The allocated scratch spaces
/// and other local variables are also needed.
struct StateProfile<'a, P: Profile> {
    query: &'a PaddedBytes,
    i: usize,
    reference: &'a P,
    j: usize,
    min_size: usize,
    max_size: usize,
    x_drop: i32
}

/// Data structure storing the settings for block aligner.
pub struct Block<const TRACE: bool, const X_DROP: bool = false, const LOCAL_START: bool = false, const FREE_QUERY_START_GAPS: bool = false, const FREE_QUERY_END_GAPS: bool = false> {
    res: AlignResult,
    allocated: Allocated
}

macro_rules! align_core_gen {
    ($fn_name:ident, $matrix_or_profile:tt, $state:tt, $place_block_right_fn:path, $place_block_down_fn:path) => {
        #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
        #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
        #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
        #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
        #[allow(non_snake_case)]
        unsafe fn $fn_name<M: $matrix_or_profile>(&mut self, mut state: $state<M>) {
            // store the best alignment ending location for x drop alignment
            let mut best_max = 0i32;
            let mut best_argmax_i = 0usize;
            let mut best_argmax_j = 0usize;

            let mut prev_dir = Direction::Grow;
            let mut dir = Direction::Grow;
            let mut prev_size = 0;
            let mut block_size = state.min_size;

            // 32-bit score offsets
            let mut off = 0i32;
            let mut prev_off;
            let mut off_max = 0i32;

            // how many steps since the latest best score was encountered
            let mut y_drop_iter = 0;

            // how many steps where the X-drop threshold is met
            let mut x_drop_iter = 0;

            let mut i_ckpt = state.i;
            let mut j_ckpt = state.j;
            let mut off_ckpt = 0i32;

            // corner value that affects the score when shifting down then right, or right then down
            let mut D_corner = simd_set1_i16(MIN);

            loop {
                #[cfg(feature = "debug")]
                {
                    println!("i: {}", state.i);
                    println!("j: {}", state.j);
                    println!("{:?}", dir);
                    println!("block size: {}", block_size);
                }

                prev_off = off;

                // grow_D_max is an auxiliary value used when growing because it requires two separate
                // place_block steps
                let mut grow_D_max = simd_set1_i16(MIN);
                let mut grow_D_argmax_i = simd_set1_i16(0);
                let mut grow_D_argmax_j = simd_set1_i16(0);
                let (D_max, D_argmax_i, D_argmax_j, mut right_max, mut down_max) = match dir {
                    Direction::Right => {
                        off = off_max;
                        #[cfg(feature = "debug")]
                        println!("off: {}", off);
                        let off_add = simd_set1_i16(clamp(prev_off - off));

                        if TRACE {
                            self.allocated.trace.add_block(state.i, state.j + block_size - STEP, STEP, block_size, true);
                        }

                        // offset previous columns with newly computed offset
                        Self::just_offset(block_size, self.allocated.D_col.as_mut_ptr(), self.allocated.C_col.as_mut_ptr(), off_add);

                        // compute new elements in the block as a result of shifting by the step size
                        // this region should be block_size x step
                        let (D_max, D_argmax_i, D_argmax_j) = $place_block_right_fn(
                            &state,
                            state.query,
                            state.reference,
                            &mut self.allocated.trace,
                            state.i,
                            state.j + block_size - STEP,
                            STEP,
                            block_size,
                            self.allocated.D_col.as_mut_ptr(),
                            self.allocated.C_col.as_mut_ptr(),
                            self.allocated.temp_buf1.as_mut_ptr(),
                            self.allocated.temp_buf2.as_mut_ptr(),
                            if prev_dir == Direction::Down { simd_adds_i16(D_corner, off_add) } else { simd_set1_i16(MIN) },
                            clamp(-off + (ZERO as i32)),
                            true
                        );

                        // sum of a couple elements on the right border
                        let right_max = Self::prefix_max(self.allocated.D_col.as_ptr());

                        // shift and offset bottom row
                        D_corner = Self::shift_and_offset(
                            block_size,
                            self.allocated.D_row.as_mut_ptr(),
                            self.allocated.R_row.as_mut_ptr(),
                            self.allocated.temp_buf1.as_mut_ptr(),
                            self.allocated.temp_buf2.as_mut_ptr(),
                            off_add
                        );
                        // sum of a couple elements on the bottom border
                        let down_max = Self::prefix_max(self.allocated.D_row.as_ptr());

                        (D_max, D_argmax_i, D_argmax_j, right_max, down_max)
                    },
                    Direction::Down => {
                        off = off_max;
                        #[cfg(feature = "debug")]
                        println!("off: {}", off);
                        let off_add = simd_set1_i16(clamp(prev_off - off));

                        if TRACE {
                            self.allocated.trace.add_block(state.i + block_size - STEP, state.j, block_size, STEP, false);
                        }

                        // offset previous rows with newly computed offset
                        Self::just_offset(block_size, self.allocated.D_row.as_mut_ptr(), self.allocated.R_row.as_mut_ptr(), off_add);

                        // compute new elements in the block as a result of shifting by the step size
                        // this region should be step x block_size
                        let (D_max, D_argmax_i, D_argmax_j) = $place_block_down_fn(
                            &state,
                            state.reference,
                            state.query,
                            &mut self.allocated.trace,
                            state.j,
                            state.i + block_size - STEP,
                            STEP,
                            block_size,
                            self.allocated.D_row.as_mut_ptr(),
                            self.allocated.R_row.as_mut_ptr(),
                            self.allocated.temp_buf1.as_mut_ptr(),
                            self.allocated.temp_buf2.as_mut_ptr(),
                            if prev_dir == Direction::Right { simd_adds_i16(D_corner, off_add) } else { simd_set1_i16(MIN) },
                            clamp(-off + (ZERO as i32)),
                            false
                        );

                        // sum of a couple elements on the bottom border
                        let down_max = Self::prefix_max(self.allocated.D_row.as_ptr());

                        // shift and offset last column
                        D_corner = Self::shift_and_offset(
                            block_size,
                            self.allocated.D_col.as_mut_ptr(),
                            self.allocated.C_col.as_mut_ptr(),
                            self.allocated.temp_buf1.as_mut_ptr(),
                            self.allocated.temp_buf2.as_mut_ptr(),
                            off_add
                        );
                        // sum of a couple elements on the right border
                        let right_max = Self::prefix_max(self.allocated.D_col.as_ptr());

                        (D_max, D_argmax_i, D_argmax_j, right_max, down_max)
                    },
                    Direction::Grow => {
                        D_corner = simd_set1_i16(MIN);
                        let grow_step = block_size - prev_size;

                        #[cfg(feature = "debug")]
                        println!("off: {}", off);
                        #[cfg(feature = "debug")]
                        println!("Grow down");

                        if TRACE {
                            self.allocated.trace.add_block(state.i + prev_size, state.j, prev_size, grow_step, false);
                        }

                        // down
                        // this region should be prev_size x prev_size
                        let (D_max1, D_argmax_i1, D_argmax_j1) = $place_block_down_fn(
                            &state,
                            state.reference,
                            state.query,
                            &mut self.allocated.trace,
                            state.j,
                            state.i + prev_size,
                            grow_step,
                            prev_size,
                            self.allocated.D_row.as_mut_ptr(),
                            self.allocated.R_row.as_mut_ptr(),
                            self.allocated.D_col.as_mut_ptr().add(prev_size),
                            self.allocated.C_col.as_mut_ptr().add(prev_size),
                            simd_set1_i16(MIN),
                            clamp(-off + (ZERO as i32)),
                            false
                        );

                        #[cfg(feature = "debug")]
                        println!("Grow right");

                        if TRACE {
                            self.allocated.trace.add_block(state.i, state.j + prev_size, grow_step, block_size, true);
                        }

                        // right
                        // this region should be block_size x prev_size
                        let (D_max2, D_argmax_i2, D_argmax_j2) = $place_block_right_fn(
                            &state,
                            state.query,
                            state.reference,
                            &mut self.allocated.trace,
                            state.i,
                            state.j + prev_size,
                            grow_step,
                            block_size,
                            self.allocated.D_col.as_mut_ptr(),
                            self.allocated.C_col.as_mut_ptr(),
                            self.allocated.D_row.as_mut_ptr().add(prev_size),
                            self.allocated.R_row.as_mut_ptr().add(prev_size),
                            simd_set1_i16(MIN),
                            clamp(-off + (ZERO as i32)),
                            true
                        );

                        let right_max = Self::prefix_max(self.allocated.D_col.as_ptr());
                        let down_max = Self::prefix_max(self.allocated.D_row.as_ptr());
                        grow_D_max = D_max1;
                        grow_D_argmax_i = D_argmax_i1;
                        grow_D_argmax_j = D_argmax_j1;

                        // must update the checkpoint saved values just in case
                        // the block must grow again from this position
                        let mut i = 0;
                        while i < block_size {
                            self.allocated.D_col_ckpt.set_vec(&self.allocated.D_col, i);
                            self.allocated.C_col_ckpt.set_vec(&self.allocated.C_col, i);
                            self.allocated.D_row_ckpt.set_vec(&self.allocated.D_row, i);
                            self.allocated.R_row_ckpt.set_vec(&self.allocated.R_row, i);
                            i += L;
                        }

                        if TRACE {
                            self.allocated.trace.save_ckpt();
                        }

                        (D_max2, D_argmax_i2, D_argmax_j2, right_max, down_max)
                    }
                };

                prev_dir = dir;
                let D_max_max = if FREE_QUERY_END_GAPS {
                    // can assume only the right region is computed when growing,
                    // since the min block size is greater than the query length
                    simd_slow_extract_i16(D_max, state.query.len() % L)
                } else {
                    simd_hmax_i16(D_max)
                };
                let grow_max = simd_hmax_i16(grow_D_max);
                // max score of the entire block
                // note that other than off_max and best_max, the other maxs are relative to the
                // offsets off and ZERO
                let max = cmp::max(D_max_max, grow_max);
                off_max = off + (max as i32) - (ZERO as i32);
                #[cfg(feature = "debug")]
                println!("down max: {}, right max: {}", down_max, right_max);

                y_drop_iter += 1;
                // if block grows but the best score does not improve, then the block must grow again
                let mut grow_no_max = dir == Direction::Grow;

                if off_max > best_max {
                    if FREE_QUERY_END_GAPS {
                        // can assume either growing (right region only) or shifting right, so
                        // can assume state.i == 0
                        let idx_j = simd_slow_extract_i16(D_argmax_j, state.query.len() % L) as usize;
                        best_argmax_i = state.query.len();
                        match dir {
                            Direction::Right => {
                                best_argmax_j = state.j + (block_size - STEP) + idx_j;
                            },
                            Direction::Grow => {
                                best_argmax_j = state.j + prev_size + idx_j;
                            },
                            _ => unreachable!(),
                        }
                    }

                    if X_DROP {
                        // TODO: move outside loop
                        // calculate location with the best score
                        let lane_idx = simd_hargmax_i16(D_max, D_max_max);
                        let idx_i = simd_slow_extract_i16(D_argmax_i, lane_idx) as usize;
                        let idx_j = simd_slow_extract_i16(D_argmax_j, lane_idx) as usize;
                        let r = idx_i + lane_idx;
                        let c = (block_size - STEP) + idx_j;

                        match dir {
                            Direction::Right => {
                                best_argmax_i = state.i + r;
                                best_argmax_j = state.j + c;
                            },
                            Direction::Down => {
                                best_argmax_i = state.i + c;
                                best_argmax_j = state.j + r;
                            },
                            Direction::Grow => {
                                // max could be in either block
                                if D_max_max >= grow_max {
                                    // grow right
                                    best_argmax_i = state.i + idx_i + lane_idx;
                                    best_argmax_j = state.j + prev_size + idx_j;
                                } else {
                                    // grow down
                                    let lane_idx = simd_hargmax_i16(grow_D_max, grow_max);
                                    let idx_i = simd_slow_extract_i16(grow_D_argmax_i, lane_idx) as usize;
                                    let idx_j = simd_slow_extract_i16(grow_D_argmax_j, lane_idx) as usize;
                                    best_argmax_i = state.i + prev_size + idx_j;
                                    best_argmax_j = state.j + idx_i + lane_idx;
                                }
                            }
                        }
                    }

                    if block_size < state.max_size {
                        // if able to grow in the future, then save the current location
                        // as a checkpoint
                        i_ckpt = state.i;
                        j_ckpt = state.j;
                        off_ckpt = off;

                        let mut i = 0;
                        while i < block_size {
                            self.allocated.D_col_ckpt.set_vec(&self.allocated.D_col, i);
                            self.allocated.C_col_ckpt.set_vec(&self.allocated.C_col, i);
                            self.allocated.D_row_ckpt.set_vec(&self.allocated.D_row, i);
                            self.allocated.R_row_ckpt.set_vec(&self.allocated.R_row, i);
                            i += L;
                        }

                        if TRACE {
                            self.allocated.trace.save_ckpt();
                        }

                        grow_no_max = false;
                    }

                    best_max = off_max;

                    y_drop_iter = 0;
                }

                if X_DROP {
                    if off_max < best_max - state.x_drop {
                        if x_drop_iter < X_DROP_ITER - 1 {
                            x_drop_iter += 1;
                        } else {
                            // x drop termination
                            break;
                        }
                    } else {
                        x_drop_iter = 0;
                    }
                }

                if state.i + block_size > state.query.len() && state.j + block_size > state.reference.len() {
                    // reached the end of the strings
                    break;
                }

                // first check if the shift direction is "forced" to avoid going out of bounds
                if state.j + block_size > state.reference.len() {
                    state.i += STEP;
                    dir = Direction::Down;
                    continue;
                }
                if state.i + block_size > state.query.len() {
                    state.j += STEP;
                    dir = Direction::Right;
                    continue;
                }

                // three decisions are made below (based on heuristics):
                // * whether to grow
                // * whether to shrink
                // * whether to shift right or down
                // TODO: better heuristics?

                // check if it is possible to grow
                let next_size = block_size * 2;
                if next_size <= state.max_size {
                    // if approximately (block_size / step) iterations has passed since the last best
                    // max, then it is time to grow
                    if y_drop_iter > (block_size / STEP) - 1 || grow_no_max {
                        // y drop grow block
                        prev_size = block_size;
                        block_size = next_size;
                        dir = Direction::Grow;

                        // return to checkpoint
                        state.i = i_ckpt;
                        state.j = j_ckpt;
                        off = off_ckpt;

                        let mut i = 0;
                        while i < prev_size {
                            self.allocated.D_col.set_vec(&self.allocated.D_col_ckpt, i);
                            self.allocated.C_col.set_vec(&self.allocated.C_col_ckpt, i);
                            self.allocated.D_row.set_vec(&self.allocated.D_row_ckpt, i);
                            self.allocated.R_row.set_vec(&self.allocated.R_row_ckpt, i);
                            i += L;
                        }

                        if TRACE {
                            self.allocated.trace.restore_ckpt();
                        }

                        y_drop_iter = 0;
                        continue;
                    }
                }

                // check if it is possible to shrink
                if SHRINK && block_size > state.min_size && y_drop_iter == 0 {
                    let shrink_max = cmp::max(
                        Self::suffix_max(self.allocated.D_row.as_ptr(), block_size),
                        Self::suffix_max(self.allocated.D_col.as_ptr(), block_size)
                    );
                    if shrink_max >= max {
                        // just to make sure it is not right or down shift so D_corner is not used
                        prev_dir = Direction::Grow;

                        block_size /= 2;
                        let mut i = 0;
                        while i < block_size {
                            self.allocated.D_col.copy_vec(i, i + block_size);
                            self.allocated.C_col.copy_vec(i, i + block_size);
                            self.allocated.D_row.copy_vec(i, i + block_size);
                            self.allocated.R_row.copy_vec(i, i + block_size);
                            i += L;
                        }

                        state.i += block_size;
                        state.j += block_size;

                        i_ckpt = state.i;
                        j_ckpt = state.j;
                        off_ckpt = off;

                        let mut i = 0;
                        while i < block_size {
                            self.allocated.D_col_ckpt.set_vec(&self.allocated.D_col, i);
                            self.allocated.C_col_ckpt.set_vec(&self.allocated.C_col, i);
                            self.allocated.D_row_ckpt.set_vec(&self.allocated.D_row, i);
                            self.allocated.R_row_ckpt.set_vec(&self.allocated.R_row, i);
                            i += L;
                        }

                        right_max = Self::prefix_max(self.allocated.D_col.as_ptr());
                        down_max = Self::prefix_max(self.allocated.D_row.as_ptr());

                        if TRACE {
                            self.allocated.trace.save_ckpt();
                        }

                        y_drop_iter = 0;
                    }
                }

                // move according to where the max is
                if down_max > right_max {
                    state.i += STEP;
                    dir = Direction::Down;
                } else {
                    state.j += STEP;
                    dir = Direction::Right;
                }
            }

            #[cfg(any(feature = "debug", feature = "debug_size"))]
            {
                println!("query size: {}, reference size: {}", state.query.len(), state.reference.len());
                println!("end block size: {}", block_size);
            }

            self.res = if X_DROP || FREE_QUERY_END_GAPS {
                AlignResult {
                    score: best_max,
                    query_idx: best_argmax_i,
                    reference_idx: best_argmax_j
                }
            } else {
                debug_assert!(state.i <= state.query.len());
                let score = off + match dir {
                    Direction::Right | Direction::Grow => {
                        let idx = state.query.len() - state.i;
                        debug_assert!(idx < block_size);
                        (self.allocated.D_col.get(idx) as i32) - (ZERO as i32)
                    },
                    Direction::Down => {
                        let idx = state.reference.len() - state.j;
                        debug_assert!(idx < block_size);
                        (self.allocated.D_row.get(idx) as i32) - (ZERO as i32)
                    }
                };
                AlignResult {
                    score,
                    query_idx: state.query.len(),
                    reference_idx: state.reference.len()
                }
            };
        }
    };
}

/// Place block right or down for sequence-profile alignment.
///
/// Although conceptually blocks are squares, this function is actually used to compute any
/// rectangular region. For example, when shifting a block right by some step
/// size, only the rectangular region with width = step size needs to be computed, since
/// the new shifted block will partially overlap with the previous block.
///
/// Assumes all inputs are already relative to the current offset.
///
/// Inside this function, everything will be treated as shifting right,
/// conceptually. The same process can be trivially used for shifting
/// down by calling this function with different parameters.
///
/// Right and down shifts must be handled separately since a sequence
/// is aligned to a profile.
macro_rules! place_block_profile_gen {
    ($fn_name:ident, $query: ident, $query_type: ty, $reference: ident, $reference_type: ty, $q: ident, $r: ident, $right: expr) => {
        #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
        #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
        #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
        #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
        #[allow(non_snake_case)]
        unsafe fn $fn_name<P: Profile>(_state: &StateProfile<P>,
                                       $query: $query_type,
                                       $reference: $reference_type,
                                       trace: &mut Trace,
                                       start_i: usize,
                                       start_j: usize,
                                       width: usize,
                                       height: usize,
                                       D_col: *mut i16,
                                       C_col: *mut i16,
                                       D_row: *mut i16,
                                       R_row: *mut i16,
                                       mut D_corner: Simd,
                                       relative_zero: i16,
                                       _right: bool) -> (Simd, Simd, Simd) {
            let gap_extend = simd_set1_i16($r.get_gap_extend() as i16);
            let (gap_extend_all, prefix_scan_consts) = get_prefix_scan_consts(gap_extend);
            let mut D_max = simd_set1_i16(MIN);
            let mut D_argmax_i = simd_set1_i16(0);
            let mut D_argmax_j = simd_set1_i16(0);

            let mut idx = 0;
            let mut gap_open_C = simd_set1_i16(MIN);
            let mut gap_close_C = simd_set1_i16(MIN);
            let mut gap_open_R = simd_set1_i16(MIN);
            let mut gap_close_R = simd_set1_i16(MIN);

            if width == 0 || height == 0 {
                return (D_max, D_argmax_i, D_argmax_j);
            }

            // hottest loop in the whole program
            for j in 0..width {
                let mut R01 = simd_set1_i16(MIN);
                let mut D11 = simd_set1_i16(MIN);
                let mut R11 = simd_set1_i16(MIN);
                let mut prev_trace_R = simd_set1_i16(0);

                if $right {
                    idx = start_j + j;
                    gap_open_C = $r.get_gap_open_right_C(idx);
                    gap_close_C = $r.get_gap_close_right_C(idx);
                    gap_open_R = $r.get_gap_open_right_R(idx);
                }

                let mut i = 0;
                while i < height {
                    let D10 = simd_load(D_col.add(i) as _);
                    let C10 = simd_load(C_col.add(i) as _);
                    let D00 = simd_sl_i16!(D10, D_corner, 1);
                    D_corner = D10;

                    if !$right {
                        idx = start_i + i;
                        gap_open_C = $r.get_gap_open_down_R(idx);
                        gap_open_R = $r.get_gap_open_down_C(idx);
                        gap_close_R = $r.get_gap_close_down_C(idx);
                    }

                    let scores = if $right {
                        $r.get_scores_pos(idx, halfsimd_loadu($q.as_ptr(start_i + i) as _), true)
                    } else {
                        $r.get_scores_aa(idx, $q.get(start_j + j), false)
                    };
                    D11 = simd_adds_i16(D00, scores);
                    if (!LOCAL_START && start_i + i == 0 && start_j + j == 0) || (FREE_QUERY_START_GAPS && $right && start_i + i == 0) {
                        D11 = simd_insert_i16!(D11, relative_zero, 0);
                    }

                    if LOCAL_START {
                        D11 = simd_max_i16(D11, simd_set1_i16(relative_zero));
                    }

                    let C11_open = simd_adds_i16(D10, simd_adds_i16(gap_open_C, gap_extend));
                    let C11 = simd_max_i16(simd_adds_i16(C10, gap_extend), C11_open);
                    let C11_end = if $right { simd_adds_i16(C11, gap_close_C) } else { C11 };
                    D11 = simd_max_i16(D11, C11_end);
                    // at this point, C11 is fully calculated and D11 is partially calculated

                    let D11_open = simd_adds_i16(D11, gap_open_R);
                    R11 = simd_prefix_scan_i16(D11_open, gap_extend, prefix_scan_consts);
                    // do prefix scan before using R01 to break up dependency chain that depends on
                    // the last element of R01 from the previous loop iteration
                    R11 = simd_max_i16(R11, simd_adds_i16(simd_broadcasthi_i16(R01), gap_extend_all));
                    // fully calculate D11 using R11
                    let R11_end = if $right { R11 } else { simd_adds_i16(R11, gap_close_R) };
                    D11 = simd_max_i16(D11, R11_end);
                    R01 = R11;

                    #[cfg(feature = "debug")]
                    {
                        print!("s:   ");
                        simd_dbg_i16(scores);
                        print!("D00: ");
                        simd_dbg_i16(simd_subs_i16(D00, simd_set1_i16(ZERO)));
                        print!("C11: ");
                        simd_dbg_i16(simd_subs_i16(C11, simd_set1_i16(ZERO)));
                        print!("R11: ");
                        simd_dbg_i16(simd_subs_i16(R11, simd_set1_i16(ZERO)));
                        print!("D11: ");
                        simd_dbg_i16(simd_subs_i16(D11, simd_set1_i16(ZERO)));
                    }

                    if TRACE {
                        let trace_D_C = simd_cmpeq_i16(D11, C11_end);
                        let trace_D_R = simd_cmpeq_i16(D11, R11_end);
                        #[cfg(feature = "debug")]
                        {
                            print!("D_C: ");
                            simd_dbg_i16(trace_D_C);
                            print!("D_R: ");
                            simd_dbg_i16(trace_D_R);
                        }
                        // compress trace with movemask to save space
                        let mask = simd_set1_i16(0xFF00u16 as i16);
                        let trace_data = simd_movemask_i8(simd_blend_i8(trace_D_C, trace_D_R, mask));
                        let temp_trace_R = simd_cmpeq_i16(R11, D11_open);
                        let trace_R = simd_sl_i16!(temp_trace_R, prev_trace_R, 1);
                        let trace_data2 = simd_movemask_i8(simd_blend_i8(simd_cmpeq_i16(C11, C11_open), trace_R, mask));
                        prev_trace_R = temp_trace_R;

                        if LOCAL_START {
                            let zero_mask = simd_cmpeq_i16(D11, simd_set1_i16(relative_zero));
                            trace.add_zero_mask(simd_movemask_i8(zero_mask) as TraceType);
                        }

                        trace.add_trace(trace_data as TraceType, trace_data2 as TraceType);
                    }

                    D_max = simd_max_i16(D_max, D11);

                    if X_DROP || (FREE_QUERY_END_GAPS && start_i + i + L > $query.len()) {
                        // keep track of the best score and its location
                        // note: can assume right = true and only the last SIMD vectors are needed for FREE_QUERY_END_GAPS,
                        // due to the limitation that the min block size must be greater than query length
                        let mask = simd_cmpeq_i16(D_max, D11);
                        D_argmax_i = simd_blend_i8(D_argmax_i, simd_set1_i16(i as i16), mask);
                        D_argmax_j = simd_blend_i8(D_argmax_j, simd_set1_i16(j as i16), mask);
                    }

                    simd_store(D_col.add(i) as _, D11);
                    simd_store(C_col.add(i) as _, C11);
                    i += L;
                }

                D_corner = simd_set1_i16(MIN);

                ptr::write(D_row.add(j), simd_extract_i16!(D11, L - 1));
                ptr::write(R_row.add(j), simd_extract_i16!(R11, L - 1));

                if !X_DROP && !FREE_QUERY_END_GAPS && start_i + height > $query.len()
                    && start_j + j >= $reference.len() {
                    if TRACE {
                        // make sure that the trace index is updated since the rest of the loop
                        // iterations are skipped
                        trace.add_trace_idx((width - 1 - j) * (height / L));
                    }
                    break;
                }
            }

            (D_max, D_argmax_i, D_argmax_j)
        }
    };
}

// increasing step size gives a bit extra speed but results in lower accuracy
// current settings are fast, at the expense of some accuracy, and step size does not grow
const STEP: usize = 8;
const X_DROP_ITER: usize = 2; // make sure that the X-drop iteration is truly met instead of just one "bad" step
const SHRINK: bool = true; // whether to allow the block size to shrink by powers of 2
const SHRINK_SUFFIX_LEN: usize = STEP / 4;
impl<const TRACE: bool, const X_DROP: bool, const LOCAL_START: bool, const FREE_QUERY_START_GAPS: bool, const FREE_QUERY_END_GAPS: bool> Block<{ TRACE }, { X_DROP }, { LOCAL_START }, { FREE_QUERY_START_GAPS }, { FREE_QUERY_END_GAPS }> {
    /// Allocate a block aligner instance with an upper bound query length,
    /// reference length, and max block size.
    ///
    /// A block aligner instance can be reused for multiple alignments as long
    /// as the aligned sequence lengths and block sizes do not exceed the specified
    /// upper bounds.
    pub fn new(query_len: usize, reference_len: usize, max_size: usize) -> Self {
        assert!(max_size.is_power_of_two(), "Block size must be a power of two!");

        Self {
            res: AlignResult { score: 0, query_idx: 0, reference_idx: 0 },
            allocated: Allocated::new(query_len, reference_len, max_size, TRACE, LOCAL_START, FREE_QUERY_START_GAPS)
        }
    }

    /// Align two sequences with block aligner.
    ///
    /// If `TRACE` is true, then information for computing the traceback will be stored.
    /// After alignment, the traceback CIGAR string can then be computed.
    /// This will slow down alignment and use a lot more memory.
    ///
    /// If `X_DROP` is true, then the alignment process will be terminated early when
    /// the max score in the current block drops by `x_drop` below the max score encountered
    /// so far. The location of the max score is stored in the alignment result.
    /// This allows the alignment to end anywhere in the DP matrix.
    /// If `X_DROP` is false, then global alignment is done.
    ///
    /// If `LOCAL_START` is true, then the alignment is allowed to start anywhere in the DP matrix.
    /// Local alignment can be accomplished by setting `LOCAL_START` and `X_DROP` to true and `x_drop`
    /// to a very large value.
    ///
    /// If `FREE_QUERY_START_GAPS` is true, then gaps before the start of the query are free.
    ///
    /// If `FREE_QUERY_END_GAPS` is true, then gaps after the end of the query are free.
    /// Note that this has a limitation: the min block size must be greater than the length of the query.
    ///
    /// Since larger scores are better, gap and mismatches penalties must be negative.
    ///
    /// The minimum and maximum sizes of the block must be powers of 2 that are greater than the
    /// number of 16-bit lanes in a SIMD vector.
    ///
    /// The block aligner algorithm will dynamically shift a block down or right and grow its size
    /// to efficiently calculate the alignment between two strings.
    /// This is fast, but it may be slightly less accurate than computing the entire the alignment
    /// dynamic programming matrix. Growing the size of the block allows larger gaps and
    /// other potentially difficult regions to be handled correctly.
    /// The algorithm also allows shrinking the block size for greater efficiency when handling
    /// regions in the sequences with no gaps.
    /// 16-bit deltas and 32-bit offsets are used to ensure that accurate scores are
    /// computed, even when the the strings are long.
    ///
    /// When aligning sequences `q` against `r`, this algorithm computes cells in the DP matrix
    /// with `|q| + 1` rows and `|r| + 1` columns.
    ///
    /// X-drop alignment with `ByteMatrix` is not supported.
    pub fn align<M: Matrix>(&mut self, query: &PaddedBytes, reference: &PaddedBytes, matrix: &M, gaps: Gaps, size: RangeInclusive<usize>, x_drop: i32) {
        // check invariants so bad stuff doesn't happen later
        assert!(gaps.open < 0 && gaps.extend < 0, "Gap costs must be negative!");
        // there are edge cases with calculating traceback that doesn't work if
        // gap open does not cost more than gap extend
        assert!(gaps.open < gaps.extend, "Gap open must cost more than gap extend!");
        let min_size = if *size.start() < L { L } else { *size.start() };
        let max_size = if *size.end() < L { L } else { *size.end() };
        assert!(min_size < (u16::MAX as usize) && max_size < (u16::MAX as usize), "Block sizes must be smaller than 2^16 - 1!");
        assert!(min_size.is_power_of_two() && max_size.is_power_of_two(), "Block sizes must be powers of two!");
        if X_DROP {
            assert!(x_drop >= 0, "X-drop threshold amount must be nonnegative!");
        }
        assert!(!LOCAL_START || !FREE_QUERY_START_GAPS, "Cannot set both LOCAL_START and FREE_QUERY_START_GAPS!");
        assert!(!X_DROP || !FREE_QUERY_END_GAPS, "Cannot set both X_DROP and FREE_QUERY_END_GAPS!");
        assert!(!FREE_QUERY_END_GAPS || min_size > query.len(), "Min block size must be larger than the query length for FREE_QUERY_END_GAPS!");

        unsafe { self.allocated.clear(query.len(), reference.len(), max_size, TRACE); }

        let s = State {
            query,
            i: 0,
            reference,
            j: 0,
            min_size,
            max_size,
            matrix,
            gaps,
            x_drop
        };
        unsafe { self.align_core(s); }
    }

    /// Align two sequences with exponential search on the min block size.
    ///
    /// This calls `align` multiple times, doubling the min block size in each iteration
    /// until either the max block size is reached or the score reaches or exceeds the target score.
    pub fn align_exp<M: Matrix>(&mut self, query: &PaddedBytes, reference: &PaddedBytes, matrix: &M, gaps: Gaps, size: RangeInclusive<usize>, x_drop: i32, target_score: i32) -> Option<usize> {
        let mut min_size = if *size.start() < L { L } else { *size.start() };
        let max_size = if *size.end() < L { L } else { *size.end() };
        assert!(min_size < (u16::MAX as usize) && max_size < (u16::MAX as usize), "Block sizes must be smaller than 2^16 - 1!");
        assert!(min_size.is_power_of_two() && max_size.is_power_of_two(), "Block sizes must be powers of two!");

        while min_size <= max_size {
            self.align(query, reference, matrix, gaps, min_size..=max_size, x_drop);
            let curr_score = self.res().score;

            if curr_score >= target_score {
                return Some(min_size);
            }

            min_size *= 2;
        }

        None
    }

    /// Align a sequence to a profile with block aligner.
    ///
    /// If `TRACE` is true, then information for computing the traceback will be stored.
    /// After alignment, the traceback CIGAR string can then be computed.
    /// This will slow down alignment and use a lot more memory.
    ///
    /// If `X_DROP` is true, then the alignment process will be terminated early when
    /// the max score in the current block drops by `x_drop` below the max score encountered
    /// so far. The location of the max score is stored in the alignment result.
    /// This allows the alignment to end anywhere in the DP matrix.
    /// If `X_DROP` is false, then global alignment is done.
    ///
    /// If `LOCAL_START` is true, then the alignment is allowed to start anywhere in the DP matrix.
    /// Local alignment can be accomplished by setting `LOCAL_START` and `X_DROP` to true and `x_drop`
    /// to a very large value.
    ///
    /// If `FREE_QUERY_START_GAPS` is true, then gaps before the start of the query are free.
    ///
    /// If `FREE_QUERY_END_GAPS` is true, then gaps after the end of the query are free.
    /// Note that this has a limitation: the min block size must be greater than the length of the query.
    ///
    /// Since larger scores are better, gap and mismatches penalties must be negative.
    ///
    /// The minimum and maximum sizes of the block must be powers of 2 that are greater than the
    /// number of 16-bit lanes in a SIMD vector.
    ///
    /// The block aligner algorithm will dynamically shift a block down or right and grow its size
    /// to efficiently calculate the alignment between two strings.
    /// This is fast, but it may be slightly less accurate than computing the entire the alignment
    /// dynamic programming matrix. Growing the size of the block allows larger gaps and
    /// other potentially difficult regions to be handled correctly.
    /// The algorithm also allows shrinking the block size for greater efficiency when handling
    /// regions in the sequences with no gaps.
    /// 16-bit deltas and 32-bit offsets are used to ensure that accurate scores are
    /// computed, even when the the strings are long.
    ///
    /// When aligning sequence `q` against profile `p`, this algorithm computes cells in the DP matrix
    /// with `|q| + 1` rows and `|p| + 1` columns.
    pub fn align_profile<P: Profile>(&mut self, query: &PaddedBytes, profile: &P, size: RangeInclusive<usize>, x_drop: i32) {
        // check invariants so bad stuff doesn't happen later
        assert!(profile.get_gap_extend() < 0, "Gap extend cost must be negative!");
        let min_size = if *size.start() < L { L } else { *size.start() };
        let max_size = if *size.end() < L { L } else { *size.end() };
        assert!(min_size < (u16::MAX as usize) && max_size < (u16::MAX as usize), "Block sizes must be smaller than 2^16 - 1!");
        assert!(min_size.is_power_of_two() && max_size.is_power_of_two(), "Block sizes must be powers of two!");
        if X_DROP {
            assert!(x_drop >= 0, "X-drop threshold amount must be nonnegative!");
        }
        assert!(!LOCAL_START || !FREE_QUERY_START_GAPS, "Cannot set both LOCAL_START and FREE_QUERY_START_GAPS!");
        assert!(!X_DROP || !FREE_QUERY_END_GAPS, "Cannot set both X_DROP and FREE_QUERY_END_GAPS!");
        assert!(!FREE_QUERY_END_GAPS || min_size > query.len(), "Min block size must be larger than the query length for FREE_QUERY_END_GAPS!");

        unsafe { self.allocated.clear(query.len(), profile.len(), max_size, TRACE); }

        let s = StateProfile {
            query,
            i: 0,
            reference: profile,
            j: 0,
            min_size,
            max_size,
            x_drop
        };
        unsafe { self.align_profile_core(s); }
    }

    /// Align a sequence to a profile with exponential search on the min block size.
    ///
    /// This calls `align_profile` multiple times, doubling the min block size in each iteration
    /// until either the max block size is reached or the score reaches or exceeds the target score.
    pub fn align_profile_exp<P: Profile>(&mut self, query: &PaddedBytes, profile: &P, size: RangeInclusive<usize>, x_drop: i32, target_score: i32) -> Option<usize> {
        let mut min_size = if *size.start() < L { L } else { *size.start() };
        let max_size = if *size.end() < L { L } else { *size.end() };
        assert!(min_size < (u16::MAX as usize) && max_size < (u16::MAX as usize), "Block sizes must be smaller than 2^16 - 1!");
        assert!(min_size.is_power_of_two() && max_size.is_power_of_two(), "Block sizes must be powers of two!");

        while min_size <= max_size {
            self.align_profile(query, profile, min_size..=max_size, x_drop);
            let curr_score = self.res().score;

            if curr_score >= target_score {
                return Some(min_size);
            }

            min_size *= 2;
        }

        None
    }

    align_core_gen!(align_core, Matrix, State, Self::place_block, Self::place_block);
    align_core_gen!(align_profile_core, Profile, StateProfile, Self::place_block_profile_right, Self::place_block_profile_down);

    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    #[allow(non_snake_case)]
    #[inline]
    unsafe fn just_offset(block_size: usize, buf1: *mut i16, buf2: *mut i16, off_add: Simd) {
        let mut i = 0;
        while i < block_size {
            let a = simd_adds_i16(simd_load(buf1.add(i) as _), off_add);
            let b = simd_adds_i16(simd_load(buf2.add(i) as _), off_add);
            simd_store(buf1.add(i) as _, a);
            simd_store(buf2.add(i) as _, b);
            i += L;
        }
    }

    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    #[allow(non_snake_case)]
    #[inline]
    unsafe fn prefix_max(buf: *const i16) -> i16 {
        simd_prefix_hmax_i16!(simd_load(buf as _), STEP)
    }

    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    #[allow(non_snake_case)]
    #[inline]
    unsafe fn suffix_max(buf: *const i16, buf_len: usize) -> i16 {
        simd_suffix_hmax_i16!(simd_load(buf.add(buf_len - L) as _), SHRINK_SUFFIX_LEN)
    }

    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    #[allow(non_snake_case)]
    #[inline]
    unsafe fn shift_and_offset(block_size: usize, buf1: *mut i16, buf2: *mut i16, temp_buf1: *mut i16, temp_buf2: *mut i16, off_add: Simd) -> Simd {
        let mut curr1 = simd_adds_i16(simd_load(buf1 as _), off_add);
        let D_corner = simd_set1_i16(simd_extract_i16!(curr1, STEP - 1));
        let mut curr2 = simd_adds_i16(simd_load(buf2 as _), off_add);

        let mut i = 0;
        while i < block_size - L {
            let next1 = simd_adds_i16(simd_load(buf1.add(i + L) as _), off_add);
            let next2 = simd_adds_i16(simd_load(buf2.add(i + L) as _), off_add);
            simd_store(buf1.add(i) as _, simd_step(next1, curr1));
            simd_store(buf2.add(i) as _, simd_step(next2, curr2));
            curr1 = next1;
            curr2 = next2;
            i += L;
        }

        let next1 = simd_load(temp_buf1 as _);
        let next2 = simd_load(temp_buf2 as _);
        simd_store(buf1.add(block_size - L) as _, simd_step(next1, curr1));
        simd_store(buf2.add(block_size - L) as _, simd_step(next2, curr2));
        D_corner
    }

    /// Place block right or down for sequence-sequence alignment.
    ///
    /// Although conceptually blocks are squares, this function is actually used to compute any
    /// rectangular region. For example, when shifting a block right by some step
    /// size, only the rectangular region with width = step size needs to be computed, since
    /// the new shifted block will partially overlap with the previous block.
    ///
    /// Assumes all inputs are already relative to the current offset.
    ///
    /// Inside this function, everything will be treated as shifting right,
    /// conceptually. The same process can be trivially used for shifting
    /// down by calling this function with different parameters.
    ///
    /// The same function can be reused for right and down shifts because
    /// sequence to sequence alignment is symmetric.
    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    #[allow(non_snake_case)]
    unsafe fn place_block<M: Matrix>(state: &State<M>,
                                     query: &PaddedBytes,
                                     reference: &PaddedBytes,
                                     trace: &mut Trace,
                                     start_i: usize,
                                     start_j: usize,
                                     width: usize,
                                     height: usize,
                                     D_col: *mut i16,
                                     C_col: *mut i16,
                                     D_row: *mut i16,
                                     R_row: *mut i16,
                                     mut D_corner: Simd,
                                     relative_zero: i16,
                                     right: bool) -> (Simd, Simd, Simd) {
        let gap_open = simd_set1_i16(state.gaps.open as i16);
        let gap_extend = simd_set1_i16(state.gaps.extend as i16);
        let (gap_extend_all, prefix_scan_consts) = get_prefix_scan_consts(gap_extend);
        let mut D_max = simd_set1_i16(MIN);
        let mut D_argmax_i = simd_set1_i16(0);
        let mut D_argmax_j = simd_set1_i16(0);

        if width == 0 || height == 0 {
            return (D_max, D_argmax_i, D_argmax_j);
        }

        // hottest loop in the whole program
        for j in 0..width {
            let mut R01 = simd_set1_i16(MIN);
            let mut D11 = simd_set1_i16(MIN);
            let mut R11 = simd_set1_i16(MIN);
            let mut prev_trace_R = simd_set1_i16(0);

            let c = reference.get(start_j + j);

            let mut i = 0;
            while i < height {
                #[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), feature = "mca"))]
                asm!("# LLVM-MCA-BEGIN place_block inner loop", options(nomem, nostack, preserves_flags));

                let D10 = simd_load(D_col.add(i) as _);
                let C10 = simd_load(C_col.add(i) as _);
                let D00 = simd_sl_i16!(D10, D_corner, 1);
                D_corner = D10;

                let scores = state.matrix.get_scores(c, halfsimd_loadu(query.as_ptr(start_i + i) as _), right);
                D11 = simd_adds_i16(D00, scores);
                if (!LOCAL_START && start_i + i == 0 && start_j + j == 0) || (FREE_QUERY_START_GAPS && right && start_i + i == 0) {
                    D11 = simd_insert_i16!(D11, relative_zero, 0);
                }

                if LOCAL_START {
                    D11 = simd_max_i16(D11, simd_set1_i16(relative_zero));
                }

                let C11_open = simd_adds_i16(D10, gap_open);
                let C11 = simd_max_i16(simd_adds_i16(C10, gap_extend), C11_open);
                D11 = simd_max_i16(D11, C11);
                // at this point, C11 is fully calculated and D11 is partially calculated

                let D11_open = simd_adds_i16(D11, simd_subs_i16(gap_open, gap_extend));
                R11 = simd_prefix_scan_i16(D11_open, gap_extend, prefix_scan_consts);
                // do prefix scan before using R01 to break up dependency chain that depends on
                // the last element of R01 from the previous loop iteration
                R11 = simd_max_i16(R11, simd_adds_i16(simd_broadcasthi_i16(R01), gap_extend_all));
                // fully calculate D11 using R11
                D11 = simd_max_i16(D11, R11);
                R01 = R11;

                #[cfg(feature = "debug")]
                {
                    print!("s:   ");
                    simd_dbg_i16(scores);
                    print!("D00: ");
                    simd_dbg_i16(simd_subs_i16(D00, simd_set1_i16(ZERO)));
                    print!("C11: ");
                    simd_dbg_i16(simd_subs_i16(C11, simd_set1_i16(ZERO)));
                    print!("R11: ");
                    simd_dbg_i16(simd_subs_i16(R11, simd_set1_i16(ZERO)));
                    print!("D11: ");
                    simd_dbg_i16(simd_subs_i16(D11, simd_set1_i16(ZERO)));
                }

                if TRACE {
                    let trace_D_C = simd_cmpeq_i16(D11, C11);
                    let trace_D_R = simd_cmpeq_i16(D11, R11);
                    #[cfg(feature = "debug")]
                    {
                        print!("D_C: ");
                        simd_dbg_i16(trace_D_C);
                        print!("D_R: ");
                        simd_dbg_i16(trace_D_R);
                    }
                    // compress trace with movemask to save space
                    let mask = simd_set1_i16(0xFF00u16 as i16);
                    let trace_data = simd_movemask_i8(simd_blend_i8(trace_D_C, trace_D_R, mask));
                    let temp_trace_R = simd_cmpeq_i16(R11, D11_open);
                    let trace_R = simd_sl_i16!(temp_trace_R, prev_trace_R, 1);
                    let trace_data2 = simd_movemask_i8(simd_blend_i8(simd_cmpeq_i16(C11, C11_open), trace_R, mask));
                    prev_trace_R = temp_trace_R;

                    if LOCAL_START {
                        let zero_mask = simd_cmpeq_i16(D11, simd_set1_i16(relative_zero));
                        trace.add_zero_mask(simd_movemask_i8(zero_mask) as TraceType);
                    }

                    trace.add_trace(trace_data as TraceType, trace_data2 as TraceType);
                }

                D_max = simd_max_i16(D_max, D11);

                if X_DROP || (FREE_QUERY_END_GAPS && start_i + i + L > query.len()) {
                    // keep track of the best score and its location
                    // note: can assume right = true and only the last SIMD vectors are needed for FREE_QUERY_END_GAPS,
                    // due to the limitation that the min block size must be greater than query length
                    let mask = simd_cmpeq_i16(D_max, D11);
                    D_argmax_i = simd_blend_i8(D_argmax_i, simd_set1_i16(i as i16), mask);
                    D_argmax_j = simd_blend_i8(D_argmax_j, simd_set1_i16(j as i16), mask);
                }

                simd_store(D_col.add(i) as _, D11);
                simd_store(C_col.add(i) as _, C11);
                i += L;

                #[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), feature = "mca"))]
                asm!("# LLVM-MCA-END", options(nomem, nostack, preserves_flags));
            }

            D_corner = simd_set1_i16(MIN);

            ptr::write(D_row.add(j), simd_extract_i16!(D11, L - 1));
            ptr::write(R_row.add(j), simd_extract_i16!(R11, L - 1));

            if !X_DROP && !FREE_QUERY_END_GAPS && start_i + height > query.len()
                && start_j + j >= reference.len() {
                if TRACE {
                    // make sure that the trace index is updated since the rest of the loop
                    // iterations are skipped
                    trace.add_trace_idx((width - 1 - j) * (height / L));
                }
                break;
            }
        }

        (D_max, D_argmax_i, D_argmax_j)
    }

    place_block_profile_gen!(place_block_profile_right, query, &PaddedBytes, reference, &P, query, reference, true);
    place_block_profile_gen!(place_block_profile_down, reference, &P, query, &PaddedBytes, query, reference, false);

    /// Get the resulting score and ending location of the alignment.
    #[inline]
    pub fn res(&self) -> AlignResult {
        self.res
    }

    /// Get the trace of the alignment, assuming `TRACE` is true.
    #[inline]
    pub fn trace(&self) -> &Trace {
        assert!(TRACE);
        &self.allocated.trace
    }
}

/// Allocated scratch spaces for alignment.
///
/// Scratch spaces can be reused for aligning strings with shorter lengths
/// and smaller block sizes.
#[allow(non_snake_case)]
struct Allocated {
    pub trace: Trace,

    // bottom and right borders of the current block
    pub D_col: Aligned,
    pub C_col: Aligned,
    pub D_row: Aligned,
    pub R_row: Aligned,

    // the state at the previous checkpoint (where latest best score was encountered)
    pub D_col_ckpt: Aligned,
    pub C_col_ckpt: Aligned,
    pub D_row_ckpt: Aligned,
    pub R_row_ckpt: Aligned,

    // reused buffers for storing values that must be shifted
    // into the other border when the block moves in one direction
    pub temp_buf1: Aligned,
    pub temp_buf2: Aligned,

    query_len: usize,
    reference_len: usize,
    max_size: usize,
    trace_flag: bool
}

impl Allocated {
    #[allow(non_snake_case)]
    fn new(query_len: usize, reference_len: usize, max_size: usize, trace_flag: bool, local_start: bool, free_query_start_gaps: bool) -> Self {
        unsafe {
            let trace = if trace_flag {
                Trace::new(query_len, reference_len, max_size, local_start, free_query_start_gaps)
            } else {
                Trace::new(0, 0, 0, false, false)
            };
            let D_col = Aligned::new(max_size);
            let C_col = Aligned::new(max_size);
            let D_row = Aligned::new(max_size);
            let R_row = Aligned::new(max_size);
            let D_col_ckpt = Aligned::new(max_size);
            let C_col_ckpt = Aligned::new(max_size);
            let D_row_ckpt = Aligned::new(max_size);
            let R_row_ckpt = Aligned::new(max_size);
            let temp_buf1 = Aligned::new(L);
            let temp_buf2 = Aligned::new(L);

            Self {
                trace,
                D_col,
                C_col,
                D_row,
                R_row,
                D_col_ckpt,
                C_col_ckpt,
                D_row_ckpt,
                R_row_ckpt,
                temp_buf1,
                temp_buf2,
                query_len,
                reference_len,
                max_size,
                trace_flag
            }
        }
    }

    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    unsafe fn clear(&mut self, query_len: usize, reference_len: usize, max_size: usize, trace_flag: bool) {
        // do not overwrite query_len, reference_len, etc. because they are upper bounds
        assert!(query_len + reference_len <= self.query_len + self.reference_len);
        assert!(max_size <= self.max_size);
        assert_eq!(trace_flag, self.trace_flag);

        self.trace.clear(query_len, reference_len);
        self.D_col.clear(max_size);
        self.C_col.clear(max_size);
        self.D_row.clear(max_size);
        self.R_row.clear(max_size);
        self.D_col_ckpt.clear(max_size);
        self.C_col_ckpt.clear(max_size);
        self.D_row_ckpt.clear(max_size);
        self.R_row_ckpt.clear(max_size);
        self.temp_buf1.clear(L);
        self.temp_buf2.clear(L);
    }
}

/// Holds the trace generated by block aligner.
#[derive(Clone)]
pub struct Trace {
    trace: Vec<TraceType>,
    trace2: Vec<TraceType>,
    right: Vec<u64>,
    block_start: Vec<u32>,
    block_size: Vec<u16>,
    zero_mask: Vec<TraceType>,
    trace_idx: usize,
    block_idx: usize,
    ckpt_trace_idx: usize,
    ckpt_block_idx: usize,
    query_len: usize,
    reference_len: usize,
    local_start: bool,
    free_query_start_gaps: bool
}

impl Trace {
    #[inline]
    fn new(query_len: usize, reference_len: usize, max_size: usize, local_start: bool, free_query_start_gaps: bool) -> Self {
        let len = query_len + reference_len + 2;
        let trace = vec![0 as TraceType; (max_size / L) * (len + max_size * 2)];
        let trace2 = vec![0 as TraceType; (max_size / L) * (len + max_size * 2)];
        let right = vec![0u64; div_ceil(len, 64)];
        let block_start = vec![0u32; len * 2];
        let block_size = vec![0u16; len * 2];
        let zero_mask = if local_start {
            vec![0 as TraceType; (max_size / L) * (len + max_size * 2)]
        } else {
            vec![]
        };

        Self {
            trace,
            trace2,
            right,
            block_start,
            block_size,
            zero_mask,
            trace_idx: 0,
            block_idx: 0,
            ckpt_trace_idx: 0,
            ckpt_block_idx: 0,
            query_len,
            reference_len,
            local_start,
            free_query_start_gaps,
        }
    }

    #[inline]
    fn clear(&mut self, query_len: usize, reference_len: usize) {
        // no need to clear trace, block_start, and block_size
        self.right.fill(0);
        self.trace_idx = 0;
        self.block_idx = 0;
        self.ckpt_trace_idx = 0;
        self.ckpt_block_idx = 0;
        self.query_len = query_len;
        self.reference_len = reference_len;
    }

    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    #[inline]
    unsafe fn add_trace(&mut self, t: TraceType, t2: TraceType) {
        debug_assert!(self.trace_idx < self.trace.len());
        store_trace(self.trace.as_mut_ptr().add(self.trace_idx), t);
        store_trace(self.trace2.as_mut_ptr().add(self.trace_idx), t2);
        self.trace_idx += 1;
    }

    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    #[inline]
    unsafe fn add_zero_mask(&mut self, mask: TraceType) {
        store_trace(self.zero_mask.as_mut_ptr().add(self.trace_idx), mask);
    }

    #[inline]
    fn add_block(&mut self, i: usize, j: usize, width: usize, height: usize, right: bool) {
        debug_assert!(self.block_idx * 2 < self.block_start.len());
        unsafe {
            *self.block_start.as_mut_ptr().add(self.block_idx * 2) = i as u32;
            *self.block_start.as_mut_ptr().add(self.block_idx * 2 + 1) = j as u32;
            *self.block_size.as_mut_ptr().add(self.block_idx * 2) = height as u16;
            *self.block_size.as_mut_ptr().add(self.block_idx * 2 + 1) = width as u16;

            let a = self.block_idx / 64;
            let b = self.block_idx % 64;
            let v = *self.right.as_ptr().add(a) & !(1 << b); // clear bit
            *self.right.as_mut_ptr().add(a) = v | ((right as u64) << b);

            self.block_idx += 1;
        }
    }

    #[inline]
    fn add_trace_idx(&mut self, add: usize) {
        self.trace_idx += add;
    }

    #[inline]
    fn save_ckpt(&mut self) {
        self.ckpt_trace_idx = self.trace_idx;
        self.ckpt_block_idx = self.block_idx;
    }

    /// The trace data structure is like a stack, so all trace values and blocks after the
    /// checkpoint is essentially popped off the stack.
    #[inline]
    fn restore_ckpt(&mut self) {
        self.trace_idx = self.ckpt_trace_idx;
        self.block_idx = self.ckpt_block_idx;
    }

    /// Create a CIGAR string that represents a single traceback path ending on the specified
    /// location.
    ///
    /// When aligning `q` against `r`, this represents the edits to go from `r` to `q`.
    /// Matches and mismatches are both represented with `M`.
    pub fn cigar(&self, i: usize, j: usize, cigar: &mut Cigar) {
        self.cigar_core::<false>(i, j, None, None, cigar);
    }

    /// Create a CIGAR string that represents a single traceback path ending on the specified
    /// location.
    ///
    /// When aligning `q` against `r`, this represents the edits to go from `r` to `q`.
    /// Matches are represented using `=` and mismatches are represented using `X`.
    pub fn cigar_eq(&self, query: &PaddedBytes, reference: &PaddedBytes, i: usize, j: usize, cigar: &mut Cigar) {
        self.cigar_core::<true>(i, j, Some(query), Some(reference), cigar);
    }

    fn cigar_core<const EQ: bool>(&self, mut i: usize, mut j: usize, q: Option<&PaddedBytes>, r: Option<&PaddedBytes>, cigar: &mut Cigar) {
        assert!(i <= self.query_len && j <= self.reference_len, "Traceback cigar end position must be in bounds!");
        if EQ {
            assert!(q.is_some() && r.is_some());
        }

        cigar.clear(i, j);

        unsafe {
            let mut block_idx = self.block_idx;
            let mut trace_idx = self.trace_idx;
            let mut block_i;
            let mut block_j;
            let mut block_width;
            let mut block_height;
            let mut right;

            #[derive(Copy, Clone, PartialEq, Debug)]
            enum Table {
                D = 0b00,
                C = 0b01,
                R = 0b10
            }

            // use lookup table instead of hard to predict branches
            // constructed at compile time!
            static OP_LUT: [[(Operation, usize, usize, Table); 64]; 2] = {
                let mut lut = [[(Operation::D, 0, 1, Table::D); 64]; 2];

                // table: the current DP table, D, C, or R (tables are standardized to right = true, C and R would be swapped for right = false)
                // trace: 2 bits, first bit is whether the max equals C table entry, second bit is
                // whether the max equals R table entry (vice versa for right = false)
                // trace2: 2 bits, first bit is whether the max in the C table is the gap beginning, second
                // bit is whether the max in the R table is the gap beginning (vice versa for right = false)
                // right: whether the current block contains vectors laid out vertically

                let mut right = 0;
                while right < 2 {
                    let mut trace = 0;
                    while trace < 4 {
                        let mut trace2 = 0;
                        while trace2 < 4 {
                            let mut table_idx = 0;
                            while table_idx < 3 {
                                let table = match table_idx {
                                    0b00 => Table::D,
                                    0b01 => Table::C,
                                    _ => Table::R
                                };

                                let res = if right == 1 {
                                    match (trace, trace2, table) {
                                        (_, 0b00 | 0b10, Table::C) => (Operation::D, 0, 1, Table::C), // C table gap extend
                                        (_, 0b01 | 0b11, Table::C) => (Operation::D, 0, 1, Table::D), // C table gap open
                                        (_, 0b00 | 0b01, Table::R) => (Operation::I, 1, 0, Table::R), // R table gap extend
                                        (_, 0b10 | 0b11, Table::R) => (Operation::I, 1, 0, Table::D), // R table gap open
                                        (0b00, _, Table::D) => (Operation::M, 1, 1, Table::D), // D table match/mismatch
                                        (0b01 | 0b11, 0b00 | 0b10, Table::D) => (Operation::D, 0, 1, Table::C), // D table C gap extend
                                        (0b01 | 0b11, 0b01 | 0b11, Table::D) => (Operation::D, 0, 1, Table::D), // D table C gap open
                                        (0b10, 0b00 | 0b01, Table::D) => (Operation::I, 1, 0, Table::R), // D table R gap extend
                                        (0b10, 0b10 | 0b11, Table::D) => (Operation::I, 1, 0, Table::D), // D table R gap open
                                        _ => (Operation::D, 0, 1, Table::D)
                                    }
                                } else {
                                    // everything is basically swapped (C/R and I/D) for down (right = false)
                                    match (trace, trace2, table) {
                                        (_, 0b00 | 0b10, Table::R) => (Operation::I, 1, 0, Table::R), // R table gap extend
                                        (_, 0b01 | 0b11, Table::R) => (Operation::I, 1, 0, Table::D), // R table gap open
                                        (_, 0b00 | 0b01, Table::C) => (Operation::D, 0, 1, Table::C), // C table gap extend
                                        (_, 0b10 | 0b11, Table::C) => (Operation::D, 0, 1, Table::D), // C table gap open
                                        (0b00, _, Table::D) => (Operation::M, 1, 1, Table::D), // D table match/mismatch
                                        (0b01 | 0b11, 0b00 | 0b10, Table::D) => (Operation::I, 1, 0, Table::R), // D table R gap extend
                                        (0b01 | 0b11, 0b01 | 0b11, Table::D) => (Operation::I, 1, 0, Table::D), // D table R gap open
                                        (0b10, 0b00 | 0b01, Table::D) => (Operation::D, 0, 1, Table::C), // D table C gap extend
                                        (0b10, 0b10 | 0b11, Table::D) => (Operation::D, 0, 1, Table::D), // D table C gap open
                                        _ => (Operation::I, 1, 0, Table::D)
                                    }
                                };

                                lut[right][(trace << 4) | (trace2 << 2) | (table as usize)] = res;
                                table_idx += 1;
                            }
                            trace2 += 1;
                        }
                        trace += 1;
                    }
                    right += 1;
                }

                lut
            };

            let mut table = Table::D;

            'outer: while i > 0 || j > 0 {
                // find the current block that contains (i, j)
                loop {
                    block_idx -= 1;
                    block_i = *self.block_start.as_ptr().add(block_idx * 2) as usize;
                    block_j = *self.block_start.as_ptr().add(block_idx * 2 + 1) as usize;
                    block_height = *self.block_size.as_ptr().add(block_idx * 2) as usize;
                    block_width = *self.block_size.as_ptr().add(block_idx * 2 + 1) as usize;
                    trace_idx -= block_width * block_height / L;

                    if i >= block_i && j >= block_j {
                        right = ((*self.right.as_ptr().add(block_idx / 64) >> (block_idx % 64)) & 0b1) as usize;
                        break;
                    }
                }

                // compute traceback within the current block
                let lut = &*OP_LUT.as_ptr().add(right);
                if right > 0 {
                    // right block
                    while i >= block_i && j >= block_j && (i > 0 || j > 0) {
                        if self.free_query_start_gaps && i == 0 {
                            // can do this because the row (i == 0) must be within right blocks
                            break 'outer;
                        }

                        let curr_i = i - block_i;
                        let curr_j = j - block_j;
                        let idx = trace_idx + curr_i / L + curr_j * (block_height / L);

                        if self.local_start && table == Table::D {
                            // terminate alignment on zero
                            let zero = ((*self.zero_mask.as_ptr().add(idx) >> ((curr_i % L) * 2)) & 0b1) > 0;
                            if zero {
                                break 'outer;
                            }
                        }

                        // build the index into the lookup table
                        let t = ((*self.trace.as_ptr().add(idx) >> ((curr_i % L) * 2)) & 0b11) as usize;
                        let t2 = ((*self.trace2.as_ptr().add(idx) >> ((curr_i % L) * 2)) & 0b11) as usize;
                        let lut_idx = (t << 4) | (t2 << 2) | (table as usize);
                        let lut_entry = &*lut.as_ptr().add(lut_idx);

                        let op = if EQ && lut_entry.0 == Operation::M {
                            if q.unwrap_unchecked().get(i) == r.unwrap_unchecked().get(j) {
                                Operation::Eq
                            } else {
                                Operation::X
                            }
                        } else {
                            lut_entry.0
                        };
                        i -= lut_entry.1;
                        j -= lut_entry.2;
                        table = lut_entry.3;
                        cigar.add(op);
                    }
                } else {
                    // down block
                    while i >= block_i && j >= block_j && (i > 0 || j > 0) {
                        let curr_i = i - block_i;
                        let curr_j = j - block_j;
                        let idx = trace_idx + curr_j / L + curr_i * (block_width / L);

                        if self.local_start && table == Table::D {
                            // terminate alignment on zero
                            let zero = ((*self.zero_mask.as_ptr().add(idx) >> ((curr_j % L) * 2)) & 0b1) > 0;
                            if zero {
                                break 'outer;
                            }
                        }

                        // build the index into the lookup table
                        let t = ((*self.trace.as_ptr().add(idx) >> ((curr_j % L) * 2)) & 0b11) as usize;
                        let t2 = ((*self.trace2.as_ptr().add(idx) >> ((curr_j % L) * 2)) & 0b11) as usize;
                        let lut_idx = (t << 4) | (t2 << 2) | (table as usize);
                        let lut_entry = &*lut.as_ptr().add(lut_idx);

                        let op = if EQ && lut_entry.0 == Operation::M {
                            if q.unwrap_unchecked().get(i) == r.unwrap_unchecked().get(j) {
                                Operation::Eq
                            } else {
                                Operation::X
                            }
                        } else {
                            lut_entry.0
                        };
                        i -= lut_entry.1;
                        j -= lut_entry.2;
                        table = lut_entry.3;
                        cigar.add(op);
                    }
                }
            }
        }
    }

    /// Return all of the rectangular regions that were calculated separately as
    /// block aligner shifts and grows.
    pub fn blocks(&self) -> Vec<Rectangle> {
        let mut res = Vec::with_capacity(self.block_idx);

        for i in 0..self.block_idx {
            unsafe {
                res.push(Rectangle {
                    row: *self.block_start.as_ptr().add(i * 2) as usize,
                    col: *self.block_start.as_ptr().add(i * 2 + 1) as usize,
                    height: *self.block_size.as_ptr().add(i * 2) as usize,
                    width: *self.block_size.as_ptr().add(i * 2 + 1) as usize
                });
            }
        }

        res
    }
}

/// A rectangular region.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct Rectangle {
    pub row: usize,
    pub col: usize,
    pub width: usize,
    pub height: usize
}

#[inline]
fn clamp(x: i32) -> i16 {
    cmp::min(cmp::max(x, i16::MIN as i32), i16::MAX as i32) as i16
}

#[inline]
fn div_ceil(n: usize, d: usize) -> usize {
    (n + d - 1) / d
}

/// Same alignment as SIMD vectors.
struct Aligned {
    layout: alloc::Layout,
    ptr: *const i16
}

impl Aligned {
    pub unsafe fn new(block_size: usize) -> Self {
        // custom alignment
        let layout = alloc::Layout::from_size_align_unchecked(block_size * 2, L_BYTES);
        let ptr = alloc::alloc_zeroed(layout) as *const i16;
        Self { layout, ptr }
    }

    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    pub unsafe fn clear(&mut self, block_size: usize) {
        let mut i = 0;
        while i < block_size {
            simd_store(self.ptr.add(i) as _, simd_set1_i16(MIN));
            i += L;
        }
    }

    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    #[inline]
    pub unsafe fn set_vec(&mut self, o: &Aligned, idx: usize) {
        simd_store(self.ptr.add(idx) as _, simd_load(o.as_ptr().add(idx) as _));
    }

    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    #[inline]
    pub unsafe fn copy_vec(&mut self, new_idx: usize, idx: usize) {
        simd_store(self.ptr.add(new_idx) as _, simd_load(self.ptr.add(idx) as _));
    }

    #[inline]
    pub fn get(&self, i: usize) -> i16 {
        unsafe { *self.ptr.add(i) }
    }

    #[allow(dead_code)]
    #[inline]
    pub fn set(&mut self, i: usize, v: i16) {
        unsafe { ptr::write(self.ptr.add(i) as _, v); }
    }

    #[inline]
    pub fn as_mut_ptr(&mut self) -> *mut i16 {
        self.ptr as _
    }

    #[inline]
    pub fn as_ptr(&self) -> *const i16 {
        self.ptr
    }
}

impl Drop for Aligned {
    fn drop(&mut self) {
        unsafe { alloc::dealloc(self.ptr as _, self.layout); }
    }
}

/// A padded string that helps avoid out of bounds access when using SIMD.
///
/// A single padding byte in inserted before the start of the string,
/// and `block_size` bytes are inserted after the end of the string.
#[derive(Clone, PartialEq, Debug)]
pub struct PaddedBytes {
    s: Vec<u8>,
    len: usize
}

impl PaddedBytes {
    /// Create an empty `PaddedBytes` instance that can hold byte strings
    /// of a specific size.
    pub fn new<M: Matrix>(len: usize, block_size: usize) -> Self {
        Self {
            s: vec![M::convert_char(M::NULL); 1 + len + block_size],
            len
        }
    }

    /// Modifies the bytes in place, filling in the rest of the memory with padding bytes.
    pub fn set_bytes<M: Matrix>(&mut self, b: &[u8], block_size: usize) {
        self.s[0] = M::convert_char(M::NULL);
        self.s[1..1 + b.len()].copy_from_slice(b);
        self.s[1..1 + b.len()].iter_mut().for_each(|c| *c = M::convert_char(*c));
        self.s[1 + b.len()..1 + b.len() + block_size].fill(M::convert_char(M::NULL));
        self.len = b.len();
    }

    /// Modifies the bytes in place in reverse, filling in the rest of the memory with padding bytes.
    pub fn set_bytes_rev<M: Matrix>(&mut self, b: &[u8], block_size: usize) {
        self.s[0] = M::convert_char(M::NULL);
        self.s[1..1 + b.len()].copy_from_slice(b);
        self.s[1..1 + b.len()].reverse();
        self.s[1..1 + b.len()].iter_mut().for_each(|c| *c = M::convert_char(*c));
        self.s[1 + b.len()..1 + b.len() + block_size].fill(M::convert_char(M::NULL));
        self.len = b.len();
    }

    /// Create from a byte slice.
    ///
    /// Make sure that `block_size` is greater than or equal to the upper bound
    /// block size used in the `Block::align` function.
    #[inline]
    pub fn from_bytes<M: Matrix>(b: &[u8], block_size: usize) -> Self {
        let mut v = b.to_owned();
        let len = v.len();
        v.insert(0, M::NULL);
        v.resize(v.len() + block_size, M::NULL);
        v.iter_mut().for_each(|c| *c = M::convert_char(*c));
        Self { s: v, len }
    }

    /// Create from the bytes in a string slice.
    ///
    /// Make sure that `block_size` is greater than or equal to the upper bound
    /// block size used in the `Block::align` function.
    #[inline]
    pub fn from_str<M: Matrix>(s: &str, block_size: usize) -> Self {
        Self::from_bytes::<M>(s.as_bytes(), block_size)
    }

    /// Create from the bytes in a string.
    ///
    /// Make sure that `block_size` is greater than or equal to the upper bound
    /// block size used in the `Block::align` function.
    #[inline]
    pub fn from_string<M: Matrix>(s: String, block_size: usize) -> Self {
        let mut v = s.into_bytes();
        let len = v.len();
        v.insert(0, M::NULL);
        v.resize(v.len() + block_size, M::NULL);
        v.iter_mut().for_each(|c| *c = M::convert_char(*c));
        Self { s: v, len }
    }

    /// Get the byte at a certain index (unchecked).
    #[inline]
    pub unsafe fn get(&self, i: usize) -> u8 {
        *self.s.as_ptr().add(i)
    }

    /// Set the byte at a certain index (unchecked).
    #[inline]
    pub unsafe fn set(&mut self, i: usize, c: u8) {
        *self.s.as_mut_ptr().add(i) = c;
    }

    /// Create a pointer to a specific index.
    #[inline]
    pub unsafe fn as_ptr(&self, i: usize) -> *const u8 {
        self.s.as_ptr().add(i)
    }

    /// Length of the original string (no padding).
    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }
}

/// Resulting score and alignment end position.
#[repr(C)]
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
    Grow
}

#[cfg(test)]
mod tests {
    use crate::scores::*;

    use super::*;

    #[test]
    fn test_no_x_drop() {
        let test_gaps = Gaps { open: -11, extend: -1 };

        let mut a = Block::<false, false>::new(100, 100, 16);

        let r = PaddedBytes::from_bytes::<AAMatrix>(b"", 16);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"", 16);
        a.align(&q, &r, &BLOSUM62, test_gaps, 16..=16, 0);
        assert_eq!(a.res().score, 0);

        let r = PaddedBytes::from_bytes::<AAMatrix>(b"AAAA", 16);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"", 16);
        a.align(&q, &r, &BLOSUM62, test_gaps, 16..=16, 0);
        assert_eq!(a.res().score, -14);

        let r = PaddedBytes::from_bytes::<AAMatrix>(b"", 16);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"AAAA", 16);
        a.align(&q, &r, &BLOSUM62, test_gaps, 16..=16, 0);
        assert_eq!(a.res().score, -14);

        let r = PaddedBytes::from_bytes::<AAMatrix>(b"AAAA", 16);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"AARA", 16);
        a.align(&q, &r, &BLOSUM62, test_gaps, 16..=16, 0);
        assert_eq!(a.res().score, 11);

        let r = PaddedBytes::from_bytes::<AAMatrix>(b"AAAAAAAA", 16);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"AARAAAA", 16);
        a.align(&q, &r, &BLOSUM62, test_gaps, 16..=16, 0);
        assert_eq!(a.res().score, 12);

        let r = PaddedBytes::from_bytes::<AAMatrix>(b"AAAA", 16);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"AAAA", 16);
        a.align(&q, &r, &BLOSUM62, test_gaps, 16..=16, 0);
        assert_eq!(a.res().score, 16);

        let r = PaddedBytes::from_bytes::<AAMatrix>(b"AAAA", 16);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"AARA", 16);
        a.align(&q, &r, &BLOSUM62, test_gaps, 16..=16, 0);
        assert_eq!(a.res().score, 11);

        let r = PaddedBytes::from_bytes::<AAMatrix>(b"AAAA", 16);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"RRRR", 16);
        a.align(&q, &r, &BLOSUM62, test_gaps, 16..=16, 0);
        assert_eq!(a.res().score, -4);

        let r = PaddedBytes::from_bytes::<AAMatrix>(b"AAAA", 16);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"AAA", 16);
        a.align(&q, &r, &BLOSUM62, test_gaps, 16..=16, 0);
        assert_eq!(a.res().score, 1);

        let test_gaps2 = Gaps { open: -2, extend: -1 };

        let r = PaddedBytes::from_bytes::<NucMatrix>(b"AAAN", 16);
        let q = PaddedBytes::from_bytes::<NucMatrix>(b"ATAA", 16);
        a.align(&q, &r, &NW1, test_gaps2, 16..=16, 0);
        assert_eq!(a.res().score, 0);

        let r = PaddedBytes::from_bytes::<NucMatrix>(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 16);
        let q = PaddedBytes::from_bytes::<NucMatrix>(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 16);
        a.align(&q, &r, &NW1, test_gaps2, 16..=16, 0);
        assert_eq!(a.res().score, 32);

        let r = PaddedBytes::from_bytes::<NucMatrix>(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 16);
        let q = PaddedBytes::from_bytes::<NucMatrix>(b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", 16);
        a.align(&q, &r, &NW1, test_gaps2, 16..=16, 0);
        assert_eq!(a.res().score, -32);

        let r = PaddedBytes::from_bytes::<NucMatrix>(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 16);
        let q = PaddedBytes::from_bytes::<NucMatrix>(b"TATATATATATATATATATATATATATATATA", 16);
        a.align(&q, &r, &NW1, test_gaps2, 16..=16, 0);
        assert_eq!(a.res().score, 0);

        let r = PaddedBytes::from_bytes::<NucMatrix>(b"TTAAAAAAATTTTTTTTTTTT", 16);
        let q = PaddedBytes::from_bytes::<NucMatrix>(b"TTTTTTTTAAAAAAATTTTTTTTT", 16);
        a.align(&q, &r, &NW1, test_gaps2, 16..=16, 0);
        assert_eq!(a.res().score, 7);

        let r = PaddedBytes::from_bytes::<NucMatrix>(b"AAAA", 16);
        let q = PaddedBytes::from_bytes::<NucMatrix>(b"C", 16);
        a.align(&q, &r, &NW1, test_gaps2, 16..=16, 0);
        assert_eq!(a.res().score, -5);
        a.align(&r, &q, &NW1, test_gaps2, 16..=16, 0);
        assert_eq!(a.res().score, -5);
    }

    #[test]
    fn test_x_drop() {
        let test_gaps = Gaps { open: -11, extend: -1 };

        let mut a = Block::<false, true>::new(100, 100, 16);

        let r = PaddedBytes::from_bytes::<AAMatrix>(b"", 16);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"", 16);
        a.align(&q, &r, &BLOSUM62, test_gaps, 16..=16, 1);
        assert_eq!(a.res(), AlignResult { score: 0, query_idx: 0, reference_idx: 0 });

        let r = PaddedBytes::from_bytes::<AAMatrix>(b"AAAA", 16);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"", 16);
        a.align(&q, &r, &BLOSUM62, test_gaps, 16..=16, 1);
        assert_eq!(a.res(), AlignResult { score: 0, query_idx: 0, reference_idx: 0 });

        let r = PaddedBytes::from_bytes::<AAMatrix>(b"", 16);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"AAAA", 16);
        a.align(&q, &r, &BLOSUM62, test_gaps, 16..=16, 1);
        assert_eq!(a.res(), AlignResult { score: 0, query_idx: 0, reference_idx: 0 });

        let r = PaddedBytes::from_bytes::<AAMatrix>(b"AAARRA", 16);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"AAAAAA", 16);
        a.align(&q, &r, &BLOSUM62, test_gaps, 16..=16, 1);
        assert_eq!(a.res(), AlignResult { score: 14, query_idx: 6, reference_idx: 6 });

        let r = PaddedBytes::from_bytes::<AAMatrix>(b"AAAAAAAAAAAAAAARRRRRRRRRRRRRRRRAAAAAAAAAAAAA", 16);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 16);
        a.align(&q, &r, &BLOSUM62, test_gaps, 16..=16, 1);
        assert_eq!(a.res(), AlignResult { score: 60, query_idx: 15, reference_idx: 15 });

        let mut a = Block::<true, true>::new(2048, 2048, 2048);
        let long_str = std::iter::repeat(b'A').take(2048).collect::<Vec<_>>();
        let r = PaddedBytes::from_bytes::<AAMatrix>(&long_str, 2048);
        let q = PaddedBytes::from_bytes::<AAMatrix>(&long_str, 2048);
        a.align(&q, &r, &BLOSUM62, test_gaps, 2048..=2048, 100);
        assert_eq!(a.res(), AlignResult { score: 8192, query_idx: 2048, reference_idx: 2048 });

        let mut a = Block::<true, true>::new(0, 0, 16);

        let r = PaddedBytes::from_bytes::<AAMatrix>(b"", 16);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"", 16);
        a.align(&q, &r, &BLOSUM62, test_gaps, 16..=16, 1);
        assert_eq!(a.res(), AlignResult { score: 0, query_idx: 0, reference_idx: 0 });

        let mut a = Block::<true, true>::new(4, 4, 16);

        let r = PaddedBytes::from_bytes::<AAMatrix>(b"AAAA", 16);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"", 16);
        a.align(&q, &r, &BLOSUM62, test_gaps, 16..=16, 1);
        assert_eq!(a.res(), AlignResult { score: 0, query_idx: 0, reference_idx: 0 });

        let r = PaddedBytes::from_bytes::<AAMatrix>(b"", 16);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"AAAA", 16);
        a.align(&q, &r, &BLOSUM62, test_gaps, 16..=16, 1);
        assert_eq!(a.res(), AlignResult { score: 0, query_idx: 0, reference_idx: 0 });
    }

    #[test]
    fn test_trace() {
        let test_gaps = Gaps { open: -11, extend: -1 };

        let mut cigar = Cigar::new(100, 100);

        let mut a = Block::<true, false>::new(100, 100, 16);

        let r = PaddedBytes::from_bytes::<AAMatrix>(b"AAARRA", 16);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"AAAAAA", 16);
        a.align(&q, &r, &BLOSUM62, test_gaps, 16..=16, 0);
        let res = a.res();
        assert_eq!(res, AlignResult { score: 14, query_idx: 6, reference_idx: 6 });
        a.trace().cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar);
        assert_eq!(cigar.to_string(), "3=2X1=");

        let r = PaddedBytes::from_bytes::<AAMatrix>(b"AAAA", 16);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"AAA", 16);
        a.align(&q, &r, &BLOSUM62, test_gaps, 16..=16, 0);
        let res = a.res();
        assert_eq!(res, AlignResult { score: 1, query_idx: 3, reference_idx: 4 });
        a.trace().cigar(res.query_idx, res.reference_idx, &mut cigar);
        assert_eq!(cigar.to_string(), "3M1D");

        let test_gaps2 = Gaps { open: -2, extend: -1 };

        let r = PaddedBytes::from_bytes::<NucMatrix>(b"TTAAAAAAATTTTTTTTTTTT", 16);
        let q = PaddedBytes::from_bytes::<NucMatrix>(b"TTTTTTTTAAAAAAATTTTTTTTT", 16);
        a.align(&q, &r, &NW1, test_gaps2, 16..=16, 0);
        let res = a.res();
        assert_eq!(res, AlignResult { score: 7, query_idx: 24, reference_idx: 21 });
        a.trace().cigar(res.query_idx, res.reference_idx, &mut cigar);
        assert_eq!(cigar.to_string(), "2M6I16M3D");

        let mut a = Block::<true, false>::new(100, 100, 32);
        let q = PaddedBytes::from_bytes::<NucMatrix>(b"AAAAAAAAATTGCGCT", 32);
        let r = PaddedBytes::from_bytes::<NucMatrix>(b"AAAAAAAAAGCGC", 32);

        a.align(&q, &r, &NW1, test_gaps2, 32..=32, 0);
        let res = a.res();
        assert_eq!(res, AlignResult { score: 8, query_idx: 16, reference_idx: 13 });
        a.trace().cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar);
        assert_eq!(cigar.to_string(), "9=2I4=1I");

        let matrix = NucMatrix::new_simple(2, -1);
        let test_gaps3 = Gaps { open: -5, extend: -2 };
        a.align(&q, &r, &matrix, test_gaps3, 32..=32, 0);
        let res = a.res();
        assert_eq!(res, AlignResult { score: 14, query_idx: 16, reference_idx: 13 });
        a.trace().cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar);
        assert_eq!(cigar.to_string(), "9=2I4=1I");
    }

    #[test]
    fn test_bytes() {
        let test_gaps = Gaps { open: -2, extend: -1 };

        let mut a = Block::<false, false>::new(100, 100, 16);

        let r = PaddedBytes::from_bytes::<ByteMatrix>(b"AAAaaA", 16);
        let q = PaddedBytes::from_bytes::<ByteMatrix>(b"AAAAAA", 16);
        a.align(&q, &r, &BYTES1, test_gaps, 16..=16, 0);
        assert_eq!(a.res().score, 2);

        let r = PaddedBytes::from_bytes::<ByteMatrix>(b"abcdefg", 16);
        let q = PaddedBytes::from_bytes::<ByteMatrix>(b"abdefg", 16);
        a.align(&q, &r, &BYTES1, test_gaps, 16..=16, 0);
        assert_eq!(a.res().score, 4);
    }

    #[test]
    fn test_profile() {
        let mut a = Block::<false, false>::new(100, 100, 16);
        let r = AAProfile::from_bytes(b"AAAA", 16, 1, -1, -1, 0, -1, -1);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"AAAA", 16);
        a.align_profile(&q, &r, 16..=16, 0);
        assert_eq!(a.res().score, 4);

        let r = AAProfile::from_bytes(b"AATTAA", 16, 1, -1, -1, 0, -1, -1);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"AAAA", 16);
        a.align_profile(&q, &r, 16..=16, 0);
        assert_eq!(a.res().score, 1);

        let r = AAProfile::from_bytes(b"AATTAA", 16, 1, -1, -1, -1, -1, -1);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"AAAA", 16);
        a.align_profile(&q, &r, 16..=16, 0);
        assert_eq!(a.res().score, 0);

        let mut a = Block::<true, false>::new(100, 100, 16);
        let mut cigar = Cigar::new(100, 100);

        let r = AAProfile::from_bytes(b"TTAAAAAAATTTTTTTTTTTT", 16, 1, -1, -1, 0, -1, -1);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"TTTTTTTTAAAAAAATTTTTTTTT", 16);
        a.align_profile(&q, &r, 16..=16, 0);
        let res = a.res();
        assert_eq!(res, AlignResult { score: 7, query_idx: 24, reference_idx: 21 });
        a.trace().cigar(res.query_idx, res.reference_idx, &mut cigar);
        assert_eq!(cigar.to_string(), "2M6I16M3D");

        let r = AAProfile::from_bytes(b"TTAAAAAAATTTTTTTTTTTT", 16, 1, -1, -1, -1, -1, -1);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"TTTTTTTTAAAAAAATTTTTTTTT", 16);
        a.align_profile(&q, &r, 16..=16, 0);
        let res = a.res();
        assert_eq!(res, AlignResult { score: 6, query_idx: 24, reference_idx: 21 });
        a.trace().cigar(res.query_idx, res.reference_idx, &mut cigar);
        assert_eq!(cigar.to_string(), "2M6I16M3D");

        let mut r = AAProfile::from_bytes(b"TTAAAAAAATTTTTTTTTTTT", 16, 1, -1, -2, -1, -1, -1);
        r.set_gap_close_C(17, -1);
        r.set_gap_close_C(19, 0);
        let q = PaddedBytes::from_bytes::<AAMatrix>(b"TTTTTTTTAAAAAAATTTTTTTTT", 16);
        a.align_profile(&q, &r, 16..=16, 0);
        let res = a.res();
        assert_eq!(res, AlignResult { score: 6, query_idx: 24, reference_idx: 21 });
        a.trace().cigar(res.query_idx, res.reference_idx, &mut cigar);
        assert_eq!(cigar.to_string(), "2M6I14M3D2M");
    }

    #[test]
    fn test_local_and_free_query_gaps() {
        let test_gaps = Gaps { open: -2, extend: -1 };

        let mut local = Block::<true, false, true, false>::new(100, 100, 32);
        let mut cigar = Cigar::new(100, 100);

        let r = PaddedBytes::from_bytes::<NucMatrix>(b"TTTTAAAAAA", 32);
        let q = PaddedBytes::from_bytes::<NucMatrix>(b"CCCCCCCCCCAAAAAA", 32);
        local.align(&q, &r, &NW1, test_gaps, 32..=32, 0);
        let res = local.res();
        assert_eq!(res, AlignResult { score: 6, query_idx: 16, reference_idx: 10 });
        local.trace().cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar);
        assert_eq!(cigar.to_string(), "6=");

        let mut local = Block::<true, true, true, false>::new(100, 100, 32);

        let r = PaddedBytes::from_bytes::<NucMatrix>(b"TTTTAAAAAATTTTTTT", 32);
        let q = PaddedBytes::from_bytes::<NucMatrix>(b"CCCCCCCCCCAAAAAACCCCCCCCCCCC", 32);
        local.align(&q, &r, &NW1, test_gaps, 32..=32, 100);
        let res = local.res();
        assert_eq!(res, AlignResult { score: 6, query_idx: 16, reference_idx: 10 });
        local.trace().cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar);
        assert_eq!(cigar.to_string(), "6=");

        let mut q_start = Block::<true, false, false, true>::new(100, 100, 32);

        let r = PaddedBytes::from_bytes::<NucMatrix>(b"CCCCCCCCCCAAAAAA", 32);
        let q = PaddedBytes::from_bytes::<NucMatrix>(b"AAAAAA", 32);
        q_start.align(&q, &r, &NW1, test_gaps, 32..=32, 0);
        let res = q_start.res();
        assert_eq!(res, AlignResult { score: 6, query_idx: 6, reference_idx: 16 });
        q_start.trace().cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar);
        assert_eq!(cigar.to_string(), "6=");

        let r = PaddedBytes::from_bytes::<NucMatrix>(b"CCCCCCCCCCAAATAA", 32);
        let q = PaddedBytes::from_bytes::<NucMatrix>(b"AAAAAA", 32);
        q_start.align(&q, &r, &NW1, test_gaps, 32..=32, 0);
        let res = q_start.res();
        assert_eq!(res, AlignResult { score: 4, query_idx: 6, reference_idx: 16 });
        q_start.trace().cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar);
        assert_eq!(cigar.to_string(), "3=1X2=");

        let mut q_end = Block::<true, false, false, false, true>::new(100, 100, 32);

        let r = PaddedBytes::from_bytes::<NucMatrix>(b"AAAAAACCCCCCCCCC", 32);
        let q = PaddedBytes::from_bytes::<NucMatrix>(b"AAAAAA", 32);
        q_end.align(&q, &r, &NW1, test_gaps, 32..=32, 0);
        let res = q_end.res();
        assert_eq!(res, AlignResult { score: 6, query_idx: 6, reference_idx: 6 });
        q_end.trace().cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar);
        assert_eq!(cigar.to_string(), "6=");

        let r = PaddedBytes::from_bytes::<NucMatrix>(b"AAATAACCCCCCCCCC", 32);
        let q = PaddedBytes::from_bytes::<NucMatrix>(b"AAAAAA", 32);
        q_end.align(&q, &r, &NW1, test_gaps, 32..=32, 0);
        let res = q_end.res();
        assert_eq!(res, AlignResult { score: 4, query_idx: 6, reference_idx: 6 });
        q_end.trace().cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar);
        assert_eq!(cigar.to_string(), "3=1X2=");
    }
}
