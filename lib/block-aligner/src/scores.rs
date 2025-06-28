//! Structs for representing match/mismatch scoring matrices.

#[cfg(feature = "simd_sse2")]
use crate::sse2::*;

#[cfg(feature = "simd_avx2")]
use crate::avx2::*;

#[cfg(feature = "simd_wasm")]
use crate::simd128::*;

#[cfg(feature = "simd_neon")]
use crate::neon::*;

use std::i8;

pub trait Matrix {
    /// Byte to use as padding.
    const NULL: u8;
    /// Create a new matrix with default (usually nonsense) values.
    ///
    /// Use `new_simple` to create a sensible scoring matrix.
    fn new() -> Self;
    /// Set the score for a pair of bytes.
    fn set(&mut self, a: u8, b: u8, score: i8);
    /// Get the score for a pair of bytes.
    fn get(&self, a: u8, b: u8) -> i8;
    /// Get the pointer for a specific index.
    fn as_ptr(&self, i: usize) -> *const i8;
    /// Get the scores for a certain byte and a certain SIMD vector of bytes.
    unsafe fn get_scores(&self, c: u8, v: HalfSimd, right: bool) -> Simd;
    /// Convert a byte to a better storage format that makes retrieving scores
    /// easier.
    fn convert_char(c: u8) -> u8;
}

/// Amino acid scoring matrix.
///
/// Supports characters `A` to `Z`. Lowercase characters are uppercased.
#[repr(C, align(32))]
#[derive(Clone, PartialEq, Debug)]
pub struct AAMatrix {
    scores: [i8; 27 * 32]
}

impl AAMatrix {
    /// Create a simple matrix with a certain match and mismatch score.
    pub const fn new_simple(match_score: i8, mismatch_score: i8) -> Self {
        let mut scores = [i8::MIN; 27 * 32];
        let mut i = b'A';
        while i <= b'Z' {
            let mut j = b'A';
            while j <= b'Z' {
                let idx = ((i - b'A') as usize) * 32 + ((j - b'A') as usize);
                scores[idx] = if i == j { match_score } else { mismatch_score };
                j += 1;
            }
            i += 1;
        }
        Self { scores }
    }

    /// Create an AAMatrix from a tab-separated table with no headers.
    ///
    /// Use `aa_order` to pass in the amino acids in order.
    pub fn from_tsv(tsv: &str, aa_order: &str) -> Self {
        let tsv = tsv.trim();
        let aa_order = aa_order.split_ascii_whitespace().map(|s| s.as_bytes()[0]).collect::<Vec<_>>();
        let mut res = Self::new();

        for (line, &a) in tsv.split("\n").zip(&aa_order) {
            for (score, &b) in line.split_ascii_whitespace().zip(&aa_order) {
                let score = score.parse::<i8>().unwrap();
                res.set(a, b, score);
            }
        }

        res
    }
}

impl Matrix for AAMatrix {
    const NULL: u8 = b'A' + 26u8;

    fn new() -> Self {
        Self { scores: [i8::MIN; 27 * 32] }
    }

    fn set(&mut self, a: u8, b: u8, score: i8) {
        let a = a.to_ascii_uppercase();
        let b = b.to_ascii_uppercase();
        assert!(b'A' <= a && a <= b'Z' + 1);
        assert!(b'A' <= b && b <= b'Z' + 1);
        let idx = ((a - b'A') as usize) * 32 + ((b - b'A') as usize);
        self.scores[idx] = score;
        let idx = ((b - b'A') as usize) * 32 + ((a - b'A') as usize);
        self.scores[idx] = score;
    }

    fn get(&self, a: u8, b: u8) -> i8 {
        let a = a.to_ascii_uppercase();
        let b = b.to_ascii_uppercase();
        assert!(b'A' <= a && a <= b'Z' + 1);
        assert!(b'A' <= b && b <= b'Z' + 1);
        let idx = ((a - b'A') as usize) * 32 + ((b - b'A') as usize);
        self.scores[idx]
    }

    #[inline]
    fn as_ptr(&self, i: usize) -> *const i8 {
        debug_assert!(i < 27);
        unsafe { self.scores.as_ptr().add(i * 32) }
    }

    // TODO: get rid of lookup for around half of the shifts by constructing position specific scoring matrix?
    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    #[inline]
    unsafe fn get_scores(&self, c: u8, v: HalfSimd, _right: bool) -> Simd {
        // efficiently lookup scores for each character in v
        let matrix_ptr = self.as_ptr(c as usize);
        let scores1 = lutsimd_load(matrix_ptr as *const LutSimd);
        let scores2 = lutsimd_load((matrix_ptr as *const LutSimd).add(1));
        halfsimd_lookup2_i16(scores1, scores2, v)
    }

    #[inline]
    fn convert_char(c: u8) -> u8 {
        let c = c.to_ascii_uppercase();
        assert!(c >= b'A' && c <= Self::NULL);
        c - b'A'
    }
}

/// Nucleotide scoring matrix.
///
/// Supports characters `A`, `C`, `G`, `N`, and `T`. Lowercase characters are uppercased.
#[repr(C, align(32))]
#[derive(Clone, PartialEq, Debug)]
pub struct NucMatrix {
    scores: [i8; 8 * 16]
}

impl NucMatrix {
    /// Create a simple matrix with a certain match and mismatch score.
    pub const fn new_simple(match_score: i8, mismatch_score: i8) -> Self {
        let mut scores = [i8::MIN; 8 * 16];
        let alpha = [b'A', b'T', b'C', b'G', b'N'];
        let mut i = 0;
        while i < alpha.len() {
            let mut j = 0;
            while j < alpha.len() {
                let idx = ((alpha[i] & 0b111) as usize) * 16 + ((alpha[j] & 0b1111) as usize);
                scores[idx] = if i == j { match_score } else { mismatch_score };
                j += 1;
            }
            i += 1;
        }
        Self { scores }
    }
}

impl Matrix for NucMatrix {
    const NULL: u8 = b'Z';

    fn new() -> Self {
        Self { scores: [i8::MIN; 8 * 16] }
    }

    fn set(&mut self, a: u8, b: u8, score: i8) {
        let a = a.to_ascii_uppercase();
        let b = b.to_ascii_uppercase();
        assert!(b'A' <= a && a <= b'Z');
        assert!(b'A' <= b && b <= b'Z');
        let idx = ((a & 0b111) as usize) * 16 + ((b & 0b1111) as usize);
        self.scores[idx] = score;
        let idx = ((b & 0b111) as usize) * 16 + ((a & 0b1111) as usize);
        self.scores[idx] = score;
    }

    fn get(&self, a: u8, b: u8) -> i8 {
        let a = a.to_ascii_uppercase();
        let b = b.to_ascii_uppercase();
        assert!(b'A' <= a && a <= b'Z');
        assert!(b'A' <= b && b <= b'Z');
        let idx = ((a & 0b111) as usize) * 16 + ((b & 0b1111) as usize);
        self.scores[idx]
    }

    #[inline]
    fn as_ptr(&self, i: usize) -> *const i8 {
        unsafe { self.scores.as_ptr().add((i & 0b111) * 16) }
    }

    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    #[inline]
    unsafe fn get_scores(&self, c: u8, v: HalfSimd, _right: bool) -> Simd {
        // efficiently lookup scores for each character in v
        let matrix_ptr = self.as_ptr(c as usize);
        let scores = lutsimd_load(matrix_ptr as *const LutSimd);
        halfsimd_lookup1_i16(scores, v)
    }

    #[inline]
    fn convert_char(c: u8) -> u8 {
        let c = c.to_ascii_uppercase();
        assert!(c >= b'A' && c <= Self::NULL);
        c
    }
}

/// Arbitrary bytes scoring matrix.
#[repr(C)]
#[derive(Clone, PartialEq, Debug)]
pub struct ByteMatrix {
    match_score: i8,
    mismatch_score: i8
}

impl ByteMatrix {
    /// Create a simple matrix with a certain match and mismatch score.
    pub const fn new_simple(match_score: i8, mismatch_score: i8) -> Self {
        Self { match_score, mismatch_score }
    }
}

impl Matrix for ByteMatrix {
    /// May lead to inaccurate results with x drop alignment,
    /// if the block reaches the ends of the strings.
    ///
    /// Avoid using `ByteMatrix` with x drop alignment.
    const NULL: u8 = b'\0';

    fn new() -> Self {
        Self { match_score: i8::MIN, mismatch_score: i8::MIN }
    }

    fn set(&mut self, _a: u8, _b: u8, _score: i8) {
        unimplemented!();
    }

    fn get(&self, a: u8, b: u8) -> i8 {
        if a == b { self.match_score } else { self.mismatch_score }
    }

    #[inline]
    fn as_ptr(&self, _i: usize) -> *const i8 {
        unimplemented!()
    }

    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    #[inline]
    unsafe fn get_scores(&self, c: u8, v: HalfSimd, _right: bool) -> Simd {
        let match_scores = halfsimd_set1_i8(self.match_score);
        let mismatch_scores = halfsimd_set1_i8(self.mismatch_score);
        halfsimd_lookup_bytes_i16(match_scores, mismatch_scores, halfsimd_set1_i8(c as i8), v)
    }

    #[inline]
    fn convert_char(c: u8) -> u8 {
        c
    }
}

/// Match = 1, mismatch = -1.
#[cfg_attr(not(target_arch = "wasm32"), no_mangle)]
pub static NW1: NucMatrix = NucMatrix::new_simple(1, -1);

#[cfg_attr(not(target_arch = "wasm32"), no_mangle)]
pub static BLOSUM45: AAMatrix = AAMatrix { scores: include!("../matrices/BLOSUM45") };

#[cfg_attr(not(target_arch = "wasm32"), no_mangle)]
pub static BLOSUM50: AAMatrix = AAMatrix { scores: include!("../matrices/BLOSUM50") };

#[cfg_attr(not(target_arch = "wasm32"), no_mangle)]
pub static BLOSUM62: AAMatrix = AAMatrix { scores: include!("../matrices/BLOSUM62") };

#[cfg_attr(not(target_arch = "wasm32"), no_mangle)]
pub static BLOSUM80: AAMatrix = AAMatrix { scores: include!("../matrices/BLOSUM80") };

#[cfg_attr(not(target_arch = "wasm32"), no_mangle)]
pub static BLOSUM90: AAMatrix = AAMatrix { scores: include!("../matrices/BLOSUM90") };

#[cfg_attr(not(target_arch = "wasm32"), no_mangle)]
pub static PAM100: AAMatrix = AAMatrix { scores: include!("../matrices/PAM100") };

#[cfg_attr(not(target_arch = "wasm32"), no_mangle)]
pub static PAM120: AAMatrix = AAMatrix { scores: include!("../matrices/PAM120") };

#[cfg_attr(not(target_arch = "wasm32"), no_mangle)]
pub static PAM160: AAMatrix = AAMatrix { scores: include!("../matrices/PAM160") };

#[cfg_attr(not(target_arch = "wasm32"), no_mangle)]
pub static PAM200: AAMatrix = AAMatrix { scores: include!("../matrices/PAM200") };

#[cfg_attr(not(target_arch = "wasm32"), no_mangle)]
pub static PAM250: AAMatrix = AAMatrix { scores: include!("../matrices/PAM250") };

/// Match = 1, mismatch = -1.
#[cfg_attr(not(target_arch = "wasm32"), no_mangle)]
pub static BYTES1: ByteMatrix = ByteMatrix::new_simple(1, -1);

/*pub trait ScoreParams {
    const GAP_OPEN: i8;
    const GAP_EXTEND: i8;
    const I: usize;
}

pub struct Params<const GAP_OPEN: i8, const GAP_EXTEND: i8, const I: usize>;

impl<const GAP_OPEN: i8, const GAP_EXTEND: i8, const I: usize> ScoreParams for Params<{ GAP_OPEN }, { GAP_EXTEND }, { I }> {
    const GAP_OPEN: i8 = GAP_OPEN;
    const GAP_EXTEND: i8 = GAP_EXTEND;
    const I: usize = I;
}

pub type GapParams<const GAP_OPEN: i8, const GAP_EXTEND: i8> = Params<{ GAP_OPEN }, { GAP_EXTEND }, 0>;*/

/// Open and extend gap costs.
///
/// Open cost must include the extend cost. For example, with `Gaps { open: -11, extend: -1 }`,
/// a gap of length 1 costs -11, and a gap of length 2 costs -12.
#[derive(Copy, Clone, PartialEq, Debug)]
#[repr(C)]
pub struct Gaps {
    pub open: i8,
    pub extend: i8
}

#[allow(non_snake_case)]
pub trait Profile {
    /// Byte to use as padding.
    const NULL: u8;

    /// Create a new profile of a specific length, with default (large negative) values.
    ///
    /// Note that internally, the created profile is longer than a conventional position-specific scoring
    /// matrix (and `str_len`) by 1, so the profile will have the same length as the number of
    /// columns in the DP matrix.
    /// The first column of scores in the profile should be large negative values (padding).
    /// This allows gap open costs to be specified for the first column of the DP matrix.
    fn new(str_len: usize, block_size: usize, gap_extend: i8) -> Self;
    /// Create a new profile from a byte string.
    fn from_bytes(b: &[u8], block_size: usize, match_score: i8, mismatch_score: i8, gap_open_C: i8, gap_close_C: i8, gap_open_R: i8, gap_extend: i8) -> Self;

    /// Get the length of the profile.
    fn len(&self) -> usize;
    /// Clear the profile so it can be reused for profile lengths less than or equal
    /// to the length this struct was created with.
    fn clear(&mut self, str_len: usize, block_size: usize);

    /// Set the score for a position and byte.
    ///
    /// The profile should be first `clear`ed before it is reused with different lengths.
    ///
    /// The first column (`i = 0`) should be padded with large negative values.
    /// Therefore, set values starting from `i = 1`.
    fn set(&mut self, i: usize, b: u8, score: i8);
    /// Set the scores for all positions in the position specific scoring matrix.
    ///
    /// The profile should be first `clear`ed before it is reused with different lengths.
    ///
    /// Use `order` to specify the order of bytes that is used in the `scores` matrix.
    /// Scores (in `scores`) should be stored in row-major order, where each row is a different position
    /// and each column is a different byte.
    ///
    /// Use `left_shift` and `right_shift` to scale all the scores.
    fn set_all(&mut self, order: &[u8], scores: &[i8], left_shift: usize, right_shift: usize);
    /// Set the scores for all positions in reverse in the position specific scoring matrix.
    ///
    /// The profile should be first `clear`ed before it is reused with different lengths.
    ///
    /// Use `order` to specify the order of bytes that is used in the `scores` matrix.
    /// Scores (in `scores`) should be stored in row-major order, where each row is a different position
    /// and each column is a different byte.
    ///
    /// Use `left_shift` and `right_shift` to scale all the scores.
    fn set_all_rev(&mut self, order: &[u8], scores: &[i8], left_shift: usize, right_shift: usize);

    /// Set the gap open cost for a column.
    ///
    /// When aligning a sequence `q` to a profile `r`, this is the gap open cost at column `i` for a
    /// column transition in the DP matrix with `|q| + 1` rows and `|r| + 1` columns.
    /// This represents starting a gap in `q`.
    fn set_gap_open_C(&mut self, i: usize, gap: i8);
    /// Set the gap close cost for a column.
    ///
    /// When aligning a sequence `q` to a profile `r`, this is the gap close cost at column `i` for
    /// ending column transitions in the DP matrix with `|q| + 1` rows and `|r| + 1` columns.
    /// This represents ending a gap in `q`.
    fn set_gap_close_C(&mut self, i: usize, gap: i8);
    /// Set the gap open cost for a row.
    ///
    /// When aligning a sequence `q` to a profile `r`, this is the gap open cost at column `i` for
    /// a row transition in the DP matrix with `|q| + 1` rows and `|r| + 1` columns.
    /// This represents starting a gap in `r`.
    fn set_gap_open_R(&mut self, i: usize, gap: i8);

    /// Set the gap open cost for all column transitions.
    fn set_all_gap_open_C(&mut self, gap: i8);
    /// Set the gap close cost for all column transitions.
    fn set_all_gap_close_C(&mut self, gap: i8);
    /// Set the gap open cost for all row transitions.
    fn set_all_gap_open_R(&mut self, gap: i8);

    /// Get the score for a position and byte.
    fn get(&self, i: usize, b: u8) -> i8;
    /// Get the gap extend cost.
    fn get_gap_extend(&self) -> i8;
    /// Get the pointer for a specific index.
    fn as_ptr_pos(&self, i: usize) -> *const i8;
    /// Get the pointer for a specific amino acid.
    fn as_ptr_aa(&self, a: usize) -> *const i16;

    /// Get the scores for a certain SIMD vector of bytes at a specific position in the profile.
    unsafe fn get_scores_pos(&self, i: usize, v: HalfSimd, right: bool) -> Simd;
    /// Get the scores for a certain byte starting at a specific position in the profile.
    unsafe fn get_scores_aa(&self, i: usize, c: u8, right: bool) -> Simd;

    /// Get the gap open cost for a column.
    unsafe fn get_gap_open_right_C(&self, i: usize) -> Simd;
    /// Get the gap close cost for a column.
    unsafe fn get_gap_close_right_C(&self, i: usize) -> Simd;
    /// Get the gap open cost for a row.
    unsafe fn get_gap_open_right_R(&self, i: usize) -> Simd;

    /// Get the gap open cost for a column.
    unsafe fn get_gap_open_down_C(&self, i: usize) -> Simd;
    /// Get the gap close cost for a column.
    unsafe fn get_gap_close_down_C(&self, i: usize) -> Simd;
    /// Get the gap open cost for a row.
    unsafe fn get_gap_open_down_R(&self, i: usize) -> Simd;

    /// Convert a byte to a better storage format that makes retrieving scores
    /// easier.
    fn convert_char(c: u8) -> u8;
}

/// Amino acid position specific scoring matrix.
///
/// Supports characters `A` to `Z`. Lowercase characters are uppercased.
#[allow(non_snake_case)]
#[derive(Clone, PartialEq, Debug)]
pub struct AAProfile {
    aa_pos: Vec<i16>,
    pos_aa: Vec<i8>,
    gap_extend: i8,
    pos_gap_open_C: Vec<i16>,
    pos_gap_close_C: Vec<i16>,
    pos_gap_open_R: Vec<i16>,
    // length used for underlying allocated vectors
    max_len: usize,
    // length used for the current padded profile
    curr_len: usize,
    // length of the profile without padding (same length as the consensus sequence of the position
    // specific scoring matrix)
    str_len: usize
}

impl Profile for AAProfile {
    const NULL: u8 = b'A' + 26u8;

    fn new(str_len: usize, block_size: usize, gap_extend: i8) -> Self {
        let max_len = str_len + block_size + 1;
        Self {
            aa_pos: vec![i8::MIN as i16; 32 * max_len],
            pos_aa: vec![i8::MIN; max_len * 32],
            gap_extend,
            pos_gap_open_C: vec![i8::MIN as i16; max_len],
            pos_gap_close_C: vec![i8::MIN as i16; max_len],
            pos_gap_open_R: vec![i8::MIN as i16; max_len],
            max_len,
            curr_len: max_len,
            str_len
        }
    }

    #[allow(non_snake_case)]
    fn from_bytes(b: &[u8], block_size: usize, match_score: i8, mismatch_score: i8, gap_open_C: i8, gap_close_C: i8, gap_open_R: i8, gap_extend: i8) -> Self {
        let mut res = Self::new(b.len(), block_size, gap_extend);

        for i in 0..b.len() {
            for c in b'A'..=b'Z' {
                res.set(i + 1, c, if c == b[i] { match_score } else { mismatch_score });
            }
        }

        for i in 0..b.len() + 1 {
            res.set_gap_open_C(i, gap_open_C);
            res.set_gap_close_C(i, gap_close_C);
            res.set_gap_open_R(i, gap_open_R);
        }

        res
    }

    fn len(&self) -> usize {
        self.str_len
    }

    fn clear(&mut self, str_len: usize, block_size: usize) {
        let curr_len = str_len + block_size + 1;
        assert!(curr_len <= self.max_len);
        self.aa_pos[..32 * curr_len].fill(i8::MIN as i16);
        self.pos_aa[..curr_len * 32].fill(i8::MIN);
        self.pos_gap_open_C[..curr_len].fill(i8::MIN as i16);
        self.pos_gap_close_C[..curr_len].fill(i8::MIN as i16);
        self.pos_gap_open_R[..curr_len].fill(i8::MIN as i16);
        self.str_len = str_len;
        self.curr_len = curr_len;
    }

    fn set(&mut self, i: usize, b: u8, score: i8) {
        let b = b.to_ascii_uppercase();
        assert!(b'A' <= b && b <= b'Z' + 1);
        let idx = i * 32 + ((b - b'A') as usize);
        self.pos_aa[idx] = score;
        let idx = ((b - b'A') as usize) * self.curr_len + i;
        self.aa_pos[idx] = score as i16;
    }

    fn set_all(&mut self, order: &[u8], scores: &[i8], left_shift: usize, right_shift: usize) {
        self.set_all_core::<false>(order, scores, left_shift, right_shift);
    }

    fn set_all_rev(&mut self, order: &[u8], scores: &[i8], left_shift: usize, right_shift: usize) {
        self.set_all_core::<true>(order, scores, left_shift, right_shift);
    }

    fn set_gap_open_C(&mut self, i: usize, gap: i8) {
        assert!(gap < 0, "Gap open cost must be negative!");
        self.pos_gap_open_C[i] = gap as i16;
    }

    fn set_gap_close_C(&mut self, i: usize, gap: i8) {
        self.pos_gap_close_C[i] = gap as i16;
    }

    fn set_gap_open_R(&mut self, i: usize, gap: i8) {
        assert!(gap < 0, "Gap open cost must be negative!");
        self.pos_gap_open_R[i] = gap as i16;
    }

    fn set_all_gap_open_C(&mut self, gap: i8) {
        assert!(gap < 0, "Gap open cost must be negative!");
        self.pos_gap_open_C[..self.str_len + 1].fill(gap as i16);
    }

    fn set_all_gap_close_C(&mut self, gap: i8) {
        self.pos_gap_close_C[..self.str_len + 1].fill(gap as i16);
    }

    fn set_all_gap_open_R(&mut self, gap: i8) {
        assert!(gap < 0, "Gap open cost must be negative!");
        self.pos_gap_open_R[..self.str_len + 1].fill(gap as i16);
    }

    fn get(&self, i: usize, b: u8) -> i8 {
        let b = b.to_ascii_uppercase();
        assert!(b'A' <= b && b <= b'Z' + 1);
        let idx = i * 32 + ((b - b'A') as usize);
        self.pos_aa[idx]
    }

    fn get_gap_extend(&self) -> i8 {
        self.gap_extend
    }

    #[inline]
    fn as_ptr_pos(&self, i: usize) -> *const i8 {
        debug_assert!(i < self.curr_len);
        unsafe { self.pos_aa.as_ptr().add(i * 32) }
    }

    #[inline]
    fn as_ptr_aa(&self, a: usize) -> *const i16 {
        debug_assert!(a < 27);
        unsafe { self.aa_pos.as_ptr().add(a * self.curr_len) }
    }

    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    #[inline]
    unsafe fn get_scores_pos(&self, i: usize, v: HalfSimd, _right: bool) -> Simd {
        // efficiently lookup scores for each character in v
        let matrix_ptr = self.as_ptr_pos(i);
        let scores1 = lutsimd_loadu(matrix_ptr as *const LutSimd);
        let scores2 = lutsimd_loadu((matrix_ptr as *const LutSimd).add(1));
        halfsimd_lookup2_i16(scores1, scores2, v)
    }

    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    #[inline]
    unsafe fn get_scores_aa(&self, i: usize, c: u8, _right: bool) -> Simd {
        let matrix_ptr = self.as_ptr_aa(c as usize);
        simd_loadu(matrix_ptr.add(i) as *const Simd)
    }

    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    #[inline]
    unsafe fn get_gap_open_right_C(&self, i: usize) -> Simd {
        simd_set1_i16(*self.pos_gap_open_C.as_ptr().add(i))
    }

    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    #[inline]
    unsafe fn get_gap_close_right_C(&self, i: usize) -> Simd {
        simd_set1_i16(*self.pos_gap_close_C.as_ptr().add(i))
    }

    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    #[inline]
    unsafe fn get_gap_open_right_R(&self, i: usize) -> Simd {
        simd_set1_i16(*self.pos_gap_open_R.as_ptr().add(i))
    }

    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    #[inline]
    unsafe fn get_gap_open_down_C(&self, i: usize) -> Simd {
        simd_loadu(self.pos_gap_open_C.as_ptr().add(i) as *const Simd)
    }

    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    #[inline]
    unsafe fn get_gap_close_down_C(&self, i: usize) -> Simd {
        simd_loadu(self.pos_gap_close_C.as_ptr().add(i) as *const Simd)
    }

    #[cfg_attr(feature = "simd_sse2", target_feature(enable = "sse2"))]
    #[cfg_attr(feature = "simd_avx2", target_feature(enable = "avx2"))]
    #[cfg_attr(feature = "simd_wasm", target_feature(enable = "simd128"))]
    #[cfg_attr(feature = "simd_neon", target_feature(enable = "neon"))]
    #[inline]
    unsafe fn get_gap_open_down_R(&self, i: usize) -> Simd {
        simd_loadu(self.pos_gap_open_R.as_ptr().add(i) as *const Simd)
    }

    #[inline]
    fn convert_char(c: u8) -> u8 {
        let c = c.to_ascii_uppercase();
        assert!(c >= b'A' && c <= Self::NULL);
        c - b'A'
    }
}

impl AAProfile {
    fn set_all_core<const REV: bool>(&mut self, order: &[u8], scores: &[i8], left_shift: usize, right_shift: usize) {
        #[repr(align(32))]
        struct A([u8; 32]);
        let mut o = A([Self::NULL - b'A'; 32]);
        assert!(order.len() <= 32);

        for (i, &b) in order.iter().enumerate() {
            let b = b.to_ascii_uppercase();
            assert!(b'A' <= b && b <= b'Z' + 1);
            o.0[i] = b - b'A';
        }
        assert_eq!(scores.len() / order.len(), self.str_len);

        let mut i = if REV { self.str_len } else { 1 };
        let mut score_idx = 0;

        while if REV { i >= 1 } else { i <= self.str_len } {
            let mut j = 0;

            while j < order.len() {
                unsafe {
                    let score = ((*scores.as_ptr().add(score_idx)) << left_shift) >> right_shift;
                    let b = *o.0.as_ptr().add(j) as usize;
                    *self.pos_aa.as_mut_ptr().add(i * 32 + b) = score;
                    *self.aa_pos.as_mut_ptr().add(b * self.curr_len + i) = score as i16;
                }

                score_idx += 1;
                j += 1;
            }

            if REV {
                i -= 1;
            } else {
                i += 1;
            }
        }
    }
}
