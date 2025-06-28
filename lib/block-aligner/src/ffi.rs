//! C bindings for block aligner.
//!
//! Generics are monomorphised manually.
//!
//! Nucleotide and arbitrary byte alignment do not have bindings yet.

use std::ffi::c_void;

use crate::scan_block::*;
use crate::scores::*;
use crate::cigar::*;

// avoid generics by using void pointer and monomorphism
/// A handle for a block in block aligner.
pub type BlockHandle = *mut c_void;

/// Represents a range that has inclusive lower and upper bounds.
#[derive(Copy, Clone, PartialEq)]
#[repr(C)]
pub struct SizeRange {
    pub min: usize,
    pub max: usize
}


// AAMatrix

/// Create a new simple AAMatrix with custom match and mismatch scores.
///
/// Note that the match score must be positive and the mismatch score must be negative.
#[no_mangle]
pub unsafe extern fn block_new_simple_aamatrix(match_score: i8, mismatch_score: i8) -> *mut AAMatrix {
    let matrix = Box::new(AAMatrix::new_simple(match_score, mismatch_score));
    Box::into_raw(matrix)
}

/// Set an entry in the AAMatrix.
#[no_mangle]
pub unsafe extern fn block_set_aamatrix(matrix: *mut AAMatrix, a: u8, b: u8, score: i8) {
    let matrix = &mut *matrix;
    matrix.set(a, b, score);
}

/// Frees an AAMatrix.
#[no_mangle]
pub unsafe extern fn block_free_aamatrix(matrix: *mut AAMatrix) {
    drop(Box::from_raw(matrix));
}


// AAProfile

/// Create a new profile of a specific length, with default (large negative) values.
///
/// Note that internally, the created profile is longer than a conventional position-specific scoring
/// matrix (and `str_len`) by 1, so the profile will have the same length as the number of
/// columns in the DP matrix.
/// The first column of scores in the profile should be large negative values (padding).
/// This allows gap open costs to be specified for the first column of the DP matrix.
#[no_mangle]
pub unsafe extern fn block_new_aaprofile(str_len: usize, block_size: usize, gap_extend: i8) -> *mut AAProfile {
    let profile = Box::new(AAProfile::new(str_len, block_size, gap_extend));
    Box::into_raw(profile)
}

/// Get the length of the profile.
#[no_mangle]
pub unsafe extern fn block_len_aaprofile(profile: *const AAProfile) -> usize {
    let profile = &*profile;
    profile.len()
}

/// Clear the profile so it can be reused for profile lengths less than or equal
/// to the length this struct was created with.
#[no_mangle]
pub unsafe extern fn block_clear_aaprofile(profile: *mut AAProfile, str_len: usize, block_size: usize) {
    let profile = &mut *profile;
    profile.clear(str_len, block_size);
}

/// Set the score for a position and byte.
///
/// The profile should be first `clear`ed before it is reused with different lengths.
///
/// The first column (`i = 0`) should be padded with large negative values.
/// Therefore, set values starting from `i = 1`.
#[no_mangle]
pub unsafe extern fn block_set_aaprofile(profile: *mut AAProfile, i: usize, b: u8, score: i8) {
    let profile = &mut *profile;
    profile.set(i, b, score);
}

/// Set the scores for all positions in the position specific scoring matrix.
///
/// The profile should be first `clear`ed before it is reused with different lengths.
///
/// Use `order` to specify the order of bytes that is used in the `scores` matrix.
/// Scores (in `scores`) should be stored in row-major order, where each row is a different position
/// and each column is a different byte.
#[no_mangle]
pub unsafe extern fn block_set_all_aaprofile(profile: *mut AAProfile, order: *const u8, order_len: usize, scores: *const i8, scores_len: usize, left_shift: usize, right_shift: usize) {
    let profile = &mut *profile;
    let order = std::slice::from_raw_parts(order, order_len);
    let scores = std::slice::from_raw_parts(scores, scores_len);
    profile.set_all(order, scores, left_shift, right_shift);
}

/// Set the scores for all positions in reverse in the position specific scoring matrix.
///
/// The profile should be first `clear`ed before it is reused with different lengths.
///
/// Use `order` to specify the order of bytes that is used in the `scores` matrix.
/// Scores (in `scores`) should be stored in row-major order, where each row is a different position
/// and each column is a different byte.
#[no_mangle]
pub unsafe extern fn block_set_all_rev_aaprofile(profile: *mut AAProfile, order: *const u8, order_len: usize, scores: *const i8, scores_len: usize, left_shift: usize, right_shift: usize) {
    let profile = &mut *profile;
    let order = std::slice::from_raw_parts(order, order_len);
    let scores = std::slice::from_raw_parts(scores, scores_len);
    profile.set_all_rev(order, scores, left_shift, right_shift);
}

/// Set the gap open cost for a column.
///
/// When aligning a sequence `q` to a profile `r`, this is the gap open cost at column `i` for a
/// column transition in the DP matrix with `|q| + 1` rows and `|r| + 1` columns.
/// This represents starting a gap in `q`.
#[no_mangle]
pub unsafe extern fn block_set_gap_open_C_aaprofile(profile: *mut AAProfile, i: usize, gap: i8) {
    let profile = &mut *profile;
    profile.set_gap_open_C(i, gap);
}

/// Set the gap close cost for a column.
///
/// When aligning a sequence `q` to a profile `r`, this is the gap close cost at column `i` for
/// ending column transitions in the DP matrix with `|q| + 1` rows and `|r| + 1` columns.
/// This represents ending a gap in `q`.
#[no_mangle]
pub unsafe extern fn block_set_gap_close_C_aaprofile(profile: *mut AAProfile, i: usize, gap: i8) {
    let profile = &mut *profile;
    profile.set_gap_close_C(i, gap);
}

/// Set the gap open cost for a row.
///
/// When aligning a sequence `q` to a profile `r`, this is the gap open cost at column `i` for
/// a row transition in the DP matrix with `|q| + 1` rows and `|r| + 1` columns.
/// This represents starting a gap in `r`.
#[no_mangle]
pub unsafe extern fn block_set_gap_open_R_aaprofile(profile: *mut AAProfile, i: usize, gap: i8) {
    let profile = &mut *profile;
    profile.set_gap_open_R(i, gap);
}

/// Set the gap open cost for all column transitions.
#[no_mangle]
pub unsafe extern fn block_set_all_gap_open_C_aaprofile(profile: *mut AAProfile, gap: i8) {
    let profile = &mut *profile;
    profile.set_all_gap_open_C(gap);
}

/// Set the gap close cost for all column transitions.
#[no_mangle]
pub unsafe extern fn block_set_all_gap_close_C_aaprofile(profile: *mut AAProfile, gap: i8) {
    let profile = &mut *profile;
    profile.set_all_gap_close_C(gap);
}

/// Set the gap open cost for all row transitions.
#[no_mangle]
pub unsafe extern fn block_set_all_gap_open_R_aaprofile(profile: *mut AAProfile, gap: i8) {
    let profile = &mut *profile;
    profile.set_all_gap_open_R(gap);
}

/// Get the score for a position and byte.
#[no_mangle]
pub unsafe extern fn block_get_aaprofile(profile: *const AAProfile, i: usize, b: u8) -> i8 {
    let profile = &*profile;
    profile.get(i, b)
}

/// Get the gap extend cost.
#[no_mangle]
pub unsafe extern fn block_get_gap_extend_aaprofile(profile: *const AAProfile) -> i8 {
    let profile = &*profile;
    profile.get_gap_extend()
}

/// Frees an AAProfile.
#[no_mangle]
pub unsafe extern fn block_free_aaprofile(profile: *mut AAProfile) {
    drop(Box::from_raw(profile));
}


// CIGAR

/// Create a new empty CIGAR string.
#[no_mangle]
pub unsafe extern fn block_new_cigar(query_len: usize, reference_len: usize) -> *mut Cigar {
    let cigar = Box::new(Cigar::new(query_len, reference_len));
    Box::into_raw(cigar)
}

/// Get the operation at a certain index in a CIGAR string.
#[no_mangle]
pub unsafe extern fn block_get_cigar(cigar: *const Cigar, i: usize) -> OpLen {
    let cigar_str = &*cigar;
    cigar_str.get(i)
}

/// Get the length of a CIGAR string.
#[no_mangle]
pub unsafe extern fn block_len_cigar(cigar: *const Cigar) -> usize {
    let cigar_str = &*cigar;
    cigar_str.len()
}

/// Frees a CIGAR string.
#[no_mangle]
pub unsafe extern fn block_free_cigar(cigar: *mut Cigar) {
    drop(Box::from_raw(cigar));
}


// PaddedBytes

/// Create a new empty padded amino acid string.
#[no_mangle]
pub unsafe extern fn block_new_padded_aa(len: usize, max_size: usize) -> *mut PaddedBytes {
    let padded_bytes = Box::new(PaddedBytes::new::<AAMatrix>(len, max_size));
    Box::into_raw(padded_bytes)
}

/// Write to a padded amino acid string.
#[no_mangle]
pub unsafe extern fn block_set_bytes_padded_aa(padded: *mut PaddedBytes, s: *const u8, len: usize, max_size: usize) {
    let bytes = std::slice::from_raw_parts(s, len);
    let padded_bytes = &mut *padded;
    padded_bytes.set_bytes::<AAMatrix>(bytes, max_size);
}

/// Write to a padded amino acid string, in reverse.
#[no_mangle]
pub unsafe extern fn block_set_bytes_rev_padded_aa(padded: *mut PaddedBytes, s: *const u8, len: usize, max_size: usize) {
    let bytes = std::slice::from_raw_parts(s, len);
    let padded_bytes = &mut *padded;
    padded_bytes.set_bytes_rev::<AAMatrix>(bytes, max_size);
}

/// Frees a padded amino acid string.
#[no_mangle]
pub unsafe extern fn block_free_padded_aa(padded: *mut PaddedBytes) {
    drop(Box::from_raw(padded));
}


// Block

macro_rules! gen_functions {
    ($new_name:ident, $new_doc:expr,
     $align_name:ident, $align_doc:expr,
     $align_profile_name:ident, $align_profile_doc:expr,
     $res_name:ident, $res_doc:expr,
     $trace_name:ident, $trace_doc:expr,
     $trace_eq_name:ident, $trace_eq_doc:expr,
     $free_name:ident, $free_doc:expr,
     $matrix:ty, $profile:ty, $trace:literal, $x_drop:literal) => {
        #[doc = $new_doc]
        #[no_mangle]
        pub unsafe extern fn $new_name(query_len: usize,
                                       reference_len: usize,
                                       max_size: usize) -> BlockHandle {
            let aligner = Box::new(Block::<$trace, $x_drop>::new(query_len, reference_len, max_size));
            Box::into_raw(aligner) as BlockHandle
        }

        #[doc = $align_doc]
        #[no_mangle]
        pub unsafe extern fn $align_name(b: BlockHandle,
                                         q: *const PaddedBytes,
                                         r: *const PaddedBytes,
                                         m: *const $matrix,
                                         g: Gaps,
                                         s: SizeRange,
                                         x: i32) {
            let aligner = &mut *(b as *mut Block<$trace, $x_drop>);
            aligner.align(&*q, &*r, &*m, g, s.min..=s.max, x);
        }

        #[doc = $align_profile_doc]
        #[no_mangle]
        pub unsafe extern fn $align_profile_name(b: BlockHandle,
                                                 q: *const PaddedBytes,
                                                 r: *const $profile,
                                                 s: SizeRange,
                                                 x: i32) {
            let aligner = &mut *(b as *mut Block<$trace, $x_drop>);
            aligner.align_profile(&*q, &*r, s.min..=s.max, x);
        }

        #[doc = $res_doc]
        #[no_mangle]
        pub unsafe extern fn $res_name(b: BlockHandle) -> AlignResult {
            let aligner = &*(b as *const Block<$trace, $x_drop>);
            aligner.res()
        }

        #[doc = $trace_doc]
        #[no_mangle]
        pub unsafe extern fn $trace_name(b: BlockHandle, query_idx: usize, reference_idx: usize, cigar: *mut Cigar) {
            let aligner = &*(b as *const Block<$trace, $x_drop>);
            aligner.trace().cigar(query_idx, reference_idx, &mut *cigar);
        }

        #[doc = $trace_eq_doc]
        #[no_mangle]
        pub unsafe extern fn $trace_eq_name(b: BlockHandle, q: *const PaddedBytes, r: *const PaddedBytes, query_idx: usize, reference_idx: usize, cigar: *mut Cigar) {
            let aligner = &*(b as *const Block<$trace, $x_drop>);
            aligner.trace().cigar_eq(&*q, &*r, query_idx, reference_idx, &mut *cigar);
        }

        #[doc = $free_doc]
        #[no_mangle]
        pub unsafe extern fn $free_name(b: BlockHandle) {
            drop(Box::from_raw(b as *mut Block<$trace, $x_drop>));
        }
    };
}

gen_functions!(
    block_new_aa,
    "Create a new block aligner instance for global alignment of amino acid strings (no traceback).",
    block_align_aa,
    "Global alignment of two amino acid strings (no traceback).",
    block_align_profile_aa,
    "Global alignment of an amino acid sequence to a profile (no traceback).",
    block_res_aa,
    "Retrieves the result of global alignment of two amino acid strings (no traceback).",
    _block_cigar_aa,
    "Don't use.",
    _block_cigar_eq_aa,
    "Don't use.",
    block_free_aa,
    "Frees the block used for global alignment of two amino acid strings (no traceback).",
    AAMatrix, AAProfile, false, false
);

gen_functions!(
    block_new_aa_xdrop,
    "Create a new block aligner instance for X-drop alignment of amino acid strings (no traceback).",
    block_align_aa_xdrop,
    "X-drop alignment of two amino acid strings (no traceback).",
    block_align_profile_aa_xdrop,
    "X-drop alignment of an amino acid sequence to a profile (no traceback).",
    block_res_aa_xdrop,
    "Retrieves the result of X-drop alignment of two amino acid strings (no traceback).",
    _block_cigar_aa_xdrop,
    "Don't use.",
    _block_cigar_eq_aa_xdrop,
    "Don't use.",
    block_free_aa_xdrop,
    "Frees the block used for X-drop alignment of two amino acid strings (no traceback).",
    AAMatrix, AAProfile, false, true
);

gen_functions!(
    block_new_aa_trace,
    "Create a new block aligner instance for global alignment of amino acid strings, with traceback.",
    block_align_aa_trace,
    "Global alignment of two amino acid strings, with traceback.",
    block_align_profile_aa_trace,
    "Global alignment of an amino acid sequence to a profile, with traceback.",
    block_res_aa_trace,
    "Retrieves the result of global alignment of two amino acid strings, with traceback.",
    block_cigar_aa_trace,
    "Retrieves the resulting CIGAR string from global alignment of two amino acid strings, with traceback.",
    block_cigar_eq_aa_trace,
    "Retrieves the resulting CIGAR string from global alignment of two amino acid strings, with traceback containing =/X.",
    block_free_aa_trace,
    "Frees the block used for global alignment of two amino acid strings, with traceback.",
    AAMatrix, AAProfile, true, false
);

gen_functions!(
    block_new_aa_trace_xdrop,
    "Create a new block aligner instance for X-drop alignment of amino acid strings, with traceback.",
    block_align_aa_trace_xdrop,
    "X-drop alignment of two amino acid strings, with traceback.",
    block_align_profile_aa_trace_xdrop,
    "X-drop alignment of an amino acid sequence to a profile, with traceback.",
    block_res_aa_trace_xdrop,
    "Retrieves the result of X-drop alignment of two amino acid strings, with traceback.",
    block_cigar_aa_trace_xdrop,
    "Retrieves the resulting CIGAR string from X-drop alignment of two amino acid strings, with traceback.",
    block_cigar_eq_aa_trace_xdrop,
    "Retrieves the resulting CIGAR string from X-drop alignment of two amino acid strings, with traceback containing =/X.",
    block_free_aa_trace_xdrop,
    "Frees the block used for X-drop alignment of two amino acid strings, with traceback.",
    AAMatrix, AAProfile, true, true
);
