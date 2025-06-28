//! SIMD-accelerated library for computing global and X-drop affine
//! gap penalty sequence-to-sequence or sequence-to-profile alignments
//! using an adaptive block-based algorithm.
//!
//! Currently, SSE2, AVX2, Neon, and WASM SIMD are supported.
//!
//! ## Example
//! ```
//! use block_aligner::{cigar::*, scan_block::*, scores::*};
//!
//! let min_block_size = 32;
//! let max_block_size = 256;
//!
//! // A gap of length n will cost: open + extend * (n - 1)
//! let gaps = Gaps { open: -2, extend: -1 };
//!
//! // Note that PaddedBytes, Block, and Cigar can be initialized with sequence length
//! // and block size upper bounds and be reused later for shorter sequences, to avoid
//! // repeated allocations.
//! let r = PaddedBytes::from_bytes::<NucMatrix>(b"TTAAAAAAATTTTTTTTTTTT", max_block_size);
//! let q = PaddedBytes::from_bytes::<NucMatrix>(b"TTTTTTTTAAAAAAATTTTTTTTT", max_block_size);
//!
//! // Align with traceback, but no X-drop threshold (global alignment).
//! let mut a = Block::<true, false>::new(q.len(), r.len(), max_block_size);
//! a.align(&q, &r, &NW1, gaps, min_block_size..=max_block_size, 0);
//! let res = a.res();
//!
//! assert_eq!(res, AlignResult { score: 7, query_idx: 24, reference_idx: 21 });
//!
//! let mut cigar = Cigar::new(res.query_idx, res.reference_idx);
//! // Compute traceback and resolve =/X (matches/mismatches).
//! a.trace().cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar);
//!
//! assert_eq!(cigar.to_string(), "2=6I16=3D");
//! ```
//!
//! ## Tuning block sizes
//!
//! For long, noisy Nanopore reads, a min block size of ~1% sequence length and a max block size
//! of ~10% sequence length performs well (tested with reads up to ~50kbps).
//! For proteins, a min block size of 32 and a max block size of 256 performs well.
//! Using a minimum block size that is at least 32 is recommended for most applications.
//! Using a maximum block size greater than `2^14 = 16384` is not recommended.
//! If the alignment scores are saturating (score too large), then use a smaller block size.
//! Let me know how block aligner performs on your data!
//!
//! When building your code that uses this library, it is important to specify the
//! correct feature flags: `simd_sse2`, `simd_avx2`, `simd_neon`, or `simd_wasm`.
//! More information on specifying different features for different platforms
//! with the same dependency [here](https://doc.rust-lang.org/cargo/reference/specifying-dependencies.html#platform-specific-dependencies).

// special SIMD instruction set modules adapted for this library
// their types and lengths are abstracted out

#[cfg(feature = "simd_sse2")]
#[macro_use]
#[doc(hidden)]
/// cbindgen:ignore
pub mod sse2;

#[cfg(feature = "simd_sse2")]
pub use sse2::L;

#[cfg(feature = "simd_avx2")]
#[macro_use]
#[doc(hidden)]
/// cbindgen:ignore
pub mod avx2;

#[cfg(feature = "simd_avx2")]
pub use avx2::L;

#[cfg(feature = "simd_wasm")]
#[macro_use]
#[doc(hidden)]
/// cbindgen:ignore
pub mod simd128;

#[cfg(feature = "simd_wasm")]
pub use simd128::L;

#[cfg(feature = "simd_neon")]
#[macro_use]
#[doc(hidden)]
/// cbindgen:ignore
pub mod neon;

#[cfg(feature = "simd_neon")]
pub use neon::L;

#[cfg(any(feature = "simd_sse2", feature = "simd_avx2", feature = "simd_wasm", feature = "simd_neon"))]
pub mod scan_block;
#[cfg(any(feature = "simd_sse2", feature = "simd_avx2", feature = "simd_wasm", feature = "simd_neon"))]
pub mod scores;
#[cfg(any(feature = "simd_sse2", feature = "simd_avx2", feature = "simd_wasm", feature = "simd_neon"))]
pub mod cigar;

#[cfg(any(feature = "simd_sse2", feature = "simd_avx2", feature = "simd_wasm", feature = "simd_neon"))]
#[doc(hidden)]
pub mod ffi;

#[cfg(not(any(feature = "no_simd", feature = "simd_sse2", feature = "simd_avx2", feature = "simd_wasm", feature = "simd_neon")))]
compile_error!("No SIMD feature flag specified! Specify \"no_simd\" to disable all SIMD features.");

/// Calculate the percentage of a length, rounded to the next power of two.
///
/// This is useful for computing the min and max block sizes for sequences of a certain
/// length by using percentages. The returned value is at least 32 and at most 16384.
pub fn percent_len(len: usize, p: f32) -> usize {
    ((p * (len as f32)).round() as usize).max(32).next_power_of_two().min(1 << 14)
}
