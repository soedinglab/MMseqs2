#![feature(test)]
#![cfg(feature = "simd_avx2")]

extern crate test;
use test::{Bencher, black_box};

use block_aligner::avx2::*;

#[repr(align(32))]
struct A([i16; L]);

#[bench]
fn bench_opt_prefix_scan(b: &mut Bencher) {
    #[target_feature(enable = "avx2")]
    unsafe fn inner(b: &mut Bencher) {
        let vec = A([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 12, 13, 14, 11]);
        let vec = simd_load(vec.0.as_ptr() as *const Simd);

        b.iter(|| {
            let gap = simd_set1_i16(-1);
            let (_, consts) = get_prefix_scan_consts(gap);
            simd_prefix_scan_i16(black_box(vec), gap, consts)
        });
    }
    unsafe { inner(b); }
}

#[bench]
fn bench_naive_prefix_scan(b: &mut Bencher) {
    #[target_feature(enable = "avx2")]
    unsafe fn inner(b: &mut Bencher) {
        let vec = A([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 12, 13, 14, 11]);
        let vec = simd_load(vec.0.as_ptr() as *const Simd);

        b.iter(|| {
            let gap = simd_set1_i16(-1);
            let (_, consts) = get_prefix_scan_consts(gap);
            simd_naive_prefix_scan_i16(black_box(vec), gap, consts)
        });
    }
    unsafe { inner(b); }
}
