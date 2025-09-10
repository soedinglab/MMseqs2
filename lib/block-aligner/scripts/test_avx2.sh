cargo test --all-targets --features simd_avx2 -- "$@"
cargo test --doc --features simd_avx2 -- "$@"
