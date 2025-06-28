CC=/usr/local/opt/llvm/bin/clang cargo run --example nanopore_bench_global --release --features simd_avx2 -- "$@"
