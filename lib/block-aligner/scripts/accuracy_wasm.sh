CARGO_TARGET_WASM32_WASI_RUNNER="wasmtime --wasm-features simd --" cargo run --target=wasm32-wasi --example accuracy --release --features simd_wasm -- "$@"
