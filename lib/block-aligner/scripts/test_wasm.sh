CARGO_TARGET_WASM32_WASI_RUNNER="wasmtime --wasm-features simd --" cargo test --target=wasm32-wasi --all-targets --features simd_wasm -- --nocapture "$@"
