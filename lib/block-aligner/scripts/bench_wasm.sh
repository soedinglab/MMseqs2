#set -e

#cargo clean

# CARGO_TARGET_WASM32_WASI_RUNNER="wasmer run --native --llvm --enable-simd --"
# CARGO_TARGET_WASM32_WASI_RUNNER="wavm run --enable simd"

CARGO_TARGET_WASM32_WASI_RUNNER="wasmtime --wasm-features simd --" cargo bench --target=wasm32-wasi --features simd_wasm -- --nocapture "$@"
#RUSTFLAGS="-C target-feature=+simd128" cargo build --release --benches --target=wasm32-wasi

# binaryen wasm-opt pass
#for f in target/wasm32-wasi/*/deps/*.wasm; do
#    wasm-opt --enable-simd --enable-sign-ext -O4 --inlining-optimizing -ifwl -ocimfs 300 -fimfs 300 -aimfs 20 -o $f.opt $f
#    echo $f.opt
#done

#for f in target/wasm32-wasi/*/deps/*.wasm.opt; do
#    $CARGO_TARGET_WASM32_WASI_RUNNER $f --bench "$@"
#done
