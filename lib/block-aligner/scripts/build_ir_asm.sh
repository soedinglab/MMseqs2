set -e

LLVM_MCA=/usr/local/opt/llvm/bin/llvm-mca

#RUSTFLAGS="-g -Z asm-comments --emit llvm-ir,asm -C llvm-args=-x86-asm-syntax=intel -C target-cpu=native" cargo build --release --example profile --features mca
RUSTFLAGS="-Z asm-comments --emit llvm-ir,asm -C llvm-args=-x86-asm-syntax=intel" cargo build --release --example profile --features mca,simd_avx2

# demangle symbols
#for f in target/release/examples/*.{s,ll}; do
#    rustfilt -i $f > $f.filt
#    echo "$f.filt"
#done

for f in target/release/examples/*.{s,ll}; do
    echo "$f"
done

# run llvm-mca
for f in target/release/examples/*.s; do
    $LLVM_MCA -output-asm-variant=1 -all-views $f > $f.mca
    echo "$f.mca"
done

# also create source/asm interleaved version with objdump
#shopt -s extglob
#for f in target/*/deps/!(*.*); do
#    objdump -drwSl -x86-asm-syntax=intel $f | rustfilt -o $f.objdump
#    echo "$f.objdump"
#done
