CC=clang
CFLAGS=-Ltarget/release -lblock_aligner_c -I.

.PHONY: all block_aligner example

all: example

block_aligner:
	cargo build --release --features simd_avx2 --offline
	cbindgen --config cbindgen.toml --crate block-aligner-c --output block_aligner.h --quiet .

example: block_aligner example.c
	$(CC) $(CFLAGS) example.c
