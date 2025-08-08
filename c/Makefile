CC=clang
CFLAGS=-Ltarget/release -lblock_aligner_c -I.

.PHONY: all block_aligner example align_prefix align_local

all: example align_prefix align_local

block_aligner:
	cargo build --release --features simd_avx2 --offline
	cbindgen --config cbindgen.toml --crate block_aligner_c --output block_aligner.h --quiet .

example: block_aligner example.c
	$(CC) $(CFLAGS) -o example example.c

align_prefix: block_aligner align_prefix.c
	$(CC) $(CFLAGS) -o align_prefix align_prefix.c

align_local: block_aligner align_local.c
	$(CC) $(CFLAGS) -o align_local align_local.c
