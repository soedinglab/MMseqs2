# C API
This directory contains an example of how to use the C API of block aligner.

Currently, only sequence to sequence and sequence to profile alignment with
amino acid scoring matrices is supported with the C API. However, this can easily be adapted
to align nucleotides by setting custom match/mismatch scores.
Other features may be added if there is demand for them.

## Running the example
1. `cd` into this directory.
2. Run `make`. This will build block aligner in release mode, use cbindgen
to generate the header file, and make sure block aligner is linked to the
example program. You only need to have [cbindgen](https://github.com/eqrion/cbindgen) installed
if you are making changes to the block aligner code and need to regenerate bindings. If you are only using
the library, a generated header file is provided in this directory.
3. Run `./a.out`. This will run the example program to perform alignment
calculations.

The generated header file, `c/block_aligner.h`, should be included in
code that calls block aligner functions. It is C++ compatible.
Like in the example `Makefile`, the `block_aligner_c` library in `c/target/release`
must be linked to any C/C++ code that calls block aligner functions.

Note that this directory has a minimal `Cargo.toml` that has no dependencies, so
block aligner can be compiled offline.
