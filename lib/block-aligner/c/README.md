# C API
This directory contains an example of how to use the C API of block aligner.

Currently, only sequence to sequence and sequence to profile alignment with
proteins is supported with the C API. Other features may be added if there
is demand for them.

## Running the example
1. `cd` into this directory.
2. Run `make`. This will build block aligner in release mode, use cbindgen
to generate the header file, and make sure block aligner is linked to the
example program. Make sure you have [cbindgen](https://github.com/eqrion/cbindgen) installed
if you are making changes to the block aligner code and need to regenerate bindings.
3. Run `./example`. This will run the example program to perform alignment
calculations.

The generated header file, `c/block_aligner.h`, should be included in
code that calls block aligner functions. It is C++ compatible.
Like in the example `Makefile`, the `block_aligner` library in `c/target/release`
must be linked to any C/C++ code that calls block aligner functions.

Note that this directory has a minimal `Cargo.toml` that has no dependencies, so
block aligner can be compiled offline.
