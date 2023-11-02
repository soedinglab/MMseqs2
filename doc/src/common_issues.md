# Commonly encountered (Non-Corrosion) Issues

## Linking Debug C/C++ libraries into Rust fails on Windows MSVC targets

`rustc` always links against the non-debug Windows runtime on `*-msvc` targets.
This is tracked [in this issue](https://github.com/rust-lang/rust/issues/39016)
and could be fixed upstream.

A typical error message for this issue is:

```
   Compiling rust_bin v0.1.0 (D:\a\corrosion\corrosion\test\cxxbridge\cxxbridge_cpp2rust\rust)
error: linking with `link.exe` failed: exit code: 1319
[ redacted ]
  = note: cxxbridge-cpp.lib(lib.cpp.obj) : error LNK2038: mismatch detected for '_ITERATOR_DEBUG_LEVEL': value '2' doesn't match value '0' in libcxx-bafec361a1a30317.rlib(cxx.o)

          cxxbridge-cpp.lib(lib.cpp.obj) : error LNK2038: mismatch detected for 'RuntimeLibrary': value 'MDd_DynamicDebug' doesn't match value 'MD_DynamicRelease' in libcxx-bafec361a1a30317.rlib(cxx.o)

          cpp_lib.lib(cpplib.cpp.obj) : error LNK2038: mismatch detected for '_ITERATOR_DEBUG_LEVEL': value '2' doesn't match value '0' in libcxx-bafec361a1a30317.rlib(cxx.o)

          cpp_lib.lib(cpplib.cpp.obj) : error LNK2038: mismatch detected for 'RuntimeLibrary': value 'MDd_DynamicDebug' doesn't match value 'MD_DynamicRelease' in libcxx-bafec361a1a30317.rlib(cxx.o)

          msvcrt.lib(initializers.obj) : warning LNK4098: defaultlib 'msvcrtd.lib' conflicts with use of other libs; use /NODEFAULTLIB:library
```

### Solutions

One solution is to also use the non-debug version when building the C/C++ libraries. 
You can set the [MSVC_RUNTIME_LIBRARY] target properties of your C/C++ libraries to the non-debug variants.
By default you will probably want to select the `MultiThreadedDLL` variant, unless you specified
[`-Ctarget-feature=+crt-static`](https://rust-lang.github.io/rfcs/1721-crt-static.html) in your
`RUSTFLAGS`.


[MSVC_RUNTIME_LIBRARY]: https://cmake.org/cmake/help/latest/prop_tgt/MSVC_RUNTIME_LIBRARY.html#prop_tgt:MSVC_RUNTIME_LIBRARY

## Linking Rust static libraries into Debug C/C++ binaries fails on Windows MSVC targets

This issue is quite similar to the previous one, except that this time it's a Rust library being linked
into a C/C++ target. If it's 100% only Rust code you likely won't even have any issues.
However, if somewhere in the dependency graph C/C++ code is built and linked into your Rust library,
you will likely encounter this issue. Please note, that using [cxx] counts as using C++ code and will
lead to this issue.

The previous solution should also work for this case, but additionally you [may also
have success](https://github.com/rust-lang/rust/issues/39016#issuecomment-853964918) by using 
`corrosion_set_env_vars(your_rust_lib "CFLAGS=-MDd" "CXXFLAGS=-MDd")` (or `-MTd` for a statically linked
runtime).
For debug builds, this is likely to be the preferable solution. It assumes that downstream C/C++ code
is built by the `cc` crate, which respects the `CFLAGS` and `CXXFLAGS` environment variables.

[cxx]: https://github.com/dtolnay/cxx
