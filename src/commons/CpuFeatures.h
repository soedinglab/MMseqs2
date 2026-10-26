#ifndef MMSEQS_CPUFEATURES_H
#define MMSEQS_CPUFEATURES_H

// Runtime CPU feature detection for the dynamically dispatched SIMD kernels.
//
// This header must stay free of simd.h (and of any intrinsics header): it is included by the
// dispatch headers, which in turn are included by translation units compiled with different
// -m flags. Everything here is plain preprocessor + builtins.
//
// MMSEQS_HAVE_AVX512_DISPATCH is defined by CMake when the toolchain was able to compile the
// AVX512 translation units. Without it no AVX512 backend exists to dispatch to, so the check
// is compiled out and the default backend is always used.

#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
#define MMSEQS_ARCH_X86 1
#endif

#if defined(MMSEQS_ARCH_X86) && (defined(__GNUC__) || defined(__clang__))
#if defined(__has_builtin)
#if __has_builtin(__builtin_cpu_supports)
#define MMSEQS_HAVE_CPU_SUPPORTS 1
#endif
#elif defined(__GNUC__) && (__GNUC__ >= 5)
// GCC has __builtin_cpu_supports since 4.8, but the avx512* names only since 5.
#define MMSEQS_HAVE_CPU_SUPPORTS 1
#endif
#endif

// True when the CPU can execute the AVX512 backends. AVX512VL is part of the requirement even
// though the gapless kernel only needs F/BW/DQ: the fwbw backend runs at 256-bit width using
// AVX512VL masks, and requiring it for both keeps one detection result for all backends. In
// practice this only excludes Knights Landing, which has no VL.
inline bool cpuHasAVX512() {
#if defined(MMSEQS_HAVE_AVX512_DISPATCH) && defined(MMSEQS_HAVE_CPU_SUPPORTS)
#if defined(__has_builtin)
#if __has_builtin(__builtin_cpu_init)
    __builtin_cpu_init();
#endif
#else
    __builtin_cpu_init();
#endif
    return __builtin_cpu_supports("avx512f")
        && __builtin_cpu_supports("avx512bw")
        && __builtin_cpu_supports("avx512dq")
        && __builtin_cpu_supports("avx512vl");
#else
    return false;
#endif
}

#endif //MMSEQS_CPUFEATURES_H
