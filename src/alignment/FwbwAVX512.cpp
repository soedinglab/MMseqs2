// Forward-backward aligner for AVX512 hosts. CMake compiles this file (and only this file) with
// -mavx512f -mavx512bw -mavx512dq -mavx512vl.
//
// MMSEQS_FORCE_SIMD256 keeps the kernel at 256-bit width while still using AVX512VL mask
// features: fwbw gains more from the halved segment count and the cheaper cross-lane shifts than
// it loses from the narrower vectors, and the k-register compares are not available in plain
// AVX2. Build with -DMMSEQS_NO_SIMD256 to A/B against the native 512-bit width.
#ifndef MMSEQS_NO_SIMD256
#define MMSEQS_FORCE_SIMD256 1
#endif
#include "simd.h"

#if !defined(MMSEQS_AVX512VL_MASKS) && !defined(AVX512)
#error "FwbwAVX512.cpp must be compiled with AVX512 enabled"
#endif

#include "Fwbw.h"

#include "StripedSmithWaterman.h"   // SmithWaterman::computeCov
#include "SubstitutionMatrix.h"
#include "Util.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <limits>

namespace FwBwAVX512 {

#include "FwbwInternal.cpp"

} // namespace FwBwAVX512
