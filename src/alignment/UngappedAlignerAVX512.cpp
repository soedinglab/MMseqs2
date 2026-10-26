// Gapless aligner at 512-bit width. CMake compiles this file (and only this file) with
// -mavx512f -mavx512bw -mavx512dq -mavx512vl; simd.h picks the AVX512 path from the resulting
// predefined macros. Entered only when cpuHasAVX512() says the host can run it.
#include "simd.h"

#ifndef AVX512
#error "UngappedAlignerAVX512.cpp must be compiled with AVX512 enabled"
#endif

#include "UngappedAligner.h"

#include "StripedSmithWaterman.h"
#include "Parameters.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "Util.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <utility>

namespace UngappedAVX512 {

#include "UngappedAlignerInternal.cpp"

} // namespace UngappedAVX512
