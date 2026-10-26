// Forward-backward aligner at the instruction set the rest of the binary was compiled for.
#include "simd.h"
#include "Fwbw.h"

#include "StripedSmithWaterman.h"   // SmithWaterman::computeCov
#include "SubstitutionMatrix.h"
#include "Util.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <limits>

namespace FwBwDefault {

#include "FwbwInternal.cpp"

} // namespace FwBwDefault
