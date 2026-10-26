// Gapless aligner at the instruction set the rest of the binary was compiled for.
#include "simd.h"
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

namespace UngappedDefault {

#include "UngappedAlignerInternal.cpp"

} // namespace UngappedDefault
