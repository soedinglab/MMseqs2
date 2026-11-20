// Kernel parameters. Buffers: 0 chars, 1 offsets, 2 lengths, 3 pssm, 4 scores, 5 params, 6 id list,
// 7 letter translation table (256 bytes, used with MARV_METAL_FLAG_TRANSLATE), 8 tile carry (multi-tile kernel only).
#ifndef MARV_METAL_SHARED_H
#define MARV_METAL_SHARED_H

#ifndef __METAL_VERSION__
#include <cstdint>
#endif

struct MarvMetalGaplessParams {
    uint32_t numItems;
    uint32_t queryLength;
    uint32_t alphabetSize;
    uint32_t flags;
    uint32_t indexOffset;
    uint32_t pssmStride;
    uint64_t charsBase;
    uint64_t offsetsBase;
    uint64_t lengthsBase;
    uint64_t pssmBase;
    uint64_t outBase;
    uint64_t idsBase;
    // multi-tile kernel: query tiles, persistent SIMD-groups and their carry slices (int32 elements)
    uint32_t numTiles;
    uint32_t numSimdgroups;
    uint64_t carryBase;
    uint64_t carryStride;
};

#define MARV_METAL_FLAG_ID_LIST      2u
#define MARV_METAL_FLAG_TRANSLATE    4u
#define MARV_METAL_FLAG_RETRY        8u
#define MARV_METAL_FLAG_ZERO_OFFSETS 16u

#define MARV_METAL_FP16_EXACT_LIMIT 2048.0f

#endif
