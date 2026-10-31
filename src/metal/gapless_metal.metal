// Gapless PSSM scan for Apple GPUs. Tiled kernel: GROUP_SIZE lanes per subject, NUM_REGS query columns per lane in
// permuted vec4 groups, PSSM cached in threadgroup memory (ScoreT or int8). Untiled kernel: one subject per SIMD-group.
#include <metal_stdlib>
#include "marv_metal_shared.h"
using namespace metal;

// short saturates so that the shared overflow correction sees every out-of-range score
inline half4  accum(half4 a, half4 b)   { return a + b; }
inline short4 accum(short4 a, short4 b) { return addsat(a, b); }
inline int4   accum(int4 a, int4 b)     { return a + b; }

template<uint GROUP_SIZE, uint NUM_REGS, class ScoreT, class TgT, bool MAY_RETRY>
kernel void gapless_tiled(
    device const char*      dbChars      [[buffer(0)]],
    device const uint64_t*  dbOffsets    [[buffer(1)]],
    device const int32_t*   dbLengths    [[buffer(2)]],
    device const int8_t*    pssm         [[buffer(3)]],
    device int32_t*         outScores    [[buffer(4)]],
    constant MarvMetalGaplessParams& params [[buffer(5)]],
    device const int32_t*   subjectIds   [[buffer(6)]],
    constant const uchar*   letterTable  [[buffer(7)]],
    threadgroup TgT*        tg_pssm      [[threadgroup(0)]],
    uint simdgroup_id        [[simdgroup_index_in_threadgroup]],
    uint thread_in_simdgroup [[thread_index_in_simdgroup]],
    uint threadgroup_id      [[threadgroup_position_in_grid]],
    uint simdgroups_per_tg   [[simdgroups_per_threadgroup]],
    uint thread_in_tg        [[thread_position_in_threadgroup]],
    uint threads_per_tg      [[threads_per_threadgroup]]
) {
    using Vec = vec<ScoreT, 4>;
    using TgVec = vec<TgT, 4>;
    constexpr uint QUARTER         = NUM_REGS / 4;
    constexpr uint PADDED_COLS     = GROUP_SIZE * NUM_REGS;
    constexpr uint SUBJECTS_PER_SG = 32u / GROUP_SIZE;

    const uint queryLength  = params.queryLength;
    const uint alphabetSize = params.alphabetSize;
    device const int8_t* pssmData = pssm + params.pssmBase;

    auto cachePssm = [&] {
        for (uint i = thread_in_tg; i < alphabetSize * PADDED_COLS; i += threads_per_tg) {
            uint letter = i / PADDED_COLS;
            uint slot   = i % PADDED_COLS;
            uint lane   = slot / NUM_REGS;
            uint r      = slot % NUM_REGS;
            uint m      = r / 4, k = r % 4;
            uint col    = lane * NUM_REGS + m + k * QUARTER;
            tg_pssm[i] = (col < queryLength) ? (TgT) pssmData[letter * params.pssmStride + col] : (TgT) 0;
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    };
    if (!MAY_RETRY) cachePssm();

    const uint lane      = thread_in_simdgroup;
    const uint groupLane = lane % GROUP_SIZE;
    const uint groupInSg = lane / GROUP_SIZE;
    const uint globalSg  = threadgroup_id * simdgroups_per_tg + simdgroup_id;
    const uint listIdx   = globalSg * SUBJECTS_PER_SG + groupInSg;
    bool active = (listIdx < params.numItems);
    threadgroup uint activeSimdgroups[32];
    if (MAY_RETRY) {
        if (params.flags & MARV_METAL_FLAG_RETRY) {
            active = active && outScores[params.outBase + listIdx] >= 2048;
            bool anyInSimdgroup = simd_any(active);
            if (lane == 0) activeSimdgroups[simdgroup_id] = anyInSimdgroup;
            threadgroup_barrier(mem_flags::mem_threadgroup);
            bool anyInThreadgroup = false;
            for (uint sg = 0; sg < simdgroups_per_tg; ++sg) anyInThreadgroup |= activeSimdgroups[sg] != 0;
            // uniform exit before the PSSM is staged
            if (!anyInThreadgroup) return;
        }
    }

    if (MAY_RETRY) cachePssm();

    uint seqIdx = 0;
    if (active) seqIdx = (params.flags & MARV_METAL_FLAG_ID_LIST) ? (uint) subjectIds[params.idsBase + listIdx] : params.indexOffset + listIdx;
    const int32_t subjectLength = active ? dbLengths[params.lengthsBase + seqIdx] : 0;
    const uint64_t baseOffset   = params.charsBase + (active ? (dbOffsets[params.offsetsBase + seqIdx] - ((params.flags & MARV_METAL_FLAG_ZERO_OFFSETS) ? dbOffsets[params.offsetsBase] : 0u)) : 0);
    const uint myStart = groupLane * NUM_REGS;
    const int32_t maxLen = simd_max(subjectLength);

    Vec penalty_here[QUARTER];
    #pragma unroll
    for (uint r = 0; r < QUARTER; r++) penalty_here[r] = Vec(0);
    ScoreT prev_border = 0;
    Vec maximum4 = Vec(0);

    // one subject residue; lanes whose subject ended skip the recurrence but still shuffle
    auto step = [&](int subjectLetter, bool valid) {
      if (valid) {
        threadgroup const TgT* pssmRow = tg_pssm + subjectLetter * PADDED_COLS;
        threadgroup const TgVec* row4 = (threadgroup const TgVec*)(pssmRow + myStart);
        Vec old_last = penalty_here[QUARTER - 1];
        Vec old_prev = penalty_here[0];
        Vec scores0 = Vec(row4[0]);
        Vec diag0 = Vec(prev_border, old_last.x, old_last.y, old_last.z);
        Vec new0 = max(accum(diag0, scores0), Vec(0));
        penalty_here[0] = new0;
        maximum4 = max(maximum4, new0);
        #pragma unroll
        for (uint m = 1; m < QUARTER; m++) {
            Vec old_cur = penalty_here[m];
            Vec new_val = max(accum(old_prev, Vec(row4[m])), Vec(0));
            penalty_here[m] = new_val;
            maximum4 = max(maximum4, new_val);
            old_prev = old_cur;
        }
      }
        if (GROUP_SIZE > 1) {
            ScoreT shuffled = simd_shuffle_up(penalty_here[QUARTER - 1].w, 1);
            prev_border = (groupLane == 0) ? (ScoreT) 0 : shuffled;
        }
    };
    const int clampLetter = (int) alphabetSize - 1;
    const bool translate = (params.flags & MARV_METAL_FLAG_TRANSLATE) != 0;
    auto code = [&](int raw) { return min(translate ? (int) letterTable[raw] : raw, clampLetter); };

    // every lane runs to the longest subject of the SIMD-group so that the shuffles stay uniform
    // packed pairs have alignment 1; the common prefix needs no masks
    const int32_t evenLen = subjectLength & ~1;
    device const packed_uchar2* two = (device const packed_uchar2*)(dbChars + baseOffset);
    int32_t j = 0;
    const int32_t commonEvenLen = simd_min(subjectLength) & ~1;
    for (; j < commonEvenLen; j += 2) {
        packed_uchar2 c = two[j >> 1];
        step(code(c.x), true);
        step(code(c.y), true);
    }
    for (; j + 1 < maxLen; j += 2) {
        int a = clampLetter, b = clampLetter;
        const bool va = j < subjectLength, vb = (j + 1) < subjectLength;
        if (j + 1 < evenLen) { packed_uchar2 c = two[j >> 1]; a = c.x; b = c.y; }
        else if (va) { a = (uchar) dbChars[baseOffset + j]; }
        step(code(a), va);
        step(code(b), vb);
    }
    if (j < maxLen) {
        const bool va = j < subjectLength;
        int a = va ? (int)(uchar) dbChars[baseOffset + j] : clampLetter;
        step(code(a), va);
    }

    ScoreT groupMax = max(max(maximum4.x, maximum4.y), max(maximum4.z, maximum4.w));
    #pragma unroll
    for (uint off = GROUP_SIZE / 2; off >= 1; off >>= 1) {
        groupMax = max(groupMax, simd_shuffle_xor(groupMax, off));
    }
    if (groupLane == 0 && active) {
        // integer sentinel: a half score can reach infinity and must still be retried
        outScores[params.outBase + listIdx] = MAY_RETRY ? (int32_t) groupMax
            : (int32_t) min(groupMax, (ScoreT) MARV_METAL_FP16_EXACT_LIMIT);
    }
}

template<uint REGS, class ScoreT>
kernel void gapless_untiled(
    device const char*      dbChars      [[buffer(0)]],
    device const uint64_t*  dbOffsets    [[buffer(1)]],
    device const int32_t*   dbLengths    [[buffer(2)]],
    device const int8_t*    pssm         [[buffer(3)]],
    device int32_t*         outScores    [[buffer(4)]],
    constant MarvMetalGaplessParams& params [[buffer(5)]],
    device const int32_t*   subjectIds   [[buffer(6)]],
    constant const uchar*   letterTable  [[buffer(7)]],
    uint simdgroup_id        [[simdgroup_index_in_threadgroup]],
    uint thread_in_simdgroup [[thread_index_in_simdgroup]],
    uint threadgroup_id      [[threadgroup_position_in_grid]],
    uint simdgroups_per_tg   [[simdgroups_per_threadgroup]]
) {
    using Vec = vec<ScoreT, 4>;
    const uint listIdx = threadgroup_id * simdgroups_per_tg + simdgroup_id;
    if (listIdx >= params.numItems) return;
    const uint seqIdx = (params.flags & MARV_METAL_FLAG_ID_LIST) ? (uint) subjectIds[params.idsBase + listIdx] : params.indexOffset + listIdx;

    const uint lane = thread_in_simdgroup;
    const int32_t subjectLength = dbLengths[params.lengthsBase + seqIdx];
    const uint64_t baseOffset = params.charsBase + (dbOffsets[params.offsetsBase + seqIdx] - ((params.flags & MARV_METAL_FLAG_ZERO_OFFSETS) ? dbOffsets[params.offsetsBase] : 0u));
    const int clampLetter = (int) params.alphabetSize - 1;
    const bool translate = (params.flags & MARV_METAL_FLAG_TRANSLATE) != 0;
    device const int8_t* pssmData = pssm + params.pssmBase;

    constexpr uint numRegsPerThread = 4 * REGS;
    constexpr uint quarter = REGS;
    const uint myStart = lane * numRegsPerThread;

    Vec penalty_here[REGS];
    ScoreT prev_border = 0;
    Vec maximum4 = Vec(0);
    #pragma unroll
    for (uint r = 0; r < quarter; r++) penalty_here[r] = Vec(0);

    for (int32_t j = 0; j < subjectLength; j++) {
        int raw = (int)(uchar) dbChars[baseOffset + j];
        int subjectLetter = min(translate ? (int) letterTable[raw] : raw, clampLetter);
        device const int8_t* pssmRow = pssmData + subjectLetter * params.pssmStride;

        Vec old_last = penalty_here[quarter - 1];
        Vec old_prev = penalty_here[0];
        uint p0 = myStart;
        Vec scores0;
        scores0.x = (ScoreT) pssmRow[p0];
        scores0.y = (ScoreT) pssmRow[p0 + quarter];
        scores0.z = (ScoreT) pssmRow[p0 + 2 * quarter];
        scores0.w = (ScoreT) pssmRow[p0 + 3 * quarter];
        Vec diag0 = Vec(prev_border, old_last.x, old_last.y, old_last.z);
        Vec new0 = max(accum(diag0, scores0), Vec(0));
        penalty_here[0] = new0;
        maximum4 = max(maximum4, new0);
        #pragma unroll
        for (uint m = 1; m < quarter; m++) {
            Vec old_cur = penalty_here[m];
            uint pm = myStart + m;
            Vec scores;
            scores.x = (ScoreT) pssmRow[pm];
            scores.y = (ScoreT) pssmRow[pm + quarter];
            scores.z = (ScoreT) pssmRow[pm + 2 * quarter];
            scores.w = (ScoreT) pssmRow[pm + 3 * quarter];
            Vec new_val = max(accum(old_prev, scores), Vec(0));
            penalty_here[m] = new_val;
            maximum4 = max(maximum4, new_val);
            old_prev = old_cur;
        }
        ScoreT shuffled = simd_shuffle_up(penalty_here[quarter - 1].w, 1);
        prev_border = (lane == 0) ? (ScoreT) 0 : shuffled;
    }
    ScoreT groupMax = max(max(maximum4.x, maximum4.y), max(maximum4.z, maximum4.w));
    groupMax = simd_max(groupMax);
    if (lane == 0) {
        outScores[params.outBase + listIdx] = (int32_t) groupMax;
    }
}

#define INSTANTIATE_TILED_T(NAME, G, N, S, T, R) \
    template [[host_name(NAME)]] kernel void gapless_tiled<G, N, S, T, R>( \
        device const char*, device const uint64_t*, device const int32_t*, device const int8_t*, device int32_t*, \
        constant MarvMetalGaplessParams&, device const int32_t*, constant const uchar*, threadgroup T*, \
        uint, uint, uint, uint, uint, uint);
#define INSTANTIATE_TILED(G, N) \
    INSTANTIATE_TILED_T("gl_tg_h_g" #G "_nr" #N, G, N, half, half, false) \
    INSTANTIATE_TILED_T("gl_tg_s_g" #G "_nr" #N, G, N, short, short, true)
#define INSTANTIATE_TILED8(G, N) \
    INSTANTIATE_TILED_T("gl_tg_h8_g" #G "_nr" #N, G, N, half, int8_t, false) \
    INSTANTIATE_TILED_T("gl_tg_s8_g" #G "_nr" #N, G, N, short, int8_t, true)

#define TILED_4_48(G) \
    INSTANTIATE_TILED(G,4)  INSTANTIATE_TILED(G,8)  INSTANTIATE_TILED(G,12) INSTANTIATE_TILED(G,16) \
    INSTANTIATE_TILED(G,20) INSTANTIATE_TILED(G,24) INSTANTIATE_TILED(G,28) INSTANTIATE_TILED(G,32) \
    INSTANTIATE_TILED(G,36) INSTANTIATE_TILED(G,40) INSTANTIATE_TILED(G,44) INSTANTIATE_TILED(G,48)
#define TILED_52_64(G) \
    INSTANTIATE_TILED(G,52) INSTANTIATE_TILED(G,56) INSTANTIATE_TILED(G,60) INSTANTIATE_TILED(G,64)
#define TILED_68_96(G) \
    INSTANTIATE_TILED(G,68) INSTANTIATE_TILED(G,72) INSTANTIATE_TILED(G,76) INSTANTIATE_TILED(G,80) \
    INSTANTIATE_TILED(G,84) INSTANTIATE_TILED(G,88) INSTANTIATE_TILED(G,92) INSTANTIATE_TILED(G,96)
#define TILED8_24_48(G) \
    INSTANTIATE_TILED8(G,24) INSTANTIATE_TILED8(G,28) INSTANTIATE_TILED8(G,32) INSTANTIATE_TILED8(G,36) \
    INSTANTIATE_TILED8(G,40) INSTANTIATE_TILED8(G,44) INSTANTIATE_TILED8(G,48)
#define TILED8_52_96(G) \
    INSTANTIATE_TILED8(G,52) INSTANTIATE_TILED8(G,56) INSTANTIATE_TILED8(G,60) INSTANTIATE_TILED8(G,64) \
    INSTANTIATE_TILED8(G,68) INSTANTIATE_TILED8(G,72) INSTANTIATE_TILED8(G,76) INSTANTIATE_TILED8(G,80) \
    INSTANTIATE_TILED8(G,84) INSTANTIATE_TILED8(G,88) INSTANTIATE_TILED8(G,92) INSTANTIATE_TILED8(G,96)

TILED_4_48(4)
TILED_4_48(8)  TILED_52_64(8)
TILED_4_48(16) TILED_52_64(16) TILED_68_96(16)
TILED_4_48(32)
INSTANTIATE_TILED8(16,48) TILED8_52_96(16)
TILED8_24_48(32)

#define UNTILED_ARGS \
    device const char*, device const uint64_t*, device const int32_t*, device const int8_t*, device int32_t*, \
    constant MarvMetalGaplessParams&, device const int32_t*, constant const uchar*, uint, uint, uint, uint
#define UNTILED(R) \
    template [[host_name("gl_untiled_s_r" #R)]] kernel void gapless_untiled<R, short>(UNTILED_ARGS); \
    template [[host_name("gl_untiled_i_r" #R)]] kernel void gapless_untiled<R, int>(UNTILED_ARGS);
UNTILED(1) UNTILED(2) UNTILED(3) UNTILED(4) UNTILED(5) UNTILED(6) UNTILED(7) UNTILED(8)
UNTILED(9) UNTILED(10) UNTILED(11) UNTILED(12) UNTILED(13) UNTILED(14) UNTILED(15) UNTILED(16)
UNTILED(20) UNTILED(24) UNTILED(32) UNTILED(40) UNTILED(48) UNTILED(64)
