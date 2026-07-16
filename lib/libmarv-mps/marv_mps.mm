// Metal/MPS implementation of the Marv GPU-prefilter interface declared in
// lib/libmarv/src/marv.h. This is the Apple Silicon counterpart to
// lib/libmarv (CUDA/HIP): it implements the same "gapless" (ungapped)
// local-alignment scan used to accelerate MMseqs2's prefilter stage, so it
// slots into the existing GpuUtil.cpp / gpuserver.cpp / ungappedprefilter.cpp
// pipeline without changes to their control flow.
//
// Scope: only AlignmentType::GAPLESS is implemented (the default prefilter
// path). The gapped Smith-Waterman kernel in libmarv (used by
// --prefilter-mode 2) has no Metal equivalent yet.
//
// Algorithm: H[i][j] = max(0, H[i-1][j-1] + score(query[i], subject[j])),
// boundary zero -- the same no-gap recurrence as libmarv's GAPLESS kernel
// (pssmkernels_gapless.cuh). Subject bytes are folded with `& 31` to
// collapse MMseqs2's soft-masked (lowercase, +32) residue encoding back to
// the base alphabet, mirroring libmarv's CaseSensitive_to_CaseInsensitive
// conversion (see convert.cuh).
//
// Parallelism (v2, simdgroup-cooperative): a first version of this kernel
// used one GPU thread per subject sequence, with its O(queryLen) DP state in
// `device` (global) memory. That was measured to be ~12x slower than
// MMseqs2's CPU prefilter on a real, length-skewed database (lengths 9 to
// 8083 residues, median 347): a handful of very long sequences dominated
// each dispatch's wall-clock while thousands of short-sequence threads sat
// idle, and every DP step round-tripped through slow global memory.
//
// This version instead mirrors libmarv's actual design: one Metal simdgroup
// (32 lanes, lockstep -- Apple's equivalent of a CUDA warp) cooperates on a
// single alignment. Each lane owns a contiguous slice of the query in fast
// `threadgroup` memory; per subject residue, one simd_shuffle_up carries the
// single boundary value from lane L's slice to lane L+1's (the only
// cross-lane dependency the recurrence has -- see GaplessPSSMState::relax /
// shuffleScores in pssmkernels_gapless.cuh for the CUDA original), then all
// 32 lanes update their slice in parallel. Unlike libmarv, the slice length
// is a runtime loop over dynamically-sized threadgroup memory rather than a
// fixed-size register array selected from a family of compiled kernel
// variants -- this trades some peak throughput for supporting arbitrary
// query lengths with a single kernel. A further simdgroup dispatched per
// database entry would still leave the same long-tail scheduling problem,
// so simdgroups instead pull work from a shared atomic counter until the
// database is exhausted: idle simdgroups immediately pick up the next
// (typically short) sequence instead of waiting on whichever simdgroup drew
// the longest one.
//
// PSSM layout (v3, coalesced): with lanes owning *contiguous* query slices,
// reading pssm[s][chunkStart+li] at loop step `li` means lane L reads
// address L*chunkSize+li -- lanes are chunkSize apart, not adjacent, so the
// 32 lanes' simultaneous reads hit scattered cache lines instead of one
// coalesced burst. submitScan() uploads each residue row pre-transposed
// into chunkSize*32 "column-major" blocks (row[li*32+lane] instead of
// row[lane*chunkSize+li]), matching libmarv's own smem layout trick (see
// SmemIndexCalculator in pssmkernels_gapless.cuh) -- the 32 lanes now read
// 32 consecutive bytes each step. The DP state in threadgroup memory (Hall)
// stays untransposed since it's private per simdgroup, not shared. Measured
// impact was modest (0-20% depending on query length, <1% in aggregate on a
// real query mix) -- Apple Silicon's unified GPU memory system apparently
// already absorbs the scattered pattern reasonably well via its cache.
//
// Register tiling -- tried, reverted: a v4 attempt moved the DP state out
// of threadgroup memory into a fixed-size per-lane register array (32
// int32 registers/lane, covering queries up to 1024 residues via 4 tiles of
// 256, falling back to the threadgroup-memory kernel above for longer
// queries), the direct Metal analogue of libmarv's register-resident,
// border_in/border_out tile-chained kernels. It was verified correct --
// byte-identical output against the threadgroup-memory kernel across 25
// boundary-length test cases (1, 2, 31-33, 63-65, ..., 1023-1025, 2000) and
// the full 500-query/20k-sequence real benchmark -- but measured ~27%
// *slower* on the queries it targeted (8.36s vs an estimated 6.59s for the
// same 457 queries under the threadgroup-memory kernel).
//
// Checked with Instruments (Metal System Trace, via `xcrun xctrace record
// --template "Metal System Trace"`) rather than guessing: the
// graphics-compiler-spill-events table showed zero spill events for either
// kernel, ruling out register spilling. metal-shader-profiler-shader-list
// showed the actual difference -- gaplessAlignRegisterTiled compiled to
// 2472 bytes of GPU code vs. gaplessAlignCooperative's 1006, 2.46x more,
// almost certainly from the compiler fully unrolling the compile-time-sized
// tile*register loops (4 tiles x 8 registers) plus the per-register
// `globalPos < qlen` bounds check and the extra cross-tile border shuffle
// bookkeeping.
//
// Follow-up: is it the unrolling (code size) or the bounds-check/shuffle
// overhead that dominates? Tested by making the inner register loop bound a
// runtime Params field instead of the compile-time constant kRegs, so the
// compiler could no longer statically unroll it. Code size did drop, to
// 1484 bytes (vs. 2472 unrolled, 1006 cooperative) -- but wall-clock got
// *worse*, not better: 15.3s vs. 11.7s (unrolled) vs. 9.9s (cooperative)
// on the same 500-query benchmark. graphics-compiler-spill-events explained
// why: this version spilled 144 bytes/thread on essentially every dispatch,
// where the unrolled version spilled nothing. GPU registers can't be
// dynamically indexed -- a loop bound the compiler can't resolve at compile
// time forces it to either unroll anyway (code bloat) or address the array
// through memory (spilling), and spilling turned out to cost more than the
// code bloat it replaced. So it's not really "unrolling vs. bounds-check
// overhead" -- both register-tiled variants lose to the threadgroup-memory
// kernel, for complementary reasons, because there's no way to keep a
// variably-sized per-lane array genuinely register-resident on this
// hardware. The threadgroup-memory kernel wins by making the memory access
// explicit and coalesced instead of fighting the register model. This
// looks like the ceiling for this kernel design; a further improvement
// would need a fundamentally different approach (e.g. libmarv's actual
// fixed per-(group_size,numRegs)-config kernel family with host-side
// config selection per query length), not another variant of this one.

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include "marv.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// Upper bound on how many simdgroups cooperate within one threadgroup.
// Actual value is chosen at dispatch time from the device's threadgroup
// memory budget; this just avoids pointlessly large threadgroups for short
// queries.
constexpr uint32_t kMaxSimdgroupsPerThreadgroup = 8;
// Upper bound on total simdgroups launched across the whole grid; with
// work-stealing this only needs to be "enough to saturate the GPU", not
// one per database entry.
constexpr uint32_t kMaxTotalSimdgroups = 4096;

const char* kKernelSource = R"METAL(
#include <metal_stdlib>
using namespace metal;

struct Params {
    uint32_t dbEntries;
    uint32_t queryLen;
};

kernel void gaplessAlignCooperative(
    constant Params& params [[buffer(0)]],
    device const int8_t* pssm [[buffer(1)]],
    device const uint8_t* dbData [[buffer(2)]],
    device const uint64_t* offsets [[buffer(3)]],
    device const uint32_t* lengths [[buffer(4)]],
    device atomic_uint* workCounter [[buffer(5)]],
    device int32_t* results [[buffer(6)]],
    threadgroup int32_t* Hall [[threadgroup(0)]],
    uint simdLane [[thread_index_in_simdgroup]],
    uint simdGroupId [[simdgroup_index_in_threadgroup]])
{
    const uint32_t qlen = params.queryLen;
    threadgroup int32_t* myH = Hall + (size_t)simdGroupId * qlen;

    // Each of the 32 lanes owns a contiguous slice of the query. `pssm` rows
    // are pre-transposed into chunkSize*32 blocks so that step `li` reads
    // row[li*32 + lane] -- 32 consecutive bytes -- instead of the scattered
    // row[lane*chunkSize + li] a plain row-major layout would give.
    const uint32_t chunkSize = (qlen + 31u) / 32u;
    const uint32_t chunkStart = min(simdLane * chunkSize, qlen);
    const uint32_t chunkEnd = min(chunkStart + chunkSize, qlen);
    const uint32_t paddedQueryLen = chunkSize * 32u;

    while (true) {
        // Lane 0 claims the next database entry for this simdgroup; the
        // claimed index is broadcast so every lane agrees on when to stop
        // (required for the simd_shuffle_up calls below to stay uniform).
        uint32_t work = 0;
        if (simdLane == 0) {
            work = atomic_fetch_add_explicit(workCounter, 1u, memory_order_relaxed);
        }
        work = simd_broadcast_first(work);
        if (work >= params.dbEntries) {
            break;
        }

        for (uint32_t i = chunkStart; i < chunkEnd; i++) {
            myH[i] = 0;
        }

        const uint64_t start = offsets[work];
        const uint32_t len = lengths[work];

        int32_t best = 0;
        for (uint32_t j = 0; j < len; j++) {
            const uint8_t s = dbData[start + j] & 31;
            device const int8_t* row = pssm + (size_t)s * paddedQueryLen;

            // The only cross-lane dependency: H[i-1][j-1] for this lane's
            // first position is H[chunkEnd-1][j-1] from the lane below,
            // i.e. that lane's last slot *before* it gets overwritten this
            // step. Read-then-shuffle happens before any lane writes to
            // `myH` for step j, so this is race-free without a barrier.
            const int32_t myLastOld = (chunkEnd > chunkStart) ? myH[chunkEnd - 1] : 0;
            int32_t borderIn = simd_shuffle_up(myLastOld, 1);
            if (simdLane == 0) {
                borderIn = 0; // query position -1 boundary
            }

            int32_t carry = borderIn;
            const uint32_t numLocal = chunkEnd - chunkStart;
            for (uint32_t li = 0; li < numLocal; li++) {
                const uint32_t i = chunkStart + li;
                const int32_t oldVal = myH[i];
                const int32_t val = max(carry + int32_t(row[li * 32u + simdLane]), 0);
                myH[i] = val;
                best = max(best, val);
                carry = oldVal;
            }
        }

        const int32_t groupBest = simd_max(best);
        if (simdLane == 0) {
            results[work] = groupBest;
        }
    }
}
)METAL";

// Score/index pair used by the bounded top-K selection in collectScan().
struct ScoreIdx {
    int32_t score;
    uint32_t id;
};

} // namespace

struct MarvDb {
    id<MTLBuffer> data = nil;
    id<MTLBuffer> offsets = nil;
    id<MTLBuffer> lengths = nil;
    size_t dbEntries = 0;
};

struct MarvContext {
    id<MTLDevice> device = nil;
    id<MTLCommandQueue> queue = nil;
    id<MTLComputePipelineState> pipeline = nil;
    size_t maxThreadgroupMemory = 0;

    id<MTLBuffer> resultsBuf = nil;      // dbEntries int32 scores, reused across scans
    id<MTLBuffer> workCounterBuf = nil;  // single atomic_uint, reset before each scan

    id<MTLBuffer> pssmBuf = nil;         // reused, grown on demand instead of realloc'd every scan
    size_t pssmCapacityBytes = 0;

    // In-flight GPU work started by submitScan(), consumed by collectScan().
    id<MTLCommandBuffer> pendingCmdBuf = nil;
    size_t pendingTotalCells = 0;
    std::chrono::high_resolution_clock::time_point pendingSubmitTime;

    std::vector<ScoreIdx> topKHeap;      // reused scratch for collectScan()'s top-K selection

    MarvDb* activeDb = nullptr;
    size_t maxResults = 0;
};

Marv::Marv(size_t dbEntries, int alphabetSize, int maxSeqLength, size_t maxSeqs, Marv::AlignmentType alignmentType)
    : dbEntries(dbEntries), alphabetSize(alphabetSize), dbmanager(NULL), alignmentType(alignmentType) {
    (void)maxSeqLength;

    if (alignmentType != AlignmentType::GAPLESS) {
        throw std::runtime_error(
            "Marv (Metal/MPS backend): only AlignmentType::GAPLESS is currently implemented. "
            "Gapped Smith-Waterman GPU alignment is not yet ported to Apple Silicon -- "
            "use the default --prefilter-mode 0 with --gpu 1.");
    }

    MarvContext* ctx = new MarvContext();
    ctx->device = MTLCreateSystemDefaultDevice();
    if (!ctx->device) {
        delete ctx;
        throw std::runtime_error("Marv (MPS): no Metal device available");
    }
    ctx->maxThreadgroupMemory = [ctx->device maxThreadgroupMemoryLength];

    MTLCompileOptions* compileOptions = [MTLCompileOptions new];
    compileOptions.languageVersion = MTLLanguageVersion3_0;

    NSError* error = nil;
    id<MTLLibrary> library = [ctx->device newLibraryWithSource:[NSString stringWithUTF8String:kKernelSource]
                                                         options:compileOptions
                                                           error:&error];
    if (!library) {
        std::string msg = "Marv (MPS): shader compile failed: ";
        msg += error ? [[error localizedDescription] UTF8String] : "unknown error";
        delete ctx;
        throw std::runtime_error(msg);
    }
    id<MTLFunction> fn = [library newFunctionWithName:@"gaplessAlignCooperative"];
    ctx->pipeline = [ctx->device newComputePipelineStateWithFunction:fn error:&error];
    if (!ctx->pipeline) {
        std::string msg = "Marv (MPS): pipeline creation failed: ";
        msg += error ? [[error localizedDescription] UTF8String] : "unknown error";
        delete ctx;
        throw std::runtime_error(msg);
    }
    ctx->queue = [ctx->device newCommandQueue];

    ctx->maxResults = std::min(maxSeqs, dbEntries > 0 ? dbEntries : maxSeqs);
    ctx->resultsBuf = [ctx->device newBufferWithLength:std::max((size_t)1, dbEntries) * sizeof(int32_t)
                                                options:MTLResourceStorageModeShared];
    ctx->workCounterBuf = [ctx->device newBufferWithLength:sizeof(uint32_t)
                                                    options:MTLResourceStorageModeShared];
    ctx->topKHeap.reserve(ctx->maxResults);

    cudasw = static_cast<void*>(ctx);
}

Marv::~Marv() {
    MarvContext* ctx = static_cast<MarvContext*>(cudasw);
    delete ctx;
}

std::vector<int> Marv::getDeviceIds() {
    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    if (!device) {
        return {};
    }
    return {0};
}

void* Marv::loadDb(char* data, size_t* offset, int32_t* length, size_t dbByteSize) {
    MarvContext* ctx = static_cast<MarvContext*>(cudasw);

    MarvDb* db = new MarvDb();
    db->dbEntries = dbEntries;
    db->data = [ctx->device newBufferWithBytes:data length:dbByteSize options:MTLResourceStorageModeShared];
    if (!db->data) {
        delete db;
        throw std::runtime_error("Marv (MPS): failed to allocate GPU buffer for database");
    }

    std::vector<uint64_t> offsets64(dbEntries);
    std::vector<uint32_t> lengths32(dbEntries);
    for (size_t i = 0; i < dbEntries; i++) {
        offsets64[i] = static_cast<uint64_t>(offset[i]);
        lengths32[i] = static_cast<uint32_t>(length[i]);
    }
    db->offsets = [ctx->device newBufferWithBytes:offsets64.data()
                                            length:offsets64.size() * sizeof(uint64_t)
                                           options:MTLResourceStorageModeShared];
    db->lengths = [ctx->device newBufferWithBytes:lengths32.data()
                                            length:lengths32.size() * sizeof(uint32_t)
                                           options:MTLResourceStorageModeShared];
    return static_cast<void*>(db);
}

void* Marv::loadDb(char* /*data*/, size_t /*dbByteSize*/, void* /*otherdb*/) {
    throw std::runtime_error("Marv (MPS): cross-process external DB sharing is not supported by the Metal backend");
}

void Marv::setDb(void* dbhandle) {
    MarvContext* ctx = static_cast<MarvContext*>(cudasw);
    ctx->activeDb = static_cast<MarvDb*>(dbhandle);
}

void Marv::setDbWithAllocation(void* /*dbhandle*/, const std::string& /*allocationinfo*/) {
    throw std::runtime_error("Marv (MPS): cross-process external DB sharing is not supported by the Metal backend");
}

std::string Marv::getDbMemoryHandle() {
    throw std::runtime_error("Marv (MPS): cross-process external DB sharing is not supported by the Metal backend");
}

void Marv::printInfo() {
    MarvContext* ctx = static_cast<MarvContext*>(cudasw);
    fprintf(stderr, "Marv (Metal/MPS): device=%s dbEntries=%zu alphabetSize=%d\n",
            [[ctx->device name] UTF8String], dbEntries, alphabetSize);
}

void Marv::prefetch() {
    // No-op: Apple Silicon uses unified memory, so there is no separate
    // device memory to migrate data into ahead of time.
}

void Marv::startTimer() {}
void Marv::stopTimer() {}

// sequence must be encoded
void Marv::submitScan(const char* /*sequence*/, size_t sequenceLength, int8_t* pssm) {
    MarvContext* ctx = static_cast<MarvContext*>(cudasw);
    if (ctx->activeDb == nullptr) {
        throw std::runtime_error("Marv (MPS): submitScan() called before setDb()");
    }
    if (ctx->pendingCmdBuf != nil) {
        throw std::runtime_error("Marv (MPS): submitScan() called while a previous scan is still pending -- call collectScan() first");
    }
    MarvDb* db = ctx->activeDb;
    const uint32_t queryLen = static_cast<uint32_t>(sequenceLength);

    const size_t bytesPerSimdgroup = (size_t)queryLen * sizeof(int32_t);
    if (bytesPerSimdgroup > ctx->maxThreadgroupMemory) {
        throw std::runtime_error(
            "Marv (MPS): query too long for the Metal backend (" + std::to_string(queryLen) +
            " residues needs " + std::to_string(bytesPerSimdgroup) + " bytes of threadgroup memory, device limit is " +
            std::to_string(ctx->maxThreadgroupMemory) + ")");
    }
    uint32_t simdgroupsPerTG = static_cast<uint32_t>(
        std::min<size_t>(kMaxSimdgroupsPerThreadgroup, ctx->maxThreadgroupMemory / bytesPerSimdgroup));
    simdgroupsPerTG = std::max(1u, simdgroupsPerTG);

    // Reuse the pssm buffer across calls instead of allocating a fresh
    // MTLBuffer every query -- it's a small copy, but resource creation has
    // real fixed overhead that adds up over hundreds/thousands of queries.
    //
    // Uploaded pre-transposed to match the kernel's per-lane striping (see
    // the file header comment): row `r`'s chunkSize*32 block holds
    // dst[li*32+lane] = pssm[r][lane*chunkSize+li], so that the 32 lanes'
    // simultaneous reads at loop step `li` land on 32 consecutive bytes
    // instead of being chunkSize apart. Padding past each lane's actual
    // chunkEnd is never read by the kernel, so it's left uninitialized.
    const uint32_t chunkSize = (queryLen + 31u) / 32u;
    const uint32_t paddedQueryLen = chunkSize * 32u;
    const size_t pssmBytesNeeded = (size_t)alphabetSize * paddedQueryLen * sizeof(int8_t);
    if (pssmBytesNeeded > ctx->pssmCapacityBytes) {
        ctx->pssmBuf = [ctx->device newBufferWithLength:pssmBytesNeeded options:MTLResourceStorageModeShared];
        if (!ctx->pssmBuf) {
            throw std::runtime_error("Marv (MPS): failed to allocate pssm buffer");
        }
        ctx->pssmCapacityBytes = pssmBytesNeeded;
    }
    int8_t* dst = static_cast<int8_t*>([ctx->pssmBuf contents]);
    for (uint32_t r = 0; r < (uint32_t)alphabetSize; r++) {
        const int8_t* srcRow = pssm + (size_t)r * queryLen;
        int8_t* dstRow = dst + (size_t)r * paddedQueryLen;
        for (uint32_t lane = 0; lane < 32u; lane++) {
            const uint32_t chunkStart = std::min(lane * chunkSize, queryLen);
            const uint32_t chunkEnd = std::min(chunkStart + chunkSize, queryLen);
            for (uint32_t li = 0; li < chunkEnd - chunkStart; li++) {
                dstRow[li * 32u + lane] = srcRow[chunkStart + li];
            }
        }
    }

    std::memset([ctx->workCounterBuf contents], 0, sizeof(uint32_t));

    struct {
        uint32_t dbEntries;
        uint32_t queryLen;
    } params{static_cast<uint32_t>(db->dbEntries), queryLen};

    ctx->pendingTotalCells = (size_t)queryLen * (size_t)[db->data length];
    ctx->pendingSubmitTime = std::chrono::high_resolution_clock::now();

    if (db->dbEntries == 0) {
        ctx->pendingCmdBuf = nil; // nothing to dispatch; collectScan() handles this
        return;
    }

    const uint32_t totalSimdgroups = static_cast<uint32_t>(std::min<size_t>(db->dbEntries, kMaxTotalSimdgroups));
    const uint32_t numThreadgroups = std::max(1u, (totalSimdgroups + simdgroupsPerTG - 1) / simdgroupsPerTG);

    id<MTLCommandBuffer> cmdBuf = [ctx->queue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [cmdBuf computeCommandEncoder];
    [encoder setComputePipelineState:ctx->pipeline];
    [encoder setBytes:&params length:sizeof(params) atIndex:0];
    [encoder setBuffer:ctx->pssmBuf offset:0 atIndex:1];
    [encoder setBuffer:db->data offset:0 atIndex:2];
    [encoder setBuffer:db->offsets offset:0 atIndex:3];
    [encoder setBuffer:db->lengths offset:0 atIndex:4];
    [encoder setBuffer:ctx->workCounterBuf offset:0 atIndex:5];
    [encoder setBuffer:ctx->resultsBuf offset:0 atIndex:6];
    [encoder setThreadgroupMemoryLength:simdgroupsPerTG * bytesPerSimdgroup atIndex:0];

    [encoder dispatchThreadgroups:MTLSizeMake(numThreadgroups, 1, 1)
             threadsPerThreadgroup:MTLSizeMake(simdgroupsPerTG * 32, 1, 1)];
    [encoder endEncoding];
    [cmdBuf commit];

    ctx->pendingCmdBuf = cmdBuf; // not waited on -- collectScan() does that
}

Marv::Stats Marv::collectScan(Result* results) {
    MarvContext* ctx = static_cast<MarvContext*>(cudasw);
    if (ctx->activeDb == nullptr) {
        throw std::runtime_error("Marv (MPS): collectScan() called before setDb()");
    }
    MarvDb* db = ctx->activeDb;

    if (ctx->pendingCmdBuf != nil) {
        [ctx->pendingCmdBuf waitUntilCompleted];
        if (ctx->pendingCmdBuf.status == MTLCommandBufferStatusError) {
            std::string msg = "Marv (MPS): command buffer failed: ";
            msg += ctx->pendingCmdBuf.error ? [[ctx->pendingCmdBuf.error localizedDescription] UTF8String] : "unknown error";
            ctx->pendingCmdBuf = nil;
            throw std::runtime_error(msg);
        }
        ctx->pendingCmdBuf = nil;
    } else if (db->dbEntries != 0) {
        throw std::runtime_error("Marv (MPS): collectScan() called without a matching submitScan()");
    }

    auto timeEnd = std::chrono::high_resolution_clock::now();
    double seconds = std::chrono::duration<double>(timeEnd - ctx->pendingSubmitTime).count();

    const int32_t* scores = static_cast<const int32_t*>([ctx->resultsBuf contents]);

    // Bounded top-K selection: keep a min-heap of at most maxResults
    // entries instead of sorting all dbEntries scores every query. Cheap
    // even for very large databases since the heap never grows past
    // maxResults (typically a few hundred).
    const size_t maxResults = std::min(ctx->maxResults, db->dbEntries);
    std::vector<ScoreIdx>& heap = ctx->topKHeap;
    heap.clear();
    auto heapCmp = [](const ScoreIdx& a, const ScoreIdx& b) { return a.score > b.score; }; // min-heap by score
    for (size_t i = 0; i < db->dbEntries; i++) {
        const int32_t sc = scores[i];
        if (heap.size() < maxResults) {
            heap.push_back(ScoreIdx{sc, static_cast<uint32_t>(i)});
            std::push_heap(heap.begin(), heap.end(), heapCmp);
        } else if (sc > heap.front().score) {
            std::pop_heap(heap.begin(), heap.end(), heapCmp);
            heap.back() = ScoreIdx{sc, static_cast<uint32_t>(i)};
            std::push_heap(heap.begin(), heap.end(), heapCmp);
        }
    }
    std::sort(heap.begin(), heap.end(), [](const ScoreIdx& a, const ScoreIdx& b) { return a.score > b.score; });

    for (size_t i = 0; i < heap.size(); i++) {
        results[i] = Result(heap[i].id, heap[i].score, -1, -1);
    }

    Stats stats;
    stats.results = heap.size();
    stats.numOverflows = 0;
    stats.seconds = seconds;
    stats.gcups = seconds > 0 ? (double)ctx->pendingTotalCells / seconds / 1e9 : 0;
    return stats;
}

Marv::Stats Marv::scan(const char* sequence, size_t sequenceLength, int8_t* pssm, Result* results) {
    submitScan(sequence, sequenceLength, pssm);
    return collectScan(results);
}
