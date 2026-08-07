/*
 * createrepdb -- builds the representative sub-database for the next linclust pass.
 *
 * Replaces `createsubdb` for this pipeline, for two reasons.
 *
 * The first is the usual one: `createsubdb` opens the database with USE_INDEX and
 * so holds the text index resident, 24 bytes per sequence -- 2.4 TB at 1e11. This
 * streams through the dense companion index instead and holds one bit per key.
 *
 * The second is specific and easy to miss: **a sub-database of representatives has
 * sparse keys**, and every stage of the distributed pipeline assumes dense ones --
 * DenseIndex treats a key as an array offset, the greedy sweep walks the key space
 * in order, and the edge and merge buckets are key ranges. Feeding a `createsubdb`
 * output into pass 2 would break all three. So the representatives are *re-keyed*
 * densely here, and a `subkey -> original key` map is written alongside so the
 * pass-2 clustering can be translated back before merging.
 *
 * Re-keying is nearly free because of how the keys were assigned in the first
 * place. Original keys are length-ranked, so any subset of them is *already* in
 * descending length order; the i-th representative in ascending original-key order
 * is exactly sub-key i. No sort is needed, and the new database is a sequential
 * copy rather than a gather.
 *
 * Memory is O(N/8 + R/CHUNK_ENTRIES), not O(R)
 * -------------------------------------------
 * The copy used to materialise, for every representative, a packed source index
 * entry (12 B), a destination offset (8 B), a length (4 B) and a key-map slot
 * (8 B). At the 37% representative rate measured on the 5B run that is ~1.29 TB at
 * 1e11 and ~12.9 TB at 1e12, on a node with 2 TB -- so the stage could not run at
 * the scale the pipeline exists for.
 *
 * None of it has to be resident. Destination offsets are a prefix sum, and a
 * prefix sum only needs a checkpoint every CHUNK_ENTRIES representatives: the
 * chunk that owns a representative recomputes the offsets inside itself from its
 * own checkpoint. Source entries and lengths come back out of the source
 * `.index.bin`, which the chunk reads sequentially anyway. The key map is written
 * straight to its file as the plan pass discovers it.
 *
 * What is left is two 8-byte checkpoints per chunk -- 195 MB at 1e12 with
 * CHUNK_ENTRIES = 8192 -- plus the representative bitmap at N/8 (125 GB at 1e12,
 * which is what the whole pipeline already budgets for a key-space bitmap).
 *
 * The cost is one extra sequential pass over the source index, since the copy
 * re-reads what the plan pass already saw. That is 12 B per source entry against
 * the 32 B per representative it replaces, and it is streaming rather than
 * resident.
 */
#include "Command.h"
#include "Debug.h"
#include "DenseIndex.h"
#include "LengthRankTable.h"
#include "FastSort.h"
#include "FileUtil.h"
#include "Parameters.h"
#include "Util.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <string>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

#ifdef OPENMP
#include <omp.h>
#endif

namespace {

// Representatives per copy chunk. Sets both the checkpoint density (two 8-byte
// values per chunk) and the size of a worker's buffers (one index entry and one
// sequence per representative), so it trades resident memory against I/O size.
// 8192 gives ~2 MB of sequence data per chunk on protein input, which is already
// a comfortable write, and 195 MB of checkpoints at 1e12.
const uint64_t CHUNK_ENTRIES = 8192;

class Bitmap {
public:
    explicit Bitmap(uint64_t bits) : words((bits + 63) / 64, 0) {}
    bool get(uint64_t i) const { return (words[i >> 6] >> (i & 63)) & 1ULL; }
    void set(uint64_t i) { words[i >> 6] |= 1ULL << (i & 63); }
    // Sets the bit and reports whether it was this call that set it. Used by the
    // threaded TSV scan below, where several threads can meet the same
    // representative and exactly one of them must count it.
    bool testAndSet(uint64_t i) {
        const uint64_t mask = 1ULL << (i & 63);
        const uint64_t previous = __sync_fetch_and_or(&words[i >> 6], mask);
        return (previous & mask) == 0;
    }
    uint64_t bytes() const { return words.size() * sizeof(uint64_t); }

    // Number of set bits below `i`. This is exactly "the sub-key of original key
    // i", because sub-keys are assigned to representatives in ascending original
    // key order.
    //
    // Backed by one checkpoint per RANK_BLOCK keys, so it costs 8 bytes per 4096
    // keys -- 2 GB at 1e12 against the 125 GB the bitmap itself already needs --
    // and a lookup is one array read plus at most 64 popcounts. Holding a
    // key -> sub-key table instead would be 8 bytes per *key*.
    static const uint64_t RANK_BLOCK = 4096;

    void buildRank() {
        rank.assign(words.size() / (RANK_BLOCK / 64) + 1, 0);
        uint64_t total = 0;
        for (size_t w = 0; w < words.size(); w++) {
            if (w % (RANK_BLOCK / 64) == 0) {
                rank[w / (RANK_BLOCK / 64)] = total;
            }
            total += static_cast<uint64_t>(__builtin_popcountll(words[w]));
        }
    }

    uint64_t rankOf(uint64_t i) const {
        const size_t word = static_cast<size_t>(i >> 6);
        const size_t block = word / (RANK_BLOCK / 64);
        uint64_t count = rank[block];
        for (size_t w = block * (RANK_BLOCK / 64); w < word; w++) {
            count += static_cast<uint64_t>(__builtin_popcountll(words[w]));
        }
        const uint64_t bit = i & 63;
        if (bit > 0) {
            count += static_cast<uint64_t>(
                __builtin_popcountll(words[word] & ((1ULL << bit) - 1)));
        }
        return count;
    }

    uint64_t rankBytes() const { return rank.size() * sizeof(uint64_t); }

private:
    std::vector<uint64_t> words;
    std::vector<uint64_t> rank;
};

int openOrDie(const std::string &path, int flags) {
    const int fd = open(path.c_str(), flags, 0666);
    if (fd < 0) {
        Debug(Debug::ERROR) << "Cannot open " << path << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    return fd;
}

void readAt(int fd, void *dst, size_t len, size_t off, const char *what) {
    char *p = static_cast<char *>(dst);
    size_t done = 0;
    while (done < len) {
        const ssize_t got = pread(fd, p + done, len - done, static_cast<off_t>(off + done));
        if (got <= 0) {
            if (got < 0 && errno == EINTR) continue;
            Debug(Debug::ERROR) << "Cannot read " << what << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        done += static_cast<size_t>(got);
    }
}

void writeAt(int fd, const void *src, size_t len, size_t off, const char *what) {
    const char *p = static_cast<const char *>(src);
    size_t done = 0;
    while (done < len) {
        const ssize_t put = pwrite(fd, p + done, len - done, static_cast<off_t>(off + done));
        if (put <= 0) {
            if (put < 0 && errno == EINTR) continue;
            Debug(Debug::ERROR) << "Cannot write " << what << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        done += static_cast<size_t>(put);
    }
}

void allocate(const std::string &path, uint64_t size) {
    const int fd = openOrDie(path, O_WRONLY | O_CREAT | O_TRUNC);
    if (size > 0 && ftruncate(fd, static_cast<off_t>(size)) != 0) {
        Debug(Debug::ERROR) << "Cannot size " << path << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    close(fd);
}

// Where each chunk of CHUNK_ENTRIES representatives begins, on both sides.
//
// srcRow[c] is the source row holding the chunk's first representative, so the
// chunk can seek straight there instead of scanning from zero. dataOffset[c] is
// the destination byte that representative starts at, so the chunk recomputes the
// offsets inside itself and nothing global has to be held.
struct CopyPlan {
    std::vector<uint64_t> srcRow;
    std::vector<uint64_t> dataOffset;
    uint64_t totalBytes;
    uint32_t maxLen;

    CopyPlan() : totalBytes(0), maxLen(0) {}
    // One entry per distinct sequence length among the representatives,
    // longest first. Written as the sub-database's .lenrank.
    std::vector<LengthRankTable::Run> lengthRuns;
};

// One streaming pass over the source index. Optionally emits the sub-key ->
// original-key map as it goes, which is what keeps that 8 B per representative
// off the heap.
CopyPlan planCopy(int srcIdx, const Bitmap &keep, uint64_t entryCount, uint64_t keptCount,
                  FILE *keyMapOut, const std::string &keyMapPath, bool lengthRanked) {
    CopyPlan plan;
    plan.srcRow.reserve(static_cast<size_t>(keptCount / CHUNK_ENTRIES + 2));
    plan.dataOffset.reserve(static_cast<size_t>(keptCount / CHUNK_ENTRIES + 2));

    const uint64_t block = 1 << 20;
    std::vector<DenseIndex::Entry> buf(block);
    std::vector<uint64_t> keyBuf;
    keyBuf.reserve(1 << 16);

    uint64_t kept = 0;
    uint64_t total = 0;
    for (uint64_t from = 0; from < entryCount; from += block) {
        const uint64_t n = std::min(block, entryCount - from);
        readAt(srcIdx, buf.data(), n * sizeof(DenseIndex::Entry), DenseIndex::entryOffset(from),
               "source index");
        for (uint64_t i = 0; i < n; i++) {
            if (keep.get(from + i) == false) continue;
            if (kept % CHUNK_ENTRIES == 0) {
                plan.srcRow.push_back(from + i);
                plan.dataOffset.push_back(total);
            }
            total += buf[i].length;
            plan.maxLen = std::max(plan.maxLen, buf[i].length);
            // The sub-database's length-rank table, accumulated in the same pass.
            //
            // It exists because pass 2 runs on this database and every stage there
            // recovers a sequence length from its key. It is free to build here:
            // original keys are length-ranked, so the representatives -- visited in
            // ascending original-key order -- already come out longest first, and
            // the runs are just the points where the length changes.
            //
            // An entry occupies its residues plus a newline and DBWriter's
            // terminating null, which is the same +2 createdbparallel plans with.
            if (lengthRanked) {
                const uint64_t residues = buf[i].length >= 2 ? buf[i].length - 2 : 0;
                if (plan.lengthRuns.empty() == false &&
                    residues > plan.lengthRuns.back().length) {
                    // Only reachable if the source database is not length-ranked,
                    // which every stage here assumes. Silently writing a
                    // non-monotone table would give pass 2 plausible wrong lengths.
                    Debug(Debug::ERROR)
                        << "Source key " << (from + i) << " is longer (" << residues
                        << ") than a representative before it (" << plan.lengthRuns.back().length
                        << "). The source database is not length-ranked.\n";
                    EXIT(EXIT_FAILURE);
                }
                if (plan.lengthRuns.empty() || plan.lengthRuns.back().length != residues) {
                    LengthRankTable::Run run;
                    run.length = residues;
                    run.firstKey = kept;
                    run.count = 0;
                    plan.lengthRuns.push_back(run);
                }
                plan.lengthRuns.back().count++;
            }
            if (keyMapOut != NULL) {
                keyBuf.push_back(from + i);
                if (keyBuf.size() == keyBuf.capacity()) {
                    if (fwrite(keyBuf.data(), sizeof(uint64_t), keyBuf.size(), keyMapOut)
                            != keyBuf.size()) {
                        Debug(Debug::ERROR) << "Cannot write " << keyMapPath << ": "
                                            << strerror(errno) << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                    keyBuf.clear();
                }
            }
            kept++;
        }
    }
    if (keyMapOut != NULL && keyBuf.empty() == false) {
        if (fwrite(keyBuf.data(), sizeof(uint64_t), keyBuf.size(), keyMapOut) != keyBuf.size()) {
            Debug(Debug::ERROR) << "Cannot write " << keyMapPath << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    if (kept != keptCount) {
        Debug(Debug::ERROR) << "Counted " << keptCount << " representatives in the clustering but "
                            << kept << " in the index. The clustering and the database do not "
                            << "match.\n";
        EXIT(EXIT_FAILURE);
    }
    plan.totalBytes = total;
    return plan;
}

// Copies the entries flagged in `keep` into a new dense database, preserving key
// order. Returns the total data bytes and the longest entry.
struct CopyResult {
    uint64_t dataBytes;
    uint32_t maxLen;
};

// lengthRanked says whether the entries are sequences ordered longest first, in
// which case the sub-database gets its own length-rank table. It is false for the
// header database: an entry's "length" there is its header text, which has no
// relation to sequence length and is not monotone in the key.
CopyResult copyFlagged(const std::string &srcDb, const std::string &dstDb, const Bitmap &keep,
                       uint64_t entryCount, uint64_t keptCount, int threads,
                       FILE *keyMapOut, const std::string &keyMapPath, bool lengthRanked) {
    const int srcIdx = openOrDie(DenseIndex::fileName(srcDb), O_RDONLY);
    const int srcData = openOrDie(srcDb, O_RDONLY);

    const CopyPlan plan =
        planCopy(srcIdx, keep, entryCount, keptCount, keyMapOut, keyMapPath, lengthRanked);

    CopyResult result;
    result.maxLen = plan.maxLen;
    result.dataBytes = plan.totalBytes;

    allocate(dstDb, plan.totalBytes);
    DenseIndex::createEmpty(dstDb, keptCount, 0, plan.totalBytes, plan.maxLen);
    // Pass 2 addresses this database by key and recovers lengths from the key, so
    // the sub-database needs its own length-rank table just as the full one does.
    if (lengthRanked) {
        LengthRankTable::write(dstDb, plan.lengthRuns, keptCount);
    }

    const int dstData = openOrDie(dstDb, O_WRONLY);
    const int dstIdx = openOrDie(DenseIndex::fileName(dstDb), O_WRONLY);

    const int64_t chunkCount = static_cast<int64_t>(plan.srcRow.size());

    // Pass 2: copy, one chunk of CHUNK_ENTRIES representatives at a time. Chunks
    // are independent -- their destination ranges are disjoint and their source
    // ranges are ascending -- so this is a threaded forward scan rather than a
    // gather, exactly as before; what changed is that a chunk re-reads its own
    // slice of the source index instead of every chunk's slice being held.
#pragma omp parallel num_threads(threads)
    {
        std::vector<char> buf;
        std::vector<char> span;
        std::vector<DenseIndex::Entry> entries;
        std::vector<DenseIndex::Entry> idxBuf;
        std::vector<DenseIndex::Entry> scanBuf(4096);
#pragma omp for schedule(dynamic, 1)
        for (int64_t c = 0; c < chunkCount; c++) {
            const uint64_t keptFrom = static_cast<uint64_t>(c) * CHUNK_ENTRIES;
            const uint64_t count = std::min(CHUNK_ENTRIES, keptCount - keptFrom);

            // Collect this chunk's source entries, starting at the row the plan
            // recorded. A run of representatives is contiguous in the *kept*
            // numbering but sparse in the source, so this reads forward until it
            // has `count` of them.
            entries.clear();
            entries.reserve(static_cast<size_t>(count));
            uint64_t row = plan.srcRow[static_cast<size_t>(c)];
            while (entries.size() < count) {
                const uint64_t n = std::min<uint64_t>(scanBuf.size(), entryCount - row);
                if (n == 0) {
                    Debug(Debug::ERROR) << "Ran out of source entries copying " << srcDb
                                        << " chunk " << c << "\n";
                    EXIT(EXIT_FAILURE);
                }
                readAt(srcIdx, scanBuf.data(), static_cast<size_t>(n) * sizeof(DenseIndex::Entry),
                       DenseIndex::entryOffset(row), "source index");
                for (uint64_t i = 0; i < n && entries.size() < count; i++) {
                    if (keep.get(row + i)) {
                        entries.push_back(scanBuf[static_cast<size_t>(i)]);
                    }
                }
                row += n;
            }

            // Destination offsets, recomputed inside the chunk from its checkpoint.
            const uint64_t base = plan.dataOffset[static_cast<size_t>(c)];
            idxBuf.resize(static_cast<size_t>(count));
            uint64_t at = base;
            for (uint64_t k = 0; k < count; k++) {
                idxBuf[static_cast<size_t>(k)].offset = at;
                idxBuf[static_cast<size_t>(k)].length = entries[static_cast<size_t>(k)].length;
                at += entries[static_cast<size_t>(k)].length;
            }
            const uint64_t bytes = at - base;
            buf.resize(static_cast<size_t>(bytes));

            // Coalesced reads. Representatives are a *subset* of the source, so
            // one pread per entry meant 353M random reads at 1e9 and made this
            // stage 4x stock's createsubdb. Consecutive representatives are only a
            // couple of sequence lengths apart, so a run of them is fetched in one
            // read and the wanted pieces copied out.
            //
            // The run is closed when the bytes read would exceed twice the bytes
            // actually wanted: coalescing across a sparse region reads the skipped
            // sequences too, and without a cap a database whose representatives
            // are thinly spread would read the whole file. That bounds the waste
            // at 2x while still collapsing dense regions into single reads.
            uint64_t i = 0;
            while (i < count) {
                const uint64_t from = entries[static_cast<size_t>(i)].offset;
                uint64_t to = from + entries[static_cast<size_t>(i)].length;
                uint64_t wanted = entries[static_cast<size_t>(i)].length;
                uint64_t j = i + 1;
                while (j < count) {
                    const uint64_t end = entries[static_cast<size_t>(j)].offset
                                         + entries[static_cast<size_t>(j)].length;
                    if (end - from > 2 * (wanted + entries[static_cast<size_t>(j)].length)) {
                        break;
                    }
                    to = end;
                    wanted += entries[static_cast<size_t>(j)].length;
                    j++;
                }
                span.resize(static_cast<size_t>(to - from));
                readAt(srcData, span.data(), static_cast<size_t>(to - from), from, "source data");
                for (uint64_t k = i; k < j; k++) {
                    memcpy(buf.data() + (idxBuf[static_cast<size_t>(k)].offset - base),
                           span.data() + (entries[static_cast<size_t>(k)].offset - from),
                           entries[static_cast<size_t>(k)].length);
                }
                i = j;
            }
            writeAt(dstData, buf.data(), static_cast<size_t>(bytes), base, "database");
            writeAt(dstIdx, idxBuf.data(), idxBuf.size() * sizeof(DenseIndex::Entry),
                    DenseIndex::entryOffset(keptFrom), "index");
        }
    }

    close(srcIdx);
    close(srcData);
    close(dstData);
    close(dstIdx);
    return result;
}

// ---------------------------------------------------------------------------
// The pass-2 filter gate, as a CSR file rather than a per-worker vector.
//
// alignparallel's pass-2 gate needs, for a candidate member sub-key t, every
// original key in t's pass-1 cluster (Align2clust.cpp:657-730, the `allpass`
// loop). Every align worker used to build that itself: the whole sub-key ->
// original-key map (8 B x R), a cluster offset per representative (8 B x R) and
// every pass-1 member key (8 B x N), by streaming the pass-1 TSV twice and
// binary-searching the key map once per line. At 1e11 with R/N = 0.4 that is
// ~1.44 TB resident *per worker* on a 2 TB node -- before the bucket's own edges
// and sequence arena -- and ~14.4 TB at 1e12. It is the reason the stage could
// not run at either target.
//
// None of it is per-worker state: it is one function of clu1.tsv and the
// representative set, both of which this stage already has in hand. Built once
// here, in sub-key order, an align worker preads only the offsets and members of
// the sub-keys its own bucket touches -- O(bucket) instead of O(N).
//
// Layout, both fixed width so a slice is addressable by pread:
//   <repDb>.gate.offsets   (R + 1) x uint64, offsets[s] .. offsets[s+1] is the
//                          member range of sub-key s
//   <repDb>.gate.members   the original keys, grouped by sub-key, ascending
//                          within a group so the file is reproducible
//   <repDb>.gate.info      the two counts, for validation on the reading side
// ---------------------------------------------------------------------------

struct __attribute__((__packed__)) GateRow {
    uint64_t sub;     // sub-key of the pass-1 representative
    uint64_t member;  // original key of one of its members
};

// Thread-private spill shards, so the scan of clu1.tsv needs no locking. Same
// shape as translatekeys' pass-2 writers: <prefix>.t<thread>.<bucket>.
class GateSpill {
public:
    GateSpill(const std::string &prefix, unsigned int buckets, size_t rowsPerBuffer)
        : prefix(prefix), buckets(buckets), rowsPerBuffer(rowsPerBuffer), closed(false) {
        opened.assign(buckets, false);
        buffers.resize(buckets);
        for (unsigned int b = 0; b < buckets; b++) {
            buffers[b].reserve(rowsPerBuffer);
        }
    }
    ~GateSpill() { close(); }

    void append(unsigned int b, uint64_t sub, uint64_t member) {
        GateRow row;
        row.sub = sub;
        row.member = member;
        buffers[b].push_back(row);
        if (buffers[b].size() >= rowsPerBuffer) {
            flush(b);
        }
    }

    void close() {
        if (closed) return;
        closed = true;
        for (unsigned int b = 0; b < buckets; b++) {
            flush(b);
        }
    }

    static std::string path(const std::string &prefix, unsigned int b) {
        return prefix + "." + SSTR(b);
    }

    static void read(const std::string &prefix, unsigned int b, std::vector<GateRow> &out) {
        const std::string p = path(prefix, b);
        if (FileUtil::fileExists(p.c_str()) == false) return;
        const size_t bytes = FileUtil::getFileSize(p);
        if (bytes == 0) return;
        if (bytes % sizeof(GateRow) != 0) {
            Debug(Debug::ERROR) << "Gate spill " << p << " is " << bytes << " bytes, not a whole "
                                << "number of rows. Remove it and re-run createrepdb.\n";
            EXIT(EXIT_FAILURE);
        }
        const size_t count = bytes / sizeof(GateRow);
        const size_t at = out.size();
        out.resize(at + count);
        FILE *f = FileUtil::openFileOrDie(p.c_str(), "rb", true);
        if (fread(out.data() + at, sizeof(GateRow), count, f) != count) {
            Debug(Debug::ERROR) << "Cannot read " << p << "\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(f);
    }

private:
    void flush(unsigned int b) {
        if (buffers[b].empty()) return;
        const std::string p = path(prefix, b);
        FILE *f = fopen(p.c_str(), opened[b] ? "ab" : "wb");
        if (f == NULL) {
            Debug(Debug::ERROR) << "Cannot open " << p << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        opened[b] = true;
        if (fwrite(buffers[b].data(), sizeof(GateRow), buffers[b].size(), f)
                != buffers[b].size()) {
            Debug(Debug::ERROR) << "Cannot write " << p << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (fclose(f) != 0) {
            Debug(Debug::ERROR) << "Cannot close " << p << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        buffers[b].clear();
    }

    std::string prefix;
    unsigned int buckets;
    size_t rowsPerBuffer;
    std::vector<std::vector<GateRow> > buffers;
    std::vector<bool> opened;
    bool closed;
};

bool byThenMember(const GateRow &a, const GateRow &b) {
    if (a.sub != b.sub) return a.sub < b.sub;
    return a.member < b.member;
}

// Offset of the first byte after the newline at or after `from`, i.e. where the
// first whole line starting at or after `from` begins. Same record-boundary rule
// createdbparallel uses to make byte ranges of a text file independent.
uint64_t nextLineStart(FILE *file, uint64_t from, uint64_t size) {
    if (from == 0 || from >= size) {
        return std::min(from, size);
    }
    if (fseeko(file, static_cast<off_t>(from - 1), SEEK_SET) != 0) {
        Debug(Debug::ERROR) << "Cannot seek the clustering: " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    int c;
    uint64_t at = from - 1;
    while ((c = getc(file)) != EOF) {
        at++;
        if (c == '\n') {
            return at;
        }
    }
    return size;
}

// Builds the CSR described above. One threaded pass over clu1.tsv into buckets of
// sub-key space, then one ascending sweep that emits each bucket's members
// contiguously and fills its slice of the offsets file.
//
// Resident memory is one bucket: its rows (16 B each) plus its offsets slice
// (8 B per sub-key in the range). The bucket count is derived from
// --split-memory-limit against a conservative upper bound on the row count, so a
// pathological clustering makes more buckets rather than a bigger one.
void buildFilterGate(const std::string &clusterTsv, const std::string &repDb, const Bitmap &isRep,
                     uint64_t entryCount, uint64_t repCount, int threads, size_t memoryLimit) {
    const std::string offsetsPath = repDb + ".gate.offsets";
    const std::string membersPath = repDb + ".gate.members";
    const std::string infoPath = repDb + ".gate.info";
    const std::string spillPrefix = repDb + ".gate.spill";

    const uint64_t tsvSize = FileUtil::getFileSize(clusterTsv);
    // A row is "0\t0\n" at the very least, so this cannot under-count.
    const uint64_t rowsUpperBound = tsvSize / 4 + 1;
    // A floor only so a tiny limit does not explode into pathologically many
    // buckets; low enough that the multi-bucket path is reachable in testing,
    // which matters because that is the path used at scale. Same rule as
    // mergeclusterparallel and translatecluster.
    const uint64_t targetBytes =
        std::max<uint64_t>(Util::computeMemory(memoryLimit) / 8, 1ULL * 1024 * 1024);
    unsigned int buckets = 1;
    while (buckets < 65536
           && ((rowsUpperBound / buckets) * sizeof(GateRow) + (repCount / buckets) * sizeof(uint64_t))
                  > targetBytes) {
        buckets *= 2;
    }
    const uint64_t span = (repCount + buckets - 1) / buckets;
    Debug(Debug::INFO) << "Building the pass-2 filter gate over " << buckets << " sub-key buckets "
                       << "of " << span << "\n";

    // Clear anything an interrupted attempt left: the sweep below reads every
    // shard it finds, so a stale one would be merged into the new gate.
    for (int t = 0; t < std::max(1, threads); t++) {
        for (unsigned int b = 0; b < buckets; b++) {
            const std::string p = GateSpill::path(spillPrefix + ".t" + SSTR(t), b);
            if (FileUtil::fileExists(p.c_str())) {
                FileUtil::remove(p.c_str());
            }
        }
    }

    const int realThreads = std::max(1, threads);
    {
        std::vector<GateSpill *> writers(realThreads, NULL);
        // Derived, not fixed. There are buckets x threads buffers alive at once,
        // so a fixed 4096 rows each is 64 KB x buckets x threads -- 8 GB at 1024
        // buckets and 128 threads, which is the configuration this exists for. The
        // budget is split across all of them, with a floor so a large bucket count
        // degrades to more frequent flushes rather than to one write per row.
        const size_t rowsPerBuffer = std::max<size_t>(
            targetBytes / (static_cast<size_t>(buckets) * static_cast<size_t>(realThreads)
                           * sizeof(GateRow)),
            64);
        for (int t = 0; t < realThreads; t++) {
            writers[t] = new GateSpill(spillPrefix + ".t" + SSTR(t), buckets, rowsPerBuffer);
        }
        const uint64_t rangeSize = (tsvSize + realThreads - 1) / realThreads;
#pragma omp parallel num_threads(realThreads)
        {
            int tid = 0;
#ifdef OPENMP
            tid = omp_get_thread_num();
#endif
            const uint64_t nominalFrom = static_cast<uint64_t>(tid) * rangeSize;
            if (nominalFrom < tsvSize) {
                const uint64_t nominalTo = std::min(nominalFrom + rangeSize, tsvSize);
                FILE *f = FileUtil::openFileOrDie(clusterTsv.c_str(), "r", true);
                const uint64_t from = nextLineStart(f, nominalFrom, tsvSize);
                if (fseeko(f, static_cast<off_t>(from), SEEK_SET) != 0) {
                    Debug(Debug::ERROR) << "Cannot seek " << clusterTsv << ": " << strerror(errno)
                                        << "\n";
                    EXIT(EXIT_FAILURE);
                }
                char *line = NULL;
                size_t cap = 0;
                ssize_t len;
                uint64_t at = from;
                while (at < nominalTo && (len = getline(&line, &cap, f)) > 0) {
                    at += static_cast<uint64_t>(len);
                    char *tab = strchr(line, '\t');
                    if (tab == NULL) continue;
                    const uint64_t rep = strtoull(line, NULL, 10);
                    const uint64_t member = strtoull(tab + 1, NULL, 10);
                    if (rep >= entryCount || member >= entryCount) {
                        Debug(Debug::ERROR) << "Clustering " << clusterTsv << " names key "
                                            << std::max(rep, member) << ", beyond the "
                                            << entryCount << " in the database\n";
                        EXIT(EXIT_FAILURE);
                    }
                    // Every first-column key is a representative by construction,
                    // so this rank is its sub-key.
                    const uint64_t sub = isRep.rankOf(rep);
                    writers[tid]->append(static_cast<unsigned int>(sub / span), sub, member);
                }
                free(line);
                fclose(f);
            }
        }
        for (int t = 0; t < realThreads; t++) {
            delete writers[t];
        }
    }

    // Sweep the buckets in ascending sub-key order, so members are emitted
    // contiguously and the running total *is* the offset.
    const std::string offsetsTmp = offsetsPath + ".tmp";
    const std::string membersTmp = membersPath + ".tmp";
    FILE *offsetsOut = FileUtil::openAndDelete(offsetsTmp.c_str(), "wb");
    FILE *membersOut = FileUtil::openAndDelete(membersTmp.c_str(), "wb");
    std::vector<GateRow> rows;
    std::vector<uint64_t> offsets;
    std::vector<uint64_t> members;
    uint64_t total = 0;
    for (unsigned int b = 0; b < buckets; b++) {
        const uint64_t lo = static_cast<uint64_t>(b) * span;
        if (lo >= repCount) break;
        const uint64_t hi = std::min(lo + span, repCount);

        rows.clear();
        for (int t = 0; t < realThreads; t++) {
            const std::string shardPrefix = spillPrefix + ".t" + SSTR(t);
            GateSpill::read(shardPrefix, b, rows);
            const std::string p = GateSpill::path(shardPrefix, b);
            if (FileUtil::fileExists(p.c_str())) {
                FileUtil::remove(p.c_str());
            }
        }
        // By (sub, member): grouping is what the CSR needs, and the second key
        // makes the file a function of the clustering alone rather than of which
        // thread happened to read which byte range.
        SORT_SERIAL(rows.begin(), rows.end(), byThenMember);

        offsets.assign(static_cast<size_t>(hi - lo), 0);
        members.clear();
        members.reserve(rows.size());
        size_t i = 0;
        for (uint64_t sub = lo; sub < hi; sub++) {
            offsets[static_cast<size_t>(sub - lo)] = total + members.size();
            while (i < rows.size() && rows[i].sub == sub) {
                members.push_back(rows[i].member);
                i++;
            }
        }
        if (i != rows.size()) {
            Debug(Debug::ERROR) << "Gate bucket " << b << " holds sub-key " << rows[i].sub
                                << ", outside its range [" << lo << ", " << hi << ")\n";
            EXIT(EXIT_FAILURE);
        }
        if (offsets.empty() == false
            && fwrite(offsets.data(), sizeof(uint64_t), offsets.size(), offsetsOut)
                   != offsets.size()) {
            Debug(Debug::ERROR) << "Cannot write " << offsetsTmp << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (members.empty() == false
            && fwrite(members.data(), sizeof(uint64_t), members.size(), membersOut)
                   != members.size()) {
            Debug(Debug::ERROR) << "Cannot write " << membersTmp << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        total += members.size();
    }
    // The closing sentinel, so offsets[s + 1] is always readable.
    if (fwrite(&total, sizeof(uint64_t), 1, offsetsOut) != 1) {
        Debug(Debug::ERROR) << "Cannot write " << offsetsTmp << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fclose(offsetsOut) != 0 || fclose(membersOut) != 0) {
        Debug(Debug::ERROR) << "Cannot close the filter gate files: " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (FileUtil::getFileSize(offsetsTmp) != (repCount + 1) * sizeof(uint64_t)) {
        Debug(Debug::ERROR) << offsetsTmp << " holds "
                            << (FileUtil::getFileSize(offsetsTmp) / sizeof(uint64_t))
                            << " offsets, not the " << (repCount + 1) << " expected\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::move(membersTmp.c_str(), membersPath.c_str());
    FileUtil::move(offsetsTmp.c_str(), offsetsPath.c_str());

    const std::string infoTmp = infoPath + ".tmp";
    FILE *infoOut = FileUtil::openAndDelete(infoTmp.c_str(), "w");
    fprintf(infoOut, "repCount\t%zu\n", (size_t)repCount);
    fprintf(infoOut, "memberCount\t%zu\n", (size_t)total);
    if (fclose(infoOut) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << infoTmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::move(infoTmp.c_str(), infoPath.c_str());

    Debug(Debug::INFO) << "Filter gate: " << repCount << " clusters, " << total << " members ("
                       << (total * sizeof(uint64_t)) / (1024 * 1024) << " MB on disk, paged by "
                       << "the align workers rather than held)\n";
}

}  // namespace

int createrepdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    const std::string seqDb = par.db1;
    const std::string clusterTsv = par.db2;
    const std::string repDb = par.db3;

    const DenseIndex::Info info = DenseIndex::readInfo(seqDb);
    if (info.entryCount == 0) {
        Debug(Debug::ERROR) << "Database " << seqDb << " is empty\n";
        EXIT(EXIT_FAILURE);
    }

    // A key is a representative exactly when it appears in the first column.
    //
    // Read by byte range across threads rather than one getline loop. At 1e12 this
    // file is ~28 TB and the scan was the serial head of a stage that is otherwise
    // threaded. Ranges are made independent the same way createdbparallel makes
    // FASTA chunks independent: a range owns the lines that *start* inside it, so
    // each thread skips to the first line boundary at or after its start and reads
    // through the boundary at or after its end.
    Bitmap isRep(info.entryCount);
    uint64_t repCount = 0;
    {
        const uint64_t tsvSize = FileUtil::getFileSize(clusterTsv);
        if (tsvSize == 0) {
            Debug(Debug::ERROR) << "Clustering " << clusterTsv << " is empty\n";
            EXIT(EXIT_FAILURE);
        }
        const int threads = std::max(1, par.threads);
        const uint64_t rangeSize = (tsvSize + threads - 1) / threads;
        uint64_t counted = 0;
#pragma omp parallel num_threads(threads) reduction(+ : counted)
        {
            int tid = 0;
#ifdef OPENMP
            tid = omp_get_thread_num();
#endif
            const uint64_t nominalFrom = static_cast<uint64_t>(tid) * rangeSize;
            if (nominalFrom < tsvSize) {
                const uint64_t nominalTo = std::min(nominalFrom + rangeSize, tsvSize);
                FILE *f = FileUtil::openFileOrDie(clusterTsv.c_str(), "r", true);
                const uint64_t from = nextLineStart(f, nominalFrom, tsvSize);
                if (fseeko(f, static_cast<off_t>(from), SEEK_SET) != 0) {
                    Debug(Debug::ERROR) << "Cannot seek " << clusterTsv << ": " << strerror(errno)
                                        << "\n";
                    EXIT(EXIT_FAILURE);
                }
                char *line = NULL;
                size_t cap = 0;
                ssize_t len;
                uint64_t at = from;
                while (at < nominalTo && (len = getline(&line, &cap, f)) > 0) {
                    at += static_cast<uint64_t>(len);
                    const uint64_t rep = strtoull(line, NULL, 10);
                    if (rep < info.entryCount && isRep.testAndSet(rep)) {
                        counted++;
                    }
                }
                free(line);
                fclose(f);
            }
        }
        repCount = counted;
    }
    Debug(Debug::INFO) << repCount << " representatives of " << info.entryCount << " sequences ("
                       << isRep.bytes() / (1024 * 1024) << " MB of flags)\n";
    if (repCount == 0) {
        Debug(Debug::ERROR) << "No representatives found in " << clusterTsv << "\n";
        EXIT(EXIT_FAILURE);
    }

    // subkey -> original key, dense in the sub-key space and in ascending order,
    // so the translation back is a sequential read rather than a lookup structure.
    // Streamed out by the sequence copy's plan pass, which is the one place that
    // already visits every representative in order.
    const std::string mapFile = repDb + ".keymap";
    const std::string keymapTmp = mapFile + ".tmp";
    FILE *m = FileUtil::openAndDelete(keymapTmp.c_str(), "wb");

    const CopyResult seq =
        copyFlagged(seqDb, repDb, isRep, info.entryCount, repCount, par.threads, m, keymapTmp,
                    true);

    if (fclose(m) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << keymapTmp << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (FileUtil::getFileSize(keymapTmp) != repCount * sizeof(uint64_t)) {
        Debug(Debug::ERROR) << keymapTmp << " holds "
                            << (FileUtil::getFileSize(keymapTmp) / sizeof(uint64_t))
                            << " keys, not the " << repCount << " representatives copied\n";
        EXIT(EXIT_FAILURE);
    }

    // Skipped when the source database has no headers (createdbparallel
    // --write-header-db 0). Nothing between here and the final TSV opens them --
    // accessions come from .lookup -- so their absence is a supported
    // configuration rather than a broken database.
    const bool haveHeaders = FileUtil::fileExists((seqDb + "_h").c_str());
    if (haveHeaders) {
        copyFlagged(seqDb + "_h", repDb + "_h", isRep, info.entryCount, repCount, par.threads, NULL,
                    "", false);
    }

    // The pass-2 filter gate, built once here instead of in every align worker.
    // Needs the rank directory, which turns an original key into its sub-key.
    isRep.buildRank();
    Debug(Debug::INFO) << "Rank directory: " << isRep.rankBytes() / (1024 * 1024) << " MB\n";
    buildFilterGate(clusterTsv, repDb, isRep, info.entryCount, repCount, par.threads,
                    par.splitMemoryLimit);

    const int dbType = FileUtil::parseDbType(seqDb.c_str());
    FileUtil::writeFile(repDb + ".dbtype", reinterpret_cast<const unsigned char *>(&dbType),
                        sizeof(int));
    if (haveHeaders) {
        const int hdrType = Parameters::DBTYPE_GENERIC_DB;
        FileUtil::writeFile(repDb + "_h.dbtype", reinterpret_cast<const unsigned char *>(&hdrType),
                            sizeof(int));
    }

    // The key map is renamed last, after the .dbtype files, because the workflow
    // guards this whole stage on its existence. Publishing it first meant a death
    // in the window before the .dbtype files were written left a marker saying the
    // stage was finished over a database that could not be opened -- and every
    // restart then skipped the stage and failed in the map, permanently. This is
    // the same last-write-wins sentinel rule createdbparallel uses.
    //
    // The rename is also what makes the map itself safe to publish: the pass-2
    // filter gate sizes itself from this file's length, so a truncated one would
    // silently be used as a short, wrong sub-key to original-key map.
    FileUtil::move(keymapTmp.c_str(), mapFile.c_str());

    Debug(Debug::INFO) << "Wrote " << repDb << ": " << repCount << " sequences, " << seq.dataBytes
                       << " data bytes, longest " << seq.maxLen << "\n";
    return EXIT_SUCCESS;
}
