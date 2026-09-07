#include "Lin8DbReader.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "Debug.h"
#include "FileUtil.h"
#include "Util.h"
#include "NodePlacement.h"
#include "SubstitutionMatrix.h"
#include "ReducedMatrix.h"
#include "FastSort.h"
#include "Timer.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <climits>
#include <fcntl.h>
#include <unistd.h>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

struct __attribute__((packed)) HashEntry {
    uint64_t hash;
    uint64_t rank;

    static const uint64_t ANCHOR = uint64_t(1) << 63;

    uint64_t rankOf() const { return rank & ~ANCHOR; }
    bool isAnchor() const { return (rank & ANCHOR) != 0; }

    static bool byHashAndRank(const HashEntry &first, const HashEntry &second) {
        if (first.hash != second.hash) {
            return first.hash < second.hash;
        }
        return first.rankOf() < second.rankOf();
    }
};

static const uint64_t LINCLUSTHASH_MAGIC = 0x4C494E4348504153ull;

struct PairFileHeader {
    uint64_t magic;
    uint64_t version;
    uint64_t keyWidth;
    uint64_t pairs;
};

struct ClusterPair {
    uint64_t member;
    uint64_t representative;

    static bool byMember(const ClusterPair &first, const ClusterPair &second) {
        return first.member < second.member;
    }

    static bool byMemberAndRepresentative(const ClusterPair &first, const ClusterPair &second) {
        if (first.member != second.member) {
            return first.member < second.member;
        }
        return first.representative < second.representative;
    }

    static bool sameMember(const ClusterPair &first, const ClusterPair &second) {
        return first.member == second.member;
    }

    static bool isSelf(const ClusterPair &pair) {
        return pair.member == pair.representative;
    }
};

static std::pair<size_t, size_t> runsInFileSlot(const SequenceLocator &runs, size_t fileSlot) {
    size_t begin = runs.size();
    size_t end = 0;
    for (size_t i = 0; i < runs.size(); i++) {
        if (runs[i].fileIdx() % runs.filesPerNode() == fileSlot) {
            begin = std::min(begin, i);
            end = std::max(end, i + 1);
        }
    }
    return std::make_pair(begin, end);
}

static const uint64_t ANCHOR_BASE = 1099511628211ull;

static double spentHashing = 0;
static double spentSorting = 0;
static double spentGrouping = 0;
static double spentComparing = 0;
static double spentSpilling = 0;
static double hashingThreadSeconds = 0;
static uint64_t lengthGroupsSeen = 0;
static uint64_t lengthGroupsOnOneThread = 0;

static uint32_t anchorLength(float identity) {
    if (identity >= 1.0f) {
        return 12;
    }
    const double fits = 1.0 / (2.0 * (1.0 - (double) identity));
    return (uint32_t) std::min(12.0, std::max(8.0, fits));
}

static const unsigned int ANCHOR_COUNT = 1;
static const unsigned int KEYS_PER_SEQUENCE = ANCHOR_COUNT + 1;
static const size_t ANCHOR_BUCKET_MAX = 256;

static uint64_t mixLengthIntoHash(uint64_t hash, uint64_t length) {
    return hash ^ (length + 0x9e3779b97f4a7c15ull + (hash << 6) + (hash >> 2));
}

static uint64_t hashSequence(const char *data, uint32_t length, const unsigned char *aa2num) {
    uint64_t hash = static_cast<uint64_t>(length);
    for (uint32_t i = 0; i < length; i++) {
        hash = (hash * 1099511628211ull) ^ aa2num[static_cast<unsigned char>(data[i])];
    }
    return mixLengthIntoHash(hash, length);
}

static void anchorHashes(const char *data, uint32_t length, uint32_t anchorK,
                         uint64_t *into) {
    if (length < anchorK) {
        for (unsigned int i = 0; i < ANCHOR_COUNT; i++) {
            into[i] = 0;
        }
        return;
    }
    uint64_t best[ANCHOR_COUNT];
    for (unsigned int i = 0; i < ANCHOR_COUNT; i++) {
        best[i] = UINT64_MAX;
    }
    uint64_t hash = 0;
    uint64_t drop = 1;
    for (uint32_t i = 0; i < anchorK; i++) {
        hash = hash * ANCHOR_BASE + static_cast<unsigned char>(data[i]);
        drop = i + 1 < anchorK ? drop * ANCHOR_BASE : drop;
    }
    for (uint32_t at = 0; at + anchorK <= length; at++) {
        if (at > 0) {
            hash -= drop * static_cast<unsigned char>(data[at - 1]);
            hash = hash * ANCHOR_BASE + static_cast<unsigned char>(data[at + anchorK - 1]);
        }
        uint64_t spread = hash * 0xff51afd7ed558ccdull;
        spread ^= spread >> 33;
        spread *= 0xc4ceb9fe1a85ec53ull;
        spread ^= spread >> 29;
        for (unsigned int slot = 0; slot < ANCHOR_COUNT; slot++) {
            if (spread < best[slot]) {
                for (unsigned int back = ANCHOR_COUNT - 1; back > slot; back--) {
                    best[back] = best[back - 1];
                }
                best[slot] = spread;
                break;
            }
            if (spread == best[slot]) {
                break;
            }
        }
    }
    for (unsigned int i = 0; i < ANCHOR_COUNT; i++) {
        into[i] = best[i];
    }
}

static bool sequencesMatch(const char *left, const char *right, uint32_t length, float identity) {
    if (identity >= 1.0f) {
        return memcmp(left, right, length) == 0;
    }
    const uint32_t need = static_cast<uint32_t>(std::ceil(identity * length - 1e-6));
    uint32_t same = 0;
    uint32_t left_over = length;
    for (uint32_t i = 0; i < length; i++) {
        if (left[i] == right[i]) {
            same++;
        }
        left_over--;
        if (same + left_over < need) {
            return false;
        }
    }
    return same >= need;
}

static void reduceOneHashBucket(const RunDbReader &reader, const HashEntry *bucket, size_t size,
                          uint32_t length, float identity,
                          std::vector<char> &claimed, std::vector<ClusterPair> &out) {
    if (size < 2) {
        return;
    }
    claimed.assign(size, 0);
    for (size_t i = 0; i < size; i++) {
        if (claimed[i]) {
            continue;
        }
        const char *query = reader.getData(bucket[i].rankOf());
        for (size_t j = i + 1; j < size; j++) {
            if (claimed[j]) {
                continue;
            }
            if (sequencesMatch(query, reader.getData(bucket[j].rankOf()), length, identity)) {
                claimed[j] = 1;
                ClusterPair pair;
                pair.member = bucket[j].rankOf();
                pair.representative = bucket[i].rankOf();
                out.push_back(pair);
            }
        }
        claimed[i] = 1;
    }
}

static const size_t MAX_HASH_PARTITIONS = 4096;

static size_t hashPartitionCount(size_t entries, size_t budget) {
    const size_t need = entries * sizeof(HashEntry) * 3;
    if (budget == 0 || need <= budget) {
        return 1;
    }
    return std::min(MAX_HASH_PARTITIONS, (need + budget - 1) / budget);
}

class HashPartitions {
public:
    HashPartitions(const std::string &prefix, size_t partitions, unsigned int threads)
        : prefix(prefix), partitions(partitions), files(partitions, NULL),
          pending(threads, std::vector<std::vector<HashEntry> >(partitions)) {
        for (size_t i = 0; i < partitions; i++) {
            files[i] = FileUtil::openAndDelete(path(i).c_str(), "wb");
        }
        depth = std::max<size_t>(64, PENDING_BYTES / (partitions * sizeof(HashEntry)));
    }

    void add(unsigned int thread, const HashEntry &entry) {
        const size_t at = (entry.hash >> HASH_SHIFT) % partitions;
        std::vector<HashEntry> &buffer = pending[thread][at];
        buffer.push_back(entry);
        if (buffer.size() >= depth) {
            flush(at, buffer);
        }
    }

    void finish() {
        for (size_t thread = 0; thread < pending.size(); thread++) {
            for (size_t at = 0; at < partitions; at++) {
                flush(at, pending[thread][at]);
            }
        }
        for (size_t i = 0; i < partitions; i++) {
            if (fclose(files[i]) != 0) {
                Debug(Debug::ERROR) << "Cannot close " << path(i) << "\n";
                EXIT(EXIT_FAILURE);
            }
            files[i] = NULL;
        }
    }

    void load(size_t at, std::vector<HashEntry> &into) {
        into.clear();
        const std::string name = path(at);
        FILE *in = fopen(name.c_str(), "rb");
        if (in == NULL) {
            Debug(Debug::ERROR) << "Cannot reopen " << name << "\n";
            EXIT(EXIT_FAILURE);
        }
        fseek(in, 0, SEEK_END);
        const size_t entries = static_cast<size_t>(ftell(in)) / sizeof(HashEntry);
        fseek(in, 0, SEEK_SET);
        into.resize(entries);
        if (entries > 0 && fread(into.data(), sizeof(HashEntry), entries, in) != entries) {
            Debug(Debug::ERROR) << "Partition " << name << " is truncated\n";
            EXIT(EXIT_FAILURE);
        }
#if defined(__linux__)
        posix_fadvise(fileno(in), 0, 0, POSIX_FADV_DONTNEED);
#endif
        fclose(in);
        FileUtil::remove(name.c_str());
    }

private:
    static const size_t PENDING_BYTES = 8 * 1024 * 1024;
    static const unsigned int HASH_SHIFT = 40;

    std::string path(size_t at) const { return prefix + "." + SSTR(at); }

    void flush(size_t at, std::vector<HashEntry> &buffer) {
        if (buffer.empty()) {
            return;
        }
#pragma omp critical(linclusthash_partition)
        {
            if (fwrite(buffer.data(), sizeof(HashEntry), buffer.size(), files[at]) != buffer.size()) {
                Debug(Debug::ERROR) << "Cannot write to " << path(at) << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        buffer.clear();
    }

    std::string prefix;
    size_t partitions;
    size_t depth;
    std::vector<FILE *> files;
    std::vector<std::vector<std::vector<HashEntry> > > pending;
};

static const size_t HASH_PER_THREAD = 256;
static const size_t READ_ARENA_BYTES = 64u << 20;

static unsigned int threadsWorthStarting(size_t work, unsigned int threads) {
    const size_t want = work / HASH_PER_THREAD;
    return (unsigned int) std::max<size_t>(1, std::min<size_t>(threads, want));
}

static void reduceEntriesToClusters(const RunDbReader &reader, std::vector<HashEntry> &entries,
                           uint32_t length, float identity,
                           unsigned int threads, std::vector<ClusterPair> &out, size_t &crowded) {
    double mark = omp_get_wtime();
    SORT_PARALLEL(entries.begin(), entries.end(), HashEntry::byHashAndRank);
    spentSorting += omp_get_wtime() - mark;
    mark = omp_get_wtime();
    std::vector<size_t> bucketStart;
    for (size_t i = 0; i < entries.size();) {
        size_t j = i + 1;
        while (j < entries.size() && entries[j].hash == entries[i].hash) {
            j++;
        }
        if (j - i > 1 && (j - i <= ANCHOR_BUCKET_MAX || entries[i].isAnchor() == false)) {
            bucketStart.push_back(i);
            bucketStart.push_back(j);
        } else if (j - i > 1) {
            crowded++;
        }
        i = j;
    }
    spentGrouping += omp_get_wtime() - mark;
    mark = omp_get_wtime();
    const unsigned int useThreads = threadsWorthStarting(bucketStart.size() / 2, threads);
    std::vector<std::vector<ClusterPair> > pairsPerThread(useThreads);
#pragma omp parallel num_threads(useThreads)
    {
        unsigned int thread = 0;
#ifdef OPENMP
        thread = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::vector<char> claimed;
#pragma omp for schedule(dynamic, 1)
        for (size_t b = 0; b < bucketStart.size() / 2; b++) {
            reduceOneHashBucket(reader, entries.data() + bucketStart[2 * b],
                              bucketStart[2 * b + 1] - bucketStart[2 * b], length, identity,
                              claimed, pairsPerThread[thread]);
        }
    }
    spentComparing += omp_get_wtime() - mark;
    for (unsigned int i = 0; i < useThreads; i++) {
        out.insert(out.end(), pairsPerThread[i].begin(), pairsPerThread[i].end());
    }
}

static void hashRankRange(const RunDbReader &reader, uint64_t from, uint64_t until, uint32_t length,
                          const unsigned char *aa2num, uint32_t anchorK, unsigned int thread,
                          std::vector<HashEntry> &into) {
    const size_t perBatch = std::max<size_t>(1, reader.batchRoomFor(length));
    std::vector<uint64_t> want[RunDbReader::LANES];
    size_t took[RunDbReader::LANES] = {0, 0};
    unsigned int lane = 0;
    uint64_t at = from;
    struct Fill {
        static size_t run(const RunDbReader &reader, uint64_t &at, uint64_t until, size_t perBatch,
                          unsigned int thread, unsigned int lane, std::vector<uint64_t> &want) {
            want.clear();
            while (at < until && want.size() < perBatch) {
                if (reader.isValid(at)) {
                    want.push_back(at);
                }
                at++;
            }
            if (want.empty()) {
                return 0;
            }
            const size_t took = 1 + reader.startBatch(want[0], want.data() + 1, want.size() - 1, thread, lane);
            if (took < want.size()) {
                at = want[took];
            }
            return took;
        }
    };
    took[lane] = Fill::run(reader, at, until, perBatch, thread, lane, want[lane]);
    while (took[lane] > 0) {
        const unsigned int next = lane ^ 1u;
        took[next] = Fill::run(reader, at, until, perBatch, thread, next, want[next]);
        reader.awaitBatch(thread, lane);
        for (size_t i = 0; i < took[lane]; i++) {
            const char *data = (i == 0) ? reader.batchQueryAt(thread, lane)
                                        : reader.batchAt(thread, lane, i - 1);
            HashEntry entry;
            entry.rank = want[lane][i];
            entry.hash = hashSequence(data, length, aa2num);
            into.push_back(entry);
            uint64_t anchor[ANCHOR_COUNT];
            anchorHashes(data, length, anchorK, anchor);
            for (unsigned int a = 0; a < ANCHOR_COUNT; a++) {
                if (anchor[a] == UINT64_MAX || anchor[a] == 0) {
                    continue;
                }
                entry.hash = anchor[a];
                entry.rank = want[lane][i] | HashEntry::ANCHOR;
                into.push_back(entry);
            }
        }
        lane = next;
    }
}

static void reduceSequencesOfOneLength(const RunDbReader &reader, uint64_t rankBegin, uint64_t rankEnd,
                             uint32_t length, const unsigned char *aa2num, float identity,
                             uint32_t anchorK, size_t budget, unsigned int threads, const std::string &tmpPrefix,
                             std::vector<ClusterPair> &out, size_t &crowded) {
    const size_t count = static_cast<size_t>(rankEnd - rankBegin);
    if (count < 2) {
        return;
    }
    const size_t partitions = hashPartitionCount(count * KEYS_PER_SEQUENCE, budget);
    std::vector<HashEntry> entries;
    if (partitions == 1) {
        const double mark = omp_get_wtime();
        const unsigned int useThreads = threadsWorthStarting(count, threads);
        std::vector<std::vector<HashEntry> > perThread(useThreads);
#pragma omp parallel num_threads(useThreads)
        {
            unsigned int thread = 0;
#ifdef OPENMP
            thread = static_cast<unsigned int>(omp_get_thread_num());
#endif
            const uint64_t span = rankEnd - rankBegin;
            const uint64_t from = rankBegin + span * thread / useThreads;
            const uint64_t until = rankBegin + span * (thread + 1) / useThreads;
            hashRankRange(reader, from, until, length, aa2num, anchorK, thread, perThread[thread]);
        }
        for (unsigned int i = 0; i < useThreads; i++) {
            entries.insert(entries.end(), perThread[i].begin(), perThread[i].end());
            std::vector<HashEntry>().swap(perThread[i]);
        }
        const double took = omp_get_wtime() - mark;
        spentHashing += took;
        hashingThreadSeconds += took * useThreads;
        lengthGroupsSeen++;
        lengthGroupsOnOneThread += (useThreads == 1);
        reduceEntriesToClusters(reader, entries, length, identity, threads, out, crowded);
        return;
    }

    double mark = omp_get_wtime();
    const unsigned int useThreads = threadsWorthStarting(count, threads);
    HashPartitions parts(tmpPrefix, partitions, useThreads);
#pragma omp parallel num_threads(useThreads)
    {
        unsigned int thread = 0;
#ifdef OPENMP
        thread = static_cast<unsigned int>(omp_get_thread_num());
#endif
        RunDbReader::Cursor cursor;
#pragma omp for schedule(static)
        for (uint64_t rank = rankBegin; rank < rankEnd; rank++) {
            if (reader.isValid(rank) == false) {
                continue;
            }
            const char *data = reader.getData(rank, cursor);
            HashEntry entry;
            entry.rank = rank;
            entry.hash = hashSequence(data, length, aa2num);
            parts.add(thread, entry);
            uint64_t anchor[ANCHOR_COUNT];
            anchorHashes(data, length, anchorK, anchor);
            for (unsigned int a = 0; a < ANCHOR_COUNT; a++) {
                if (anchor[a] == UINT64_MAX || anchor[a] == 0) {
                    continue;
                }
                entry.hash = anchor[a];
                entry.rank = rank | HashEntry::ANCHOR;
                parts.add(thread, entry);
            }
        }
    }
    parts.finish();
    {
        const double took = omp_get_wtime() - mark;
        spentHashing += took;
        hashingThreadSeconds += took * useThreads;
        lengthGroupsSeen++;
        lengthGroupsOnOneThread += (useThreads == 1);
    }
    for (size_t at = 0; at < partitions; at++) {
        mark = omp_get_wtime();
        parts.load(at, entries);
        spentSpilling += omp_get_wtime() - mark;
        reduceEntriesToClusters(reader, entries, length, identity, threads, out, crowded);
    }
}

static uint64_t mergeNodeShards(const std::string &base, const std::string &keptBitmap,
                            unsigned int nodes, const RunDbReader &reader, const std::string &tag) {
    std::vector<uint64_t> valid;
    if (reader.hasValid()) {
        valid.assign(reader.validWords(), reader.validWords() + reader.validWordCount());
    } else {
        valid.assign(reader.getSize() / 64 + (reader.getSize() % 64 != 0), UINT64_MAX);
        if (reader.getSize() % 64 != 0) {
            valid.back() = (uint64_t(1) << (reader.getSize() % 64)) - 1;
        }
    }
    const std::string pairsTmp = base + ".tmp" + tag;
    FILE *out = FileUtil::openAndDelete(pairsTmp.c_str(), "wb");
    PairFileHeader header;
    header.magic = LINCLUSTHASH_MAGIC;
    header.version = 1;
    header.keyWidth = sizeof(uint64_t);
    header.pairs = 0;
    if (fwrite(&header, sizeof(PairFileHeader), 1, out) != 1) {
        Debug(Debug::ERROR) << "Cannot write to " << pairsTmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    std::vector<ClusterPair> batch(1024 * 1024);
    for (unsigned int node = 0; node < nodes; node++) {
        const std::string shard = base + "." + SSTR(node);
        FILE *in = fopen(shard.c_str(), "rb");
        if (in == NULL) {
            Debug(Debug::ERROR) << "Cannot open the shard " << shard << "\n";
            EXIT(EXIT_FAILURE);
        }
        PairFileHeader part;
        if (fread(&part, sizeof(PairFileHeader), 1, in) != 1 || part.magic != LINCLUSTHASH_MAGIC) {
            Debug(Debug::ERROR) << "File " << shard << " is not a cluster pair file\n";
            EXIT(EXIT_FAILURE);
        }
        uint64_t left = part.pairs;
        while (left > 0) {
            const size_t now = static_cast<size_t>(std::min<uint64_t>(left, batch.size()));
            if (fread(batch.data(), sizeof(ClusterPair), now, in) != now
                || fwrite(batch.data(), sizeof(ClusterPair), now, out) != now) {
                Debug(Debug::ERROR) << "Shard " << shard << " is truncated\n";
                EXIT(EXIT_FAILURE);
            }
            for (size_t i = 0; i < now; i++) {
                valid[batch[i].member >> 6] &= ~(uint64_t(1) << (batch[i].member & 63));
            }
            left -= now;
        }
        fclose(in);
        header.pairs += part.pairs;
    }
    rewind(out);
    if (fwrite(&header, sizeof(PairFileHeader), 1, out) != 1 || fclose(out) != 0) {
        Debug(Debug::ERROR) << "Cannot finish " << pairsTmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::publishAtomically(pairsTmp, base);

    const std::string keptTmp = keptBitmap + ".tmp" + tag;
    FILE *bits = FileUtil::openAndDelete(keptTmp.c_str(), "wb");
    const uint64_t validHeader[3] = {RunDbReader::VALID_MAGIC, reader.getSize(), 1};
    if (fwrite(validHeader, sizeof(uint64_t), 3, bits) != 3) {
        Debug(Debug::ERROR) << "Cannot write the kept bitmap header to " << keptTmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (valid.empty() == false
        && fwrite(valid.data(), sizeof(uint64_t), valid.size(), bits) != valid.size()) {
        Debug(Debug::ERROR) << "Cannot write the kept bitmap to " << keptTmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fclose(bits) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << keptTmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::publishAtomically(keptTmp, keptBitmap);

    uint64_t alive = 0;
    for (size_t i = 0; i < valid.size(); i++) {
        alive += static_cast<uint64_t>(__builtin_popcountll(valid[i]));
    }
    return alive;
}

int lin8clusthash(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.alphabetSize = MultiParam<NuclAA<int> >(NuclAA<int>(Parameters::CLUST_HASH_DEFAULT_ALPH_SIZE, 5));
    par.seqIdThr = static_cast<float>(Parameters::CLUST_HASH_DEFAULT_MIN_SEQ_ID) / 100.0f;
    par.parseParameters(argc, argv, command, true, 0, 0);
    FileUtil::fixRlimitNoFile();

    const NodePlacement node = NodePlacement::resolve(par);
    // a retry only waits for the other nodes; skip the reader, the matrix and the header it prints
    const bool alreadyFolded = FileUtil::fileExists(nodeDonePath(par.db2, node.index).c_str());
    RunDbReader reader(par.db1);
    // node 0 merges even on a retry, and the merge needs the reader's size and valid bitmap
    if (alreadyFolded == false || node.index == 0) {
        reader.open();
    }
    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);
    BaseMatrix *hashMat = &subMat;
    ReducedMatrix *reduced = NULL;
    size_t budget = 0;
    float identity = 1.0f;
    const std::string shard = par.db2 + "." + SSTR(node.index);
    if (alreadyFolded == false) {
        if (par.alphabetSize.values.aminoacid() < subMat.alphabetSize) {
            reduced = new ReducedMatrix(subMat.probMatrix, subMat.subMatrixPseudoCounts, subMat.aa2num,
                                        subMat.num2aa, subMat.alphabetSize,
                                        par.alphabetSize.values.aminoacid(), 2.0);
            hashMat = reduced;
        }
        budget = static_cast<size_t>(Util::computeMemory(par.splitMemoryLimit) * 0.95);
        identity = reduced != NULL ? (float) par.seqIdThr : 1.0f;
        Debug(Debug::INFO) << "Reducing to a " << hashMat->alphabetSize << "-letter alphabet at identity " << identity << "\n";
        Debug(Debug::INFO) << "Database size: " << reader.getSize() << " sequences in "
                           << reader.getSequenceLocator().size() << " length groups, budget "
                           << (budget / (1024 * 1024)) << " MB\n";
    }

    if (alreadyFolded == false) {
        const std::string pairsTmp = shard + ".tmp" + uniqueTmpSuffix();
        FILE *out = FileUtil::openAndDelete(pairsTmp.c_str(), "wb");
        PairFileHeader header;
        header.magic = LINCLUSTHASH_MAGIC;
        header.version = 1;
        header.keyWidth = sizeof(uint64_t);
        header.pairs = 0;
        if (fwrite(&header, sizeof(PairFileHeader), 1, out) != 1) {
            Debug(Debug::ERROR) << "Cannot write to " << pairsTmp << "\n";
            EXIT(EXIT_FAILURE);
        }

        Timer timer;
        reader.openBatch(par.threads, READ_ARENA_BYTES, budget, RunDbReader::READ_AGAIN);
        const SequenceLocator &runs = reader.getSequenceLocator();
        const std::vector<size_t> mine = nodeFileSlots(runs, node);
        Debug(Debug::INFO) << "Node " << node.index << " of " << node.count << " takes " << mine.size()
                           << " of " << runs.filesPerNode() << " length groups\n";

        std::vector<ClusterPair> pairs;
        uint64_t written = 0;
        size_t oversized = 0;
        size_t crowded = 0;
        size_t moved = 0;
        size_t stayed = 0;
        size_t lengths = 0;
        for (size_t at = 0; at < mine.size(); at++) {
            const std::pair<size_t, size_t> span = runsInFileSlot(runs, mine[at]);
            for (size_t segment = span.first; segment < span.second; segment++) {
                lengths += (segment == span.first || runs[segment].seqLen() != runs[segment - 1].seqLen());
            }
        }
        Debug(Debug::INFO) << "Hashing and comparing " << lengths << " length groups\n";
        Debug::Progress progress(lengths);
        for (size_t at = 0; at < mine.size(); at++) {
            const std::pair<size_t, size_t> span = runsInFileSlot(runs, mine[at]);
            for (size_t segment = span.first; segment < span.second;) {
                const uint32_t length = runs[segment].seqLen();
                const uint64_t rankBegin = runs[segment].rankBase();
                size_t next = segment + 1;
                while (next < span.second && runs[next].seqLen() == length) {
                    next++;
                }
                const uint64_t rankEnd = runs.rankEnd(next - 1);
                if (hashPartitionCount(static_cast<size_t>(rankEnd - rankBegin), budget) > 1) {
                    oversized++;
                }
                reduceSequencesOfOneLength(reader, rankBegin, rankEnd, length, hashMat->aa2num, identity,
                                 anchorLength(identity), budget,
                                 par.threads, par.db2 + ".part" + uniqueTmpSuffix(), pairs, crowded);
                SORT_PARALLEL(pairs.begin(), pairs.end(), ClusterPair::byMemberAndRepresentative);
                pairs.erase(std::unique(pairs.begin(), pairs.end(), ClusterPair::sameMember),
                            pairs.end());
                for (size_t i = 0; i < pairs.size(); i++) {
                    ClusterPair want;
                    want.member = pairs[i].representative;
                    want.representative = 0;
                    std::vector<ClusterPair>::iterator at =
                        std::lower_bound(pairs.begin(), pairs.end(), want, ClusterPair::byMember);
                    if (at == pairs.end() || at->member != pairs[i].representative
                        || ClusterPair::isSelf(*at)) {
                        continue;
                    }
                    moved++;
                    const uint64_t root = at->representative;
                    ClusterPair above;
                    above.member = root;
                    above.representative = 0;
                    std::vector<ClusterPair>::iterator up =
                        std::lower_bound(pairs.begin(), pairs.end(), above, ClusterPair::byMember);
                    const bool rootWasRemoved = up != pairs.end() && up->member == root
                                                && ClusterPair::isSelf(*up) == false;
                    if (rootWasRemoved == false
                        && sequencesMatch(reader.getData(pairs[i].member), reader.getData(root),
                                          length, identity)) {
                        pairs[i].representative = root;
                    } else {
                        pairs[i].representative = pairs[i].member;
                        stayed++;
                    }
                }
                pairs.erase(std::remove_if(pairs.begin(), pairs.end(), ClusterPair::isSelf),
                            pairs.end());
                if (pairs.empty() == false) {
                    if (fwrite(pairs.data(), sizeof(ClusterPair), pairs.size(), out) != pairs.size()) {
                        Debug(Debug::ERROR) << "Cannot write cluster pairs to " << pairsTmp << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                    written += pairs.size();
                    pairs.clear();
                }
                progress.updateProgress();
                segment = next;
            }
            reader.releaseFileSlot(mine[at]);
        }
        if (oversized > 0) {
            Debug(Debug::INFO) << "Split " << oversized << " length groups by hash to fit the memory budget\n";
        }
        if (moved > 0) {
            Debug(Debug::INFO) << "Chained members: " << moved << ", " << stayed
                               << " did not match their representative and stayed on their own\n";
        }
        if (crowded > 0) {
            Debug(Debug::INFO) << "Skipped " << crowded << " k-mers shared by more than " << ANCHOR_BUCKET_MAX
                               << " sequences of one length\n";
        }
        header.pairs = written;
        rewind(out);
        if (fwrite(&header, sizeof(PairFileHeader), 1, out) != 1 || fclose(out) != 0) {
            Debug(Debug::ERROR) << "Cannot finish " << pairsTmp << "\n";
            EXIT(EXIT_FAILURE);
        }
        FileUtil::publishAtomically(pairsTmp, shard);
        markNodeDone(par.db2, node.index);

        Debug(Debug::INFO) << "Found " << written << " redundant sequences in " << timer.lap() << "\n";
    }

    if (node.index != 0) {
        delete reduced;
        if (alreadyFolded == false) {
            reader.close();
        }
        return EXIT_SUCCESS;
    }
    Debug(Debug::INFO) << "Folding done, merging the other nodes\n";
    waitEveryNodeDone(par.db2, node.count, 86400);

    const uint64_t all = mergeNodeShards(par.db2, par.db1 + RunDbReader::KEPT_BITMAP_SUFFIX,
                                         node.count, reader, uniqueTmpSuffix());
    if (lengthGroupsSeen > 0)
    Debug(Debug::INFO) << "Time for hashing: " << (uint64_t) spentHashing
                       << "s sorting: " << (uint64_t) spentSorting << "s grouping: "
                       << (uint64_t) spentGrouping << "s comparing: " << (uint64_t) spentComparing
                       << "s writing: " << (uint64_t) spentSpilling << "s\n";
    if (lengthGroupsSeen > 0)
    Debug(Debug::INFO) << "Hashing used " << (spentHashing > 0 ? hashingThreadSeconds / spentHashing : 0)
                       << " threads on average over " << lengthGroupsSeen << " length groups, "
                       << lengthGroupsOnOneThread << " on one thread\n";
    Debug(Debug::INFO) << "Kept " << all << " of " << reader.getSize() << " sequences ("
                       << (100.0 * static_cast<double>(reader.getSize() - all)
                           / static_cast<double>(std::max<uint64_t>(reader.getSize(), 1)))
                       << "% removed)\n";
    reader.close();
    return EXIT_SUCCESS;
}
