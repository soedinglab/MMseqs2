/*
 * greedycluster -- the final clustering step of the distributed linclust.
 *
 * Reads the surviving edges, runs linclust's greedy, and writes the clustering.
 *
 * The whole stage rests on length-ranked keys. Stock's greedy walks sequences
 * longest-first and makes each unassigned one a representative; key 0 is the
 * longest sequence, so that order is exactly ascending key order. The greedy
 * therefore collapses to a single left-to-right sweep -- no connected components,
 * no iteration to a fixpoint, no cross-worker negotiation.
 *
 * The sweep must visit **every key**, not only the keys that appear as
 * representatives in some edge. When it reaches key k unclaimed, k becomes a
 * representative *there and then*, which is what stops a later, higher-keyed
 * representative from claiming it. Skipping edge-less keys and sweeping them up
 * as singletons afterwards is subtly wrong: with `--cov-mode 1` assignGroup emits
 * reversed edges (kmermatcher.cpp, COV_MODE_TARGET branch), so a representative's
 * key can exceed its member's, and those late claims would steal keys that stock
 * had already made representatives. Measured, that mis-assigned 54 of 1,000,000
 * sequences.
 *
 * The second thing that falls out is the memory. Stock keeps
 * `assignedCluster[dbSize]` at 8 bytes per sequence (Align2clust.cpp:445) --
 * 800 GB at 1e11 and 8 TB at 1e12, on a node that has 2 TB. Here the sweep needs
 * only two bits per key:
 *
 *   claimed  -- some earlier representative took this sequence
 *   hasMembers -- this key is a representative that claimed at least one member
 *
 * which is 25 GB at 1e11 and 250 GB at 1e12, both resident on one node. The
 * actual assignment never has to be stored, because the sweep visits
 * representatives in increasing key order and can therefore emit each cluster
 * exactly when it is finished.
 *
 * The sweep is sequential by necessity -- that is what makes it exact -- but the
 * work per bucket (reading and sorting edges) is threaded, and the sweep itself
 * is a bit test per edge.
 */
#include "CandidateEdge.h"
#include "Command.h"
#include "Debug.h"
#include "DenseIndex.h"
#include "FastSort.h"
#include "FileUtil.h"
#include "ParallelCoordination.h"
#include "Parameters.h"
#include "Util.h"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstring>
#include <string>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

#ifdef OPENMP
#include <omp.h>
#endif

// fwrite that fails loudly. Every caller here is building a final result file that
// the workflow renames into place on success; a short write that goes unnoticed
// becomes a truncated clustering the next restart treats as finished.
//
// Outside the OPENMP guard: the callers are not, so a -DREQUIRE_OPENMP=0 build
// (the "Old compilers" CI task) would not compile with it inside.
static void writeAllOrDie(const void *data, size_t bytes, FILE *file, const std::string &path) {
    if (bytes == 0) {
        return;
    }
    if (fwrite(data, 1, bytes, file) != bytes) {
        Debug(Debug::ERROR) << "Cannot write " << bytes << " bytes to " << path << ": "
                            << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
}

namespace {

// One bit per key -- "decided": either claimed by a representative or made one.
// The only structure in this stage that scales with the database, at 12.5 GB for
// 1e11 sequences and 125 GB for 1e12, against stock's 8 bytes per sequence.
class KeyFlags {
public:
    explicit KeyFlags(uint64_t keyCount) : words((keyCount + 63) / 64, 0) {}

    bool isDecided(uint64_t key) const { return (words[key >> 6] >> (key & 63)) & 1ULL; }
    void decide(uint64_t key) { words[key >> 6] |= 1ULL << (key & 63); }

    uint64_t bytes() const { return words.size() * sizeof(uint64_t); }

private:
    std::vector<uint64_t> words;
};

// Appends "a\tb\n". snprintf costs more than the rest of the sweep at scale.
inline void appendPair(std::string &out, uint64_t a, uint64_t b) {
    char tmp[48];
    int n = 0;
    char digits[24];
    int d = 0;
    do { digits[d++] = static_cast<char>('0' + (a % 10)); a /= 10; } while (a);
    while (d) tmp[n++] = digits[--d];
    tmp[n++] = '\t';
    do { digits[d++] = static_cast<char>('0' + (b % 10)); b /= 10; } while (b);
    while (d) tmp[n++] = digits[--d];
    tmp[n++] = '\n';
    out.append(tmp, n);
}

bool compareByRepThenMember(const CandidateEdge &a, const CandidateEdge &b) {
    const uint64_t ra = a.getRep(), rb = b.getRep();
    if (ra != rb) return ra < rb;
    return a.getMember() < b.getMember();
}

// Reads one bucket's edges with parallel pread.
//
// The sweep that consumes them has to be sequential, but this does not: at 1e11
// sequences the surviving edges are ~7.6 TB and reading them is the dominant cost
// of the stage, so it is split across threads into a preallocated buffer.
size_t readBucket(const std::string &alnDir, unsigned int bucket, std::vector<CandidateEdge> &out,
                  int threads) {
    out.clear();
    const std::string path = EdgeWriter::partitionPath(alnDir, bucket);
    if (FileUtil::fileExists(path.c_str()) == false) {
        return 0;
    }
    const size_t bytes = FileUtil::getFileSize(path);
    if (bytes % sizeof(CandidateEdge) != 0) {
        Debug(Debug::ERROR) << "Edge file " << path << " is " << bytes
                            << " bytes, not a whole number of edges\n";
        EXIT(EXIT_FAILURE);
    }
    if (bytes == 0) {
        return 0;
    }
    out.resize(bytes / sizeof(CandidateEdge));

    const int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        Debug(Debug::ERROR) << "Cannot open " << path << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    // Whole records per chunk, so no edge straddles two reads.
    const size_t chunkBytes = 64 * 1024 * 1024;
    const size_t chunks = (bytes + chunkBytes - 1) / chunkBytes;
    char *base = reinterpret_cast<char *>(out.data());
    // Both atomic. Threads used to write a plain `bool failed` concurrently, which
    // is a data race and so undefined behaviour -- in practice the main thread
    // could miss the flag entirely and treat a short or failed read as a complete
    // bucket, dropping edges with no diagnostic. The error code is captured by the
    // failing thread too, because reporting the *calling* thread's errno after the
    // region names whatever happened to fail last in that thread, which is usually
    // nothing.
    std::atomic<bool> failed(false);
    std::atomic<int> failedErrno(0);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (size_t c = 0; c < chunks; c++) {
        const size_t from = c * chunkBytes;
        const size_t want = std::min(chunkBytes, bytes - from);
        size_t done = 0;
        while (done < want) {
            const ssize_t got = pread(fd, base + from + done, want - done,
                                      static_cast<off_t>(from + done));
            if (got <= 0) {
                if (got < 0 && errno == EINTR) {
                    continue;
                }
                int expected = 0;
                failedErrno.compare_exchange_strong(expected, got < 0 ? errno : EIO);
                failed.store(true);
                break;
            }
            done += static_cast<size_t>(got);
        }
    }
    close(fd);
    if (failed.load()) {
        Debug(Debug::ERROR) << "Cannot read " << path << ": " << strerror(failedErrno.load())
                            << "\n";
        EXIT(EXIT_FAILURE);
    }
    return out.size();
}

}  // namespace

int greedycluster(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_CLUST);

    const std::string seqDb = par.db1;
    const std::string alnDir = par.db2;
    const std::string outFile = par.db3;

    const DenseIndex::Info info = DenseIndex::readInfo(seqDb);
    par.printParameters(command.cmd, argc, argv, *command.params);

    // Bucket count and span come from the layout alignparallel recorded, and the
    // align queue is checked complete before anything is read.
    //
    // Counting consecutive p<n>.edges files instead -- which is what this used to
    // do -- treats a *prefix* of the buckets as the whole layout. An align stage
    // that had produced only the first few buckets therefore looked complete, and
    // bucketSpan was recomputed from the smaller count, so the sweep walked key
    // ranges that did not correspond to any bucket. The result is a successful,
    // silently wrong clustering, and a restart preserves it because the greedy
    // output now exists.
    const std::string layoutPath = alnDir + "/coord/align.info";
    if (FileUtil::fileExists(layoutPath.c_str()) == false) {
        Debug(Debug::ERROR) << "No align layout at " << layoutPath
                            << ". Run alignparallel first (with this build: older ones did not "
                            << "record the layout).\n";
        EXIT(EXIT_FAILURE);
    }
    unsigned int bucketCount = 0;
    uint64_t bucketSpan = 0;
    uint64_t layoutEntryCount = 0;
    {
        FILE *f = FileUtil::openFileOrDie(layoutPath.c_str(), "r", true);
        char name[64];
        size_t value;
        while (fscanf(f, "%63s\t%zu\n", name, &value) == 2) {
            const std::string key = name;
            if (key == "bucketCount") bucketCount = static_cast<unsigned int>(value);
            else if (key == "bucketSpan") bucketSpan = value;
            else if (key == "entryCount") layoutEntryCount = value;
        }
        fclose(f);
    }
    if (bucketCount == 0 || bucketSpan == 0) {
        Debug(Debug::ERROR) << "Align layout " << layoutPath << " is incomplete\n";
        EXIT(EXIT_FAILURE);
    }
    if (layoutEntryCount != info.entryCount) {
        Debug(Debug::ERROR) << "The alignments in " << alnDir << " were produced for "
                            << layoutEntryCount << " sequences, but " << seqDb << " holds "
                            << info.entryCount << ". They are from different runs.\n";
        EXIT(EXIT_FAILURE);
    }
    {
        std::vector<int64_t> workers;
        if (WorkQueue::readCompletedWorkers(alnDir + "/coord/align.queue", workers) == false) {
            Debug(Debug::ERROR) << "No align work queue in " << alnDir
                                << "/coord. Run alignparallel first.\n";
            EXIT(EXIT_FAILURE);
        }
        if (workers.size() != bucketCount) {
            Debug(Debug::ERROR) << "The align queue covers " << workers.size()
                                << " buckets but the layout records " << bucketCount << "\n";
            EXIT(EXIT_FAILURE);
        }
        for (size_t i = 0; i < workers.size(); i++) {
            if (workers[i] < 0) {
                Debug(Debug::ERROR) << "Bucket " << i << " was never aligned. Re-run "
                                    << "alignparallel before clustering.\n";
                EXIT(EXIT_FAILURE);
            }
        }
    }
    for (unsigned int b = 0; b < bucketCount; b++) {
        if (FileUtil::fileExists(EdgeWriter::partitionPath(alnDir, b).c_str()) == false) {
            Debug(Debug::ERROR) << "Bucket " << b << " is recorded aligned but "
                                << EdgeWriter::partitionPath(alnDir, b) << " is missing\n";
            EXIT(EXIT_FAILURE);
        }
    }

    KeyFlags flags(info.entryCount);
    Debug(Debug::INFO) << "Sweeping " << info.entryCount << " keys over " << bucketCount
                       << " edge buckets, " << flags.bytes() / (1024 * 1024) << " MB of flags\n";

    FILE *out = FileUtil::openAndDelete(outFile.c_str(), "w");
    std::string buffer;
    buffer.reserve(64 * 1024 * 1024);

    uint64_t clusterCount = 0;
    uint64_t assignedCount = 0;
    std::vector<CandidateEdge> edges;

    // Edge buckets are contiguous ascending representative-key ranges, the same
    // ranges kmerreduceparallel derived, so walking buckets in order walks the key
    // space in order. bucketSpan comes from the layout above rather than being
    // re-derived, so it cannot disagree with the ranges the edges were bucketed on.

    for (unsigned int bucket = 0; bucket < bucketCount; bucket++) {
        const uint64_t lo = bucket * bucketSpan;
        if (lo >= info.entryCount) {
            break;
        }
        const uint64_t hi = std::min(lo + bucketSpan, info.entryCount);
        readBucket(alnDir, bucket, edges, par.threads);
        // Verified, not sorted. alignparallel produces this file already in
        // (rep, member) order: mergePairCopies sorts by (rep, member, diagonal,
        // strand), compacts to one record per pair in place, and the survivors are
        // then appended in index order into one file per bucket, published by
        // atomic rename. Re-sorting up to ~74 GB of edges per bucket at 1e12 was
        // therefore pure cost -- but the sweep below is silently wrong if the order
        // ever does not hold, so the assumption is checked rather than assumed. A
        // linear scan against an O(n log n) sort is a trade worth making in both
        // directions.
        if (edges.empty() == false
                && std::is_sorted(edges.begin(), edges.end(), compareByRepThenMember) == false) {
            Debug(Debug::ERROR) << "Edge bucket " << EdgeWriter::partitionPath(alnDir, bucket)
                                << " is not ordered by (representative, member). The greedy sweep "
                                << "consumes it in one forward pass and would silently drop the "
                                << "out-of-order edges. It was not written by this build's "
                                << "alignparallel.\n";
            EXIT(EXIT_FAILURE);
        }

        size_t p = 0;
        for (uint64_t key = lo; key < hi; key++) {
            while (p < edges.size() && edges[p].getRep() < key) {
                p++;
            }
            const size_t from = p;
            while (p < edges.size() && edges[p].getRep() == key) {
                p++;
            }
            if (flags.isDecided(key)) {
                continue;
            }

            // Stock forms a cluster only when at least two of its members are still
            // unassigned -- the representative itself plus one other
            // (Align2clust.cpp:362, `validMemberIds.size() <= 1`). A key whose
            // partners have all been taken is therefore *not* made its own cluster
            // here; it stays undecided and a later, higher-keyed representative may
            // still claim it. Deciding it early costs 27 wrong assignments per
            // million; never deciding it costs 54 the other way.
            bool hasFreePartner = false;
            for (size_t e = from; e < p; e++) {
                if (flags.isDecided(edges[e].getMember()) == false) {
                    hasFreePartner = true;
                    break;
                }
            }
            if (hasFreePartner == false) {
                continue;
            }

            flags.decide(key);
            appendPair(buffer, key, key);
            clusterCount++;
            for (size_t e = from; e < p; e++) {
                const uint64_t member = edges[e].getMember();
                if (flags.isDecided(member)) {
                    continue;
                }
                flags.decide(member);
                appendPair(buffer, key, member);
                assignedCount++;
            }
            if (buffer.size() > 32 * 1024 * 1024) {
                writeAllOrDie(buffer.data(), buffer.size(), out, outFile);
                buffer.clear();
            }
        }
    }
    if (buffer.empty() == false) {
        writeAllOrDie(buffer.data(), buffer.size(), out, outFile);
        buffer.clear();
    }

    // Whatever the sweep never decided is a singleton. Each key is independent, so
    // this is threaded; blocks keep the buffered text bounded, and `schedule(static)`
    // hands out ascending contiguous ranges so thread order is key order.
    uint64_t singletonCount = 0;
    {
        const uint64_t blockKeys = 64 * 1024 * 1024;
        std::vector<std::string> threadBuffers(par.threads);
        for (uint64_t blockStart = 0; blockStart < info.entryCount; blockStart += blockKeys) {
            const uint64_t blockEnd = std::min(blockStart + blockKeys, info.entryCount);
            uint64_t blockSingletons = 0;
#pragma omp parallel num_threads(par.threads) reduction(+ : blockSingletons)
            {
                int tid = 0;
#ifdef OPENMP
                tid = omp_get_thread_num();
#endif
                std::string &mine = threadBuffers[tid];
                mine.clear();
#pragma omp for schedule(static)
                for (uint64_t key = blockStart; key < blockEnd; key++) {
                    if (flags.isDecided(key) == false) {
                        appendPair(mine, key, key);
                        blockSingletons++;
                    }
                }
            }
            for (int t = 0; t < par.threads; t++) {
                if (threadBuffers[t].empty() == false) {
                    writeAllOrDie(threadBuffers[t].data(), threadBuffers[t].size(), out, outFile);
                }
            }
            singletonCount += blockSingletons;
        }
    }

    if (fclose(out) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << outFile << "\n";
        EXIT(EXIT_FAILURE);
    }

    Debug(Debug::INFO) << clusterCount + singletonCount << " clusters (" << singletonCount
                       << " singletons), " << assignedCount
                       << " sequences assigned to another representative\n";
    return EXIT_SUCCESS;
}
