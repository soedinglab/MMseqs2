/*
 * translatecluster -- rewrites a clustering from sub-key space into original keys.
 *
 * `createrepdb` re-keys the representatives densely so the next linclust pass can
 * run on them, which means that pass's clustering comes back in sub-key space and
 * has to be translated before it can be merged with the first pass.
 *
 * The map is dense and ascending in sub-key space, so translating one column is a
 * sequential read. The difficulty is that a clustering has *two* key columns and
 * only one of them can be in order at a time. Holding the whole map resident would
 * be 8 bytes per representative -- ~296 GB at 1e11 given the 37% representative
 * fraction measured on the 5B run, and ~3 TB at 1e12, which does not fit.
 *
 * So each column is translated in its own bucketed pass: partition by the sub-key
 * being translated, then per bucket read exactly that contiguous slice of the map.
 * Peak memory is one slice, and every read of the map is sequential.
 *
 * Multi-worker
 * ------------
 * Every phase runs across the whole allocation, coordinating only through files as
 * the map/reduce/align stages do:
 *
 *   scatter    one item per chunk of the input, bucketed by member sub-key.
 *   remap      one item per bucket: translate the member against that bucket's
 *              slice of the map, then re-scatter by representative sub-key.
 *   emit       one item per bucket: translate the representative and write the
 *              bucket's text to its own part.
 *   assemble   one item per part, pwritten into the output at a prefix-sum offset.
 *
 * Output is byte-identical to a single-worker run, and stays so however the work
 * was distributed. A shard is named by (work item, worker): the item fixes the
 * order it is read back in -- which is the order a single sequential pass would
 * have written it -- and the worker keeps two attempts at one item off a single
 * path, since the writers append and a second attempt truncating the first
 * mid-flight leaves a file that is still a whole number of records, still in
 * range, and silently short. Which attempt counts is recorded by the queue, not
 * guessed.
 *
 *
 * Peak scratch, which --scratch-budget does not see. The workflow only measures
 * the k-mer shuffle against the budget, and these stages hold their spill shards
 * until the whole stage finishes rather than dropping each bucket as it is
 * consumed -- a worker whose lease lapsed may still be reading them, and letting
 * it publish a truncated part over a complete one is worse than the disk. So the
 * peak is the spill *plus* the parts, roughly twice what the single-writer loop
 * needed: about 11 TB at 1e11 against 5 TB. Size the filesystem for that; the
 * budget will not warn.
 * The remap's work item is a *band* of source buckets rather than one bucket.
 * That phase writes into every destination bucket, so one item per source bucket
 * would mean buckets^2 shards -- 6.7e7 at 1e11 on 4096 workers, which no parallel
 * filesystem will walk or unlink in reasonable time. Bands bound it to
 * buckets x bands, and a band is still a contiguous, deterministic range.
 */
#include "Command.h"
#include "Debug.h"
#include "FileUtil.h"
#include "ParallelCoordination.h"
#include "Parameters.h"
#include "StageParts.h"
#include "Util.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <string>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

// fwrite that fails loudly. Every caller here is building a final result file that
// the workflow renames into place on success; a short write that goes unnoticed
// becomes a truncated clustering the next restart treats as finished.
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

struct __attribute__((__packed__)) Pair {
    uint64_t key;    // the sub-key being translated in this pass
    uint64_t other;  // the column carried through untouched
};

class Buckets {
public:
    // Named by (item, worker). The item fixes the order the shards are read back
    // in; the worker keeps two attempts at one item off a single path.
    Buckets(const std::string &prefix, unsigned int count, uint64_t item, int64_t worker,
            size_t bufferPairs = 1 << 16)
        : prefix(prefix), count(count), item(item), worker(worker), bufferPairs(bufferPairs),
          closed(false) {
        opened.assign(count, false);
        buffers.resize(count);
    }
    ~Buckets() { close(); }

    void append(unsigned int b, uint64_t key, uint64_t other) {
        Pair p;
        p.key = key;
        p.other = other;
        buffers[b].push_back(p);
        if (buffers[b].size() >= bufferPairs) flush(b);
    }

    void close() {
        if (closed) return;
        closed = true;
        for (unsigned int b = 0; b < count; b++) {
            flush(b);
        }
    }

    // A bucket's shards, concatenated in work-item order, taking only the attempt
    // at each item that the queue recorded as finished.
    static std::vector<Pair> read(const std::string &prefix, unsigned int b,
                                  const std::vector<int64_t> &authority) {
        std::vector<Pair> out;
        for (size_t item = 0; item < authority.size(); item++) {
        if (authority[item] < 0) {
            Debug(Debug::ERROR) << "Work item " << item << " of " << prefix
                                << " has no completion record although its phase reported "
                                << "done. Remove the working directory and re-run the stage.\n";
            EXIT(EXIT_FAILURE);
        }
        const std::string p = StageParts::shardPath(prefix, b, item, authority[item]);
        if (FileUtil::fileExists(p.c_str()) == false) continue;
        const size_t bytes = FileUtil::getFileSize(p);
        // An empty bucket is legitimate -- a key range no assignment fell into.
        // A *torn* one is not, and returning it as empty silently dropped every
        // assignment in that key range, leaving a smaller clustering that still
        // looks valid. mergeclusterparallel already makes the identical situation
        // fatal, with a comment saying exactly why; the two now agree.
        if (bytes % sizeof(Pair) != 0) {
            Debug(Debug::ERROR) << "Spill bucket " << p << " is " << bytes << " bytes, not a whole "
                                << "number of " << sizeof(Pair) << "-byte pairs. It was left by an "
                                << "interrupted run; remove the spill files and re-run the stage.\n";
            EXIT(EXIT_FAILURE);
        }
        if (bytes == 0) continue;
        const size_t had = out.size();
        out.resize(had + bytes / sizeof(Pair));
        FILE *f = FileUtil::openFileOrDie(p.c_str(), "rb", true);
        if (fread(out.data() + had, sizeof(Pair), out.size() - had, f) != out.size() - had) {
            Debug(Debug::ERROR) << "Cannot read " << p << "\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(f);
        }
        return out;
    }

private:
    void flush(unsigned int b) {
        if (buffers[b].empty()) return;
        // Append-and-close per flush: one descriptor per bucket can reach the
        // derived bucket count, and the close is where a buffered or remote
        // filesystem first reports ENOSPC or a quota. Losing whole records keeps
        // the file record-aligned, so downstream never notices the truncation.
        FILE *file = fopen(StageParts::shardPath(prefix, b, item, worker).c_str(),
                           opened[b] ? "ab" : "wb");
        if (file == NULL) {
            Debug(Debug::ERROR) << "Cannot open "
                                << StageParts::shardPath(prefix, b, item, worker) << ": "
                                << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        opened[b] = true;
        if (fwrite(buffers[b].data(), sizeof(Pair), buffers[b].size(), file) !=
            buffers[b].size()) {
            Debug(Debug::ERROR) << "Cannot write bucket " << b << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (fclose(file) != 0) {
            Debug(Debug::ERROR) << "Cannot close "
                                << StageParts::shardPath(prefix, b, item, worker) << ": " << strerror(errno)
                                << "\n";
            EXIT(EXIT_FAILURE);
        }
        buffers[b].clear();
    }

    std::string prefix;
    unsigned int count;
    uint64_t item;
    int64_t worker;
    size_t bufferPairs;
    std::vector<std::vector<Pair> > buffers;
    std::vector<bool> opened;
    bool closed;
};

// Reads map[lo, hi) -- one contiguous slice, never the whole file.
std::vector<uint64_t> mapSlice(int fd, uint64_t lo, uint64_t hi) {
    std::vector<uint64_t> out(static_cast<size_t>(hi - lo));
    size_t want = out.size() * sizeof(uint64_t), done = 0;
    char *p = reinterpret_cast<char *>(out.data());
    while (done < want) {
        const ssize_t got = pread(fd, p + done, want - done,
                                  static_cast<off_t>(lo * sizeof(uint64_t) + done));
        if (got <= 0) {
            if (got < 0 && errno == EINTR) continue;
            Debug(Debug::ERROR) << "Cannot read key map: " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        done += static_cast<size_t>(got);
    }
    return out;
}

inline void appendPair(std::string &out, uint64_t a, uint64_t b) {
    char tmp[48], digits[24];
    int n = 0, d = 0;
    do { digits[d++] = static_cast<char>('0' + (a % 10)); a /= 10; } while (a);
    while (d) tmp[n++] = digits[--d];
    tmp[n++] = '\t';
    do { digits[d++] = static_cast<char>('0' + (b % 10)); b /= 10; } while (b);
    while (d) tmp[n++] = digits[--d];
    tmp[n++] = '\n';
    out.append(tmp, n);
}

}  // namespace

int translatecluster(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    const std::string inTsv = par.db1;
    const std::string mapFile = par.db2;
    const std::string outTsv = par.db3;

    const size_t mapBytes = FileUtil::getFileSize(mapFile);
    if (mapBytes % sizeof(uint64_t) != 0) {
        Debug(Debug::ERROR) << mapFile << " is not a whole number of keys\n";
        EXIT(EXIT_FAILURE);
    }
    const uint64_t subCount = mapBytes / sizeof(uint64_t);
    // Refused rather than divided by. span would be 0 for an empty map, and the
    // `key / span` that buckets every row is the first thing this does.
    if (subCount == 0) {
        Debug(Debug::ERROR) << "Key map " << mapFile << " is empty, so there is nothing to "
                            << "translate into. createrepdb did not finish.\n";
        EXIT(EXIT_FAILURE);
    }

    // Sized against the *rows* as well as the keys. A bucket holds the 8 B/key map
    // slice over its key range and, alongside it, every row whose key falls in that
    // range at 16 B each. Sizing on keys alone bounded the slice and left the rows
    // unbounded, so the stage overran the limit it was told to honour -- cluster
    // sizes are skewed, so a bucket of `span` keys can hold far more than `span`
    // rows. translatekeys already sizes this way and says so in the same terms.
    const uint64_t targetBytes = std::max<uint64_t>(Util::computeMemory(par.splitMemoryLimit) / 8,
                                                    1ULL * 1024 * 1024);
    const uint64_t inBytes = FileUtil::getFileSize(inTsv);
    const uint64_t rowsUpperBound = inBytes / 4 + 1;
    unsigned int buckets = 1;
    while (buckets < 65536
           && ((subCount / buckets) * sizeof(uint64_t) > targetBytes
               || (rowsUpperBound / buckets) * sizeof(Pair) > targetBytes)) {
        buckets *= 2;
    }
    buckets = StageParts::raiseBucketsForWorkers(buckets, subCount, par.workerCount);

    unsigned int chunkCount = StageParts::deriveChunkCount(par.workerCount, inBytes, buckets);

    const std::string coordDir = outTsv + ".coord";
    if (FileUtil::directoryExists(coordDir.c_str()) == false) {
        FileUtil::makeDir(coordDir.c_str());
    }
    // Decided once and reused by every later worker, so a resume with a different
    // --workers reads the shards that are actually on disk.
    const StageParts::Layout layout = StageParts::publishLayout(coordDir, buckets, chunkCount);
    buckets = layout.buckets;
    chunkCount = layout.chunks;
    const uint64_t span = (subCount + buckets - 1) / buckets;

    // The remap writes into every destination bucket, so its work item is a band of
    // source buckets rather than a single one: one item per source bucket would put
    // buckets^2 shards on the filesystem.
    //
    // Derived from the recorded layout, not from --workers. Both counts it uses
    // come back from publishLayout, so a run resumed with a different --workers
    // bands the buckets exactly as the run that wrote the shards did; deriving it
    // afresh would have read a band's shards against a different band count, which
    // is the silent-wrong-output case the layout exists to prevent.
    const unsigned int bandCount = std::min(buckets, chunkCount);
    const unsigned int bucketsPerBand = (buckets + bandCount - 1) / bandCount;

    SharedCounter workerCounter(coordDir + "/worker.counter");
    const int64_t workerId = workerCounter.fetchAdd();

    Debug(Debug::INFO) << "Worker " << workerId << " joined; translating " << subCount
                       << " sub-keys over " << buckets << " buckets of " << span
                       << ", scattering in " << chunkCount << " chunks, remapping in "
                       << bandCount << " bands\n";

    const std::string tmpA = outTsv + ".bymember";
    const std::string tmpB = outTsv + ".byrep";
    StageParts::ensureShardDirs(tmpA, buckets);
    StageParts::ensureShardDirs(tmpB, buckets);

    // --- scatter: bucket the input by the member sub-key -----------------------
    {
        WorkQueue queue(coordDir + "/scatter.queue", static_cast<int64_t>(chunkCount));
        const bool finished = queue.drain(workerId, [&](size_t item) {
            const unsigned int chunk = static_cast<unsigned int>(item);
            const StageParts::ChunkRange range =
                StageParts::chunkOf(inTsv, inBytes, chunk, chunkCount);
            Buckets byMember(tmpA, buckets, chunk, workerId);
            if (range.end <= range.begin) {
                return;
            }
            FILE *f = FileUtil::openFileOrDie(inTsv.c_str(), "r", true);
            if (fseeko(f, static_cast<off_t>(range.begin), SEEK_SET) != 0) {
                Debug(Debug::ERROR) << "Cannot seek " << inTsv << "\n";
                EXIT(EXIT_FAILURE);
            }
            uint64_t at = range.begin;
            char *line = NULL;
            size_t cap = 0;
            ssize_t len;
            while (at < range.end && (len = getline(&line, &cap, f)) > 0) {
                at += static_cast<uint64_t>(len);
                char *tab = strchr(line, '\t');
                if (tab == NULL) continue;
                const uint64_t subRep = strtoull(line, NULL, 10);
                const uint64_t subMember = strtoull(tab + 1, NULL, 10);
                // Validated before it is divided into a bucket index, not after.
                // Both columns index `buffers[b]` immediately, so an out-of-range key
                // -- a clustering paired with the wrong key map -- wrote past the end
                // of the bucket vector before the range check further down could
                // report it.
                if (subRep >= subCount || subMember >= subCount) {
                    Debug(Debug::ERROR) << "Clustering " << inTsv << " names sub-key "
                                        << std::max(subRep, subMember) << ", beyond the "
                                        << subCount << " keys in " << mapFile
                                        << ". The clustering and the key map are from "
                                        << "different runs.\n";
                    EXIT(EXIT_FAILURE);
                }
                byMember.append(static_cast<unsigned int>(subMember / span), subMember, subRep);
            }
            free(line);
            fclose(f);
            byMember.close();
        }, StageParts::POLL_MILLIS);
        if (finished == false) {
            Debug(Debug::ERROR) << "The scatter phase stalled with work outstanding\n";
            EXIT(EXIT_FAILURE);
        }
    }

    // --- remap: translate the member, then re-scatter by representative --------
    // The shard index here is the *source* bucket, so byRep bucket d receives its
    // rows in ascending source-bucket order however the work was distributed --
    // which is the order the single-threaded loop produced them in.
    const std::vector<int64_t> scatterAuthority =
        StageParts::authoritativeWorkers(coordDir + "/scatter.queue", chunkCount);
    {
        WorkQueue queue(coordDir + "/remap.queue", static_cast<int64_t>(bandCount));
        const bool finished = queue.drain(workerId, [&](size_t item) {
            const unsigned int band = static_cast<unsigned int>(item);
            Buckets byRep(tmpB, buckets, band, workerId);
            const unsigned int firstBucket = band * bucketsPerBand;
            const unsigned int lastBucket = std::min(firstBucket + bucketsPerBand, buckets);
            for (unsigned int b = firstBucket; b < lastBucket; b++) {
            const uint64_t lo = static_cast<uint64_t>(b) * span;
            if (lo >= subCount) {
                continue;
            }
            const uint64_t hi = std::min(lo + span, subCount);
            const int mapFd = open(mapFile.c_str(), O_RDONLY);
            if (mapFd < 0) {
                Debug(Debug::ERROR) << "Cannot open " << mapFile << ": " << strerror(errno) << "\n";
                EXIT(EXIT_FAILURE);
            }
            const std::vector<uint64_t> slice = mapSlice(mapFd, lo, hi);
            close(mapFd);
            const std::vector<Pair> pairs = Buckets::read(tmpA, b, scatterAuthority);
            for (size_t i = 0; i < pairs.size(); i++) {
                if (pairs[i].key < lo || pairs[i].key - lo >= slice.size()) {
                    Debug(Debug::ERROR) << "Sub-key " << pairs[i].key << " is outside the key map "
                                        << mapFile << ", which holds " << subCount << " keys\n";
                    EXIT(EXIT_FAILURE);
                }
                const uint64_t origMember = slice[static_cast<size_t>(pairs[i].key - lo)];
                // Checked again here: `other` is the column this pass carries
                // through untouched, so it was validated on the way in but has
                // been through a spill file since, and it is about to index
                // byRep's bucket vector.
                if (pairs[i].other >= subCount) {
                    Debug(Debug::ERROR) << "Spill bucket " << b << " holds sub-key "
                                        << pairs[i].other << ", beyond the " << subCount
                                        << " keys in " << mapFile << "\n";
                    EXIT(EXIT_FAILURE);
                }
                byRep.append(static_cast<unsigned int>(pairs[i].other / span), pairs[i].other,
                             origMember);
            }
            }
            byRep.close();
        }, StageParts::POLL_MILLIS);
        if (finished == false) {
            Debug(Debug::ERROR) << "The remap phase stalled with work outstanding\n";
            EXIT(EXIT_FAILURE);
        }
    }

    // --- emit: translate the representative into this bucket's part ------------
    const std::vector<int64_t> remapAuthority =
        StageParts::authoritativeWorkers(coordDir + "/remap.queue", bandCount);
    {
        WorkQueue queue(coordDir + "/emit.queue", static_cast<int64_t>(buckets));
        // Said out loud when a worker cannot possibly get work. The unit counts
        // come from --split-memory-limit and --workers, not from how many workers
        // turn up, so an allocation larger than the unit count leaves the surplus
        // with nothing to claim -- they sleep in the backoff and exit 0, and the
        // stage looks like it succeeded on every node while running on a few.
        // alignparallel and kmerreduceparallel both say this; these did not.
        if (workerId >= queue.getItemCount()) {
            Debug(Debug::WARNING)
                << "Worker " << workerId << " has no work: this phase was sized for "
                << queue.getItemCount() << " work unit(s) and cannot use more workers than "
                << "that. Pass --workers <n> so the unit count is sized for the allocation.\n";
        }
        const bool finished = queue.drain(workerId, [&](size_t item) {
            const unsigned int b = static_cast<unsigned int>(item);
            const std::string finalPath = StageParts::partPath(outTsv, b);
            // The worker id, not getpid(): pids are unique per node and this path
            // is on a shared filesystem, so two workers on two nodes drawing the
            // same pid would truncate each other's temporary.
            const std::string tmpPath = finalPath + ".tmp.w" + SSTR(workerId);
            FILE *out = FileUtil::openAndDelete(tmpPath.c_str(), "w");
            const uint64_t lo = static_cast<uint64_t>(b) * span;
            if (lo < subCount) {
                const uint64_t hi = std::min(lo + span, subCount);
                const int mapFd = open(mapFile.c_str(), O_RDONLY);
                if (mapFd < 0) {
                    Debug(Debug::ERROR) << "Cannot open " << mapFile << ": " << strerror(errno)
                                        << "\n";
                    EXIT(EXIT_FAILURE);
                }
                const std::vector<uint64_t> slice = mapSlice(mapFd, lo, hi);
                close(mapFd);
                // Read against the *band* count: the remap phase produced one shard
                // per band, not one per input chunk.
                const std::vector<Pair> pairs = Buckets::read(tmpB, b, remapAuthority);
                std::string buffer;
                buffer.reserve(64 * 1024 * 1024);
                for (size_t i = 0; i < pairs.size(); i++) {
                    if (pairs[i].key < lo || pairs[i].key - lo >= slice.size()) {
                        Debug(Debug::ERROR) << "Sub-key " << pairs[i].key
                                            << " is outside the key map " << mapFile
                                            << ", which holds " << subCount << " keys\n";
                        EXIT(EXIT_FAILURE);
                    }
                    appendPair(buffer, slice[static_cast<size_t>(pairs[i].key - lo)],
                               pairs[i].other);
                    if (buffer.size() > 32 * 1024 * 1024) {
                        writeAllOrDie(buffer.data(), buffer.size(), out, tmpPath);
                        buffer.clear();
                    }
                }
                if (buffer.empty() == false) {
                    writeAllOrDie(buffer.data(), buffer.size(), out, tmpPath);
                }
            }
            if (fclose(out) != 0 || rename(tmpPath.c_str(), finalPath.c_str()) != 0) {
                Debug(Debug::ERROR) << "Cannot publish " << finalPath << ": " << strerror(errno)
                                    << "\n";
                EXIT(EXIT_FAILURE);
            }
        }, StageParts::POLL_MILLIS);
        if (finished == false) {
            Debug(Debug::ERROR) << "The emit phase stalled with work outstanding\n";
            EXIT(EXIT_FAILURE);
        }
    }

    // --- assemble: parts into one file -----------------------------------------
    // Opened before the offsets are read, because whether the parts still have to
    // exist is exactly what this queue answers. Once assemble has drained, the
    // cleanup phase deletes them; a run resumed past that point finds no parts and
    // must not treat it as damage.
    {
        WorkQueue queue(coordDir + "/assemble.queue", static_cast<int64_t>(buckets));
        std::vector<uint64_t> offset;
        bool haveParts = false;
        if (queue.allDone() == false) {
            haveParts = StageParts::partOffsets(outTsv, buckets, offset);
            // Re-checked rather than believed: the allDone() above is a snapshot,
            // and the scan takes long enough at scale that the last assemble item
            // can finish inside it, after which cleanup starts deleting the parts
            // being scanned. Only a queue still unfinished makes that damage.
            if (haveParts == false && queue.allDone() == false) {
                Debug(Debug::ERROR) << "A part is missing beside " << outTsv
                                    << " while the assemble queue still has work outstanding. "
                                    << "Remove the working directory and re-run the stage.\n";
                EXIT(EXIT_FAILURE);
            }
        }
        if (haveParts) {
            StageParts::createAssembled(outTsv, offset[buckets], coordDir);
            const bool finished = queue.drain(workerId, [&](size_t item) {
                const unsigned int b = static_cast<unsigned int>(item);
                StageParts::assemblePart(outTsv, outTsv, b, offset[b]);
            }, StageParts::POLL_MILLIS);
            if (finished == false) {
                Debug(Debug::ERROR) << "The assemble phase stalled with work outstanding\n";
                EXIT(EXIT_FAILURE);
            }
        }
    }

    // --- cleanup ---------------------------------------------------------------
    // Its own phase, after assemble has drained: dropping a bucket's spill as soon
    // as it was consumed would be cheaper on scratch, but a worker whose lease had
    // expired could still be reading those shards and would then publish a
    // truncated part over a complete one. unlink direct and its result ignored,
    // because a re-claimed item legitimately finds the files already gone.
    {
        WorkQueue queue(coordDir + "/cleanup.queue", static_cast<int64_t>(buckets));
        const bool finished = queue.drain(workerId, [&](size_t item) {
            const unsigned int b = static_cast<unsigned int>(item);
            // The whole directory, walked: an attempt abandoned mid-item leaves a
            // shard under its own worker id that no derived name reaches.
            StageParts::removeShardDir(tmpA, b);
            StageParts::removeShardDir(tmpB, b);
            unlink(StageParts::partPath(outTsv, b).c_str());
        }, StageParts::POLL_MILLIS);
        if (finished == false) {
            Debug(Debug::ERROR) << "The cleanup phase stalled with work outstanding\n";
            EXIT(EXIT_FAILURE);
        }
    }
    StageParts::removeShardRoot(tmpA);
    StageParts::removeShardRoot(tmpB);

    Debug(Debug::INFO) << "Worker " << workerId << " finished " << outTsv << "\n";
    return EXIT_SUCCESS;
}
