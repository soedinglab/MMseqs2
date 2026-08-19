/*
 * mergeclusterparallel -- composes the clusterings of successive linclust passes.
 *
 * linclust clusters, then re-clusters the representatives, then folds the second
 * result into the first: a sequence assigned to representative r in pass 1 ends
 * up under whatever representative r was given in pass 2. Stock does this in
 * `mergeclusters` with `std::list<size_t>[dbSize]` (mergeclusters.cpp:28) and
 * splices lists between entries. That array is 24 bytes per sequence of *empty
 * list headers* before a single member is stored -- 2.4 TB at 1e11 and 24 TB at
 * 1e12 -- and every member is held in memory besides.
 *
 * The composition is really just a join, `final(x) = later(earlier(x))`, so with
 * dense keys it needs no resident per-key state at all:
 *
 *   1. radix-partition both clusterings into key-range buckets -- the later one by
 *      the member it reassigns, the earlier one by its representative, so the two
 *      sides of the join land in the same bucket;
 *   2. per bucket, the later clustering's mapping is a *dense array over that key
 *      range*, because the keys in a range are contiguous integers;
 *   3. stream the earlier clustering's bucket through that array and emit.
 *
 * Peak memory is one bucket's remap array plus the two sides of its join, all
 * three bounded by the bucket count, rather than anything proportional to the
 * database.
 *
 * Multi-worker
 * ------------
 * All three phases run across every worker of the allocation, coordinating only
 * through files exactly as the map/reduce/align stages do -- same byte-identical
 * command line, same WorkQueue, same crash tolerance:
 *
 *   scatter    one work item per chunk of each input TSV; a chunk is a byte range
 *              whose boundaries are resolved to line starts, so chunks are
 *              independent and cover the file exactly once.
 *   compose    one work item per key bucket, which is what the join was already
 *              partitioned into. Each bucket's result goes to its own part file.
 *   assemble   one work item per part: the parts are concatenated into the output
 *              by pwrite at a prefix-sum offset, so this phase is parallel too
 *              rather than a single-node copy of the whole clustering.
 *
 * Output is byte-identical to a single-worker run, and stays so however the work
 * was distributed. A shard is named by (work item, worker): the item fixes the
 * order readBucket concatenates in -- item order is input order -- and the worker
 * keeps two attempts at one item off a single path, since the writers append and
 * a second attempt truncating the first mid-flight would leave a file that is
 * still a whole number of records, still in range, and silently short. Which of
 * the two attempts counts is not guessed: the queue records it, and readBucket
 * reads only that one.
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
 * Shards live one directory per bucket. A flat directory is where a parallel
 * filesystem's metadata server stops being usable at these counts, and the counts
 * themselves are bounded -- see StageParts::deriveChunkCount -- because buckets
 * and chunks both rise with the allocation and their product is a file count.
 */
#include "Command.h"
#include "Debug.h"
#include "DenseIndex.h"
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

// (key, value) as written to a bucket: for the later clustering the key is the
// member being reassigned, for the earlier one it is the representative.
struct __attribute__((__packed__)) Pair {
    uint64_t key;
    uint64_t value;
};

const uint64_t INVALID = UINT64_MAX;

class BucketWriter {
public:
    // Named by (item, worker), which is what makes a redone item safe.
    //
    // The item keeps the concatenation order -- and so the output -- independent
    // of how the work was distributed. The worker keeps two attempts at one item
    // off the same path: the flushes below append, so a second attempt opening
    // "wb" over a first one still running truncates it mid-flight and the file
    // ends up a whole number of Pairs, in range, silently short.
    BucketWriter(const std::string &prefix, unsigned int buckets, uint64_t item,
                 int64_t worker, size_t bufferPairs = 1 << 16)
        : prefix(prefix), buckets(buckets), item(item), worker(worker),
          bufferPairs(bufferPairs), closed(false) {
        opened.assign(buckets, false);
        this->buffers.resize(buckets);
    }
    ~BucketWriter() { close(); }

    void append(unsigned int bucket, uint64_t key, uint64_t value) {
        Pair p;
        p.key = key;
        p.value = value;
        buffers[bucket].push_back(p);
        if (buffers[bucket].size() >= bufferPairs) {
            flush(bucket);
        }
    }

    void close() {
        if (closed) {
            return;
        }
        closed = true;
        for (unsigned int b = 0; b < buckets; b++) {
            flush(b);
        }
    }

private:
    void flush(unsigned int bucket) {
        if (buffers[bucket].empty()) {
            return;
        }
        // Append-and-close per flush: one descriptor per bucket can reach the
        // derived bucket count, and the close is where a buffered or remote
        // filesystem first reports ENOSPC or a quota. Losing whole records keeps
        // the file record-aligned, so downstream never notices the truncation.
        FILE *file = fopen(StageParts::shardPath(prefix, bucket, item, worker).c_str(),
                            opened[bucket] ? "ab" : "wb");
        if (file == NULL) {
            Debug(Debug::ERROR) << "Cannot open "
                                << StageParts::shardPath(prefix, bucket, item, worker) << ": "
                                << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        opened[bucket] = true;
        if (fwrite(buffers[bucket].data(), sizeof(Pair), buffers[bucket].size(), file) !=
            buffers[bucket].size()) {
            Debug(Debug::ERROR) << "Cannot write bucket " << bucket << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (fclose(file) != 0) {
            Debug(Debug::ERROR) << "Cannot close "
                                << StageParts::shardPath(prefix, bucket, item, worker) << ": "
                                << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        buffers[bucket].clear();
    }

    std::string prefix;
    unsigned int buckets;
    uint64_t item;
    int64_t worker;
    size_t bufferPairs;
    std::vector<std::vector<Pair> > buffers;
    std::vector<bool> opened;
    bool closed;
};

// Streams one chunk of a "rep<TAB>member" TSV, bucketing on whichever column the
// join needs.
//
// Both columns are checked against the key space *before* either is divided into a
// bucket index. `key / bucketSpan` indexes the writer's bucket vector immediately,
// and the value ends up indexing `remap[key - lo]` in compose(), so a key from a
// clustering built over a different database wrote out of bounds long before
// anything downstream could report it. Dense keys make the valid range exactly
// [0, entryCount).
void partitionRange(const std::string &tsv, const StageParts::ChunkRange &range, BucketWriter &writer,
                    uint64_t bucketSpan, bool byRep, uint64_t entryCount) {
    if (range.end <= range.begin) {
        return;
    }
    FILE *file = FileUtil::openFileOrDie(tsv.c_str(), "r", true);
    if (fseeko(file, static_cast<off_t>(range.begin), SEEK_SET) != 0) {
        Debug(Debug::ERROR) << "Cannot seek " << tsv << " to " << range.begin << "\n";
        EXIT(EXIT_FAILURE);
    }
    uint64_t at = range.begin;
    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    while (at < range.end && (len = getline(&line, &cap, file)) > 0) {
        at += static_cast<uint64_t>(len);
        char *tab = strchr(line, '\t');
        if (tab == NULL) {
            continue;
        }
        const uint64_t rep = strtoull(line, NULL, 10);
        const uint64_t member = strtoull(tab + 1, NULL, 10);
        if (rep >= entryCount || member >= entryCount) {
            Debug(Debug::ERROR) << "Clustering " << tsv << " names key "
                                << std::max(rep, member) << ", beyond the " << entryCount
                                << " keys of the database. The clusterings being merged and the "
                                << "database are from different runs.\n";
            EXIT(EXIT_FAILURE);
        }
        const uint64_t key = byRep ? rep : member;
        const uint64_t value = byRep ? member : rep;
        writer.append(static_cast<unsigned int>(key / bucketSpan), key, value);
    }
    free(line);
    fclose(file);
}

// A bucket's shards, concatenated in work-item order.
//
// Item order is input order, so this reproduces exactly the sequence a single scan
// of the input would have written into the bucket -- which is what makes the
// merged output independent of how many workers ran and which items they took.
//
// `authority[i]` is the worker whose attempt at item i the queue recorded as
// finished. Reading only that one is what makes a redone item harmless: an
// abandoned attempt leaves its own shard behind, under its own worker id, and
// nothing reads it.
std::vector<Pair> readBucket(const std::string &prefix, unsigned int bucket,
                             const std::vector<int64_t> &authority) {
    std::vector<Pair> out;
    for (size_t item = 0; item < authority.size(); item++) {
        if (authority[item] < 0) {
            // No completion recorded. Only reachable if the phase was not drained,
            // which the caller has already established it was.
            Debug(Debug::ERROR) << "Work item " << item << " of " << prefix
                                << " has no completion record although its phase reported done. "
                                << "Remove the working directory and re-run the merge.\n";
            EXIT(EXIT_FAILURE);
        }
        const std::string p = StageParts::shardPath(prefix, bucket, item, authority[item]);
        if (FileUtil::fileExists(p.c_str()) == false) {
            // An item that produced nothing for this bucket writes no file.
            continue;
        }
        const size_t bytes = FileUtil::getFileSize(p);
        if (bytes == 0) {
            continue;
        }
        if (bytes % sizeof(Pair) != 0) {
            // Not silently skipped. These shards are written by this same command a
            // moment earlier, in whole Pair units, so a partial one means a failed
            // write or a truncated filesystem -- not something the input can cause.
            // Returning empty dropped every cluster in the bucket's key range from
            // the merged output, which is a smaller clustering that still looks
            // valid.
            Debug(Debug::ERROR) << "Shard " << p << " is " << bytes
                                << " bytes, not a whole number of " << sizeof(Pair)
                                << "-byte pairs. It was truncated after this stage wrote it; "
                                << "remove the working directory and re-run the merge.\n";
            EXIT(EXIT_FAILURE);
        }
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

inline void appendPair(std::string &out, uint64_t a, uint64_t b) {
    char tmp[48];
    char digits[24];
    int n = 0, d = 0;
    do { digits[d++] = static_cast<char>('0' + (a % 10)); a /= 10; } while (a);
    while (d) tmp[n++] = digits[--d];
    tmp[n++] = '\t';
    do { digits[d++] = static_cast<char>('0' + (b % 10)); b /= 10; } while (b);
    while (d) tmp[n++] = digits[--d];
    tmp[n++] = '\n';
    out.append(tmp, n);
}


// earlier: member -> representative. later: that representative -> its new one.
// Composes one bucket: build the dense remap over the bucket's key range from the
// later clustering, then stream the earlier one through it.
//
// Written to a per-process temporary and renamed, like every other output in this
// pipeline. Two workers can legitimately compose the same bucket -- a lease
// expires while its holder is merely slow -- and both produce the same bytes, so
// the rename makes the duplicate harmless instead of a torn file.
void composeBucket(const std::string &tmpPrefix, unsigned int b, uint64_t bucketSpan,
                   uint64_t entryCount, const std::vector<int64_t> &laterAuthority,
                   const std::vector<int64_t> &earlierAuthority, int64_t workerId) {
    const std::string finalPath = StageParts::partPath(tmpPrefix, b);
    // The worker id, not getpid(): pids are unique per node and this path is on a
    // shared filesystem, so two workers on two nodes drawing the same pid -- routine
    // across a homogeneous allocation -- would open the same "w" path and truncate
    // each other. CandidateEdge.cpp:123 rejects the same construct for the same
    // reason.
    const std::string tmpPath = finalPath + ".tmp.w" + SSTR(workerId);
    FILE *result = FileUtil::openAndDelete(tmpPath.c_str(), "w");

    const uint64_t lo = static_cast<uint64_t>(b) * bucketSpan;
    // Buckets past the end of the key space still publish an empty part: assemble
    // takes the part set as the authority on what the output holds, and a missing
    // file there cannot be told apart from one nobody has written yet.
    if (lo < entryCount) {
        const uint64_t hi = std::min(lo + bucketSpan, entryCount);

        // Dense over the key range, which is what dense keys buy: no hash table,
        // no per-key state outside this bucket.
        std::vector<uint64_t> remap(static_cast<size_t>(hi - lo), INVALID);
        const std::vector<Pair> laterPairs = readBucket(tmpPrefix + ".later", b, laterAuthority);
        for (size_t i = 0; i < laterPairs.size(); i++) {
            // Re-checked on the way out of the spill file, not only on the way in:
            // this is the index into `remap`, and a bucket file written by an
            // earlier attempt with a different bucketSpan would land here holding
            // keys for a different range.
            if (laterPairs[i].key < lo || laterPairs[i].key >= hi) {
                Debug(Debug::ERROR) << "Bucket " << b << " of " << tmpPrefix << ".later holds key "
                                    << laterPairs[i].key << ", outside its range [" << lo << ", "
                                    << hi << "). Remove the working directory and re-run the "
                                    << "merge.\n";
                EXIT(EXIT_FAILURE);
            }
            remap[static_cast<size_t>(laterPairs[i].key - lo)] = laterPairs[i].value;
        }

        std::string buffer;
        buffer.reserve(64 * 1024 * 1024);
        const std::vector<Pair> earlierPairs =
            readBucket(tmpPrefix + ".earlier", b, earlierAuthority);
        for (size_t i = 0; i < earlierPairs.size(); i++) {
            const uint64_t rep = earlierPairs[i].key;
            if (rep < lo || rep >= hi) {
                Debug(Debug::ERROR) << "Bucket " << b << " of " << tmpPrefix
                                    << ".earlier holds key " << rep << ", outside its range ["
                                    << lo << ", " << hi << "). Remove the working directory and "
                                    << "re-run the merge.\n";
                EXIT(EXIT_FAILURE);
            }
            const uint64_t mapped = remap[static_cast<size_t>(rep - lo)];
            appendPair(buffer, mapped == INVALID ? rep : mapped, earlierPairs[i].value);
            if (buffer.size() > 32 * 1024 * 1024) {
                writeAllOrDie(buffer.data(), buffer.size(), result, tmpPath);
                buffer.clear();
            }
        }
        if (buffer.empty() == false) {
            writeAllOrDie(buffer.data(), buffer.size(), result, tmpPath);
        }
    }

    if (fclose(result) != 0 || rename(tmpPath.c_str(), finalPath.c_str()) != 0) {
        Debug(Debug::ERROR) << "Cannot publish " << finalPath << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
}


// earlier: member -> representative. later: that representative -> its new one.
void compose(const std::string &earlier, const std::string &later, const std::string &out,
             uint64_t entryCount, const std::string &tmpPrefix, uint64_t bucketSpan,
             unsigned int buckets, unsigned int chunkCount, int64_t workerId,
             const std::string &coordDir) {
    const uint64_t laterSize = FileUtil::getFileSize(later);
    const uint64_t earlierSize = FileUtil::getFileSize(earlier);

    // --- scatter: bucket both sides of the join -------------------------------
    // One queue over both inputs, so a worker that finishes the smaller one keeps
    // pulling from the larger instead of idling at a barrier between them.
    StageParts::ensureShardDirs(tmpPrefix + ".later", buckets);
    StageParts::ensureShardDirs(tmpPrefix + ".earlier", buckets);
    {
        WorkQueue queue(coordDir + "/scatter.queue", 2 * static_cast<int64_t>(chunkCount));
        const bool finished = queue.drain(workerId, [&](size_t item) {
            const bool isLater = item < chunkCount;
            const unsigned int chunk = static_cast<unsigned int>(item % chunkCount);
            const std::string &src = isLater ? later : earlier;
            const uint64_t size = isLater ? laterSize : earlierSize;
            BucketWriter writer(tmpPrefix + (isLater ? ".later" : ".earlier"), buckets, chunk,
                                workerId);
            partitionRange(src, StageParts::chunkOf(src, size, chunk, chunkCount), writer, bucketSpan,
                           isLater == false, entryCount);
            writer.close();
        }, StageParts::POLL_MILLIS);
        if (finished == false) {
            Debug(Debug::ERROR) << "The scatter phase stalled with work outstanding\n";
            EXIT(EXIT_FAILURE);
        }
    }

    // --- compose: one bucket at a time ----------------------------------------
    // Split per side: the scatter queue numbers the later clustering's chunks
    // [0, chunkCount) and the earlier one's [chunkCount, 2*chunkCount), so each
    // side's authority is its own slice of the record.
    std::vector<int64_t> laterAuthority, earlierAuthority;
    {
        const std::vector<int64_t> all = StageParts::authoritativeWorkers(
            coordDir + "/scatter.queue", 2 * static_cast<uint64_t>(chunkCount));
        laterAuthority.assign(all.begin(), all.begin() + chunkCount);
        earlierAuthority.assign(all.begin() + chunkCount, all.begin() + 2 * chunkCount);
    }
    {
        WorkQueue queue(coordDir + "/compose.queue", static_cast<int64_t>(buckets));
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
        const bool finished = queue.drain(workerId, [&](size_t bucket) {
            composeBucket(tmpPrefix, static_cast<unsigned int>(bucket), bucketSpan, entryCount,
                          laterAuthority, earlierAuthority, workerId);
        }, StageParts::POLL_MILLIS);
        if (finished == false) {
            Debug(Debug::ERROR) << "The compose phase stalled with work outstanding\n";
            EXIT(EXIT_FAILURE);
        }
    }

    // --- assemble: parts into one file ----------------------------------------
    // Opened before the offsets are read, because whether the parts still have to
    // exist is exactly what this queue answers. Once assemble has drained, the
    // cleanup phase deletes them; a run resumed past that point finds no parts and
    // must not treat it as damage.
    {
        WorkQueue queue(coordDir + "/assemble.queue", static_cast<int64_t>(buckets));
        if (queue.allDone() == false) {
            std::vector<uint64_t> offset;
            if (StageParts::partOffsets(tmpPrefix, buckets, offset) == false) {
                // Re-checked rather than believed. The allDone() above is a
                // snapshot, and scanning the parts takes long enough at scale that
                // the last assemble item can finish inside it -- after which the
                // cleanup phase starts deleting the very parts being scanned. Only
                // a queue that is still unfinished makes a missing part damage.
                if (queue.allDone()) {
                    return;
                }
                Debug(Debug::ERROR) << "A part is missing under " << tmpPrefix
                                    << " while the assemble queue still has work outstanding. "
                                    << "Remove the working directory and re-run the merge.\n";
                EXIT(EXIT_FAILURE);
            }
            StageParts::createAssembled(out, offset[buckets], coordDir);
            const bool finished = queue.drain(workerId, [&](size_t bucket) {
                StageParts::assemblePart(out, tmpPrefix, static_cast<unsigned int>(bucket),
                                         offset[bucket]);
            }, StageParts::POLL_MILLIS);
            if (finished == false) {
                Debug(Debug::ERROR) << "The assemble phase stalled with work outstanding\n";
                EXIT(EXIT_FAILURE);
            }
        }
    }

    // --- cleanup: drop the spill and the parts --------------------------------
    // Its own phase, after assemble has drained, rather than inside compose or
    // assemble. Dropping a bucket's spill as soon as it was composed would be
    // cheaper on scratch, but a worker whose lease had expired could still be
    // reading those shards, and it would then publish a truncated part over a
    // complete one. A queue rather than a loop on one worker so the unlinks --
    // buckets x chunks x 2 of them, which is a lot of metadata traffic at scale --
    // spread over the allocation and are each issued once.
    //
    // unlink direct, and its result ignored: a re-claimed item legitimately finds
    // the files already gone, and a leftover here costs scratch, never
    // correctness.
    {
        WorkQueue queue(coordDir + "/cleanup.queue", static_cast<int64_t>(buckets));
        const bool finished = queue.drain(workerId, [&](size_t bucket) {
            const unsigned int b = static_cast<unsigned int>(bucket);
            // The whole directory, walked: an attempt that was abandoned mid-item
            // leaves a shard under its own worker id that no derived name reaches.
            StageParts::removeShardDir(tmpPrefix + ".later", b);
            StageParts::removeShardDir(tmpPrefix + ".earlier", b);
            unlink(StageParts::partPath(tmpPrefix, b).c_str());
        }, StageParts::POLL_MILLIS);
        if (finished == false) {
            Debug(Debug::ERROR) << "The cleanup phase stalled with work outstanding\n";
            EXIT(EXIT_FAILURE);
        }
    }
    StageParts::removeShardRoot(tmpPrefix + ".later");
    StageParts::removeShardRoot(tmpPrefix + ".earlier");
}

}  // namespace

int mergeclusterparallel(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_VARIADIC, 0);

    std::vector<std::string> files(par.filenames);
    const std::string seqDb = files.front();
    files.erase(files.begin());
    const std::string out = files.back();
    files.pop_back();
    if (files.size() < 2) {
        Debug(Debug::ERROR) << "Need at least two clusterings to merge\n";
        EXIT(EXIT_FAILURE);
    }
    par.printParameters(command.cmd, argc, argv, *command.params);

    const DenseIndex::Info info = DenseIndex::readInfo(seqDb);
    // Refused rather than divided by: bucketSpan would be 0 for an empty database
    // and `key / bucketSpan` is the first thing partition() does with every row.
    if (info.entryCount == 0) {
        Debug(Debug::ERROR) << "Database " << seqDb << " is empty\n";
        EXIT(EXIT_FAILURE);
    }

    // Sized against the *rows* as well as the keys.
    //
    // A bucket holds three things at once: the 8 B/key remap array over its key
    // range, and both sides of the join loaded whole at 16 B per row. Sizing on
    // keys alone -- which is what this did -- bounded the remap array and left the
    // two row vectors unbounded, so the stage overran the very limit it was told
    // to honour. At 1e11 entries on a 512 GB node that settled on 16 buckets: a
    // budgeted 50 GB of remap alongside 100 GB of rows on a perfectly uniform
    // input, and several times that on a bucket carrying hub representatives.
    //
    // Cluster sizes are skewed, so a bucket of `span` keys can hold far more than
    // `span` rows; the row count is bounded from the input files instead, where a
    // row is "0\t0\n" at the very least. translatekeys already sizes this way and
    // says so in the same terms.
    const uint64_t targetBytes = std::max<uint64_t>(Util::computeMemory(par.splitMemoryLimit) / 8,
                                                    1ULL * 1024 * 1024);
    const uint64_t rowsUpperBound =
        (FileUtil::getFileSize(files[0]) + FileUtil::getFileSize(files[1])) / 4 + 1;
    const uint64_t bytesPerRow = sizeof(Pair);
    unsigned int buckets = 1;
    while (buckets < 65536
           && ((info.entryCount / buckets) * sizeof(uint64_t) > targetBytes
               || (rowsUpperBound / buckets) * bytesPerRow > targetBytes)) {
        buckets *= 2;
    }

    // A bucket is also the unit of work for the compose, assemble and cleanup
    // phases, so the count is a ceiling on how many workers those can use, and the
    // memory limit alone can leave it at 1. Raised for the allocation by the same
    // rule the other converted stages use -- shared rather than repeated, because
    // this file having its own copy is precisely how the two came to disagree:
    // the cost-based floor landed in StageParts and the merge kept over-splitting
    // against the old one.
    buckets = StageParts::raiseBucketsForWorkers(buckets, info.entryCount, par.workerCount);
    const uint64_t bucketSpan = (info.entryCount + buckets - 1) / buckets;

    // Bounded against the bucket count, not derived from the worker count alone:
    // the scatter writes one shard per (bucket, chunk), so letting both rise with
    // the allocation makes the file count rise with its square.
    const unsigned int chunkCount = StageParts::deriveChunkCount(
        par.workerCount, std::max(FileUtil::getFileSize(files[0]),
                                  FileUtil::getFileSize(files[1])), buckets);

    const std::string coordRoot = out + ".coord";
    if (FileUtil::directoryExists(coordRoot.c_str()) == false) {
        FileUtil::makeDir(coordRoot.c_str());
    }
    SharedCounter workerCounter(coordRoot + "/worker.counter");
    const int64_t workerId = workerCounter.fetchAdd();

    Debug(Debug::INFO) << "Worker " << workerId << " joined; merging " << files.size()
                       << " clusterings over " << buckets << " key buckets of " << bucketSpan
                       << " keys, scattering in " << chunkCount << " chunks\n";

    // Fold left: the accumulated clustering is always the "earlier" side.
    std::string current = files[0];
    std::string scratch = out + ".tmp";
    for (size_t i = 1; i < files.size(); i++) {
        const std::string target = (i + 1 == files.size()) ? out : (scratch + ".step" + SSTR(i));
        // Per step, so a fold of more than two clusterings does not have its steps
        // claim each other's work items.
        const std::string coordDir = coordRoot + "/step" + SSTR(i);
        if (FileUtil::directoryExists(coordDir.c_str()) == false) {
            FileUtil::makeDir(coordDir.c_str());
        }
        // Decided once per step and reused by every later worker, so a resume with
        // a different --workers reads the shards that are actually on disk.
        const StageParts::Layout layout =
            StageParts::publishLayout(coordDir, buckets, chunkCount);
        const uint64_t stepSpan = (info.entryCount + layout.buckets - 1) / layout.buckets;
        compose(current, files[i], target, info.entryCount, out + ".part", stepSpan,
                layout.buckets, layout.chunks, workerId, coordDir);
        // Safe for any worker: compose() does not return until its cleanup queue
        // has drained, which every worker of the step has therefore left behind.
        // unlink rather than FileUtil::remove because whoever gets there second
        // finds it already gone.
        if (i > 1) {
            unlink(current.c_str());
        }
        current = target;
    }

    Debug(Debug::INFO) << "Worker " << workerId << " finished " << out << "\n";
    return EXIT_SUCCESS;
}
