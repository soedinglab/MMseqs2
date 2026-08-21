#ifndef MMSEQS_STAGEPARTS_H
#define MMSEQS_STAGEPARTS_H

#include <cstdint>
#include <string>
#include <vector>

// Shared machinery for the single-node stages that were made multi-worker.
//
// mergeclusterparallel and translatecluster are the same shape: scatter a
// line-oriented input into key buckets, do something per bucket, then emit one
// text file. Making them multi-worker needs the same pieces, and they are here
// rather than copied because getting any of them subtly different between stages
// is exactly the kind of defect that shows up as a clustering quietly missing a
// key range.
//
//   chunking    a line-oriented input is cut into byte ranges whose boundaries
//               are resolved to line starts, so chunks are independent and cover
//               the input exactly once.
//   shards      a scatter writes one file per (bucket, work item, worker). The
//               work item is in the name so a reader can take only the attempt
//               that the queue recorded as finished; the worker is in the name so
//               two attempts at one item never share a path.
//   part sets   each bucket's output goes to its own file; the files are then
//               concatenated at prefix-sum offsets by pwrite, so the assembly is
//               parallel rather than one worker copying the whole result.
//   sizing      work-unit counts are raised for the allocation but bounded, so
//               the file count does not grow with the square of the allocation.
namespace StageParts {

// First idle wait in these stages' queues, in milliseconds.
//
// Not lower, and this is the reason: WorkQueue::drain resets its backoff whenever
// *any* worker completes *any* item, so while a phase is making progress the
// doubling never engages and every idle worker polls at exactly this rate. Each
// pass takes the one cluster-wide fcntl lock three times and fsyncs once, so the
// cost is (idle workers / interval) lock round trips per second against a single
// file -- and it peaks at a phase tail, precisely when the last holders need that
// lock to renew their leases. A lapsed lease there means an item runs twice.
//
// 100 ms was measured on four workers and is wrong by two orders of magnitude at
// a thousand. One second keeps the straggler wake-up well inside the noise for a
// phase that runs for minutes, which is the case that matters at scale.
extern const unsigned int POLL_MILLIS;

// Contiguous byte range of one chunk, boundaries already resolved to line starts.
struct ChunkRange {
    uint64_t begin;
    uint64_t end;
};

// The byte offset of the first line that starts at or after `from`.
//
// Offset 0 is already a line start. Anywhere else, the bytes up to and including
// the next newline belong to the *previous* chunk, whose owner reads past its own
// end to finish that line. That is what makes a chunk boundary unambiguous: every
// line is owned by exactly the chunk its first byte falls in, so no line is split
// between chunks or claimed by two of them.
uint64_t resolveLineStart(const std::string &path, uint64_t from, uint64_t size);

// The chunk'th of chunkCount equal byte spans of `path`, snapped to line starts.
ChunkRange chunkOf(const std::string &path, uint64_t size, unsigned int chunk,
                   unsigned int chunkCount);

// How many work items to cut a scatter input into.
//
// Two floors and a ceiling. A work item costs two fsyncs of queue state, so an
// item covering less than TARGET bytes spends more on bookkeeping than on work;
// an allocation with fewer items than workers leaves workers idle; and -- the one
// that decides this at scale -- a scatter writes one file per (bucket, item), so
// the item count multiplied by the bucket count is a file count on a shared
// filesystem. Both counts otherwise rise with the worker count, which makes the
// file count rise with its square: at 1e11 on 4096 workers the unbounded form
// derived 8192 buckets and 16384 chunks, i.e. 1.3e8 files per side, which no
// parallel filesystem will create, walk or unlink in reasonable time.
//
// So the product is capped explicitly and the item count gives way, because
// buckets are bounded by memory and cannot.
unsigned int deriveChunkCount(int workerCount, uint64_t inputBytes, unsigned int buckets);

// Raises a memory-derived bucket count so the allocation has work to claim.
//
// A bucket is the unit of work for everything after the scatter, so the count is a
// ceiling on how many workers those phases can use. Deriving it from the memory
// limit alone is what the k-mer partitioning already learned not to do: at 1e6
// keys with a 400 GB limit it settles on a single bucket and every worker but one
// sits polling an empty queue.
//
// Only ever raises the count, so the memory sizing still bounds what one worker
// holds, and never so far that a bucket is not worth the two fsyncs it costs.
unsigned int raiseBucketsForWorkers(unsigned int buckets, uint64_t keyCount, int workerCount);

// The work-unit layout a stage runs against.
struct Layout {
    unsigned int buckets;
    unsigned int chunks;
};

// Records the layout the first worker derived, and returns whatever is recorded.
//
// The counts are derived from --split-memory-limit and --workers, so a run
// resumed with either changed would derive different ones -- and then read shard
// 5 of 8 for a spill written as 1 of 4, or size an output against bucket
// boundaries the parts were never written to. Both are silent: the stage reports
// success over a clustering quietly missing a key range.
//
// So the layout is decided once, published under a lock next to the queues, and
// every later worker takes it from there. This is the same guard, for the same
// reason, that alignparallel writes align.info for.
Layout publishLayout(const std::string &coordDir, unsigned int buckets, unsigned int chunks);

// One directory per bucket, holding that bucket's shards.
//
// Not one flat directory: at the scales this pipeline exists for a scatter writes
// millions of shards, and a single directory holding them is where a parallel
// filesystem's metadata server stops being usable. createdbparallel splits its
// coordination files the same way and for the same reason.
std::string shardDir(const std::string &prefix, unsigned int bucket);

// <prefix>/b<bucket>/i<item>.w<worker>
//
// The item is in the name because the queue records which attempt at an item
// finished, so a reader can take that one and ignore an abandoned attempt's
// leftovers. The worker is in the name because two attempts must never share a
// path: the writers append across flushes, so a second attempt truncating the
// first one's file mid-flight leaves a whole number of records that passes every
// downstream check while having lost or duplicated rows.
std::string shardPath(const std::string &prefix, unsigned int bucket, uint64_t item,
                      int64_t worker);

// Creates the per-bucket directories. Safe to call from every worker.
void ensureShardDirs(const std::string &prefix, unsigned int buckets);

// Which worker's attempt at each item counts, from the queue that produced them.
//
// Exits if the queue is unreadable: it is the only record of which shards are
// real, and guessing would mean reading an abandoned attempt's output.
std::vector<int64_t> authoritativeWorkers(const std::string &queuePath, uint64_t itemCount);

// Removes a bucket's whole shard directory, including any abandoned attempt's
// leftovers, which is why this walks the directory rather than unlinking names it
// derives.
void removeShardDir(const std::string &prefix, unsigned int bucket);

// Removes the directory the per-bucket ones live in, once they are gone. Silent
// if it is not empty: whichever worker gets there last succeeds, and a directory
// left behind costs an inode, not correctness.
void removeShardRoot(const std::string &prefix);

// Where a bucket's output is staged before assembly.
std::string partPath(const std::string &prefix, unsigned int bucket);

// Byte offset of each part in the assembled output, plus the total at [buckets].
//
// Returns false, rather than exiting, when a part is gone or cannot be stat'd.
// A missing part is not necessarily corruption: the cleanup phase deletes them
// once they are in the output, so a run resumed after the stage had finished
// legitimately finds none. The caller distinguishes the two by asking whether the
// assemble queue is complete.
bool partOffsets(const std::string &prefix, unsigned int buckets,
                 std::vector<uint64_t> &offsets);

// Creates `out`, or resizes an existing one, to exactly totalBytes.
//
// The truncation is the point, not the creation. A previous attempt can leave a
// longer file -- a different layout, or a run that failed after assembling -- and
// the workers only pwrite the prefix they know about, so without this the stale
// tail beyond totalBytes is published as part of the result. The single-writer
// loop this replaced got truncate-on-open for free.
void createAssembled(const std::string &out, uint64_t totalBytes, const std::string &coordDir);

// Copies one part into its place in the output.
//
// pwrite rather than an append: the offset is known from the prefix sum, so parts
// go in concurrently and the concatenation costs no single-node pass over a result
// that is hundreds of gigabytes at the scale this pipeline exists for.
void assemblePart(const std::string &out, const std::string &prefix, unsigned int bucket,
                  uint64_t offset);

}  // namespace StageParts

#endif
