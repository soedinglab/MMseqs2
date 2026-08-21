#include "StageParts.h"

#include "Debug.h"
#include "FileUtil.h"
#include "ParallelCoordination.h"
#include "Util.h"

#include <algorithm>
#include <cerrno>
#include <cstring>

#include <dirent.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

namespace StageParts {

const unsigned int POLL_MILLIS = 1000;

uint64_t resolveLineStart(const std::string &path, uint64_t from, uint64_t size) {
    if (from == 0 || from >= size) {
        return std::min(from, size);
    }
    const int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        Debug(Debug::ERROR) << "Cannot open " << path << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    // Read forward in blocks rather than a byte at a time: at 1e10 a scatter input
    // is hundreds of gigabytes and every worker resolves two boundaries per chunk
    // it claims, so a getc loop over a long line is real time on a shared
    // filesystem.
    char buf[64 * 1024];
    uint64_t at = from;
    while (at < size) {
        const size_t want = static_cast<size_t>(std::min<uint64_t>(sizeof(buf), size - at));
        const ssize_t got = pread(fd, buf, want, static_cast<off_t>(at));
        if (got < 0) {
            if (errno == EINTR) {
                continue;
            }
            Debug(Debug::ERROR) << "Cannot read " << path << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (got == 0) {
            break;
        }
        const char *nl = static_cast<const char *>(memchr(buf, '\n', static_cast<size_t>(got)));
        if (nl != NULL) {
            at += static_cast<uint64_t>(nl - buf) + 1;
            close(fd);
            return std::min(at, size);
        }
        at += static_cast<uint64_t>(got);
    }
    close(fd);
    return std::min(at, size);
}

ChunkRange chunkOf(const std::string &path, uint64_t size, unsigned int chunk,
                   unsigned int chunkCount) {
    const uint64_t span = (size + chunkCount - 1) / chunkCount;
    ChunkRange r;
    r.begin = resolveLineStart(path, std::min(span * chunk, size), size);
    r.end = resolveLineStart(path, std::min(span * (chunk + 1), size), size);
    if (r.end < r.begin) {
        r.end = r.begin;
    }
    return r;
}

unsigned int deriveChunkCount(int workerCount, uint64_t inputBytes, unsigned int buckets) {
    // A chunk is a work item, and an item costs two fsyncs of queue state, so it
    // should cover enough input to make that irrelevant.
    const uint64_t TARGET_CHUNK_BYTES = 1024ULL * 1024 * 1024;
    // Shards a scatter may create in total, across every bucket. Chosen so that
    // even at the largest bucket counts the per-bucket directories stay in the
    // tens of thousands of entries, which a parallel filesystem handles; the
    // unbounded form reached 1e8 and does not.
    const uint64_t MAX_SHARD_FILES = 4ULL * 1000 * 1000;

    uint64_t chunks = inputBytes / TARGET_CHUNK_BYTES;
    // At least a few per worker, so no worker sits idle through the scatter...
    const uint64_t perWorker = static_cast<uint64_t>(std::max(1, workerCount)) * 4;
    chunks = std::max(chunks, perWorker);
    // ...but never so many that buckets x chunks becomes unmanageable. Buckets are
    // bounded by memory and cannot give way, so the chunk count does.
    const uint64_t byFiles = std::max<uint64_t>(1, MAX_SHARD_FILES / std::max(1u, buckets));
    chunks = std::min(chunks, byFiles);
    chunks = std::max<uint64_t>(chunks, 1);
    chunks = std::min<uint64_t>(chunks, 65536);
    return static_cast<unsigned int>(chunks);
}

unsigned int raiseBucketsForWorkers(unsigned int buckets, uint64_t keyCount, int workerCount) {
    if (workerCount <= 0) {
        return buckets;
    }
    // Two per worker so the makespan is not set by whichever worker draws the
    // densest bucket.
    const uint64_t BUCKETS_PER_WORKER = 2;
    // ...but never so many that a bucket is not worth claiming. A work item costs
    // two fsyncs of queue state, and these stages run at roughly a million keys a
    // second per worker, so a bucket of a million keys spends ~1 s working against
    // ~10 ms of bookkeeping. 65536 was the first value here and inverted that: at
    // 1e6 keys on four workers it turned a 0.26 s stage into 1.0-3.2 s.
    const uint64_t MIN_USEFUL_BUCKET_KEYS = 1000000;
    const uint64_t wanted = static_cast<uint64_t>(workerCount) * BUCKETS_PER_WORKER;
    unsigned int floorB = 1;
    while (floorB < wanted && floorB < 65536) {
        floorB *= 2;
    }
    while (floorB > 1 && keyCount / floorB < MIN_USEFUL_BUCKET_KEYS) {
        floorB /= 2;
    }
    return std::max(buckets, floorB);
}

Layout publishLayout(const std::string &coordDir, unsigned int buckets, unsigned int chunks) {
    Layout out;
    out.buckets = buckets;
    out.chunks = chunks;
    const std::string path = coordDir + "/layout.info";

    FileLock lock(coordDir + "/layout.lock");
    lock.lock();
    if (FileUtil::fileExists(path.c_str())) {
        FILE *f = fopen(path.c_str(), "r");
        if (f == NULL) {
            Debug(Debug::ERROR) << "Cannot read " << path << ": " << strerror(errno) << "\n";
            lock.unlock();
            EXIT(EXIT_FAILURE);
        }
        unsigned int b = 0, c = 0;
        if (fscanf(f, "buckets %u chunks %u", &b, &c) != 2 || b == 0 || c == 0) {
            Debug(Debug::ERROR) << path << " is unreadable. It records the work-unit layout this "
                                << "run's spill files were written against; remove the working "
                                << "directory and start the stage again.\n";
            fclose(f);
            lock.unlock();
            EXIT(EXIT_FAILURE);
        }
        fclose(f);
        if (b != buckets || c != chunks) {
            // Said out loud rather than silently obeyed: the difference means the
            // resume was given a different --workers or --split-memory-limit, and
            // the recorded layout is the one that matches what is on disk.
            Debug(Debug::WARNING)
                << "Resuming against the recorded layout of " << b << " buckets and " << c
                << " chunks, not the " << buckets << " and " << chunks
                << " this invocation derived. The spill files were written against the "
                << "recorded one.\n";
        }
        out.buckets = b;
        out.chunks = c;
    } else {
        // Written to a temporary and renamed, so a reader never sees a half-written
        // layout even though the lock already makes that unreachable here. The
        // temporary is keyed on the path, which the lock already serialises.
        const std::string tmp = path + ".tmp";
        FILE *f = FileUtil::openAndDelete(tmp.c_str(), "w");
        fprintf(f, "buckets %u chunks %u\n", buckets, chunks);
        if (fclose(f) != 0 || rename(tmp.c_str(), path.c_str()) != 0) {
            Debug(Debug::ERROR) << "Cannot publish " << path << ": " << strerror(errno) << "\n";
            lock.unlock();
            EXIT(EXIT_FAILURE);
        }
    }
    lock.unlock();
    return out;
}

std::string shardDir(const std::string &prefix, unsigned int bucket) {
    return prefix + "/b" + SSTR(bucket);
}

std::string shardPath(const std::string &prefix, unsigned int bucket, uint64_t item,
                      int64_t worker) {
    return shardDir(prefix, bucket) + "/i" + SSTR(item) + ".w" + SSTR(worker);
}

void ensureShardDirs(const std::string &prefix, unsigned int buckets) {
    // EEXIST is the normal case: every worker calls this.
    if (mkdir(prefix.c_str(), 0755) != 0 && errno != EEXIST) {
        Debug(Debug::ERROR) << "Cannot create " << prefix << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    for (unsigned int b = 0; b < buckets; b++) {
        const std::string d = shardDir(prefix, b);
        if (mkdir(d.c_str(), 0755) != 0 && errno != EEXIST) {
            Debug(Debug::ERROR) << "Cannot create " << d << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
}

std::vector<int64_t> authoritativeWorkers(const std::string &queuePath, uint64_t itemCount) {
    std::vector<int64_t> workers;
    if (WorkQueue::readCompletedWorkers(queuePath, workers) == false) {
        Debug(Debug::ERROR) << "Cannot read " << queuePath << ", which records which attempt at "
                            << "each work item finished. Without it there is no way to tell a "
                            << "finished shard from one an abandoned attempt left behind.\n";
        EXIT(EXIT_FAILURE);
    }
    if (workers.size() < itemCount) {
        Debug(Debug::ERROR) << queuePath << " records " << workers.size() << " items, not the "
                            << itemCount << " this stage was partitioned into.\n";
        EXIT(EXIT_FAILURE);
    }
    return workers;
}

void removeShardDir(const std::string &prefix, unsigned int bucket) {
    const std::string d = shardDir(prefix, bucket);
    DIR *dir = opendir(d.c_str());
    if (dir == NULL) {
        return;
    }
    // Walked rather than unlinked by derived name: an attempt that was abandoned
    // mid-item leaves a shard nobody reads, and only the directory knows it is
    // there. Errors ignored throughout -- a re-claimed cleanup item legitimately
    // finds the files already gone, and a leftover costs scratch, never
    // correctness.
    struct dirent *e;
    while ((e = readdir(dir)) != NULL) {
        if (e->d_name[0] == '.') {
            continue;
        }
        unlink((d + "/" + e->d_name).c_str());
    }
    closedir(dir);
    rmdir(d.c_str());
}

void removeShardRoot(const std::string &prefix) {
    rmdir(prefix.c_str());
}

std::string partPath(const std::string &prefix, unsigned int bucket) {
    return prefix + ".part." + SSTR(bucket);
}

bool partOffsets(const std::string &prefix, unsigned int buckets,
                 std::vector<uint64_t> &offsets) {
    offsets.assign(static_cast<size_t>(buckets) + 1, 0);
    for (unsigned int b = 0; b < buckets; b++) {
        const std::string p = partPath(prefix, b);
        // stat once and use that result, rather than testing existence and then
        // asking for the size: the two calls can disagree, and FileUtil's size
        // helper reports failure as (size_t)-1, which folded straight into the
        // prefix sum would run every later offset backwards and overlap the parts
        // already written.
        struct stat st;
        if (stat(p.c_str(), &st) != 0) {
            return false;
        }
        if (st.st_size < 0) {
            return false;
        }
        offsets[b + 1] = offsets[b] + static_cast<uint64_t>(st.st_size);
    }
    return true;
}

void createAssembled(const std::string &out, uint64_t totalBytes, const std::string &coordDir) {
    FileLock createLock(coordDir + "/assemble.lock");
    createLock.lock();
    const int fd = open(out.c_str(), O_WRONLY | O_CREAT, 0666);
    if (fd < 0) {
        Debug(Debug::ERROR) << "Cannot create " << out << ": " << strerror(errno) << "\n";
        createLock.unlock();
        EXIT(EXIT_FAILURE);
    }
    // Unconditionally, not only when creating. An earlier attempt can leave a
    // longer file and the workers only pwrite the prefix this layout describes,
    // so the tail beyond totalBytes would otherwise be published as result.
    if (ftruncate(fd, static_cast<off_t>(totalBytes)) != 0) {
        Debug(Debug::ERROR) << "Cannot size " << out << " to " << totalBytes << " bytes: "
                            << strerror(errno) << "\n";
        close(fd);
        createLock.unlock();
        EXIT(EXIT_FAILURE);
    }
    close(fd);
    createLock.unlock();
}

void assemblePart(const std::string &out, const std::string &prefix, unsigned int bucket,
                  uint64_t offset) {
    const std::string p = partPath(prefix, bucket);
    struct stat st;
    if (stat(p.c_str(), &st) != 0) {
        Debug(Debug::ERROR) << "Cannot stat " << p << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    const uint64_t bytes = static_cast<uint64_t>(st.st_size);
    if (bytes == 0) {
        return;
    }
    const int in = open(p.c_str(), O_RDONLY);
    if (in < 0) {
        Debug(Debug::ERROR) << "Cannot open " << p << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    const int fd = open(out.c_str(), O_WRONLY);
    if (fd < 0) {
        Debug(Debug::ERROR) << "Cannot open " << out << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    std::vector<char> buf(8 * 1024 * 1024);
    uint64_t done = 0;
    while (done < bytes) {
        const size_t want = static_cast<size_t>(std::min<uint64_t>(buf.size(), bytes - done));
        const ssize_t got = read(in, buf.data(), want);
        if (got <= 0) {
            if (got < 0 && errno == EINTR) {
                continue;
            }
            Debug(Debug::ERROR) << "Cannot read " << p << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        ssize_t put = 0;
        while (put < got) {
            const ssize_t w = pwrite(fd, buf.data() + put, static_cast<size_t>(got - put),
                                     static_cast<off_t>(offset + done + static_cast<uint64_t>(put)));
            if (w <= 0) {
                if (w < 0 && errno == EINTR) {
                    continue;
                }
                Debug(Debug::ERROR) << "Cannot write " << out << ": " << strerror(errno) << "\n";
                EXIT(EXIT_FAILURE);
            }
            put += w;
        }
        done += static_cast<uint64_t>(got);
    }
    close(in);
    // Durable before the work queue records this item DONE, which is the next
    // thing that happens: assemblePart runs as a drain() body, and drain marks the
    // item complete -- fsyncing that record -- once the body returns.
    //
    // Without this the two are the wrong way round. WorkQueue::completeLocked
    // fsyncs the DONE record while these bytes are still only in the page cache,
    // so a node lost in that window comes back with the queue saying assembly
    // finished and the bytes never written. The resume then takes the
    // `queue.allDone()` branch, skips the whole assemble block, and cleanup
    // unlinks the parts -- and since createAssembled ftruncate'd the output to its
    // full length, the missing range reads as NUL bytes. Simulated by punching the
    // hole a lost writeback leaves: the resume exits 0, logs nothing, and leaves
    // 72,594 NUL bytes of 290,376 in the published result, with the parts gone and
    // nothing left to rebuild from.
    //
    // createdbparallel.cpp:789-806 already orders these correctly, fsyncing the
    // database before its finalize sentinel and saying why; this stage was the one
    // that did not.
    //
    // One fsync per bucket, not per write, so the cost is a few thousand at 1e10
    // against a stage that moves terabytes.
    if (fsync(fd) != 0) {
        Debug(Debug::ERROR) << "Cannot flush " << out << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (close(fd) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << out << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
}

}  // namespace StageParts
