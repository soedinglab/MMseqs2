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
 * Peak memory is one bucket's remap array, which the bucket count sets, rather
 * than anything proportional to the database.
 */
#include "Command.h"
#include "Debug.h"
#include "DenseIndex.h"
#include "FileUtil.h"
#include "Parameters.h"
#include "Util.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <string>
#include <vector>

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
    BucketWriter(const std::string &prefix, unsigned int buckets, size_t bufferPairs = 1 << 16)
        : prefix(prefix), buckets(buckets), bufferPairs(bufferPairs), closed(false) {
        files.assign(buckets, NULL);
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
            if (files[b] != NULL) {
                fclose(files[b]);
                files[b] = NULL;
            }
        }
    }

    static std::string path(const std::string &prefix, unsigned int bucket) {
        return prefix + "." + SSTR(bucket);
    }

private:
    void flush(unsigned int bucket) {
        if (buffers[bucket].empty()) {
            return;
        }
        if (files[bucket] == NULL) {
            files[bucket] = fopen(path(prefix, bucket).c_str(), "wb");
            if (files[bucket] == NULL) {
                Debug(Debug::ERROR) << "Cannot open " << path(prefix, bucket) << ": "
                                    << strerror(errno) << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        if (fwrite(buffers[bucket].data(), sizeof(Pair), buffers[bucket].size(), files[bucket]) !=
            buffers[bucket].size()) {
            Debug(Debug::ERROR) << "Cannot write bucket " << bucket << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        buffers[bucket].clear();
    }

    std::string prefix;
    unsigned int buckets;
    size_t bufferPairs;
    std::vector<std::vector<Pair> > buffers;
    std::vector<FILE *> files;
    bool closed;
};

// Streams a "rep<TAB>member" TSV, bucketing on whichever column the join needs.
void partition(const std::string &tsv, BucketWriter &writer, uint64_t bucketSpan, bool byRep) {
    FILE *file = FileUtil::openFileOrDie(tsv.c_str(), "r", true);
    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    while ((len = getline(&line, &cap, file)) > 0) {
        char *tab = strchr(line, '\t');
        if (tab == NULL) {
            continue;
        }
        const uint64_t rep = strtoull(line, NULL, 10);
        const uint64_t member = strtoull(tab + 1, NULL, 10);
        const uint64_t key = byRep ? rep : member;
        const uint64_t value = byRep ? member : rep;
        writer.append(static_cast<unsigned int>(key / bucketSpan), key, value);
    }
    free(line);
    fclose(file);
}

std::vector<Pair> readBucket(const std::string &prefix, unsigned int bucket) {
    std::vector<Pair> out;
    const std::string p = BucketWriter::path(prefix, bucket);
    if (FileUtil::fileExists(p.c_str()) == false) {
        return out;
    }
    const size_t bytes = FileUtil::getFileSize(p);
    if (bytes == 0 || bytes % sizeof(Pair) != 0) {
        return out;
    }
    out.resize(bytes / sizeof(Pair));
    FILE *f = FileUtil::openFileOrDie(p.c_str(), "rb", true);
    if (fread(out.data(), sizeof(Pair), out.size(), f) != out.size()) {
        Debug(Debug::ERROR) << "Cannot read " << p << "\n";
        EXIT(EXIT_FAILURE);
    }
    fclose(f);
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
void compose(const std::string &earlier, const std::string &later, const std::string &out,
             uint64_t entryCount, const std::string &tmpPrefix, uint64_t bucketSpan,
             unsigned int buckets) {
    {
        BucketWriter laterW(tmpPrefix + ".later", buckets);
        BucketWriter earlierW(tmpPrefix + ".earlier", buckets);
        // Bucket the later clustering by the member it reassigns, and the earlier
        // one by its representative: those are the two sides of the join, so they
        // meet in the same bucket.
        partition(later, laterW, bucketSpan, false);
        partition(earlier, earlierW, bucketSpan, true);
    }

    FILE *result = FileUtil::openAndDelete(out.c_str(), "w");
    std::string buffer;
    buffer.reserve(64 * 1024 * 1024);
    std::vector<uint64_t> remap;

    for (unsigned int b = 0; b < buckets; b++) {
        const uint64_t lo = b * bucketSpan;
        if (lo >= entryCount) {
            break;
        }
        const uint64_t hi = std::min(lo + bucketSpan, entryCount);

        // Dense over the key range, which is what dense keys buy: no hash table,
        // no per-key state outside this bucket.
        remap.assign(static_cast<size_t>(hi - lo), INVALID);
        const std::vector<Pair> laterPairs = readBucket(tmpPrefix + ".later", b);
        for (size_t i = 0; i < laterPairs.size(); i++) {
            remap[static_cast<size_t>(laterPairs[i].key - lo)] = laterPairs[i].value;
        }

        const std::vector<Pair> earlierPairs = readBucket(tmpPrefix + ".earlier", b);
        for (size_t i = 0; i < earlierPairs.size(); i++) {
            const uint64_t rep = earlierPairs[i].key;
            const uint64_t mapped = remap[static_cast<size_t>(rep - lo)];
            appendPair(buffer, mapped == INVALID ? rep : mapped, earlierPairs[i].value);
            if (buffer.size() > 32 * 1024 * 1024) {
                fwrite(buffer.data(), 1, buffer.size(), result);
                buffer.clear();
            }
        }
        FileUtil::remove(BucketWriter::path(tmpPrefix + ".later", b).c_str());
        FileUtil::remove(BucketWriter::path(tmpPrefix + ".earlier", b).c_str());
    }
    if (buffer.empty() == false) {
        fwrite(buffer.data(), 1, buffer.size(), result);
    }
    if (fclose(result) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << out << "\n";
        EXIT(EXIT_FAILURE);
    }
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

    // One bucket's remap array is 8 bytes per key in its range; size the buckets so
    // that stays a modest slice of the memory limit.
    // A floor only so a tiny limit does not explode into pathologically many
    // buckets; low enough that the multi-bucket path is reachable in testing,
    // which matters because that is the path used at scale.
    const uint64_t targetBytes = std::max<uint64_t>(Util::computeMemory(par.splitMemoryLimit) / 8,
                                                    1ULL * 1024 * 1024);
    unsigned int buckets = 1;
    while (buckets < 65536 && (info.entryCount / buckets) * sizeof(uint64_t) > targetBytes) {
        buckets *= 2;
    }
    const uint64_t bucketSpan = (info.entryCount + buckets - 1) / buckets;
    Debug(Debug::INFO) << "Merging " << files.size() << " clusterings over " << buckets
                       << " key buckets of " << bucketSpan << " keys\n";

    // Fold left: the accumulated clustering is always the "earlier" side.
    std::string current = files[0];
    std::string scratch = out + ".tmp";
    for (size_t i = 1; i < files.size(); i++) {
        const std::string target = (i + 1 == files.size()) ? out : (scratch + ".step" + SSTR(i));
        compose(current, files[i], target, info.entryCount, out + ".part", bucketSpan, buckets);
        if (i > 1) {
            FileUtil::remove(current.c_str());
        }
        current = target;
    }

    Debug(Debug::INFO) << "Wrote " << out << "\n";
    return EXIT_SUCCESS;
}
