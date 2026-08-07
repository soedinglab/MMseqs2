/*
 * translatekeys -- rewrites a clustering from database keys into accessions.
 *
 * The pipeline works in dense keys throughout, so its final `clusters.tsv` reads
 * `0<TAB>27`. Stock reaches accessions with `createtsv`, which needs a result
 * *database* and builds a resident id->name table; both are exactly the per-key
 * state that cannot exist at 1e11 sequences. This does the same join streaming.
 *
 * The map is `<db>.lookup`, which `createdbparallel` writes in ascending key
 * order. Unlike the fixed-width key map `translatecluster` reads, its lines vary
 * in width, so a key's accession is at no computable offset and cannot be fetched
 * by pread. What replaces that is ordering: buckets are visited in ascending key
 * order and tile the key space, so one *sequential* pass over the lookup serves a
 * whole column -- the cursor only ever moves forward.
 *
 * A clustering has two key columns and only one can be in order at a time, so
 * each gets its own bucketed pass, as in translatecluster:
 *
 *   1. bucket the input by member key                      (fixed-width spill)
 *   2. stream the lookup, translate the member, re-bucket
 *      by representative key                               (variable-width spill)
 *   3. stream the lookup again, translate the representative, emit
 *
 * Two sequential reads of the lookup, two of the input, and peak memory of one
 * bucket -- never the whole map.
 */
#include "Command.h"
#include "Debug.h"
#include "FileUtil.h"
#include "Parameters.h"
#include "Util.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <string>
#include <vector>

#include <dirent.h>

#ifdef OPENMP
#include <omp.h>
#endif

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

// Spill file holding (key, payload) records whose payload is a key too.
struct __attribute__((__packed__)) KeyPair {
    uint64_t key;    // the column being translated in this pass
    uint64_t other;  // the column carried through untouched
};

// One row after its member has been translated: the representative still a key,
// the member already an accession. The member's key is carried alongside purely
// to order the output (see the sort in pass 3).
struct TranslatedRow {
    uint64_t repKey;
    uint64_t memberKey;
    std::string memberAccession;
};

// Removes every spill file belonging to a prefix, regardless of the thread and
// bucket counts the previous attempt used.
//
// Needed because a rerun only truncates the shards it happens to touch: pass 2
// spills under <prefix>.t<thread>.<bucket> and assigns source buckets with
// schedule(dynamic), so a different schedule -- or a different --threads -- leaves
// shards from the earlier attempt on disk, and pass 3 reads every shard it finds.
// Both attempts' rows would then be emitted.
//
// Matched against the exact generated grammar, not merely "starts with the
// prefix". The spill prefix is derived from the *output* path, so with
// clusters.tsv.tmp as the output the prefix is clusters.tsv.tmp.bymember -- and a
// plain prefix compare also removed an unrelated neighbouring file called
// clusters.tsv.tmp.bymember-notes. The three shapes actually generated are
//
//   <prefix>.<bucket>            pass 1's fixed-width spill
//   <prefix>.t<thread>.<bucket>  pass 2's thread-private spill
//   <prefix><bucket>             pass 3's output pieces, whose prefix ends in "p"
//
// so the suffix is a run of segments, each `<digits>` or `t<digits>`, separated by
// dots, with the leading dot optional. Requiring that shape makes the deletion
// exactly as wide as what was created.
bool isGeneratedSpillName(const std::string &name, const std::string &base) {
    if (name.size() <= base.size() || name.compare(0, base.size(), base) != 0) {
        return false;
    }
    size_t at = base.size();
    bool first = true;
    while (at < name.size()) {
        if (name[at] == '.') {
            at++;
        } else if (first == false) {
            return false;  // segments after the first must be dot-separated
        }
        first = false;
        if (at < name.size() && name[at] == 't') {
            at++;
        }
        const size_t digitsFrom = at;
        while (at < name.size() && name[at] >= '0' && name[at] <= '9') {
            at++;
        }
        if (at == digitsFrom) {
            return false;
        }
    }
    return true;
}

void removeSpillFiles(const std::string &prefix) {
    const size_t slash = prefix.find_last_of('/');
    const std::string dir = slash == std::string::npos ? std::string(".") : prefix.substr(0, slash);
    const std::string base = slash == std::string::npos ? prefix : prefix.substr(slash + 1);
    DIR *handle = opendir(dir.c_str());
    if (handle == NULL) {
        return;  // nothing written yet
    }
    struct dirent *entry;
    errno = 0;
    while ((entry = readdir(handle)) != NULL) {
        const std::string name = entry->d_name;
        if (isGeneratedSpillName(name, base)) {
            FileUtil::remove((dir + "/" + name).c_str());
        }
    }
    const int readErr = errno;
    if (readErr != 0) {
        Debug(Debug::ERROR) << "Cannot read " << dir << ": " << strerror(readErr) << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (closedir(handle) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << dir << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
}


bool byRepThenMember(const TranslatedRow &a, const TranslatedRow &b) {
    if (a.repKey != b.repKey) return a.repKey < b.repKey;
    return a.memberKey < b.memberKey;
}

// Spill file holding (key, key, text) records. Needed for the second pass, where
// the column carried through has already become an accession and so varies in
// width.
class TextBuckets {
public:
    TextBuckets(const std::string &prefix, unsigned int count, size_t flushBytes = 8 << 20)
        : prefix(prefix), count(count), flushBytes(flushBytes), closed(false) {
        opened.assign(count, false);
        buffers.resize(count);
    }
    ~TextBuckets() { close(); }

    void append(unsigned int b, uint64_t key, uint64_t other, const char *text, size_t length) {
        std::string &buf = buffers[b];
        const uint32_t len = static_cast<uint32_t>(length);
        buf.append(reinterpret_cast<const char *>(&key), sizeof(key));
        buf.append(reinterpret_cast<const char *>(&other), sizeof(other));
        buf.append(reinterpret_cast<const char *>(&len), sizeof(len));
        buf.append(text, length);
        if (buf.size() >= flushBytes) {
            flush(b);
        }
    }

    void close() {
        if (closed) return;
        closed = true;
        for (unsigned int b = 0; b < count; b++) {
            flush(b);
        }
    }

    static std::string path(const std::string &prefix, unsigned int b) {
        return prefix + "." + SSTR(b);
    }

    // Appends one bucket's rows to out, which the caller accumulates across shards.
    static void read(const std::string &prefix, unsigned int b, std::vector<TranslatedRow> &out) {
        const std::string p = path(prefix, b);
        if (FileUtil::fileExists(p.c_str()) == false) return;
        const size_t bytes = FileUtil::getFileSize(p);
        if (bytes == 0) return;
        std::string blob(bytes, '\0');
        FILE *f = FileUtil::openFileOrDie(p.c_str(), "rb", true);
        if (fread(&blob[0], 1, bytes, f) != bytes) {
            Debug(Debug::ERROR) << "Cannot read " << p << "\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(f);
        size_t at = 0;
        while (at + 2 * sizeof(uint64_t) + sizeof(uint32_t) <= bytes) {
            TranslatedRow row;
            uint32_t len;
            memcpy(&row.repKey, blob.data() + at, sizeof(row.repKey));
            at += sizeof(row.repKey);
            memcpy(&row.memberKey, blob.data() + at, sizeof(row.memberKey));
            at += sizeof(row.memberKey);
            memcpy(&len, blob.data() + at, sizeof(len));
            at += sizeof(len);
            row.memberAccession.assign(blob.data() + at, len);
            at += len;
            out.push_back(row);
        }
    }

private:
    // Append-and-close per flush, as the k-mer and edge bucket writers do. Holding
    // one descriptor per bucket until close needed up to 65536 of them against the
    // 8192 FileUtil::fixRlimitNoFile raises to, so translation hit EMFILE at
    // exactly the scale it exists for. The first flush truncates, later ones
    // append.
    void flush(unsigned int b) {
        if (buffers[b].empty()) return;
        FILE *file = fopen(path(prefix, b).c_str(), opened[b] ? "ab" : "wb");
        if (file == NULL) {
            Debug(Debug::ERROR) << "Cannot open " << path(prefix, b) << ": " << strerror(errno)
                                << "\n";
            EXIT(EXIT_FAILURE);
        }
        opened[b] = true;
        if (fwrite(buffers[b].data(), 1, buffers[b].size(), file) != buffers[b].size()) {
            Debug(Debug::ERROR) << "Cannot write bucket " << b << "\n";
            EXIT(EXIT_FAILURE);
        }
        // Checked: on a buffered or remote filesystem this is where ENOSPC and
        // quota failures first surface, and a lost record here is invisible
        // downstream because whole rows keep the file self-consistent.
        if (fclose(file) != 0) {
            Debug(Debug::ERROR) << "Cannot close " << path(prefix, b) << ": " << strerror(errno)
                                << "\n";
            EXIT(EXIT_FAILURE);
        }
        buffers[b].clear();
    }

    std::string prefix;
    unsigned int count;
    size_t flushBytes;
    std::vector<std::string> buffers;
    std::vector<bool> opened;
    bool closed;
};

// Fixed-width spill, for the first pass where both columns are still keys.
class KeyBuckets {
public:
    KeyBuckets(const std::string &prefix, unsigned int count, size_t bufferPairs = 1 << 16)
        : prefix(prefix), count(count), bufferPairs(bufferPairs), closed(false) {
        opened.assign(count, false);
        buffers.resize(count);
    }
    ~KeyBuckets() { close(); }

    void append(unsigned int b, uint64_t key, uint64_t other) {
        KeyPair p;
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

    static std::string path(const std::string &prefix, unsigned int b) {
        return prefix + "." + SSTR(b);
    }

    static std::vector<KeyPair> read(const std::string &prefix, unsigned int b) {
        std::vector<KeyPair> out;
        const std::string p = path(prefix, b);
        if (FileUtil::fileExists(p.c_str()) == false) return out;
        const size_t bytes = FileUtil::getFileSize(p);
        if (bytes == 0 || bytes % sizeof(KeyPair) != 0) return out;
        out.resize(bytes / sizeof(KeyPair));
        FILE *f = FileUtil::openFileOrDie(p.c_str(), "rb", true);
        if (fread(out.data(), sizeof(KeyPair), out.size(), f) != out.size()) {
            Debug(Debug::ERROR) << "Cannot read " << p << "\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(f);
        return out;
    }

private:
    // Append-and-close per flush; see the note on TextBuckets::flush.
    void flush(unsigned int b) {
        if (buffers[b].empty()) return;
        FILE *file = fopen(path(prefix, b).c_str(), opened[b] ? "ab" : "wb");
        if (file == NULL) {
            Debug(Debug::ERROR) << "Cannot open " << path(prefix, b) << ": " << strerror(errno)
                                << "\n";
            EXIT(EXIT_FAILURE);
        }
        opened[b] = true;
        if (fwrite(buffers[b].data(), sizeof(KeyPair), buffers[b].size(), file) !=
            buffers[b].size()) {
            Debug(Debug::ERROR) << "Cannot write bucket " << b << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (fclose(file) != 0) {
            Debug(Debug::ERROR) << "Cannot close " << path(prefix, b) << ": " << strerror(errno)
                                << "\n";
            EXIT(EXIT_FAILURE);
        }
        buffers[b].clear();
    }

    std::string prefix;
    unsigned int count;
    size_t bufferPairs;
    std::vector<std::vector<KeyPair> > buffers;
    std::vector<bool> opened;
    bool closed;
};

// Forward-only cursor over a `.lookup`.
//
// Keys ascend, and every caller asks for a higher range than the last, so the
// file is read once from start to end however many buckets there are. Gaps are
// tolerated (a key map from a sparse database leaves holes); a hole is only an
// error if a clustering actually refers to that key, which slice()'s caller
// detects by finding an empty accession.
class LookupCursor {
public:
    explicit LookupCursor(const std::string &path, uint64_t startOffset = 0)
        : line(NULL), cap(0), pending(false), pendingKey(0) {
        file = FileUtil::openFileOrDie(path.c_str(), "r", true);
        if (startOffset > 0) {
            fseeko(file, static_cast<off_t>(startOffset), SEEK_SET);
        }
    }
    ~LookupCursor() {
        free(line);
        if (file != NULL) fclose(file);
    }

    // Fills out[0 .. hi-lo) with the accessions of keys lo .. hi-1.
    void slice(uint64_t lo, uint64_t hi, std::vector<std::string> &out) {
        out.clear();
        out.resize(static_cast<size_t>(hi - lo));
        while (true) {
            if (pending == false) {
                if (getline(&line, &cap, file) <= 0) return;
                char *tab = strchr(line, '\t');
                if (tab == NULL) continue;
                pendingKey = strtoull(line, NULL, 10);
                char *end = strchr(tab + 1, '\t');
                pendingAccession.assign(tab + 1, end != NULL ? (size_t)(end - tab - 1)
                                                             : strlen(tab + 1));
                pending = true;
            }
            if (pendingKey >= hi) return;  // belongs to a later bucket, keep it
            if (pendingKey >= lo) {
                out[static_cast<size_t>(pendingKey - lo)].swap(pendingAccession);
            }
            pending = false;  // below lo: already consumed by an earlier bucket
        }
    }

private:
    FILE *file;
    char *line;
    size_t cap;
    bool pending;
    uint64_t pendingKey;
    std::string pendingAccession;
};

// Byte offset in the lookup where each bucket's key range begins.
//
// The lookup is variable-width text, so a bucket cannot seek to its own keys --
// which is why the three passes below used to share one forward cursor and run
// strictly one bucket at a time. One sequential scan recording `buckets` offsets
// (a few KB) makes every bucket independent, and the passes become parallel.
static std::vector<uint64_t> indexLookup(const std::string &path, uint64_t span,
                                         unsigned int buckets) {
    std::vector<uint64_t> at(buckets, UINT64_MAX);
    FILE *f = FileUtil::openFileOrDie(path.c_str(), "r", true);
    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    uint64_t offset = 0;
    while ((len = getline(&line, &cap, f)) > 0) {
        const uint64_t key = strtoull(line, NULL, 10);
        const unsigned int b = static_cast<unsigned int>(key / span);
        if (b < buckets && at[b] == UINT64_MAX) {
            at[b] = offset;
        }
        offset += static_cast<uint64_t>(len);
    }
    free(line);
    fclose(f);
    return at;
}

const std::string &accessionOf(const std::vector<std::string> &slice, uint64_t key, uint64_t lo,
                               const std::string &lookupFile) {
    // Checked before indexing, not after. A clustering naming a key the lookup
    // does not have is exactly the case the message below describes, and reading
    // slice[key - lo] first would go out of bounds before it could be printed.
    if (key < lo || key - lo >= slice.size()) {
        Debug(Debug::ERROR) << "Key " << key << " of the clustering is outside the range the "
                            << "lookup " << lookupFile << " covers. The clustering and the lookup "
                            << "are from different databases.\n";
        EXIT(EXIT_FAILURE);
    }
    const std::string &name = slice[static_cast<size_t>(key - lo)];
    if (name.empty()) {
        Debug(Debug::ERROR) << "Key " << key << " of the clustering has no entry in " << lookupFile
                            << ". The clustering and the lookup are from different databases.\n";
        EXIT(EXIT_FAILURE);
    }
    return name;
}

}  // namespace

int translatekeys(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    const std::string inTsv = par.db1;
    const std::string lookupFile = par.db2;
    const std::string outTsv = par.db3;
    const int threads = std::max(1, par.threads);

    // Sized so all threads' resident state together fits the memory budget.
    //
    // Every thread works a bucket of its own and holds that bucket's accession
    // slice, so the budget has to be divided by the thread count -- sizing one
    // slice against the whole budget, as this used to, overran it by exactly the
    // thread count and turned a 64-thread run into ~18x its --split-memory-limit.
    //
    // 192 bytes per key, not 64: an accession is ~30 bytes plus the std::string
    // that holds it, and pass 3 keeps three of these alive at once -- the slice,
    // the sorted TranslatedRow vector (which carries the member accession), and
    // the output buffer being assembled. Being wrong only changes how many buckets
    // are used.
    uint64_t keyCount = 0;
    {
        // Counting lines is one extra sequential pass, but the alternative is
        // guessing the key space from the clustering, which need not mention the
        // highest key at all.
        FILE *f = FileUtil::openFileOrDie(lookupFile.c_str(), "r", true);
        char *line = NULL;
        size_t cap = 0;
        while (getline(&line, &cap, f) > 0) keyCount++;
        free(line);
        fclose(f);
    }
    if (keyCount == 0) {
        Debug(Debug::ERROR) << lookupFile << " is empty\n";
        EXIT(EXIT_FAILURE);
    }
    const uint64_t targetBytes =
        std::max<uint64_t>(Util::computeMemory(par.splitMemoryLimit) / 8, 1ULL * 1024 * 1024);
    const uint64_t bytesPerKey = 192;
    const uint64_t perThreadBudget =
        std::max<uint64_t>(targetBytes / static_cast<uint64_t>(threads), 1ULL * 1024 * 1024);
    // Sized against the *rows* as well as the keys.
    //
    // A bucket is a range of representative key, but what pass 3 holds is one
    // TranslatedRow per *member* assigned into that range -- and cluster sizes are
    // skewed, so a bucket of `span` keys can hold far more than `span` rows. Sizing
    // on keys alone therefore bounded the accession slice and left the row vector
    // unbounded, which is the one that carries a std::string each.
    //
    // The row count is bounded from the input file: a row is "0\t0\n" at the very
    // least. ~128 B per row is the pair of keys plus a short accession and the
    // std::string that holds it; being wrong only changes how many buckets are used.
    const uint64_t rowsUpperBound = FileUtil::getFileSize(inTsv) / 4 + 1;
    const uint64_t bytesPerRow = 128;
    unsigned int buckets = 1;
    while (buckets < 65536
           && ((keyCount / buckets) * bytesPerKey > perThreadBudget
               || (rowsUpperBound / buckets) * bytesPerRow > perThreadBudget)) {
        buckets *= 2;
    }
    const uint64_t span = (keyCount + buckets - 1) / buckets;
    Debug(Debug::INFO) << "Translating " << keyCount << " keys over " << buckets << " buckets of "
                       << span << "\n";

    // Spill files go where --spill-prefix says, which the workflow points inside
    // its tmp directory. Deriving them from the output path -- which is what
    // happens when the option is not given -- puts every intermediate this command
    // writes on the *output* filesystem: at 1e11 the first fixed-width spill alone
    // is on the order of 1.6 TB, none of it counted by --scratch-budget, and a
    // small output filesystem fills even when scratch has room.
    const std::string spillBase = par.spillPrefix.empty() ? outTsv : par.spillPrefix;
    const std::string tmpA = spillBase + ".bymember";
    const std::string tmpB = spillBase + ".byrep";
    const std::string piecePrefix = spillBase + ".p";
    // Clear anything an interrupted earlier attempt left behind. A rerun only
    // truncates the shards it happens to touch, and pass 2's thread/bucket
    // assignment is dynamic, so shards from the previous attempt would otherwise
    // survive and pass 3 would emit both attempts' rows.
    removeSpillFiles(tmpA);
    removeSpillFiles(tmpB);
    // The output pieces too. Pass 3 writes <out>.p<bucket> and then concatenates
    // every one it finds, but it `continue`s past a bucket whose key range the
    // lookup never covers -- skipping the overwrite, not the concatenation. A
    // <out>.p<b> left by an interrupted earlier attempt was therefore appended to
    // the final result of the next one.
    removeSpillFiles(piecePrefix);

    // Pass 1: bucket by member key.
    {
        KeyBuckets byMember(tmpA, buckets);
        FILE *f = FileUtil::openFileOrDie(inTsv.c_str(), "r", true);
        char *line = NULL;
        size_t cap = 0;
        while (getline(&line, &cap, f) > 0) {
            char *tab = strchr(line, '\t');
            if (tab == NULL) continue;
            const uint64_t rep = strtoull(line, NULL, 10);
            const uint64_t member = strtoull(tab + 1, NULL, 10);
            if (member / span >= buckets || rep / span >= buckets) {
                Debug(Debug::ERROR) << "Clustering names key " << std::max(member, rep)
                                    << ", beyond the " << keyCount << " keys in " << lookupFile
                                    << ". The clustering and the lookup are from different "
                                    << "databases.\n";
                EXIT(EXIT_FAILURE);
            }
            byMember.append(static_cast<unsigned int>(member / span), member, rep);
        }
        free(line);
        fclose(f);
    }

    // Pass 2: translate the member, re-bucket on the representative.
    //
    // Parallel over buckets, each thread seeking straight to its own key range.
    // Its output goes to a thread-private set of representative buckets, so no
    // two threads share a writer; pass 3 reads all shards of a bucket.
    const std::vector<uint64_t> lookupAt = indexLookup(lookupFile, span, buckets);
    {
        std::vector<TextBuckets *> writers(threads, NULL);
        for (int t = 0; t < threads; t++) {
            writers[t] = new TextBuckets(tmpB + ".t" + SSTR(t), buckets);
        }
#pragma omp parallel num_threads(threads)
        {
            int tid = 0;
#ifdef OPENMP
            tid = omp_get_thread_num();
#endif
            std::vector<std::string> slice;
#pragma omp for schedule(dynamic, 1)
            for (int64_t bi = 0; bi < static_cast<int64_t>(buckets); bi++) {
                const unsigned int b = static_cast<unsigned int>(bi);
                const uint64_t lo = b * span;
                if (lo >= keyCount || lookupAt[b] == UINT64_MAX) {
                    continue;
                }
                const uint64_t hi = std::min(lo + span, keyCount);
                LookupCursor cursor(lookupFile, lookupAt[b]);
                cursor.slice(lo, hi, slice);
                const std::vector<KeyPair> pairs = KeyBuckets::read(tmpA, b);
                for (size_t i = 0; i < pairs.size(); i++) {
                    const std::string &name = accessionOf(slice, pairs[i].key, lo, lookupFile);
                    writers[tid]->append(static_cast<unsigned int>(pairs[i].other / span),
                                         pairs[i].other, pairs[i].key, name.c_str(), name.size());
                }
                FileUtil::remove(KeyBuckets::path(tmpA, b).c_str());
            }
        }
        for (int t = 0; t < threads; t++) {
            delete writers[t];
        }
    }

    // Pass 3: translate the representative and emit.
    //
    // Also parallel over buckets, each writing its own piece; the pieces are then
    // concatenated in bucket order.
    //
    // Each bucket's rows are sorted by (representative key, member key) before
    // they are written, which is what makes the output reproducible. Pass 2 spills
    // to thread-private shards, so which shard a row lands in depends on which
    // thread happened to take its source bucket -- and with `schedule(dynamic)`
    // that varies between two runs of the same binary on the same input, not just
    // between thread counts. Sorting on the keys, which are a property of the data
    // alone, removes the schedule from the result entirely. Keys rather than
    // accessions because the comparison is then an integer one, and because it
    // makes the output order match the key-space clustering this translates.
    uint64_t written = 0;
    {
#pragma omp parallel num_threads(threads) reduction(+ : written)
        {
            std::vector<std::string> slice;
            std::vector<TranslatedRow> entries;
            std::string buffer;
#pragma omp for schedule(dynamic, 1)
            for (int64_t bi = 0; bi < static_cast<int64_t>(buckets); bi++) {
                const unsigned int b = static_cast<unsigned int>(bi);
                const uint64_t lo = b * span;
                if (lo >= keyCount || lookupAt[b] == UINT64_MAX) {
                    continue;
                }
                const uint64_t hi = std::min(lo + span, keyCount);
                LookupCursor cursor(lookupFile, lookupAt[b]);
                cursor.slice(lo, hi, slice);
                buffer.clear();
                entries.clear();
                for (int t = 0; t < threads; t++) {
                    TextBuckets::read(tmpB + ".t" + SSTR(t), b, entries);
                    // A (thread, bucket) shard exists only if that thread wrote
                    // to that bucket, which is sparse.
                    const std::string shard = TextBuckets::path(tmpB + ".t" + SSTR(t), b);
                    if (FileUtil::fileExists(shard.c_str())) {
                        FileUtil::remove(shard.c_str());
                    }
                }
                // std::sort, not SORT_PARALLEL: this already runs inside the
                // parallel region, one bucket per thread.
                std::sort(entries.begin(), entries.end(), byRepThenMember);
                for (size_t i = 0; i < entries.size(); i++) {
                    const std::string &repName =
                        accessionOf(slice, entries[i].repKey, lo, lookupFile);
                    buffer.append(repName);
                    buffer.push_back('\t');
                    buffer.append(entries[i].memberAccession);
                    buffer.push_back('\n');
                    written++;
                }
                const std::string piecePath = piecePrefix + SSTR(b);
                FILE *piece = FileUtil::openAndDelete(piecePath.c_str(), "w");
                writeAllOrDie(buffer.data(), buffer.size(), piece, piecePath);
                if (fclose(piece) != 0) {
                    Debug(Debug::ERROR) << "Cannot close " << piecePath << "\n";
                    EXIT(EXIT_FAILURE);
                }
            }
        }
        FILE *out = FileUtil::openAndDelete(outTsv.c_str(), "w");
        std::vector<char> copy(8 << 20);
        for (unsigned int b = 0; b < buckets; b++) {
            const std::string piecePath = piecePrefix + SSTR(b);
            if (FileUtil::fileExists(piecePath.c_str()) == false) {
                continue;
            }
            FILE *piece = FileUtil::openFileOrDie(piecePath.c_str(), "rb", true);
            size_t got;
            while ((got = fread(copy.data(), 1, copy.size(), piece)) > 0) {
                writeAllOrDie(copy.data(), got, out, outTsv);
            }
            fclose(piece);
            FileUtil::remove(piecePath.c_str());
        }
        if (fclose(out) != 0) {
            Debug(Debug::ERROR) << "Cannot close " << outTsv << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

    Debug(Debug::INFO) << "Translated " << written << " assignments into " << outTsv << "\n";
    return EXIT_SUCCESS;
}
