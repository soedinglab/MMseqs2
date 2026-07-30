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

namespace {

// Spill file holding (key, payload) records whose payload is a key too.
struct __attribute__((__packed__)) KeyPair {
    uint64_t key;    // the column being translated in this pass
    uint64_t other;  // the column carried through untouched
};

// Spill file holding (key, text) records. Needed for the second pass, where the
// column carried through has already become an accession and so varies in width.
class TextBuckets {
public:
    TextBuckets(const std::string &prefix, unsigned int count, size_t flushBytes = 8 << 20)
        : prefix(prefix), count(count), flushBytes(flushBytes), closed(false) {
        files.assign(count, NULL);
        buffers.resize(count);
    }
    ~TextBuckets() { close(); }

    void append(unsigned int b, uint64_t key, const char *text, size_t length) {
        std::string &buf = buffers[b];
        const uint32_t len = static_cast<uint32_t>(length);
        buf.append(reinterpret_cast<const char *>(&key), sizeof(key));
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
            if (files[b] != NULL) { fclose(files[b]); files[b] = NULL; }
        }
    }

    static std::string path(const std::string &prefix, unsigned int b) {
        return prefix + "." + SSTR(b);
    }

    // Reads one bucket back as (key, accession) pairs.
    static void read(const std::string &prefix, unsigned int b,
                     std::vector<std::pair<uint64_t, std::string> > &out) {
        out.clear();
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
        while (at + sizeof(uint64_t) + sizeof(uint32_t) <= bytes) {
            uint64_t key;
            uint32_t len;
            memcpy(&key, blob.data() + at, sizeof(key));
            at += sizeof(key);
            memcpy(&len, blob.data() + at, sizeof(len));
            at += sizeof(len);
            out.push_back(std::make_pair(key, blob.substr(at, len)));
            at += len;
        }
    }

private:
    void flush(unsigned int b) {
        if (buffers[b].empty()) return;
        if (files[b] == NULL) {
            files[b] = fopen(path(prefix, b).c_str(), "wb");
            if (files[b] == NULL) {
                Debug(Debug::ERROR) << "Cannot open " << path(prefix, b) << ": " << strerror(errno)
                                    << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        if (fwrite(buffers[b].data(), 1, buffers[b].size(), files[b]) != buffers[b].size()) {
            Debug(Debug::ERROR) << "Cannot write bucket " << b << "\n";
            EXIT(EXIT_FAILURE);
        }
        buffers[b].clear();
    }

    std::string prefix;
    unsigned int count;
    size_t flushBytes;
    std::vector<std::string> buffers;
    std::vector<FILE *> files;
    bool closed;
};

// Fixed-width spill, for the first pass where both columns are still keys.
class KeyBuckets {
public:
    KeyBuckets(const std::string &prefix, unsigned int count, size_t bufferPairs = 1 << 16)
        : prefix(prefix), count(count), bufferPairs(bufferPairs), closed(false) {
        files.assign(count, NULL);
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
            if (files[b] != NULL) { fclose(files[b]); files[b] = NULL; }
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
    void flush(unsigned int b) {
        if (buffers[b].empty()) return;
        if (files[b] == NULL) {
            files[b] = fopen(path(prefix, b).c_str(), "wb");
            if (files[b] == NULL) {
                Debug(Debug::ERROR) << "Cannot open " << path(prefix, b) << ": " << strerror(errno)
                                    << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        if (fwrite(buffers[b].data(), sizeof(KeyPair), buffers[b].size(), files[b]) !=
            buffers[b].size()) {
            Debug(Debug::ERROR) << "Cannot write bucket " << b << "\n";
            EXIT(EXIT_FAILURE);
        }
        buffers[b].clear();
    }

    std::string prefix;
    unsigned int count;
    size_t bufferPairs;
    std::vector<std::vector<KeyPair> > buffers;
    std::vector<FILE *> files;
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
    explicit LookupCursor(const std::string &path) : line(NULL), cap(0), pending(false), pendingKey(0) {
        file = FileUtil::openFileOrDie(path.c_str(), "r", true);
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

const std::string &accessionOf(const std::vector<std::string> &slice, uint64_t key, uint64_t lo,
                               const std::string &lookupFile) {
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

    // Sized so one slice of accessions fits the memory budget. An accession is
    // ~30 bytes plus the std::string that holds it, so 64 bytes per key is a safe
    // working figure; being wrong only changes how many buckets are used.
    uint64_t keyCount = 0;
    {
        LookupCursor probe(lookupFile);
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
    unsigned int buckets = 1;
    while (buckets < 65536 && (keyCount / buckets) * 64 > targetBytes) {
        buckets *= 2;
    }
    const uint64_t span = (keyCount + buckets - 1) / buckets;
    Debug(Debug::INFO) << "Translating " << keyCount << " keys over " << buckets << " buckets of "
                       << span << "\n";

    const std::string tmpA = outTsv + ".bymember";
    const std::string tmpB = outTsv + ".byrep";

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
            byMember.append(static_cast<unsigned int>(member / span), member, rep);
        }
        free(line);
        fclose(f);
    }

    // Pass 2: translate the member, re-bucket on the representative.
    {
        LookupCursor cursor(lookupFile);
        TextBuckets byRep(tmpB, buckets);
        std::vector<std::string> slice;
        for (unsigned int b = 0; b < buckets; b++) {
            const uint64_t lo = b * span;
            if (lo >= keyCount) break;
            const uint64_t hi = std::min(lo + span, keyCount);
            cursor.slice(lo, hi, slice);
            const std::vector<KeyPair> pairs = KeyBuckets::read(tmpA, b);
            for (size_t i = 0; i < pairs.size(); i++) {
                const std::string &name = accessionOf(slice, pairs[i].key, lo, lookupFile);
                byRep.append(static_cast<unsigned int>(pairs[i].other / span), pairs[i].other,
                             name.c_str(), name.size());
            }
            FileUtil::remove(KeyBuckets::path(tmpA, b).c_str());
        }
    }

    // Pass 3: translate the representative and emit.
    LookupCursor cursor(lookupFile);
    FILE *out = FileUtil::openAndDelete(outTsv.c_str(), "w");
    std::string buffer;
    buffer.reserve(64 * 1024 * 1024);
    std::vector<std::string> slice;
    std::vector<std::pair<uint64_t, std::string> > entries;
    uint64_t written = 0;
    for (unsigned int b = 0; b < buckets; b++) {
        const uint64_t lo = b * span;
        if (lo >= keyCount) break;
        const uint64_t hi = std::min(lo + span, keyCount);
        cursor.slice(lo, hi, slice);
        TextBuckets::read(tmpB, b, entries);
        for (size_t i = 0; i < entries.size(); i++) {
            const std::string &repName = accessionOf(slice, entries[i].first, lo, lookupFile);
            buffer.append(repName);
            buffer.push_back('\t');
            buffer.append(entries[i].second);
            buffer.push_back('\n');
            written++;
            if (buffer.size() > 32 * 1024 * 1024) {
                fwrite(buffer.data(), 1, buffer.size(), out);
                buffer.clear();
            }
        }
        FileUtil::remove(TextBuckets::path(tmpB, b).c_str());
    }
    if (buffer.empty() == false) fwrite(buffer.data(), 1, buffer.size(), out);
    if (fclose(out) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << outTsv << "\n";
        EXIT(EXIT_FAILURE);
    }

    Debug(Debug::INFO) << "Translated " << written << " assignments into " << outTsv << "\n";
    return EXIT_SUCCESS;
}
