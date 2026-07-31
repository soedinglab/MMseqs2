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
    Buckets(const std::string &prefix, unsigned int count, size_t bufferPairs = 1 << 16)
        : prefix(prefix), count(count), bufferPairs(bufferPairs), closed(false) {
        files.assign(count, NULL);
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
            if (files[b] != NULL) { fclose(files[b]); files[b] = NULL; }
        }
    }

    static std::string path(const std::string &prefix, unsigned int b) {
        return prefix + "." + SSTR(b);
    }

    static std::vector<Pair> read(const std::string &prefix, unsigned int b) {
        std::vector<Pair> out;
        const std::string p = path(prefix, b);
        if (FileUtil::fileExists(p.c_str()) == false) return out;
        const size_t bytes = FileUtil::getFileSize(p);
        if (bytes == 0 || bytes % sizeof(Pair) != 0) return out;
        out.resize(bytes / sizeof(Pair));
        FILE *f = FileUtil::openFileOrDie(p.c_str(), "rb", true);
        if (fread(out.data(), sizeof(Pair), out.size(), f) != out.size()) {
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
                Debug(Debug::ERROR) << "Cannot open " << path(prefix, b) << ": "
                                    << strerror(errno) << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        if (fwrite(buffers[b].data(), sizeof(Pair), buffers[b].size(), files[b]) !=
            buffers[b].size()) {
            Debug(Debug::ERROR) << "Cannot write bucket " << b << "\n";
            EXIT(EXIT_FAILURE);
        }
        buffers[b].clear();
    }

    std::string prefix;
    unsigned int count;
    size_t bufferPairs;
    std::vector<std::vector<Pair> > buffers;
    std::vector<FILE *> files;
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

    const uint64_t targetBytes = std::max<uint64_t>(Util::computeMemory(par.splitMemoryLimit) / 8,
                                                    1ULL * 1024 * 1024);
    unsigned int buckets = 1;
    while (buckets < 65536 && (subCount / buckets) * sizeof(uint64_t) > targetBytes) {
        buckets *= 2;
    }
    const uint64_t span = (subCount + buckets - 1) / buckets;
    Debug(Debug::INFO) << "Translating " << subCount << " sub-keys over " << buckets
                       << " buckets of " << span << "\n";

    const int mapFd = open(mapFile.c_str(), O_RDONLY);
    if (mapFd < 0) {
        Debug(Debug::ERROR) << "Cannot open " << mapFile << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }

    const std::string tmpA = outTsv + ".bymember";
    const std::string tmpB = outTsv + ".byrep";

    // Pass 1: bucket by the member sub-key, then translate it.
    {
        Buckets byMember(tmpA, buckets);
        FILE *f = FileUtil::openFileOrDie(inTsv.c_str(), "r", true);
        char *line = NULL;
        size_t cap = 0;
        while (getline(&line, &cap, f) > 0) {
            char *tab = strchr(line, '\t');
            if (tab == NULL) continue;
            const uint64_t subRep = strtoull(line, NULL, 10);
            const uint64_t subMember = strtoull(tab + 1, NULL, 10);
            byMember.append(static_cast<unsigned int>(subMember / span), subMember, subRep);
        }
        free(line);
        fclose(f);
    }
    {
        Buckets byRep(tmpB, buckets);
        for (unsigned int b = 0; b < buckets; b++) {
            const uint64_t lo = b * span;
            if (lo >= subCount) break;
            const uint64_t hi = std::min(lo + span, subCount);
            const std::vector<uint64_t> slice = mapSlice(mapFd, lo, hi);
            const std::vector<Pair> pairs = Buckets::read(tmpA, b);
            for (size_t i = 0; i < pairs.size(); i++) {
                if (pairs[i].key < lo || pairs[i].key - lo >= slice.size()) {
                    Debug(Debug::ERROR) << "Sub-key " << pairs[i].key << " is outside the key map "
                                        << mapFile << ", which holds " << subCount << " keys\n";
                    EXIT(EXIT_FAILURE);
                }
                const uint64_t origMember = slice[static_cast<size_t>(pairs[i].key - lo)];
                // Re-bucket on the representative sub-key for the second pass.
                byRep.append(static_cast<unsigned int>(pairs[i].other / span), pairs[i].other,
                             origMember);
            }
            FileUtil::remove(Buckets::path(tmpA, b).c_str());
        }
    }

    // Pass 2: translate the representative sub-key.
    FILE *out = FileUtil::openAndDelete(outTsv.c_str(), "w");
    std::string buffer;
    buffer.reserve(64 * 1024 * 1024);
    uint64_t written = 0;
    for (unsigned int b = 0; b < buckets; b++) {
        const uint64_t lo = b * span;
        if (lo >= subCount) break;
        const uint64_t hi = std::min(lo + span, subCount);
        const std::vector<uint64_t> slice = mapSlice(mapFd, lo, hi);
        const std::vector<Pair> pairs = Buckets::read(tmpB, b);
        for (size_t i = 0; i < pairs.size(); i++) {
            if (pairs[i].key < lo || pairs[i].key - lo >= slice.size()) {
                Debug(Debug::ERROR) << "Sub-key " << pairs[i].key << " is outside the key map "
                                    << mapFile << ", which holds " << subCount << " keys\n";
                EXIT(EXIT_FAILURE);
            }
            appendPair(buffer, slice[static_cast<size_t>(pairs[i].key - lo)], pairs[i].other);
            written++;
            if (buffer.size() > 32 * 1024 * 1024) {
                writeAllOrDie(buffer.data(), buffer.size(), out, outTsv);
                buffer.clear();
            }
        }
        FileUtil::remove(Buckets::path(tmpB, b).c_str());
    }
    if (buffer.empty() == false) writeAllOrDie(buffer.data(), buffer.size(), out, outTsv);
    if (fclose(out) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << outTsv << "\n";
        EXIT(EXIT_FAILURE);
    }
    close(mapFd);

    Debug(Debug::INFO) << "Translated " << written << " assignments into " << outTsv << "\n";
    return EXIT_SUCCESS;
}
