/*
 * createrepdb -- builds the representative sub-database for the next linclust pass.
 *
 * Replaces `createsubdb` for this pipeline, for two reasons.
 *
 * The first is the usual one: `createsubdb` opens the database with USE_INDEX and
 * so holds the text index resident, 24 bytes per sequence -- 2.4 TB at 1e11. This
 * streams through the dense companion index instead and holds one bit per key.
 *
 * The second is specific and easy to miss: **a sub-database of representatives has
 * sparse keys**, and every stage of the distributed pipeline assumes dense ones --
 * DenseIndex treats a key as an array offset, the greedy sweep walks the key space
 * in order, and the edge and merge buckets are key ranges. Feeding a `createsubdb`
 * output into pass 2 would break all three. So the representatives are *re-keyed*
 * densely here, and a `subkey -> original key` map is written alongside so the
 * pass-2 clustering can be translated back before merging.
 *
 * Re-keying is nearly free because of how the keys were assigned in the first
 * place. Original keys are length-ranked, so any subset of them is *already* in
 * descending length order; the i-th representative in ascending original-key order
 * is exactly sub-key i. No sort is needed, and the new database is a sequential
 * copy rather than a gather.
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

#include <fcntl.h>
#include <unistd.h>

#ifdef OPENMP
#include <omp.h>
#endif

namespace {

class Bitmap {
public:
    explicit Bitmap(uint64_t bits) : words((bits + 63) / 64, 0) {}
    bool get(uint64_t i) const { return (words[i >> 6] >> (i & 63)) & 1ULL; }
    void set(uint64_t i) { words[i >> 6] |= 1ULL << (i & 63); }
    uint64_t bytes() const { return words.size() * sizeof(uint64_t); }

private:
    std::vector<uint64_t> words;
};

int openOrDie(const std::string &path, int flags) {
    const int fd = open(path.c_str(), flags, 0666);
    if (fd < 0) {
        Debug(Debug::ERROR) << "Cannot open " << path << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    return fd;
}

void readAt(int fd, void *dst, size_t len, size_t off, const char *what) {
    char *p = static_cast<char *>(dst);
    size_t done = 0;
    while (done < len) {
        const ssize_t got = pread(fd, p + done, len - done, static_cast<off_t>(off + done));
        if (got <= 0) {
            if (got < 0 && errno == EINTR) continue;
            Debug(Debug::ERROR) << "Cannot read " << what << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        done += static_cast<size_t>(got);
    }
}

void writeAt(int fd, const void *src, size_t len, size_t off, const char *what) {
    const char *p = static_cast<const char *>(src);
    size_t done = 0;
    while (done < len) {
        const ssize_t put = pwrite(fd, p + done, len - done, static_cast<off_t>(off + done));
        if (put <= 0) {
            if (put < 0 && errno == EINTR) continue;
            Debug(Debug::ERROR) << "Cannot write " << what << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        done += static_cast<size_t>(put);
    }
}

void allocate(const std::string &path, uint64_t size) {
    const int fd = openOrDie(path, O_WRONLY | O_CREAT | O_TRUNC);
    if (size > 0 && ftruncate(fd, static_cast<off_t>(size)) != 0) {
        Debug(Debug::ERROR) << "Cannot size " << path << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    close(fd);
}

// Copies the entries flagged in `keep` into a new dense database, preserving key
// order. Returns the total data bytes and the longest entry.
struct CopyResult {
    uint64_t dataBytes;
    uint32_t maxLen;
};

CopyResult copyFlagged(const std::string &srcDb, const std::string &dstDb, const Bitmap &keep,
                       uint64_t entryCount, uint64_t keptCount, int threads,
                       std::vector<uint64_t> *keyMap) {
    const int srcIdx = openOrDie(DenseIndex::fileName(srcDb), O_RDONLY);
    const int srcData = openOrDie(srcDb, O_RDONLY);

    // Pass 1: per-key kept lengths, so the destination offsets are a prefix sum.
    // The index is read sequentially in blocks rather than randomly per key.
    std::vector<uint64_t> offsets(keptCount + 1, 0);
    std::vector<uint32_t> lengths(keptCount, 0);
    std::vector<DenseIndex::Entry> srcEntries(keptCount);
    {
        const uint64_t block = 1 << 20;
        std::vector<DenseIndex::Entry> buf(block);
        uint64_t out = 0;
        for (uint64_t from = 0; from < entryCount; from += block) {
            const uint64_t n = std::min(block, entryCount - from);
            readAt(srcIdx, buf.data(), n * sizeof(DenseIndex::Entry),
                   DenseIndex::entryOffset(from), "source index");
            for (uint64_t i = 0; i < n; i++) {
                if (keep.get(from + i) == false) continue;
                srcEntries[out] = buf[i];
                lengths[out] = buf[i].length;
                if (keyMap != NULL) (*keyMap)[out] = from + i;
                out++;
            }
        }
    }
    CopyResult result;
    result.maxLen = 0;
    uint64_t total = 0;
    for (uint64_t i = 0; i < keptCount; i++) {
        offsets[i] = total;
        total += lengths[i];
        result.maxLen = std::max(result.maxLen, lengths[i]);
    }
    offsets[keptCount] = total;
    result.dataBytes = total;

    allocate(dstDb, total);
    DenseIndex::createEmpty(dstDb, keptCount, 0, total, result.maxLen);

    const int dstData = openOrDie(dstDb, O_WRONLY);
    const int dstIdx = openOrDie(DenseIndex::fileName(dstDb), O_WRONLY);

    // Pass 2: copy. Entries are independent and both sides are in ascending offset
    // order, so this is a threaded forward scan rather than a gather.
#pragma omp parallel num_threads(threads)
    {
        std::vector<char> buf;
        std::vector<char> span;
        std::vector<DenseIndex::Entry> idxBuf;
        const uint64_t chunk = 4096;
#pragma omp for schedule(dynamic, 1)
        for (uint64_t start = 0; start < keptCount; start += chunk) {
            const uint64_t stop = std::min(start + chunk, keptCount);
            const uint64_t bytes = offsets[stop] - offsets[start];
            buf.resize(static_cast<size_t>(bytes));
            idxBuf.resize(static_cast<size_t>(stop - start));
            for (uint64_t i = start; i < stop; i++) {
                idxBuf[i - start].offset = offsets[i];
                idxBuf[i - start].length = lengths[i];
            }
            // Coalesced reads. Representatives are a *subset* of the source, so
            // one pread per entry meant 353M random reads at 1e9 and made this
            // stage 4x stock's createsubdb. Consecutive representatives are only a
            // couple of sequence lengths apart, so a run of them is fetched in one
            // read and the wanted pieces copied out.
            //
            // The run is closed when the bytes read would exceed twice the bytes
            // actually wanted: coalescing across a sparse region reads the skipped
            // sequences too, and without a cap a database whose representatives
            // are thinly spread would read the whole file. That bounds the waste
            // at 2x while still collapsing dense regions into single reads.
            uint64_t i = start;
            while (i < stop) {
                const uint64_t from = srcEntries[i].offset;
                uint64_t to = from + lengths[i];
                uint64_t wanted = lengths[i];
                uint64_t j = i + 1;
                while (j < stop) {
                    const uint64_t end = srcEntries[j].offset + lengths[j];
                    if (end - from > 2 * (wanted + lengths[j])) {
                        break;
                    }
                    to = end;
                    wanted += lengths[j];
                    j++;
                }
                span.resize(static_cast<size_t>(to - from));
                readAt(srcData, span.data(), static_cast<size_t>(to - from), from, "source data");
                for (uint64_t k = i; k < j; k++) {
                    memcpy(buf.data() + (offsets[k] - offsets[start]),
                           span.data() + (srcEntries[k].offset - from), lengths[k]);
                }
                i = j;
            }
            writeAt(dstData, buf.data(), static_cast<size_t>(bytes), offsets[start], "database");
            writeAt(dstIdx, idxBuf.data(), idxBuf.size() * sizeof(DenseIndex::Entry),
                    DenseIndex::entryOffset(start), "index");
        }
    }

    close(srcIdx);
    close(srcData);
    close(dstData);
    close(dstIdx);
    return result;
}

}  // namespace

int createrepdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    const std::string seqDb = par.db1;
    const std::string clusterTsv = par.db2;
    const std::string repDb = par.db3;

    const DenseIndex::Info info = DenseIndex::readInfo(seqDb);

    // A key is a representative exactly when it appears in the first column.
    Bitmap isRep(info.entryCount);
    uint64_t repCount = 0;
    {
        FILE *f = FileUtil::openFileOrDie(clusterTsv.c_str(), "r", true);
        char *line = NULL;
        size_t cap = 0;
        while (getline(&line, &cap, f) > 0) {
            const uint64_t rep = strtoull(line, NULL, 10);
            if (rep < info.entryCount && isRep.get(rep) == false) {
                isRep.set(rep);
                repCount++;
            }
        }
        free(line);
        fclose(f);
    }
    Debug(Debug::INFO) << repCount << " representatives of " << info.entryCount << " sequences ("
                       << isRep.bytes() / (1024 * 1024) << " MB of flags)\n";
    if (repCount == 0) {
        Debug(Debug::ERROR) << "No representatives found in " << clusterTsv << "\n";
        EXIT(EXIT_FAILURE);
    }

    std::vector<uint64_t> keyMap(repCount, 0);
    const CopyResult seq = copyFlagged(seqDb, repDb, isRep, info.entryCount, repCount, par.threads,
                                       &keyMap);
    copyFlagged(seqDb + "_h", repDb + "_h", isRep, info.entryCount, repCount, par.threads, NULL);

    // subkey -> original key, dense in the sub-key space and in ascending order,
    // so the translation back is a sequential read rather than a lookup structure.
    const std::string mapFile = repDb + ".keymap";
    const std::string keymapTmp = mapFile + ".tmp";
    FILE *m = FileUtil::openAndDelete(keymapTmp.c_str(), "wb");
    if (fwrite(keyMap.data(), sizeof(uint64_t), keyMap.size(), m) != keyMap.size()) {
        Debug(Debug::ERROR) << "Cannot write " << keymapTmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fclose(m) != 0) {

        Debug(Debug::ERROR) << "Cannot close " << keymapTmp << ": " << strerror(errno) << "\n";

        EXIT(EXIT_FAILURE);

    }

    // Renamed only once complete: the pass-2 filter gate sizes itself from this

    // file's length, so a truncated one would silently be used as a short,

    // wrong sub-key to original-key map.

    FileUtil::move(keymapTmp.c_str(), mapFile.c_str());

    const int dbType = FileUtil::parseDbType(seqDb.c_str());
    FileUtil::writeFile(repDb + ".dbtype", reinterpret_cast<const unsigned char *>(&dbType),
                        sizeof(int));
    const int hdrType = Parameters::DBTYPE_GENERIC_DB;
    FileUtil::writeFile(repDb + "_h.dbtype", reinterpret_cast<const unsigned char *>(&hdrType),
                        sizeof(int));

    Debug(Debug::INFO) << "Wrote " << repDb << ": " << repCount << " sequences, " << seq.dataBytes
                       << " data bytes, longest " << seq.maxLen << "\n";
    return EXIT_SUCCESS;
}
