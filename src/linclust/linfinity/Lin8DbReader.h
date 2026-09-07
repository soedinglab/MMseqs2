#ifndef MMSEQS_READLINDB_H
#define MMSEQS_READLINDB_H

#include "Lin8Db.h"
#include <cstddef>
#include <cstdint>
#include <vector>
#include <string>
#include "Debug.h"
#include "FileUtil.h"
#include <cstdio>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <ctime>

struct IoRing {
    struct Read {
        void *into;
        int fd;
        uint64_t offset;
        size_t length;
        size_t required;
    };

    IoRing();
    ~IoRing();

    bool open(unsigned depth);
    bool isOpen() const { return ready; }

    std::vector<Read> &list() { return reads; }
    void submit(const char *what);
    void await(const char *what);

private:
    IoRing(const IoRing &);
    IoRing &operator=(const IoRing &);
    void preadAll(const char *what);
    void pump(const char *what, bool untilDone);

    bool ready;
    void *state;
    std::vector<Read> reads;
    size_t queued;
    size_t done;
    unsigned inflight;
};

class Lin8DbReader {
public:
    static const uint64_t KEPT_BITMAP_MAGIC;

    class Cursor {
    public:
        Cursor() : at(0) {}
        size_t at;
    };

    static const char *KEPT_BITMAP_SUFFIX;

    Lin8DbReader(const std::string &db, bool withHeaders = false);
    ~Lin8DbReader();

    void open();
    void close();

    uint64_t getSize() const { return index.getSize(); }
    uint64_t getDataSize() const { return index.getDataSize(); }
    const Lin8DbIndex &getIndex() const { return index; }

    uint32_t getSeqLen(uint64_t rank) const { return index.getSeqLen(rank); }
    const char *getData(uint64_t rank) const;

    uint32_t getSeqLen(uint64_t rank, Cursor &cursor) const;
    const char *getData(uint64_t rank, Cursor &cursor) const;

    static const int READ_ONCE = 0;
    static const int READ_AGAIN = 1;
    void openBatch(unsigned int threads, size_t arenaBytes, size_t memoryBudget, int revisit = READ_ONCE);

    size_t batchRoomFor(uint32_t seqLen) const;

    size_t startBatch(uint64_t queryRank, const uint64_t *members, size_t n, unsigned int thread,
                      unsigned int lane) const;
    void awaitBatch(unsigned int thread, unsigned int lane) const;
    const char *batchQueryAt(unsigned int thread, unsigned int lane) const;
    const char *batchAt(unsigned int thread, unsigned int lane, size_t member) const;

    static const unsigned int LANES = 2;

    bool isKept(uint64_t rank) const;
    uint64_t countKept() const;
    bool hasKeptBitmap() const { return keptLoaded; }
    const uint64_t *keptWords() const { return kept; }
    size_t keptWordCount() const { return keptCount; }

    void releaseFileSlot(size_t suffix);

    uint64_t rankAtByte(uint64_t globalByte) const { return index.rankAtByte(globalByte); }

    // --compressed header files: zstd frames plus a table keyed by uncompressed offset
    struct HeaderFrame {
        uint64_t rawAt;
        uint64_t packedAt;
        uint64_t rawSize;
        uint64_t packedSize;
    };
    static const uint64_t HEADER_FRAMES_MAGIC;

    class HeaderStream {
    public:
        HeaderStream(const Lin8DbReader &owner);
        // begin stays kept only until the next call
        bool next(const char *&begin, size_t &length);
    private:
        const char *frameText(uint32_t file, size_t &avail);
        const Lin8DbReader &owner;
        size_t range;
        uint64_t left;
        size_t at;
        uint32_t frameFile;
        size_t frame;
        std::vector<char> raw;
    };

private:
    struct BatchLane {
        std::vector<char> arena;
        char *aligned;
        const char *queryAt;
        std::vector<const char *> memberAt;
        IoRing ring;
        BatchLane() : aligned(NULL), queryAt(NULL) {}
    };
    struct BatchWorker {
        BatchLane lane[LANES];
    };
    mutable std::vector<BatchWorker *> batch;

    bool appendBatchRead(BatchLane &lane, uint64_t rank, Cursor &cursor, const char *&at) const;
    int directOf(uint32_t file) const;
    const char *fileData(uint32_t file, uint64_t offset) const;
    void mapFile(uint32_t file) const;
    void mapHeader(uint32_t file) const;

    std::string db;
    bool withHeaders;
    Lin8DbIndex index;
    mutable std::vector<char *> data;
    mutable std::vector<int> dataFd;
    mutable std::vector<int> directFd;
    mutable std::vector<int> headerFd;
    std::vector<size_t> dataSize;
    mutable std::vector<char *> headers;
    std::vector<size_t> headerSize;
    mutable std::vector<std::vector<HeaderFrame> > headerFrames;
    mutable std::vector<uint64_t> headerRawSize;
    void loadHeaderFrames(uint32_t file, const char *mapped) const;
    const uint64_t *kept;
    void *keptMap;
    size_t keptSize;
    size_t keptCount;
    bool keptLoaded;
    mutable bool wantDirect;
};

class ClusterAssignmentBitmap {
public:
    ClusterAssignmentBitmap() : words(NULL), ranks(0), wordCount(0), appliedRepRankBlocks(0) {}

    ~ClusterAssignmentBitmap() {
        if (words != NULL) {
            munmap(words, bytes());
        }
    }

    void open(const std::string &path, uint64_t ranks) {
        this->ranks = ranks;
        this->cache = path;
        wordCount = (ranks + 63) / 64;
        void *at = mmap(NULL, bytes(), PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        if (at == MAP_FAILED) {
            Debug(Debug::ERROR) << "Cannot reserve " << bytes() << " byte for the cluster assignment bitmap\n";
            EXIT(EXIT_FAILURE);
        }
        words = static_cast<uint64_t *>(at);
        if (path.empty() || FileUtil::fileExists(path.c_str()) == false) {
            return;
        }
        FILE *in = FileUtil::openFileOrDie(path.c_str(), "r", true);
        uint64_t header[HEADER_WORDS] = {0, 0, 0};
        if (fread(header, sizeof(uint64_t), HEADER_WORDS, in) != HEADER_WORDS) {
            Debug(Debug::ERROR) << "Cannot read the header of " << path << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (header[0] != MAGIC || header[1] != ranks) {
            Debug(Debug::ERROR) << path << " covers " << header[1] << " sequences and this database "
                                << "holds " << ranks << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (wordCount > 0 && fread(words, sizeof(uint64_t), wordCount, in) != wordCount) {
            Debug(Debug::ERROR) << "Cannot read " << path << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (fclose(in) != 0) {
            Debug(Debug::ERROR) << "Cannot close " << path << "\n";
            EXIT(EXIT_FAILURE);
        }
        appliedRepRankBlocks = header[2];
    }

    bool isAssigned(uint64_t rank) const { return (words[rank >> 6] >> (rank & 63) & 1) != 0; }
    void assign(uint64_t rank) { words[rank >> 6] |= uint64_t(1) << (rank & 63); }
    uint64_t size() const { return ranks; }
    size_t applied() const { return appliedRepRankBlocks; }

    size_t bytesHeld() const { return bytes(); }

    void catchUpTo(const std::string &decided, size_t repRankBlock) {
        std::vector<PairRecord> buffer(1u << 16);
        for (size_t at = applied(); at < repRankBlock; at++) {
            const std::string path = decided + ".0." + SSTR(at);
            FILE *in = fopen(path.c_str(), "r");
            if (in == NULL) {
                Debug(Debug::ERROR) << "Cannot open " << path
                                    << ", which repRankBlock " << at << " should have decided\n";
                EXIT(EXIT_FAILURE);
            }
            size_t read = 0;
            while ((read = readRecords(buffer.data(), buffer.size(), in)) > 0) {
                for (size_t k = 0; k < read; k++) {
                    assign(buffer[k].rep());
                    assign(buffer[k].member());
                }
            }
            if (ferror(in) != 0) {
                Debug(Debug::ERROR) << "Cannot read " << path << "\n";
                EXIT(EXIT_FAILURE);
            }
            fclose(in);
        }
        if (repRankBlock > appliedRepRankBlocks) {
            appliedRepRankBlocks = repRankBlock;
        }
    }

    // past the floor a missing decided block is one to skip, because deciding reads it again anyway
    void catchUpToAvailable(const std::string &decided, size_t repRankBlock,
                            size_t floorRepRankBlock) {
        std::vector<PairRecord> buffer;
        size_t at = applied();
        unsigned int waited = 0;
        while (at < repRankBlock) {
            const std::string path = decided + ".0." + SSTR(at);
            FILE *in = fopen(path.c_str(), "r");
            if (in == NULL) {
                if (at >= floorRepRankBlock) {
                    break;
                }
                static time_t lastWaitLog = 0;
                const time_t nowSec = time(NULL);
                if (waited >= 30 && nowSec - lastWaitLog >= 60) {
                    Debug(Debug::INFO) << "Still waiting for the decider to publish block " << at
                                       << " (" << waited << "s)\n";
                    lastWaitLog = nowSec;
                }
                if (waited >= 3600) {
                    Debug(Debug::ERROR) << "Waited " << waited << "s for " << path
                                        << ", which repRankBlock " << at << " should have decided\n";
                    EXIT(EXIT_FAILURE);
                }
                sleep(1);
                waited++;
                continue;
            }
            waited = 0;
            buffer.resize(1u << 16);
            size_t read = 0;
            while ((read = readRecords(buffer.data(), buffer.size(), in)) > 0) {
                for (size_t k = 0; k < read; k++) {
                    assign(buffer[k].rep());
                    assign(buffer[k].member());
                }
            }
            if (ferror(in) != 0) {
                Debug(Debug::ERROR) << "Cannot read " << path << "\n";
                EXIT(EXIT_FAILURE);
            }
            fclose(in);
            at++;
        }
        if (at > appliedRepRankBlocks) {
            appliedRepRankBlocks = at;
        }
    }

    void save(size_t repRankBlock) {
        if (cache.empty()) {
            return;
        }
        if (repRankBlock > appliedRepRankBlocks) {
            appliedRepRankBlocks = repRankBlock;
        }
        const std::string tmp = cache + ".tmp";
        FILE *out = FileUtil::openAndDelete(tmp.c_str(), "w");
        const uint64_t header[HEADER_WORDS] = {MAGIC, ranks, appliedRepRankBlocks};
        if (fwrite(header, sizeof(uint64_t), HEADER_WORDS, out) != HEADER_WORDS
            || (wordCount > 0 && fwrite(words, sizeof(uint64_t), wordCount, out) != wordCount)) {
            Debug(Debug::ERROR) << "Cannot write " << tmp << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (fclose(out) != 0) {
            Debug(Debug::ERROR) << "Cannot close " << tmp << "\n";
            EXIT(EXIT_FAILURE);
        }
        FileUtil::publishAtomically(tmp, cache);
    }

private:
    ClusterAssignmentBitmap(const ClusterAssignmentBitmap &);
    ClusterAssignmentBitmap &operator=(const ClusterAssignmentBitmap &);

    size_t bytes() const { return wordCount * sizeof(uint64_t); }

    static const size_t HEADER_WORDS = 3;
    static const uint64_t MAGIC = 0x4C494E4354414B4Eull;
    uint64_t *words;
    uint64_t ranks;
    size_t wordCount;
    size_t appliedRepRankBlocks;
    std::string cache;
};

inline bool assignCluster(uint64_t rep, const uint64_t *members, size_t count,
                          ClusterAssignmentBitmap &assignedCluster, std::vector<PairRecord> &out,
                          uint64_t &assigned) {
    if (assignedCluster.isAssigned(rep)) {
        return false;
    }
    assignedCluster.assign(rep);
    PairRecord line;
    line.set(rep, rep, 0);
    out.push_back(line);
    for (size_t i = 0; i < count; i++) {
        const uint64_t member = members[i];
        if (member == rep || assignedCluster.isAssigned(member)) {
            continue;
        }
        assignedCluster.assign(member);
        line.set(rep, member, 0);
        out.push_back(line);
        assigned++;
    }
    return true;
}

#endif
