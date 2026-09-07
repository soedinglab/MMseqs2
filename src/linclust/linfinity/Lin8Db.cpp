#include "Lin8Db.h"
#include <functional>
#include "NodePlacement.h"

#include "Debug.h"
#include "ByteParser.h"
#include "FileUtil.h"
#include "Util.h"
#include "IndexTypes.h"
#include <climits>
#include <ctime>
#include <unistd.h>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include "KSeqWrapper.h"
#define ZSTD_STATIC_LINKING_ONLY
#include <zstd.h>
#include <fcntl.h>

namespace {
const char RUN_TABLE_MAGIC[8] = {'M', 'M', 'R', 'U', 'N', 'S', '\0', '\0'};
const uint32_t RUN_TABLE_VERSION = 4;

struct Header {
    char magic[8];
    uint32_t version;
    uint32_t keyWidth;
    uint32_t nodeCount;
    uint32_t filesPerNode;
    uint64_t segmentCount;
    uint64_t entryCount;
    uint64_t byteCount;
};
}

SequenceLocator::SequenceLocator() : entries(0), bytes(0), nodes(1), perNodeFiles(1) {}

void SequenceLocator::reserve(size_t count) {
    runs.reserve(count);
    byteStarts.reserve(count + 1);
}

void SequenceLocator::append(uint64_t rankBase, uint32_t len, uint64_t byteBase, uint32_t file,
                      uint64_t hdrBase) {
    if (rankBase > MAX_RANK || byteBase > MAX_BYTE || file > MAX_FILE || len > MAX_ENTRY_LEN) {
        Debug(Debug::ERROR) << "Run table segment out of repRankBlock: rank " << rankBase << " byte "
                            << byteBase << " file " << file << " length " << len << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (runs.empty() == false && rankBase <= runs.back().rankBase()) {
        Debug(Debug::ERROR) << "Run table runs must ascend by rank, got " << rankBase
                            << " after " << runs.back().rankBase() << "\n";
        EXIT(EXIT_FAILURE);
    }
    LengthRun segment;
    segment.rankAndLen = rankBase | (static_cast<uint64_t>(len) << RANK_BITS);
    segment.byteAndFile = byteBase | (static_cast<uint64_t>(file) << 48);
    segment.hdrByte = hdrBase;
    runs.push_back(segment);
    byteStarts.push_back(bytes);
    if (runs.size() > 1) {
        const LengthRun &previous = runs[runs.size() - 2];
        bytes += (rankBase - previous.rankBase()) * previous.seqLen();
        byteStarts.back() = bytes;
    }
    entries = rankBase;
}

void SequenceLocator::checkLengthsDescend() const {
    for (size_t i = 1; i < runs.size(); i++) {
        if (runs[i].seqLen() > runs[i - 1].seqLen()) {
            Debug(Debug::ERROR) << "Length run " << i << " holds length " << runs[i].seqLen()
                                << " after " << runs[i - 1].seqLen()
                                << ", the database is not sorted by length descending\n";
            EXIT(EXIT_FAILURE);
        }
    }
}

void SequenceLocator::finish(uint64_t totalEntries) {
    if (totalEntries > MAX_RANK + 1) {
        Debug(Debug::ERROR) << "Run table holds " << totalEntries << " entries, the rank field fits "
                            << (MAX_RANK + 1) << "\n";
        EXIT(EXIT_FAILURE);
    }
    entries = totalEntries;
    rebuildByteStarts();
}

void SequenceLocator::rebuildByteStarts() {
    byteStarts.assign(runs.size(), 0);
    bytes = 0;
    for (size_t i = 0; i < runs.size(); i++) {
        byteStarts[i] = bytes;
        const uint64_t next = (i + 1 < runs.size()) ? runs[i + 1].rankBase() : entries;
        bytes += (next - runs[i].rankBase()) * runs[i].seqLen();
    }
}

size_t SequenceLocator::runOf(uint64_t rank) const {
    if (runs.empty() || rank >= entries) {
        Debug(Debug::ERROR) << "Run table lookup for rank " << rank << " of " << entries << "\n";
        EXIT(EXIT_FAILURE);
    }
    size_t low = 0;
    size_t high = runs.size() - 1;
    while (low < high) {
        const size_t mid = low + (high - low + 1) / 2;
        low = (runs[mid].rankBase() <= rank) ? mid : low;
        high = (runs[mid].rankBase() <= rank) ? high : mid - 1;
    }
    return low;
}

size_t SequenceLocator::runOfFrom(uint64_t rank, size_t cursor) const {
    if (rank < runs[cursor].rankBase()) {
        return runOf(rank);
    }
    while (cursor + 1 < runs.size() && runs[cursor + 1].rankBase() <= rank) {
        cursor++;
    }
    return cursor;
}

uint64_t SequenceLocator::offsetIn(size_t segment, uint64_t rank) const {
    const LengthRun &at = runs[segment];
    return at.byteBase() + (rank - at.rankBase()) * at.seqLen();
}

uint64_t SequenceLocator::byteAtRank(uint64_t rank) const {
    if (rank >= entries) {
        return bytes;
    }
    const size_t segment = runOf(rank);
    return byteStarts[segment] + (rank - runs[segment].rankBase()) * runs[segment].seqLen();
}

uint64_t SequenceLocator::rankAtByte(uint64_t globalByte) const {
    if (runs.empty() || globalByte >= bytes) {
        return entries;
    }
    const size_t segment = static_cast<size_t>(
        std::upper_bound(byteStarts.begin(), byteStarts.end(), globalByte) - byteStarts.begin() - 1);
    const uint64_t inRun = (globalByte - byteStarts[segment]) / runs[segment].seqLen();
    return runs[segment].rankBase() + inRun;
}

void SequenceLocator::write(const std::string &path) const {
    char host[HOST_NAME_MAX + 1];
    memset(host, 0, sizeof(host));
    gethostname(host, HOST_NAME_MAX);
    const std::string tmp = path + ".tmp" + host + "." + SSTR(getpid());
    FILE *out = FileUtil::openAndDelete(tmp.c_str(), "wb");
    Header header;
    memcpy(header.magic, RUN_TABLE_MAGIC, sizeof(header.magic));
    header.version = RUN_TABLE_VERSION;
    header.keyWidth = static_cast<uint32_t>(sizeof(DBKeyType));
    header.nodeCount = nodes;
    header.filesPerNode = perNodeFiles;
    header.segmentCount = runs.size();
    header.entryCount = entries;
    header.byteCount = bytes;
    if (fwrite(&header, sizeof(Header), 1, out) != 1) {
        Debug(Debug::ERROR) << "Cannot write the sequence locator header to " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (runs.empty() == false
        && fwrite(runs.data(), sizeof(LengthRun), runs.size(), out) != runs.size()) {
        Debug(Debug::ERROR) << "Cannot write the length runs to " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fclose(out) != 0) {
        Debug(Debug::ERROR) << "Cannot close the sequence locator " << tmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::publishAtomically(tmp, path);
}

void SequenceLocator::read(const std::string &path) {
    FILE *in = fopen(path.c_str(), "rb");
    if (in == NULL) {
        Debug(Debug::ERROR) << "Cannot open the sequence locator " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    Header header;
    if (fread(&header, sizeof(Header), 1, in) != 1
        || memcmp(header.magic, RUN_TABLE_MAGIC, sizeof(header.magic)) != 0) {
        Debug(Debug::ERROR) << "File " << path << " is not a sequence locator\n";
        EXIT(EXIT_FAILURE);
    }
    if (header.version != RUN_TABLE_VERSION) {
        Debug(Debug::ERROR) << "Run table " << path << " has version " << header.version
                            << ", expected " << RUN_TABLE_VERSION << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (header.keyWidth != sizeof(DBKeyType)) {
        Debug(Debug::ERROR) << "Run table " << path << " was written with a " << header.keyWidth
                            << " byte key, this build uses " << sizeof(DBKeyType) << "\n";
        EXIT(EXIT_FAILURE);
    }
    runs.resize(header.segmentCount);
    if (header.segmentCount > 0
        && fread(runs.data(), sizeof(LengthRun), runs.size(), in) != runs.size()) {
        Debug(Debug::ERROR) << "Run table " << path << " is truncated\n";
        EXIT(EXIT_FAILURE);
    }
    fclose(in);
    entries = header.entryCount;
    nodes = (header.nodeCount > 0) ? header.nodeCount : 1;
    perNodeFiles = (header.filesPerNode > 0) ? header.filesPerNode : 1;
    rebuildByteStarts();
    if (bytes != header.byteCount) {
        Debug(Debug::ERROR) << "Run table " << path << " covers " << bytes << " byte, header says "
                            << header.byteCount << "\n";
        EXIT(EXIT_FAILURE);
    }
}

static std::string fileNameExtension(const std::string &name) {
    const size_t dot = name.find_last_of('.');
    return dot == std::string::npos ? std::string() : name.substr(dot);
}

static bool isPlain(const std::string &name) {
    const std::string suffix = fileNameExtension(name);
    return suffix != ".gz" && suffix != ".zst" && suffix != ".zstd" && suffix != ".bz2"
           && suffix != ".xz";
}

static bool isZstd(const std::string &name) {
    const std::string suffix = fileNameExtension(name);
    return suffix == ".zst" || suffix == ".zstd";
}

static uint32_t littleEndian32(const unsigned char *at) {
    return at[0] | (at[1] << 8) | (at[2] << 16) | ((uint32_t) at[3] << 24);
}

static const size_t ZSTD_BLOCK_BUDGET = 1024;

static uint64_t zstdFrameSize(FILE *in, uint64_t at, uint64_t fileSize, size_t &budget) {
    unsigned char head[ZSTD_FRAMEHEADERSIZE_MAX];
    const size_t want = (size_t) std::min<uint64_t>(sizeof(head), fileSize - at);
    if (want < 8 || fseeko(in, (off_t) at, SEEK_SET) != 0 || fread(head, 1, want, in) != want) {
        return 0;
    }
    if ((littleEndian32(head) & 0xFFFFFFF0u) == 0x184D2A50u) {
        return 8 + littleEndian32(head + 4);
    }
    ZSTD_frameHeader header;
    if (ZSTD_getFrameHeader(&header, head, want) != 0) {
        return 0;
    }
    uint64_t span = header.headerSize;
    while (true) {
        unsigned char block[3];
        if (budget == 0) {
            return 0;
        }
        budget--;
        if (fseeko(in, (off_t) (at + span), SEEK_SET) != 0 || fread(block, 1, 3, in) != 3) {
            return 0;
        }
        const uint32_t described = block[0] | (block[1] << 8) | ((uint32_t) block[2] << 16);
        const bool last = (described & 1) != 0;
        const uint32_t type = (described >> 1) & 3;
        const uint32_t size = described >> 3;
        if (type == 3) {
            return 0;
        }
        span += 3 + (type == 1 ? 1 : size);
        if (last) {
            break;
        }
    }
    if (header.checksumFlag) {
        span += 4;
    }
    return (at + span <= fileSize) ? span : 0;
}

static std::vector<uint64_t> zstdFrames(const std::string &name) {
    std::vector<uint64_t> frames;
    FILE *in = fopen(name.c_str(), "rb");
    if (in == NULL) {
        Debug(Debug::ERROR) << "Cannot open " << name << "\n";
        EXIT(EXIT_FAILURE);
    }
    const uint64_t size = FileUtil::getFileSize(name);
    size_t budget = ZSTD_BLOCK_BUDGET;
    uint64_t at = 0;
    while (at < size) {
        const uint64_t span = zstdFrameSize(in, at, size, budget);
        if (span == 0) {
            break;
        }
        frames.push_back(at);
        at += span;
    }
    fclose(in);
    if (frames.empty()) {
        frames.push_back(0);
    }
    return frames;
}

std::vector<InputSplit> planInputSplits(const std::vector<std::string> &filenames, size_t want) {
    std::vector<uint64_t> sizes(filenames.size(), 0);
    uint64_t total = 0;
    for (size_t i = 0; i < filenames.size(); i++) {
        sizes[i] = std::max<uint64_t>(FileUtil::getFileSize(filenames[i]), 1);
        total += sizes[i];
    }
    const uint64_t per = std::max<uint64_t>(1, (total + std::max<size_t>(want, 1) - 1) / std::max<size_t>(want, 1));

    std::vector<InputSplit> chunks;
    for (size_t i = 0; i < filenames.size(); i++) {
        if (isZstd(filenames[i])) {
            const std::vector<uint64_t> frames = zstdFrames(filenames[i]);
            size_t begin = 0;
            while (begin < frames.size()) {
                size_t end = begin + 1;
                while (end < frames.size() && frames[end] - frames[begin] < per) {
                    end++;
                }
                InputSplit chunk;
                chunk.file = i;
                chunk.from = frames[begin];
                chunk.until = (end < frames.size()) ? frames[end] : sizes[i];
                chunk.compressed = true;
                chunks.push_back(chunk);
                begin = end;
            }
            continue;
        }
        if (isPlain(filenames[i]) == false || sizes[i] <= per) {
            InputSplit whole;
            whole.file = i;
            whole.from = 0;
            whole.until = sizes[i];
            whole.compressed = false;
            chunks.push_back(whole);
            continue;
        }
        for (uint64_t at = 0; at < sizes[i]; at += per) {
            InputSplit chunk;
            chunk.file = i;
            chunk.from = at;
            chunk.until = std::min<uint64_t>(at + per, sizes[i]);
            chunk.compressed = false;
            chunks.push_back(chunk);
        }
    }
    return chunks;
}

InputSplitReader::InputSplitReader(const std::string &filename, const InputSplit &chunk)
    : chunk(chunk), name(filename), buffer(0), at(0), filled(0), readTo(chunk.from),
      started(chunk.from == 0), whole(NULL), stream(NULL), packedAt(0), packedFilled(0),
      packedTo(chunk.from), fileSize(0), endsAt(chunk.until) {
    if (isPlain(filename) == false && chunk.compressed == false) {
        whole = KSeqFactory(filename.c_str());
        fd = -1;
        return;
    }
    buffer.resize(1u << 20);
    fd = ::open(filename.c_str(), O_RDONLY);
    if (fd < 0) {
        Debug(Debug::ERROR) << "Cannot open " << filename << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (chunk.compressed) {
        ZSTD_DStream *ds = ZSTD_createDStream();
        if (ds == NULL || ZSTD_isError(ZSTD_initDStream(ds))) {
            Debug(Debug::ERROR) << "Cannot start decompressing " << filename << "\n";
            EXIT(EXIT_FAILURE);
        }
        stream = ds;
        packed.resize(std::max<size_t>(ZSTD_DStreamInSize(), PACKED_READ_BYTES));
        fileSize = FileUtil::getFileSize(filename);
        endsAt = UINT64_MAX;
        readTo = 0;
#ifdef HAVE_POSIX_FADVISE
        posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
        return;
    }
#ifdef HAVE_POSIX_FADVISE
    posix_fadvise(fd, static_cast<off_t>(chunk.from), static_cast<off_t>(chunk.bytes()),
                  POSIX_FADV_SEQUENTIAL);
#endif
}

InputSplitReader::~InputSplitReader() {
    if (fd >= 0) {
        ::close(fd);
    }
    if (stream != NULL) {
        ZSTD_freeDStream(static_cast<ZSTD_DStream *>(stream));
    }
    delete whole;
}

bool InputSplitReader::fill() {
    if (at < filled) {
        return true;
    }
    if (stream != NULL) {
        ZSTD_outBuffer out = {buffer.data(), buffer.size(), 0};
        while (out.pos == 0) {
            if (packedAt == packedFilled) {
                const uint64_t left = fileSize - packedTo;
                if (left == 0) {
                    break;
                }
                const size_t want = std::min<uint64_t>(packed.size(), left);
                const ssize_t got = pread(fd, packed.data(), want, static_cast<off_t>(packedTo));
                if (got <= 0) {
                    break;
                }
                packedAt = 0;
                packedFilled = static_cast<size_t>(got);
                packedTo += packedFilled;
            }
            ZSTD_inBuffer in = {packed.data() + packedAt, packedFilled - packedAt, 0};
            const size_t code = ZSTD_decompressStream(static_cast<ZSTD_DStream *>(stream), &out, &in);
            if (ZSTD_isError(code)) {
                Debug(Debug::ERROR) << "Cannot decompress " << name << ": "
                                    << ZSTD_getErrorName(code) << "\n";
                EXIT(EXIT_FAILURE);
            }
            packedAt += in.pos;
            if (code == 0 && endsAt == UINT64_MAX
                && packedTo - (packedFilled - packedAt) >= chunk.until) {
                endsAt = readTo + out.pos;
            }
        }
        at = 0;
        filled = out.pos;
        readTo += filled;
        return filled > 0;
    }
    const ssize_t got = pread(fd, buffer.data(), buffer.size(), static_cast<off_t>(readTo));
    if (got < 0) {
        Debug(Debug::ERROR) << "Cannot read " << name << " at " << readTo << "\n";
        EXIT(EXIT_FAILURE);
    }
    at = 0;
    filled = static_cast<size_t>(got);
    readTo += filled;
    return filled > 0;
}

bool InputSplitReader::next(const char *&headerOut, size_t &headerLength, const char *&sequenceOut,
                            size_t &length) {
    if (whole != NULL) {
        if (whole->ReadEntry() == false) {
            return false;
        }
        headerOut = whole->entry.name.s;
        headerLength = whole->entry.name.l;
        sequenceOut = whole->entry.sequence.s;
        length = whole->entry.sequence.l;
        return true;
    }
    if (started == false) {
        started = true;
        char previous = 0;
        bool found = false;
        while (found == false && fill()) {
            while (at < filled) {
                if (buffer[at] == '>' && previous == '\n') {
                    found = true;
                    break;
                }
                previous = buffer[at];
                at++;
            }
        }
        if (found == false) {
            return false;
        }
    }

    const uint64_t beginsAt = readTo - filled + at;
    if (fill() == false || beginsAt > endsAt) {
        return false;
    }
    if (buffer[at] != '>') {
        return false;
    }

    header.clear();
    sequence.clear();
    at++;
    bool inHeader = true;
    while (fill()) {
        const size_t begin = at;
        while (at < filled && buffer[at] != '\n' && buffer[at] != '\r'
               && (inHeader || buffer[at] != '>')) {
            at++;
        }
        std::string &into = inHeader ? header : sequence;
        into.append(buffer.data() + begin, at - begin);
        if (at == filled) {
            continue;
        }
        if (inHeader == false && buffer[at] == '>') {
            break;
        }
        at++;
        if (inHeader) {
            inHeader = false;
        }
    }
    const size_t blank = header.find_first_of(" \t");
    headerOut = header.c_str();
    headerLength = blank == std::string::npos ? header.size() : blank;
    sequenceOut = sequence.c_str();
    length = sequence.size();
    return true;
}

void publishAllAtomically(std::vector<std::pair<std::string, std::string> > &pending,
                          unsigned int threads) {
#pragma omp parallel for schedule(dynamic, 4) num_threads(threads)
    for (size_t i = 0; i < pending.size(); i++) {
        const int fd = open(pending[i].first.c_str(), O_RDONLY);
        if (fd < 0 || fsync(fd) != 0 || close(fd) != 0) {
            Debug(Debug::ERROR) << "Cannot flush " << pending[i].first << " to storage\n";
            EXIT(EXIT_FAILURE);
        }
    }
    for (size_t i = 0; i < pending.size(); i++) {
        if (rename(pending[i].first.c_str(), pending[i].second.c_str()) != 0) {
            Debug(Debug::ERROR) << "Cannot publish " << pending[i].first << " as "
                                << pending[i].second << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    pending.clear();
}

void requireArena(const std::string &what, size_t bytes, size_t budget,
                  const std::string &narrower) {
    if (bytes <= budget) {
        return;
    }
    Debug(Debug::ERROR) << what << " holds " << ByteParser::format(bytes, 'a', 'h')
                        << " and the memory limit leaves " << ByteParser::format(budget, 'a', 'h')
                        << ". Raise --split-memory-limit past " << ByteParser::format(bytes, 'a', 'h')
                        << ", or " << narrower << ". Every machine reads the whole of one of these, "
                        << "so running on more of them does not make it smaller\n";
    EXIT(EXIT_FAILURE);
}

void writeBucketManifest(const std::string &path, const std::vector<uint64_t> &counts,
                         const std::vector<uint64_t> &bytes, const std::vector<uint64_t> &indexBytes,
                         const std::string &spanKey, uint64_t spanBegin, uint64_t spanEnd) {
    const std::string tmp = path + ".tmp";
    FILE *out = FileUtil::openAndDelete(tmp.c_str(), "w");
    fprintf(out, "#%sBegin\t%zu\n#%sEnd\t%zu\n#buckets\t%zu\n", spanKey.c_str(),
            static_cast<size_t>(spanBegin), spanKey.c_str(), static_cast<size_t>(spanEnd),
            counts.size());
    for (size_t i = 0; i < counts.size(); i++) {
        if (counts[i] > 0) {
            fprintf(out, "%zu\t%zu\t%zu\t%zu\n", i, static_cast<size_t>(counts[i]), static_cast<size_t>(bytes[i]),
                    static_cast<size_t>(indexBytes[i]));
        }
    }
    if (fclose(out) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << tmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::publishAtomically(tmp, path);
}

size_t readBucketManifests(const std::string &prefix, size_t chunks,
                           std::vector<uint64_t> &bytesInto, std::vector<uint64_t> &indexBytesInto, uint64_t *resumeAt) {
    size_t done = 0;
    for (; done < chunks; done++) {
        const std::string path = prefix + "." + SSTR(done) + ".manifest";
        FILE *in = fopen(path.c_str(), "r");
        if (in == NULL) {
            break;
        }
        char line[256];
        while (fgets(line, sizeof(line), in) != NULL) {
            size_t spanEnd = 0;
            if (resumeAt != NULL && sscanf(line, "#%*[^E]End\t%zu", &spanEnd) == 1) {
                *resumeAt = spanEnd;
            }
            size_t width = 0;
            if (sscanf(line, "#buckets\t%zu", &width) == 1 && width != bytesInto.size()) {
                Debug(Debug::ERROR) << path << " was written with " << width << " of them and this "
                                    << "run has " << bytesInto.size() << ". Start over rather than "
                                    << "resume: the counts would name different files\n";
                EXIT(EXIT_FAILURE);
            }
            if (line[0] == '#') {
                continue;
            }
            size_t bucket = 0;
            size_t count = 0;
            size_t bytes = 0;
            size_t indexBytes = 0;
            if (sscanf(line, "%zu\t%zu\t%zu\t%zu", &bucket, &count, &bytes, &indexBytes) != 4 || bucket >= bytesInto.size()) {
                Debug(Debug::ERROR) << path << " names " << bucket << " of " << bytesInto.size() << "\n";
                EXIT(EXIT_FAILURE);
            }
            bytesInto[bucket] += bytes;
            indexBytesInto[bucket] += indexBytes;
        }
        fclose(in);
    }
    return done;
}

static const uint64_t COUNTS_MAGIC = 0x4C494E43434E5453ull;

void writeSubBucketCounts(const std::string &path, const std::vector<uint64_t> &base,
                          const std::vector<std::vector<uint64_t> > &perThread, size_t chunks) {
    std::vector<uint64_t> total(base);
    for (size_t thread = 0; thread < perThread.size(); thread++) {
        for (size_t i = 0; i < total.size(); i++) {
            total[i] += perThread[thread][i];
        }
    }
    const uint64_t header[3] = {COUNTS_MAGIC, chunks, total.size()};
    const std::string tmp = path + ".tmp";
    FILE *out = FileUtil::openAndDelete(tmp.c_str(), "w");
    if (fwrite(header, sizeof(uint64_t), 3, out) != 3
        || fwrite(total.data(), sizeof(uint64_t), total.size(), out) != total.size()) {
        Debug(Debug::ERROR) << "Cannot write " << tmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fclose(out) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << tmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::publishAtomically(tmp, path);
}

std::vector<uint64_t> readSubBucketCounts(const std::string &path, size_t entries, size_t &chunks) {
    std::vector<uint64_t> counts(entries, 0);
    chunks = 0;
    FILE *in = fopen(path.c_str(), "r");
    if (in == NULL) {
        return counts;
    }
    uint64_t header[3] = {0, 0, 0};
    if (fread(header, sizeof(uint64_t), 3, in) != 3 || header[0] != COUNTS_MAGIC
        || header[2] != entries) {
        Debug(Debug::ERROR) << path << " is not a count table of " << entries << " entries\n";
        EXIT(EXIT_FAILURE);
    }
    if (fread(counts.data(), sizeof(uint64_t), entries, in) != entries) {
        Debug(Debug::ERROR) << path << " is shorter than the " << entries << " counts it holds\n";
        EXIT(EXIT_FAILURE);
    }
    fclose(in);
    chunks = header[1];
    return counts;
}

BucketCounts::BucketCounts(const std::string &prefix, unsigned int nodes, size_t subBuckets,
                           size_t buckets)
    : files(nodes, NULL), paths(nodes), subBuckets(subBuckets) {
    for (unsigned int node = 0; node < nodes; node++) {
        paths[node] = prefix + "." + SSTR(node) + ".counts";
        files[node] = fopen(paths[node].c_str(), "r");
        if (files[node] == NULL) {
            Debug(Debug::ERROR) << "Cannot open " << paths[node] << ", which node " << node
                                << " should have written\n";
            EXIT(EXIT_FAILURE);
        }
        uint64_t header[3] = {0, 0, 0};
        if (fread(header, sizeof(uint64_t), 3, files[node]) != 3 || header[0] != COUNTS_MAGIC
            || header[2] != subBuckets * buckets) {
            Debug(Debug::ERROR) << paths[node] << " is not a count table for this pass\n";
            EXIT(EXIT_FAILURE);
        }
    }
}

BucketCounts::~BucketCounts() {
    for (size_t node = 0; node < files.size(); node++) {
        if (files[node] != NULL) {
            fclose(files[node]);
        }
    }
}

std::vector<uint64_t> BucketCounts::of(size_t bucket, size_t node) const {
    std::vector<uint64_t> counts(subBuckets, 0);
    const long at = static_cast<long>((3 + bucket * subBuckets) * sizeof(uint64_t));
    if (fseek(files[node], at, SEEK_SET) != 0
        || fread(counts.data(), sizeof(uint64_t), subBuckets, files[node]) != subBuckets) {
        Debug(Debug::ERROR) << "Cannot read bucket " << bucket << " counts from " << paths[node]
                            << "\n";
        EXIT(EXIT_FAILURE);
    }
    return counts;
}

std::vector<uint64_t> BucketCounts::of(size_t bucket) const {
    std::vector<uint64_t> counts(subBuckets, 0);
    std::vector<uint64_t> one(subBuckets, 0);
    const long at = static_cast<long>((3 + bucket * subBuckets) * sizeof(uint64_t));
    for (size_t node = 0; node < files.size(); node++) {
        if (fseek(files[node], at, SEEK_SET) != 0
            || fread(one.data(), sizeof(uint64_t), subBuckets, files[node]) != subBuckets) {
            Debug(Debug::ERROR) << "Cannot read bucket " << bucket << " counts from " << paths[node]
                                << "\n";
            EXIT(EXIT_FAILURE);
        }
        for (size_t i = 0; i < counts.size(); i++) {
            counts[i] += one[i];
        }
    }
    return counts;
}

void markNodeDone(const std::string &path, unsigned int node) {
    const std::string done = nodeDonePath(path, node);
    const std::string tmp = done + ".tmp";
    FILE *out = FileUtil::openAndDelete(tmp.c_str(), "w");
    fprintf(out, "%u\n", node);
    if (fclose(out) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << tmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::publishAtomically(tmp, done);
}

void publishProgress(const std::string &path, uint64_t value) {
    const size_t bytes = 21;
    char line[bytes + 1];
    snprintf(line, sizeof(line), "%020llu\n", (unsigned long long) value);
    const int fd = open(path.c_str(), O_WRONLY | O_CREAT | O_SYNC, 0644);
    if (fd < 0) {
        Debug(Debug::ERROR) << "Cannot open " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (pwrite(fd, line, bytes, 0) != (ssize_t) bytes) {
        Debug(Debug::ERROR) << "Cannot write " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (close(fd) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
}

void dropConsumed(const std::string &prefix, unsigned int nodes, size_t from, size_t until,
                  size_t step) {
    for (size_t at = from; at < until; at += step) {
        for (unsigned int node = 0; node < nodes; node++) {
            const std::string path = prefix + "." + SSTR(node) + "." + SSTR(at);
            if (FileUtil::fileExists(path.c_str())) {
                FileUtil::remove(path.c_str());
            }
            if (FileUtil::fileExists((path + ".idx").c_str())) {
                FileUtil::remove((path + ".idx").c_str());
            }
        }
    }
}

void requireEveryNodeDone(const std::string &path, unsigned int nodes) {
    for (unsigned int node = 0; node < nodes; node++) {
        const std::string done = path + "." + SSTR(node) + ".done";
        if (FileUtil::fileExists(done.c_str()) == false) {
            Debug(Debug::ERROR) << "Node " << node << " of " << nodes
                                << " has not finished: " << done << " is missing\n";
            EXIT(EXIT_FAILURE);
        }
    }
}

void waitNodeDone(const std::string &path, unsigned int node, unsigned int limitSeconds) {
    const std::string done = nodeDonePath(path, node);
    unsigned int waited = 0;
    while (FileUtil::fileExists(done.c_str()) == false) {
        if (waited >= limitSeconds) {
            Debug(Debug::ERROR) << done << " did not appear within " << limitSeconds << "s\n";
            EXIT(EXIT_FAILURE);
        }
        static time_t lastWaitLog = 0;
        const time_t nowSec = time(NULL);
        if (waited >= 30 && nowSec - lastWaitLog >= 60) {
            Debug(Debug::INFO) << "Still waiting for " << done << " (" << waited << "s)\n";
            lastWaitLog = nowSec;
        }
        sleep(1);
        waited++;
    }
}
void waitEveryNodeDone(const std::string &path, unsigned int nodes, unsigned int limitSeconds) {
    for (unsigned int node = 0; node < nodes; node++) {
        waitNodeDone(path, node, limitSeconds);
    }
}

std::vector<size_t> nodeFileSlots(const SequenceLocator &runs, const NodePlacement &node,
                                  const std::function<uint64_t(uint32_t)> &costOfLength) {
    std::vector<uint64_t> weight(runs.filesPerNode(), 0);
    uint64_t total = 0;
    for (size_t i = 0; i < runs.size(); i++) {
        const uint64_t ranks = runs.rankEnd(i) - runs[i].rankBase();
        const uint64_t span = costOfLength
                                  ? ranks * costOfLength(runs[i].seqLen())
                                  : ranks * runs[i].seqLen();
        weight[runs[i].fileIdx() % runs.filesPerNode()] += span;
        total += span;
    }
    std::vector<size_t> mine;
    uint64_t at = 0;
    for (size_t fileSlot = 0; fileSlot < weight.size(); fileSlot++) {
        const uint64_t middle = at + weight[fileSlot] / 2;
        const unsigned int owner = static_cast<unsigned int>(
            std::min<uint64_t>(total == 0 ? 0 : middle * node.count / total, node.count - 1));
        if (owner == node.index) {
            mine.push_back(fileSlot);
        }
        at += weight[fileSlot];
    }
    return mine;
}

void preadFully(int fd, void *into, size_t bytes, uint64_t at, const std::string &what) {
    size_t got = 0;
    while (got < bytes) {
        const ssize_t read = pread(fd, static_cast<char *>(into) + got, bytes - got, static_cast<off_t>(at + got));
        if (read <= 0) {
            Debug(Debug::ERROR) << "Cannot read " << what << "\n";
            EXIT(EXIT_FAILURE);
        }
        got += static_cast<size_t>(read);
    }
}

unsigned int adjacentBitsFor(unsigned int classCount) {
    uint64_t span = 1;
    for (unsigned int slot = 0; slot < KmerRecord::ADJACENT_COUNT; slot++) {
        span *= classCount;
    }
    return bitsFor(span - 1);
}

static uint64_t packAdjacentClasses(const KmerRecord &record, unsigned int classCount) {
    uint64_t packed = 0;
    for (unsigned int slot = KmerRecord::ADJACENT_COUNT; slot-- > 0;) {
        packed = packed * classCount + record.adjacentAt(slot);
    }
    return packed;
}

static uint64_t unpackAdjacentClasses(uint64_t packed, unsigned int classCount) {
    uint64_t adjacent = 0;
    for (unsigned int slot = 0; slot < KmerRecord::ADJACENT_COUNT; slot++) {
        adjacent |= (packed % classCount) << (slot * KmerRecord::ADJACENT_BITS);
        packed /= classCount;
    }
    return adjacent;
}


void KmerSlabCodec::encode(const std::vector<KmerRecord> &records, std::vector<unsigned char> &out,
                           SlabHeader &header) const {
    memset(&header, 0, sizeof(header));
    header.magic = MAGIC;
    header.records = records.size();
    header.rankBits = rankBits;
    header.classCount = classCount;
    const unsigned int adjacentBits = adjacentBitsFor(classCount);
    out.clear();
    unsigned int sub = 0;
    uint64_t keys[SlabHeader::BLOCK_RECORDS];
    for (size_t at = 0; at < records.size();) {
        const unsigned int here = records[at].subBucket();
        markSubStart(header, sub, here, out.size());
        size_t end = at + 1;
        while (end < records.size() && end - at < SlabHeader::BLOCK_RECORDS && records[end].subBucket() == here) {
            end++;
        }
        uint64_t widestPos = 0;
        for (size_t i = at; i < end; i++) {
            keys[i - at] = records[i].key();
            widestPos = std::max(widestPos, records[i].pos());
        }
        BlockHeader block;
        countKeys(keys, end - at, block);
        block.widthBits = static_cast<uint8_t>(bitsFor(widestPos));
        putBlockHeader(out, block);
        BitWriter bits(out);
        putKeys(bits, keys, end - at, KmerRecord::KEY_BITS, block.keyBits);
        for (size_t i = at; i < end; i++) {
            bits.put(records[i].rank(), rankBits);
        }
        for (size_t i = at; i < end; i++) {
            bits.put(records[i].pos(), block.widthBits);
        }
        for (size_t i = at; i < end; i++) {
            bits.put(packAdjacentClasses(records[i], classCount), adjacentBits);
        }
        bits.finish();
        at = end;
    }
    markSubStart(header, sub, SlabHeader::SUBS, out.size());
}

size_t KmerSlabCodec::decode(const unsigned char *from, const unsigned char *to, const SlabHeader &header,
                             KmerRecord *into, size_t capacity) {
    const unsigned int adjacentBits = adjacentBitsFor(header.classCount);
    size_t count = 0;
    while (from < to) {
        BlockHeader block;
        memcpy(&block, from, sizeof(block));
        from += sizeof(block);
        BitReader bits(from);
        from += block.payloadBytes(KmerRecord::KEY_BITS, header.rankBits + block.widthBits + adjacentBits);
        const size_t n = block.records;
        if (count + n > capacity) {
            Debug(Debug::ERROR) << "A k-mer slab holds more records than its bucket counts allow\n";
            EXIT(EXIT_FAILURE);
        }
        uint64_t key[SlabHeader::BLOCK_RECORDS];
        uint64_t rank[SlabHeader::BLOCK_RECORDS];
        uint64_t pos[SlabHeader::BLOCK_RECORDS];
        getKeys(bits, key, n, KmerRecord::KEY_BITS, block.keyBits);
        for (size_t i = 0; i < n; i++) {
            rank[i] = bits.get(header.rankBits);
        }
        for (size_t i = 0; i < n; i++) {
            pos[i] = bits.get(block.widthBits);
        }
        for (size_t i = 0; i < n; i++) {
            into[count++].set(key[i], rank[i], pos[i], unpackAdjacentClasses(bits.get(adjacentBits), header.classCount));
        }
    }
    return count;
}

void PairSlabCodec::encode(const std::vector<PairRecord> &records, std::vector<unsigned char> &out,
                           SlabHeader &header) const {
    memset(&header, 0, sizeof(header));
    header.magic = MAGIC;
    header.records = records.size();
    header.rankBits = rankBits;
    out.clear();
    unsigned int sub = 0;
    uint64_t reps[SlabHeader::BLOCK_RECORDS];
    for (size_t at = 0; at < records.size();) {
        const unsigned int here = PairRecord::repRankSubBlockOf(records[at].rep(), ranks, repRankBlocks);
        markSubStart(header, sub, here, out.size());
        size_t end = at + 1;
        while (end < records.size() && end - at < SlabHeader::BLOCK_RECORDS
               && PairRecord::repRankSubBlockOf(records[end].rep(), ranks, repRankBlocks) == here) {
            end++;
        }
        uint64_t widestDiagonal = 0;
        for (size_t i = at; i < end; i++) {
            reps[i - at] = records[i].rep();
            widestDiagonal = std::max(widestDiagonal, zigzag(records[i].diagonal()));
        }
        BlockHeader block;
        countKeys(reps, end - at, block);
        block.widthBits = static_cast<uint8_t>(bitsFor(widestDiagonal));
        putBlockHeader(out, block);
        BitWriter bits(out);
        putKeys(bits, reps, end - at, rankBits, block.keyBits);
        for (size_t i = at; i < end; i++) {
            bits.put(records[i].member(), rankBits);
        }
        for (size_t i = at; i < end; i++) {
            bits.put(zigzag(records[i].diagonal()), block.widthBits);
        }
        bits.finish();
        at = end;
    }
    markSubStart(header, sub, SlabHeader::SUBS, out.size());
}

size_t PairSlabCodec::decode(const unsigned char *from, const unsigned char *to, const SlabHeader &header,
                             PairRecord *into, size_t capacity) {
    size_t count = 0;
    while (from < to) {
        BlockHeader block;
        memcpy(&block, from, sizeof(block));
        from += sizeof(block);
        BitReader bits(from);
        from += block.payloadBytes(header.rankBits, header.rankBits + block.widthBits);
        const size_t n = block.records;
        if (count + n > capacity) {
            Debug(Debug::ERROR) << "A pair slab holds more records than its bucket counts allow\n";
            EXIT(EXIT_FAILURE);
        }
        uint64_t rep[SlabHeader::BLOCK_RECORDS];
        uint64_t member[SlabHeader::BLOCK_RECORDS];
        getKeys(bits, rep, n, header.rankBits, block.keyBits);
        for (size_t i = 0; i < n; i++) {
            member[i] = bits.get(header.rankBits);
        }
        for (size_t i = 0; i < n; i++) {
            into[count++].set(rep[i], member[i], static_cast<int>(unzigzag(bits.get(block.widthBits))));
        }
    }
    return count;
}
