#include "Lin8Db.h"
#include "Lin8DbReader.h"
#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "FileUtil.h"
#include "Util.h"
#include "KSeqWrapper.h"
#include "NodePlacement.h"
#include "Timer.h"
#include "FastSort.h"

#include <algorithm>
#include <cctype>
#include <climits>
#include <cstring>
#include <fcntl.h>
#include <unistd.h>
#include <vector>
#include <zstd.h>

#ifdef OPENMP
#include <omp.h>
#endif

static const unsigned int LINCLUSTERDB_MAX_SEQ_LEN = SequenceLocator::MAX_ENTRY_LEN;
static const size_t LINCLUSTERDB_HISTOGRAM_SIZE = LINCLUSTERDB_MAX_SEQ_LEN + 1;
__extension__ typedef unsigned __int128 WideRoom;
static const uint64_t LINCLUSTERDB_HISTOGRAM_MAGIC = 0x4C494E4348495354ull;

static const uint64_t PROGRESS_STEP = 1000000;

static std::vector<InputSplit> inputSplitsForNode(const std::vector<std::string> &filenames,
                                             const NodePlacement &nodePlacement, unsigned int threads) {
    const size_t targetSplitCount =
        std::max<size_t>(1, (size_t) nodePlacement.count * std::max(threads, 1u) * 2);
    const std::vector<InputSplit> everySplit = planInputSplits(filenames, targetSplitCount);
    std::vector<uint64_t> splitByteOffsets(everySplit.size() + 1, 0);
    for (size_t i = 0; i < everySplit.size(); i++) {
        splitByteOffsets[i + 1] = splitByteOffsets[i] + std::max<uint64_t>(everySplit[i].bytes(), 1);
    }
    const uint64_t totalBytes = splitByteOffsets.back();
    std::vector<InputSplit> nodeInputSplits;
    for (size_t i = 0; i < everySplit.size(); i++) {
        const uint64_t splitMiddle = (splitByteOffsets[i] + splitByteOffsets[i + 1]) / 2;
        const unsigned int ownerNode = static_cast<unsigned int>(std::min<uint64_t>(
            totalBytes == 0 ? 0 : splitMiddle * nodePlacement.count / totalBytes,
            nodePlacement.count - 1));
        if (ownerNode == nodePlacement.index) {
            nodeInputSplits.push_back(everySplit[i]);
        }
    }
    return nodeInputSplits;
}

static std::string nodePartName(const std::string &db, const char *suffix, unsigned int node) {
    return db + "." + suffix + "." + SSTR(node);
}

struct NodeSequenceDistribution {
    std::vector<uint64_t> sequencesByLength;
    std::vector<uint64_t> sequencesByLengthAndSplit;
    std::vector<uint64_t> headerBytesByLengthAndSplit;
    unsigned int splitCount;
    uint64_t sequenceCount;
    uint64_t residueCount;
    uint64_t headerBytes;
    uint64_t rejectedSequenceCount;

    NodeSequenceDistribution()
        : splitCount(0), sequenceCount(0), residueCount(0), headerBytes(0),
          rejectedSequenceCount(0) {}

    uint64_t &sequencesAt(size_t length, unsigned int split) {
        return sequencesByLengthAndSplit[length * splitCount + split];
    }
    uint64_t sequencesAt(size_t length, unsigned int split) const {
        return sequencesByLengthAndSplit[length * splitCount + split];
    }
    uint64_t &headerBytesAt(size_t length, unsigned int split) {
        return headerBytesByLengthAndSplit[length * splitCount + split];
    }
    uint64_t headerBytesAt(size_t length, unsigned int split) const {
        return headerBytesByLengthAndSplit[length * splitCount + split];
    }
};

static NodeSequenceDistribution countSequenceLengths(const std::vector<std::string> &filenames,
                                   const std::vector<InputSplit> &nodeInputSplits, unsigned int threads,
                                   uint32_t maxSeqLen) {
    NodeSequenceDistribution nodeDistribution;
    nodeDistribution.splitCount = static_cast<unsigned int>(std::max<size_t>(nodeInputSplits.size(), 1));
    nodeDistribution.sequencesByLength.assign(LINCLUSTERDB_HISTOGRAM_SIZE, 0);
    nodeDistribution.sequencesByLengthAndSplit.assign(LINCLUSTERDB_HISTOGRAM_SIZE * nodeDistribution.splitCount, 0);
    nodeDistribution.headerBytesByLengthAndSplit.assign(LINCLUSTERDB_HISTOGRAM_SIZE * nodeDistribution.splitCount, 0);
    std::vector<std::pair<std::string, size_t> > skippedAll;
    std::vector<std::pair<std::string, size_t> > aloneAll;
    Debug::Progress progress(nodeInputSplits.size());
#pragma omp parallel num_threads(threads)
    {
        uint64_t sequenceCount = 0;
        uint64_t residues = 0;
        uint64_t headerBytes = 0;
        uint64_t rejected = 0;
        std::vector<std::pair<std::string, size_t> > skipped;
        std::vector<std::pair<std::string, size_t> > alone;
#pragma omp for schedule(dynamic, 1) nowait
        for (size_t split = 0; split < nodeInputSplits.size(); split++) {
            InputSplitReader reader(filenames[nodeInputSplits[split].file],
                                    nodeInputSplits[split]);
            const char *name = NULL;
            const char *seq = NULL;
            size_t nameLength = 0;
            size_t length = 0;
            while (reader.next(name, nameLength, seq, length)) {
                if (length == 0 || length > maxSeqLen) {
                    rejected++;
                    skipped.push_back(std::make_pair(std::string(name, nameLength), length));
                    continue;
                }
                if (length > SequenceLocator::MAX_SEQ_LEN) {
                    alone.push_back(std::make_pair(std::string(name, nameLength), length));
                }
                const size_t header = nameLength + 1;
                nodeDistribution.sequencesAt(length, split)++;
                nodeDistribution.headerBytesAt(length, split) += header;
                sequenceCount++;
                residues += length;
                headerBytes += header;
            }
            progress.updateProgress();
        }
#pragma omp critical
        {
            nodeDistribution.sequenceCount += sequenceCount;
            nodeDistribution.residueCount += residues;
            nodeDistribution.headerBytes += headerBytes;
            nodeDistribution.rejectedSequenceCount += rejected;
            skippedAll.insert(skippedAll.end(), skipped.begin(), skipped.end());
            aloneAll.insert(aloneAll.end(), alone.begin(), alone.end());
        }
    }
    std::sort(aloneAll.begin(), aloneAll.end());
    std::sort(skippedAll.begin(), skippedAll.end());
    for (size_t i = 0; i < aloneAll.size() && i < 10; i++) {
        Debug(Debug::WARNING) << "Clustered alone, too long for k-mers: " << aloneAll[i].first
                              << " (" << aloneAll[i].second << " residues)\n";
    }
    if (aloneAll.size() > 10) {
        Debug(Debug::WARNING) << "and " << (aloneAll.size() - 10) << " more too long for k-mers\n";
    }
    for (size_t i = 0; i < skippedAll.size() && i < 10; i++) {
        Debug(Debug::WARNING) << "Outside 1 to " << maxSeqLen << " residues, left out of the database: "
                            << skippedAll[i].first << " (" << skippedAll[i].second << " residues)\n";
    }
    if (skippedAll.size() > 10) {
        Debug(Debug::WARNING) << "and " << (skippedAll.size() - 10) << " more outside 1 to " << maxSeqLen << " residues\n";
    }
    for (size_t length = 0; length < LINCLUSTERDB_HISTOGRAM_SIZE; length++) {
        for (unsigned int split = 0; split < nodeDistribution.splitCount; split++) {
            nodeDistribution.sequencesByLength[length] += nodeDistribution.sequencesAt(length, split);
        }
    }
    return nodeDistribution;
}

static void writeLengthDistribution(const std::string &path, const NodeSequenceDistribution &nodeDistribution, unsigned int threads) {
    const std::string tmp = path + ".tmp" + uniqueTmpSuffix();
    FILE *out = FileUtil::openAndDelete(tmp.c_str(), "wb");
    const uint64_t fields[6] = {LINCLUSTERDB_HISTOGRAM_MAGIC, threads, nodeDistribution.sequenceCount,
                                nodeDistribution.residueCount, nodeDistribution.headerBytes, nodeDistribution.rejectedSequenceCount};
    if (fwrite(fields, sizeof(uint64_t), 6, out) != 6
        || fwrite(nodeDistribution.sequencesByLength.data(), sizeof(uint64_t), LINCLUSTERDB_HISTOGRAM_SIZE, out)
               != LINCLUSTERDB_HISTOGRAM_SIZE) {
        Debug(Debug::ERROR) << "Cannot write the counts of node to " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fclose(out) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << tmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::publishAtomically(tmp, path);
}

static NodeSequenceDistribution readLengthDistribution(const std::string &path, unsigned int threads) {
    FILE *in = fopen(path.c_str(), "rb");
    if (in == NULL) {
        Debug(Debug::ERROR) << "Cannot open the counts of node " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    NodeSequenceDistribution nodeDistribution;
    nodeDistribution.splitCount = 0;
    uint64_t fields[6];
    nodeDistribution.sequencesByLength.assign(LINCLUSTERDB_HISTOGRAM_SIZE, 0);
    if (fread(fields, sizeof(uint64_t), 6, in) != 6
        || fread(nodeDistribution.sequencesByLength.data(), sizeof(uint64_t), LINCLUSTERDB_HISTOGRAM_SIZE, in)
               != LINCLUSTERDB_HISTOGRAM_SIZE) {
        Debug(Debug::ERROR) << "The counts of node " << path << " are truncated\n";
        EXIT(EXIT_FAILURE);
    }
    fclose(in);
    if (fields[0] != LINCLUSTERDB_HISTOGRAM_MAGIC) {
        Debug(Debug::ERROR) << "File " << path << " is not a count table for a node\n";
        EXIT(EXIT_FAILURE);
    }
    if (fields[1] != threads) {
        Debug(Debug::ERROR) << "Histogram " << path << " was written with " << fields[1]
                            << " threads, this node uses " << threads
                            << ", every node must use the same value\n";
        EXIT(EXIT_FAILURE);
    }
    nodeDistribution.sequenceCount = fields[2];
    nodeDistribution.residueCount = fields[3];
    nodeDistribution.headerBytes = fields[4];
    nodeDistribution.rejectedSequenceCount = fields[5];
    return nodeDistribution;
}

struct RankLayout {
    std::vector<uint64_t> rankBase;
    uint64_t sequenceCount;
    uint64_t bytes;
    uint64_t headerBytes;
};

static RankLayout computeGlobalRankLayout(const std::vector<NodeSequenceDistribution> &nodeDistributions, unsigned int self) {
    RankLayout ranks;
    ranks.rankBase.assign(LINCLUSTERDB_HISTOGRAM_SIZE, 0);
    ranks.sequenceCount = 0;
    ranks.bytes = 0;
    ranks.headerBytes = 0;
    uint64_t rank = 0;
    for (size_t length = LINCLUSTERDB_MAX_SEQ_LEN; length >= 1; length--) {
        for (size_t node = 0; node < nodeDistributions.size(); node++) {
            if (node == self) {
                ranks.rankBase[length] = rank;
            }
            rank += nodeDistributions[node].sequencesByLength[length];
        }
    }
    ranks.sequenceCount = rank;
    for (size_t node = 0; node < nodeDistributions.size(); node++) {
        ranks.headerBytes += nodeDistributions[node].headerBytes;
    }
    for (size_t length = 1; length < LINCLUSTERDB_HISTOGRAM_SIZE; length++) {
        uint64_t total = 0;
        for (size_t node = 0; node < nodeDistributions.size(); node++) {
            total += nodeDistributions[node].sequencesByLength[length];
        }
        ranks.bytes += total * length;
    }
    return ranks;
}

struct BytePlan {
    std::vector<uint64_t> sequenceByteAt;
    std::vector<uint64_t> headerByteAt;
    std::vector<uint64_t> sequenceFileStart;
    std::vector<uint64_t> headerFileStart;
    std::vector<unsigned int> fileOfLength;
    unsigned int splitCount;
    uint64_t sequenceBytes;
    uint64_t headerBytesTotal;
};

static void planDatabaseLayout(const std::vector<NodeSequenceDistribution> &nodeDistributions, const NodeSequenceDistribution &nodeDistribution,
                       const RankLayout &ranks, unsigned int files, unsigned int fileBase,
                       BytePlan &bytePlan, SequenceLocator &nodeLocator) {
    std::vector<uint64_t> bytesAbove(LINCLUSTERDB_HISTOGRAM_SIZE, 0);
    uint64_t total = 0;
    for (size_t length = LINCLUSTERDB_MAX_SEQ_LEN; length >= 1; length--) {
        bytesAbove[length] = total;
        uint64_t count = 0;
        for (size_t node = 0; node < nodeDistributions.size(); node++) {
            count += nodeDistributions[node].sequencesByLength[length];
        }
        total += count * length;
    }
    const uint64_t perFile = std::max<uint64_t>(1, (total + files - 1) / files);

    const unsigned int splitCount = nodeDistribution.splitCount;
    bytePlan.splitCount = splitCount;
    bytePlan.sequenceByteAt.assign(LINCLUSTERDB_HISTOGRAM_SIZE * splitCount, 0);
    bytePlan.headerByteAt.assign(LINCLUSTERDB_HISTOGRAM_SIZE * splitCount, 0);
    bytePlan.sequenceFileStart.assign(files + 1, 0);
    bytePlan.headerFileStart.assign(files + 1, 0);
    bytePlan.fileOfLength.assign(LINCLUSTERDB_HISTOGRAM_SIZE, 0);

    uint64_t sequenceByteAt = 0;
    uint64_t headerByteAt = 0;
    unsigned int open = 0;
    for (size_t length = LINCLUSTERDB_MAX_SEQ_LEN; length >= 1; length--) {
        const unsigned int file = static_cast<unsigned int>(
            std::min<uint64_t>(bytesAbove[length] / perFile, files - 1));
        bytePlan.fileOfLength[length] = file;
        while (open < file) {
            open++;
            bytePlan.sequenceFileStart[open] = sequenceByteAt;
            bytePlan.headerFileStart[open] = headerByteAt;
        }
        const uint64_t seqFirst = sequenceByteAt;
        const uint64_t hdrFirst = headerByteAt;
        for (unsigned int split = 0; split < splitCount; split++) {
            bytePlan.sequenceByteAt[length * splitCount + split] = sequenceByteAt;
            bytePlan.headerByteAt[length * splitCount + split] = headerByteAt;
            sequenceByteAt += nodeDistribution.sequencesAt(length, split)
                     * length;
            headerByteAt += nodeDistribution.headerBytesAt(length, split);
        }
        if (nodeDistribution.sequencesByLength[length] > 0) {
            nodeLocator.append(ranks.rankBase[length], static_cast<uint32_t>(length),
                           seqFirst - bytePlan.sequenceFileStart[file], fileBase + file,
                           hdrFirst - bytePlan.headerFileStart[file]);
        }
    }
    while (open < files) {
        open++;
        bytePlan.sequenceFileStart[open] = sequenceByteAt;
        bytePlan.headerFileStart[open] = headerByteAt;
    }
    bytePlan.sequenceBytes = sequenceByteAt;
    bytePlan.headerBytesTotal = headerByteAt;
}

static void pwriteAll(int fd, const unsigned char *data, size_t bytes, off_t at, const char *what) {
    static const size_t WRITE_CHUNK = 1ull << 30;
    size_t done = 0;
    while (done < bytes) {
        const size_t chunk = std::min<size_t>(bytes - done, WRITE_CHUNK);
        const ssize_t wrote = pwrite(fd, data + done, chunk, at + static_cast<off_t>(done));
        if (wrote <= 0) {
            Debug(Debug::ERROR) << "Cannot write " << chunk << " byte to " << what << "\n";
            EXIT(EXIT_FAILURE);
        }
        done += static_cast<size_t>(wrote);
    }
}

// --compressed: zstd frames in finish order, a table keyed by uncompressed offset last
class HeaderPacker {
public:
    static const size_t FRAME_BYTES = 8ull << 20;
    static const int LEVEL = 3;

    HeaderPacker(const std::vector<int> &fd) : fd(fd), packedEnd(fd.size(), 0), frames(fd.size()) {}

    void append(unsigned int file, uint64_t rawAt, const unsigned char *packed, size_t packedSize, size_t rawSize) {
        const RunDbReader::HeaderFrame frame = {rawAt, __sync_fetch_and_add(&packedEnd[file], packedSize), rawSize, packedSize};
#pragma omp critical(header_frames)
        frames[file].push_back(frame);
        pwriteAll(fd[file], packed, packedSize, static_cast<off_t>(frame.packedAt), "header file");
    }

    void finish() {
        for (size_t file = 0; file < fd.size(); file++) {
            std::vector<RunDbReader::HeaderFrame> &table = frames[file];
            SORT_SERIAL(table.begin(), table.end(),
                        [](const RunDbReader::HeaderFrame &a, const RunDbReader::HeaderFrame &b) { return a.rawAt < b.rawAt; });
            const uint64_t tail[2] = {table.size(), RunDbReader::HEADER_FRAMES_MAGIC};
            const size_t bytes = table.size() * sizeof(RunDbReader::HeaderFrame);
            pwriteAll(fd[file], reinterpret_cast<const unsigned char *>(table.data()), bytes,
                      static_cast<off_t>(packedEnd[file]), "header file");
            pwriteAll(fd[file], reinterpret_cast<const unsigned char *>(tail), sizeof(tail),
                      static_cast<off_t>(packedEnd[file] + bytes), "header file");
        }
    }

private:
    const std::vector<int> &fd;
    std::vector<uint64_t> packedEnd;
    std::vector<std::vector<RunDbReader::HeaderFrame> > frames;
};

class SequenceWriter {
public:
    SequenceWriter(const std::vector<int> &seqFd, const std::vector<int> &hdrFd, HeaderPacker *packer,
               const NodeSequenceDistribution &nodeDistribution, const BytePlan &bytePlan, unsigned int file, size_t budget)
        : seqFd(seqFd), hdrFd(hdrFd), packer(packer), cctx(packer == NULL ? NULL : ZSTD_createCCtx()),
          bytePlan(bytePlan), file(file) {
        seqBuf.resize(LINCLUSTERDB_HISTOGRAM_SIZE);
        hdrBuf.resize(LINCLUSTERDB_HISTOGRAM_SIZE);
        seqRoom.assign(LINCLUSTERDB_HISTOGRAM_SIZE, 0);
        hdrRoom.assign(LINCLUSTERDB_HISTOGRAM_SIZE, 0);
        sequenceByteAt.assign(LINCLUSTERDB_HISTOGRAM_SIZE, 0);
        headerByteAt.assign(LINCLUSTERDB_HISTOGRAM_SIZE, 0);
        uint64_t want = 0;
        for (size_t length = 1; length < LINCLUSTERDB_HISTOGRAM_SIZE; length++) {
            sequenceByteAt[length] = bytePlan.sequenceByteAt[length * bytePlan.splitCount + file];
            headerByteAt[length] = bytePlan.headerByteAt[length * bytePlan.splitCount + file];
            seqRoom[length] = nodeDistribution.sequencesAt(length, file)
                              * length;
            hdrRoom[length] = nodeDistribution.headerBytesAt(length, file);
            want += seqRoom[length] + hdrRoom[length];
        }
        if (want > budget && want > 0) {
            const uint64_t floorRoom = (budget / 2) / LINCLUSTERDB_HISTOGRAM_SIZE;
            for (size_t length = 1; length < LINCLUSTERDB_HISTOGRAM_SIZE; length++) {
                const uint64_t share = budget / 2;
                seqRoom[length] = std::min<uint64_t>(
                    seqRoom[length],
                    std::max<uint64_t>((uint64_t) ((WideRoom) seqRoom[length] * share / want), floorRoom));
                hdrRoom[length] = std::min<uint64_t>(
                    hdrRoom[length],
                    std::max<uint64_t>((uint64_t) ((WideRoom) hdrRoom[length] * share / want), floorRoom));
            }
        }
    }

    ~SequenceWriter() {
        if (cctx != NULL) {
            ZSTD_freeCCtx(cctx);
        }
    }

    void add(size_t length, const unsigned char *sequence, const char *header, size_t headerLen) {
        seqBuf[length].insert(seqBuf[length].end(), sequence, sequence + length);
        hdrBuf[length].insert(hdrBuf[length].end(), header, header + headerLen);
        if (seqBuf[length].size() >= seqRoom[length] || hdrBuf[length].size() >= hdrRoom[length]) {
            flushLength(length);
        }
    }

    void flush() {
        for (size_t length = LINCLUSTERDB_MAX_SEQ_LEN; length >= 1; length--) {
            flushLength(length);
        }
    }

    uint64_t bytesWritten() const { return written; }

private:
    void flushLength(size_t length) {
        if (seqBuf[length].empty() == false) {
            writeLength(seqFd, bytePlan.sequenceFileStart, bytePlan.fileOfLength[length],
                        sequenceByteAt[length], seqBuf[length]);
            sequenceByteAt[length] += seqBuf[length].size();
            seqBuf[length].clear();
        }
        if (hdrBuf[length].empty() == false) {
            if (packer != NULL) {
                packHeaders(bytePlan.fileOfLength[length], headerByteAt[length], hdrBuf[length]);
            } else {
                writeLength(hdrFd, bytePlan.headerFileStart, bytePlan.fileOfLength[length],
                            headerByteAt[length], hdrBuf[length]);
            }
            headerByteAt[length] += hdrBuf[length].size();
            hdrBuf[length].clear();
        }
    }

    // frames end on a header boundary so a reader never needs two frames for one name
    void packHeaders(unsigned int file, uint64_t at, const std::vector<unsigned char> &data) {
        const std::vector<uint64_t> &fileStart = bytePlan.headerFileStart;
        if (at < fileStart[file] || at + data.size() > fileStart[file + 1]) {
            Debug(Debug::ERROR) << "Write of " << data.size() << " byte at " << at
                                << " does not fit header file " << file << ", which holds "
                                << fileStart[file] << " to " << fileStart[file + 1] << "\n";
            EXIT(EXIT_FAILURE);
        }
        size_t done = 0;
        while (done < data.size()) {
            size_t take = std::min<size_t>(data.size() - done, HeaderPacker::FRAME_BYTES);
            if (done + take < data.size()) {
                const void *cut = memrchr(data.data() + done, '\n', take);
                if (cut == NULL) {
                    cut = memchr(data.data() + done + take, '\n', data.size() - done - take);
                }
                take = static_cast<const unsigned char *>(cut) - (data.data() + done) + 1;
            }
            packed.resize(ZSTD_compressBound(take));
            const size_t got = ZSTD_compressCCtx(cctx, packed.data(), packed.size(), data.data() + done, take, HeaderPacker::LEVEL);
            if (ZSTD_isError(got)) {
                Debug(Debug::ERROR) << "Cannot compress headers: " << ZSTD_getErrorName(got) << "\n";
                EXIT(EXIT_FAILURE);
            }
            packer->append(file, at + done - bytePlan.headerFileStart[file], packed.data(), got, take);
            done += take;
        }
        written += data.size();
    }

    void writeLength(const std::vector<int> &fd, const std::vector<uint64_t> &fileStart,
                     unsigned int file, uint64_t at, const std::vector<unsigned char> &data) {
        if (at < fileStart[file] || at + data.size() > fileStart[file + 1]) {
            Debug(Debug::ERROR) << "Write of " << data.size() << " byte at " << at
                                << " does not fit data file " << file << ", which holds "
                                << fileStart[file] << " to " << fileStart[file + 1] << "\n";
            EXIT(EXIT_FAILURE);
        }
        pwriteAll(fd[file], data.data(), data.size(), static_cast<off_t>(at - fileStart[file]), "data file");
        written += data.size();
    }

    const std::vector<int> &seqFd;
    const std::vector<int> &hdrFd;
    HeaderPacker *packer;
    ZSTD_CCtx *cctx;
    std::vector<unsigned char> packed;
    const BytePlan &bytePlan;
    unsigned int file;
    uint64_t written = 0;
    std::vector<std::vector<unsigned char> > seqBuf;
    std::vector<std::vector<unsigned char> > hdrBuf;
    std::vector<uint64_t> seqRoom;
    std::vector<uint64_t> hdrRoom;
    std::vector<uint64_t> sequenceByteAt;
    std::vector<uint64_t> headerByteAt;
};

static void writeSequencesAndHeaders(const std::vector<std::string> &filenames, const std::vector<InputSplit> &nodeInputSplits,
                      const NodeSequenceDistribution &nodeDistribution, const BytePlan &bytePlan,
                      const std::vector<int> &seqFd, const std::vector<int> &hdrFd, HeaderPacker *packer,
                      unsigned int threads, size_t budget, uint32_t maxSeqLen) {
    uint64_t written = 0;
    Debug::Progress progress(nodeDistribution.sequenceCount / PROGRESS_STEP + 1);
#pragma omp parallel num_threads(threads)
    {
        std::string header;
#pragma omp for schedule(dynamic, 1)
        for (size_t split = 0; split < nodeInputSplits.size(); split++) {
            const size_t room = std::max<size_t>(budget / threads, 1u << 20);
            SequenceWriter writer(seqFd, hdrFd, packer, nodeDistribution, bytePlan,
                              static_cast<unsigned int>(split),
                              packer == NULL ? room : room - std::min(room / 2, HeaderPacker::FRAME_BYTES));
            InputSplitReader reader(filenames[nodeInputSplits[split].file],
                                    nodeInputSplits[split]);
            const char *name = NULL;
            const char *seq = NULL;
            size_t nameLength = 0;
            size_t length = 0;
            uint64_t placed = 0;
            while (reader.next(name, nameLength, seq, length)) {
                if (length == 0 || length > maxSeqLen) {
                    continue;
                }
                if (++placed % PROGRESS_STEP == 0) {
                    progress.updateProgress();
                }
                header.assign(name, nameLength);
                header.append(1, '\n');
                writer.add(length, reinterpret_cast<const unsigned char *>(seq), header.c_str(),
                           header.size());
            }
            if (placed % PROGRESS_STEP != 0) {
                progress.updateProgress();
            }
            writer.flush();
#pragma omp atomic
            written += writer.bytesWritten();
        }
    }
    const uint64_t expect = bytePlan.sequenceBytes + bytePlan.headerBytesTotal;
    if (written != expect) {
        Debug(Debug::ERROR) << "Placed " << written << " byte but the ranks reserved " << expect
                            << ", the output would hold unwritten gaps\n";
        EXIT(EXIT_FAILURE);
    }
}

static std::vector<int> openDataFiles(const std::string &prefix, unsigned int fileBase,
                                      unsigned int files,
                                      const std::vector<uint64_t> &fileStart, bool sized) {
    std::vector<int> fd(files, -1);
    for (unsigned int i = 0; i < files; i++) {
        const std::string path = prefix + "." + SSTR(fileBase + i);
        fd[i] = open(path.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0666);
        if (fd[i] < 0) {
            Debug(Debug::ERROR) << "Cannot open " << path << " for writing\n";
            EXIT(EXIT_FAILURE);
        }
        if (sized && ftruncate(fd[i], static_cast<off_t>(fileStart[i + 1] - fileStart[i])) != 0) {
            Debug(Debug::ERROR) << "Cannot size " << path << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    return fd;
}

static void dropCacheAndCloseFiles(std::vector<int> &fd) {
    for (size_t i = 0; i < fd.size(); i++) {
        if (fd[i] < 0) {
            continue;
        }
#if defined(__linux__)
        sync_file_range(fd[i], 0, 0, SYNC_FILE_RANGE_WAIT_BEFORE | SYNC_FILE_RANGE_WRITE
                                         | SYNC_FILE_RANGE_WAIT_AFTER);
        posix_fadvise(fd[i], 0, 0, POSIX_FADV_DONTNEED);
#endif
        if (close(fd[i]) != 0) {
            Debug(Debug::ERROR) << "Cannot close data file " << i << "\n";
            EXIT(EXIT_FAILURE);
        }
        fd[i] = -1;
    }
}
static void writeFileManifest(const std::string &db, const SequenceLocator &runs, unsigned int filesPerNode) {
    const size_t files = runs.fileCount();
    std::vector<uint32_t> maxLen(files, 0);
    std::vector<uint32_t> minLen(files, 0);
    std::vector<uint64_t> entries(files, 0);
    for (size_t i = 0; i < runs.size(); i++) {
        const uint32_t file = runs[i].fileIdx();
        if (entries[file] == 0) {
            maxLen[file] = runs[i].seqLen();
        }
        minLen[file] = runs[i].seqLen();
        entries[file] += runs.rankEnd(i) - runs[i].rankBase();
    }
    const std::string tmp = db + ".files.tmp" + uniqueTmpSuffix();
    FILE *out = FileUtil::openAndDelete(tmp.c_str(), "w");
    fprintf(out, "#file\tnode\tsuffix\tmaxLen\tminLen\tentries\n");
    for (size_t file = 0; file < files; file++) {
        fprintf(out, "%zu\t%zu\t%zu\t%u\t%u\t%zu\n", file, file / filesPerNode, file % filesPerNode,
                maxLen[file], minLen[file], static_cast<size_t>(entries[file]));
    }
    if (fclose(out) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << tmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::publishAtomically(tmp, db + ".files");
}

static SequenceLocator mergeSequenceLocators(const std::string &db, unsigned int nodeCount,
                               unsigned int filesPerNode, uint64_t sequenceCount) {
    std::vector<SequenceLocator::LengthRun> all;
    for (unsigned int node = 0; node < nodeCount; node++) {
        SequenceLocator part;
        part.read(nodePartName(db, "runs", node));
        all.insert(all.end(), part.data(), part.data() + part.size());
    }
    SORT_SERIAL(all.begin(), all.end(),
                  [](const SequenceLocator::LengthRun &first, const SequenceLocator::LengthRun &second) {
                      return first.rankBase() < second.rankBase();
                  });
    SequenceLocator merged;
    merged.setLayout(nodeCount, filesPerNode);
    merged.reserve(all.size());
    for (size_t i = 0; i < all.size(); i++) {
        merged.append(all[i].rankBase(), all[i].seqLen(), all[i].byteBase(), all[i].fileIdx(),
                      all[i].hdrBase());
    }
    merged.finish(sequenceCount);
    merged.write(db + ".runs");
    return merged;
}

int lin8createdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    const size_t inputCount = par.filenames.size() - 1;
    const std::string db = par.filenames.back();
    const std::string headerDb = db + "_h";
    std::vector<std::string> filenames(par.filenames.begin(), par.filenames.begin() + inputCount);
    SORT_SERIAL(filenames.begin(), filenames.end());

    const NodePlacement node = NodePlacement::resolve(par);
    const uint32_t maxSeqLen =
        static_cast<uint32_t>(std::min<size_t>(par.maxSeqLen, LINCLUSTERDB_MAX_SEQ_LEN));
    const bool isNodeDone = FileUtil::fileExists(nodeDonePath(db, node.index).c_str());
    std::vector<InputSplit> nodeInputSplits;
    size_t budget = 0;
    Timer timer;
    NodeSequenceDistribution nodeDistribution;
    if (isNodeDone) {
        nodeDistribution = readLengthDistribution(nodePartName(db, "hist", node.index), par.threads);
        Debug(Debug::INFO) << "Node " << node.index << " of " << node.count << " placed "
                           << nodeDistribution.sequenceCount << " sequences, merging\n";
    } else {
        nodeInputSplits = inputSplitsForNode(filenames, node, par.threads);
        budget = Util::computeMemory(par.splitMemoryLimit);
        Debug(Debug::INFO) << "Node " << node.index << " of " << node.count << " takes " << nodeInputSplits.size() << " split" << (nodeInputSplits.size() == 1 ? "" : "s")
                           << " of " << filenames.size() << " input file"
                           << (filenames.size() == 1 ? "" : "s") << ", buffer budget "
                           << (budget / (1024 * 1024)) << " MB\n";
        Debug(Debug::INFO) << "Counting the sequences and their lengths\n";
        nodeDistribution = countSequenceLengths(filenames, nodeInputSplits, par.threads, maxSeqLen);
        Debug(Debug::INFO) << "Counted " << nodeDistribution.sequenceCount << " sequences, " << nodeDistribution.residueCount
                           << " residues in " << timer.lap() << "\n";
        if (nodeDistribution.rejectedSequenceCount > 0) {
            Debug(Debug::WARNING) << "Skipped " << nodeDistribution.rejectedSequenceCount << " sequences outside 1 to "
                                  << maxSeqLen << " residues\n";
        }
        writeLengthDistribution(nodePartName(db, "hist", node.index), nodeDistribution, par.threads);
    }

    std::vector<NodeSequenceDistribution> nodeDistributions;
    if (node.count == 1) {
        nodeDistributions.push_back(nodeDistribution);
    } else {
        for (unsigned int other = 0; other < node.count; other++) {
            const std::string path = nodePartName(db, "hist", other);
            if (FileUtil::fileExists(path.c_str()) == false) {
                Debug(Debug::INFO) << "Wrote " << nodePartName(db, "hist", node.index)
                                   << ", waiting for the other nodes\n";
                return EXIT_SUCCESS;
            }
            nodeDistributions.push_back(readLengthDistribution(path, par.threads));
        }
    }
    const RankLayout ranks = computeGlobalRankLayout(nodeDistributions, node.index);
    if (ranks.sequenceCount == 0) {
        Debug(Debug::ERROR) << "The input files have no usable entry\n";
        EXIT(EXIT_FAILURE);
    }

    if (isNodeDone == false) {
        for (int file = 0; file < par.threads; file++) {
            const std::string data = db + "." + SSTR(node.index * par.threads + file);
            if (FileUtil::fileExists(data.c_str()) && FileUtil::getFileSize(data) > 0) {
                Debug(Debug::ERROR) << data << " holds data but node " << node.index
                                    << " left no marker. Another run wrote it and did not finish;"
                                    << " remove this database and build it again\n";
                EXIT(EXIT_FAILURE);
            }
        }
        BytePlan bytePlan;
        SequenceLocator nodeLocator;
        planDatabaseLayout(nodeDistributions, nodeDistribution, ranks, par.threads,
                           node.index * par.threads, bytePlan, nodeLocator);
        nodeLocator.setLayout(node.count, par.threads);
        nodeLocator.finish(ranks.sequenceCount);
        Debug(Debug::INFO) << "This node writes " << (bytePlan.sequenceBytes >> 30)
                           << " GB of sequences\n";

        timer.reset();
        const bool pack = par.compressed != 0;
        std::vector<int> seqFd = openDataFiles(db, node.index * par.threads, par.threads,
                                               bytePlan.sequenceFileStart, true);
        std::vector<int> hdrFd = openDataFiles(headerDb, node.index * par.threads, par.threads,
                                               bytePlan.headerFileStart, pack == false);
        HeaderPacker packer(hdrFd);
        Debug(Debug::INFO) << "Placing " << nodeDistribution.sequenceCount << " sequences by length\n";
        writeSequencesAndHeaders(filenames, nodeInputSplits, nodeDistribution, bytePlan, seqFd, hdrFd,
                                 pack ? &packer : NULL, par.threads, budget, maxSeqLen);
        if (pack) {
            packer.finish();
        }
        dropCacheAndCloseFiles(seqFd);
        dropCacheAndCloseFiles(hdrFd);
        nodeLocator.write(nodePartName(db, "runs", node.index));
        Debug(Debug::INFO) << "Placed " << nodeDistribution.sequenceCount << " sequences in " << timer.lap() << "\n";

        markNodeDone(db, node.index);
    }

    bool complete = true;
    for (unsigned int other = 0; other < node.count; other++) {
        complete = complete && FileUtil::fileExists(nodePartName(db, "runs", other).c_str());
    }
    if (complete && node.index == 0) {
        const SequenceLocator merged = mergeSequenceLocators(db, node.count, par.threads, ranks.sequenceCount);
        writeFileManifest(db, merged, par.threads);
        const int dbtype = DBReader<DBKeyType>::setExtendedDbtype(Parameters::DBTYPE_AMINO_ACIDS,
                                                                  Parameters::DBTYPE_EXTENDED_RUNS);
        const std::string dbtypeBase = db + ".new" + uniqueTmpSuffix();
        DBWriter::writeDbtypeFile(dbtypeBase.c_str(), dbtype, false);
        FileUtil::publishAtomically(dbtypeBase + ".dbtype", db + ".dbtype");
        Debug(Debug::INFO) << "Database holds " << ranks.sequenceCount << " sequences, "
                           << (ranks.bytes >> 30) << " GB\n";
    } else {
        Debug(Debug::INFO) << "Wrote the parts of node " << node.index << ", waiting for the other nodes\n";
    }
    return EXIT_SUCCESS;
}
