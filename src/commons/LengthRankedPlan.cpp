#include "LengthRankedPlan.h"

#include "Debug.h"
#include "FileUtil.h"

#include <algorithm>
#include <cstdio>

namespace {

// Both files are written and read as a fixed header followed by a POD array, so
// a worker reads back exactly what the planner wrote with no parsing.
const uint64_t HIST_MAGIC = 0x4d4d4c5248495331ULL;  // "MMLRHIS1"
const uint64_t PLAN_MAGIC = 0x4d4d4c52504c4e31ULL;  // "MMLRPLN1"

struct FileHeader {
    uint64_t magic;
    uint64_t chunkIdx;
    uint64_t fileIdx;
    uint64_t entryCount;
    uint64_t seqCount;
    uint64_t nuclVotes;
    uint64_t sampleCount;
    uint64_t reserved;
};

void writeBlock(const std::string &path, const FileHeader &header,
                const void *entries, size_t entryBytes) {
    // Write to a temporary and rename, so a worker that dies mid-write leaves no
    // truncated file that a later reader would mistake for a complete one.
    std::string tmp = path + ".tmp";
    FILE *file = FileUtil::openAndDelete(tmp.c_str(), "wb");
    if (fwrite(&header, sizeof(FileHeader), 1, file) != 1) {
        Debug(Debug::ERROR) << "Cannot write header to " << tmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (entryBytes > 0 && fwrite(entries, entryBytes, 1, file) != 1) {
        Debug(Debug::ERROR) << "Cannot write entries to " << tmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fclose(file) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << tmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::move(tmp.c_str(), path.c_str());
}

FileHeader readHeader(FILE *file, const std::string &path, uint64_t magic) {
    FileHeader header;
    if (fread(&header, sizeof(FileHeader), 1, file) != 1) {
        Debug(Debug::ERROR) << "Cannot read header from " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (header.magic != magic) {
        Debug(Debug::ERROR) << "File " << path << " is not a length-ranked plan file\n";
        EXIT(EXIT_FAILURE);
    }
    return header;
}

}  // namespace

void ChunkHistogram::write(const std::string &path) const {
    FileHeader header;
    header.magic = HIST_MAGIC;
    header.chunkIdx = chunkIdx;
    header.fileIdx = fileIdx;
    header.entryCount = buckets.size();
    header.seqCount = seqCount;
    header.nuclVotes = nuclVotes;
    header.sampleCount = sampleCount;
    header.reserved = 0;
    writeBlock(path, header, buckets.data(), buckets.size() * sizeof(Bucket));
}

ChunkHistogram ChunkHistogram::read(const std::string &path) {
    FILE *file = FileUtil::openFileOrDie(path.c_str(), "rb", true);
    FileHeader header = readHeader(file, path, HIST_MAGIC);

    ChunkHistogram histogram;
    histogram.chunkIdx = header.chunkIdx;
    histogram.fileIdx = header.fileIdx;
    histogram.seqCount = header.seqCount;
    histogram.nuclVotes = header.nuclVotes;
    histogram.sampleCount = header.sampleCount;
    histogram.buckets.resize(header.entryCount);
    if (header.entryCount > 0 &&
        fread(histogram.buckets.data(), sizeof(Bucket), header.entryCount, file) != header.entryCount) {
        Debug(Debug::ERROR) << "Cannot read buckets from " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    fclose(file);
    return histogram;
}

void ChunkPlan::write(const std::string &path) const {
    FileHeader header;
    header.magic = PLAN_MAGIC;
    header.chunkIdx = chunkIdx;
    header.fileIdx = fileIdx;
    header.entryCount = entries.size();
    header.seqCount = 0;
    header.nuclVotes = 0;
    header.sampleCount = 0;
    header.reserved = 0;
    writeBlock(path, header, entries.data(), entries.size() * sizeof(Entry));
}

ChunkPlan ChunkPlan::read(const std::string &path) {
    FILE *file = FileUtil::openFileOrDie(path.c_str(), "rb", true);
    FileHeader header = readHeader(file, path, PLAN_MAGIC);

    ChunkPlan plan;
    plan.chunkIdx = header.chunkIdx;
    plan.fileIdx = header.fileIdx;
    plan.entries.resize(header.entryCount);
    if (header.entryCount > 0 &&
        fread(plan.entries.data(), sizeof(Entry), header.entryCount, file) != header.entryCount) {
        Debug(Debug::ERROR) << "Cannot read entries from " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    fclose(file);
    return plan;
}

LengthRankedTotals buildLengthRankedPlan(std::vector<ChunkHistogram> &histograms,
                                         std::vector<ChunkPlan> &plans) {
    // Order by chunk index so the plan never depends on the order the histograms
    // happened to be collected in.
    std::sort(histograms.begin(), histograms.end(),
              [](const ChunkHistogram &a, const ChunkHistogram &b) { return a.chunkIdx < b.chunkIdx; });

    plans.clear();
    plans.resize(histograms.size());
    for (size_t i = 0; i < histograms.size(); i++) {
        plans[i].chunkIdx = histograms[i].chunkIdx;
        plans[i].fileIdx = histograms[i].fileIdx;
    }

    // Every distinct length in the input, longest first: that order is the key
    // order of the finished database.
    std::vector<uint64_t> lengths;
    for (size_t i = 0; i < histograms.size(); i++) {
        for (size_t b = 0; b < histograms[i].buckets.size(); b++) {
            lengths.push_back(histograms[i].buckets[b].length);
        }
    }
    std::sort(lengths.begin(), lengths.end(), std::greater<uint64_t>());
    lengths.erase(std::unique(lengths.begin(), lengths.end()), lengths.end());

    // Cursor into each chunk's bucket list. The buckets are sorted ascending by
    // length and we walk lengths descending, so each cursor moves backwards and
    // the whole sweep stays linear in the number of buckets.
    std::vector<size_t> cursor(histograms.size());
    for (size_t i = 0; i < histograms.size(); i++) {
        cursor[i] = histograms[i].buckets.size();
    }

    LengthRankedTotals totals;
    uint64_t nextKey = 0;
    uint64_t nextDataOffset = 0;
    uint64_t nextHdrOffset = 0;

    for (size_t l = 0; l < lengths.size(); l++) {
        const uint64_t length = lengths[l];
        for (size_t i = 0; i < histograms.size(); i++) {
            if (cursor[i] == 0 || histograms[i].buckets[cursor[i] - 1].length != length) {
                continue;
            }
            const ChunkHistogram::Bucket &bucket = histograms[i].buckets[cursor[i] - 1];
            cursor[i]--;

            ChunkPlan::Entry entry;
            entry.length = length;
            entry.count = bucket.count;
            entry.keyStart = nextKey;
            entry.dataOffset = nextDataOffset;
            entry.hdrOffset = nextHdrOffset;
            plans[i].entries.push_back(entry);

            nextKey += bucket.count;
            // A sequence of length L always occupies L + 2 data bytes, so this
            // stays exact without looking at the sequences themselves.
            nextDataOffset += bucket.count * (length + 2);
            nextHdrOffset += bucket.headerBytes;
        }
    }

    for (size_t i = 0; i < histograms.size(); i++) {
        if (cursor[i] != 0) {
            Debug(Debug::ERROR) << "Chunk " << histograms[i].chunkIdx
                                << " has buckets that are not sorted by length\n";
            EXIT(EXIT_FAILURE);
        }
        // Written ascending by length; pass 2 looks entries up by length, and a
        // sorted list keeps that a binary search.
        std::reverse(plans[i].entries.begin(), plans[i].entries.end());
        totals.nuclVotes += histograms[i].nuclVotes;
        totals.sampleCount += histograms[i].sampleCount;
    }

    totals.seqCount = nextKey;
    totals.dataBytes = nextDataOffset;
    totals.headerBytes = nextHdrOffset;
    totals.maxSeqLen = lengths.empty() ? 0 : lengths[0];
    return totals;
}
