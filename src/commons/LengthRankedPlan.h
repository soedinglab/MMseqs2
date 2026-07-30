#ifndef MMSEQS_LENGTHRANKEDPLAN_H
#define MMSEQS_LENGTHRANKEDPLAN_H

#include <cstdint>
#include <string>
#include <vector>

// Placement plan for a length-ranked, densely-keyed sequence database that many
// independent workers build into one set of output files without ever merging.
//
// The distributed createdb runs two passes over the input FASTA. Pass 1 has each
// worker scan a byte range (a "chunk") and report only a histogram: for every
// sequence length in that chunk, how many sequences have it and how many bytes
// their headers occupy. Pass 2 has each worker rescan its chunk and write the
// sequences straight into their final position in the output files.
//
// What makes pass 2 possible with no merge step is that the histograms alone
// determine every byte offset in the finished database:
//
//   - keys are assigned longest sequence first, so all sequences of one length
//     occupy a contiguous key range whose start is simply the number of longer
//     sequences in the whole input;
//   - a sequence of length L occupies exactly L + 2 bytes of the data file (the
//     residues, a newline, and the NUL terminator DBWriter appends), so a length
//     bucket's data offset is a weighted sum over the longer buckets;
//   - header lengths vary per sequence, but pass 1 already totalled them per
//     (chunk, length) -- exactly the granularity at which offsets are needed.
//
// Ties are broken by chunk index and then by position within the chunk, which is
// input order. The key assignment therefore depends only on the input and the
// chunk size, never on how many workers ran or in what order they finished. That
// is what makes the build both restartable and topology-invariant.
//
// The whole scheme rests on keys being dense and length-ordered, which is also
// what lets DenseIndex address an entry by key and what makes the distributed
// greedy in the reduce stage exact.

// What one worker reports after scanning one chunk in pass 1.
class ChunkHistogram {
public:
    // One sequence length present in the chunk. Sorted ascending by length.
    struct Bucket {
        uint64_t length;
        uint64_t count;
        // Total bytes the headers of those sequences occupy in the header
        // database, i.e. the header text including its newline plus the NUL
        // that DBWriter appends. Summed here because header size cannot be
        // derived from sequence length the way data size can.
        uint64_t headerBytes;
        // Total bytes of the accessions of those sequences, as
        // Util::parseFastaHeader extracts them. Only this part of a `.lookup`
        // line is unknowable without reading the input; the key and the file
        // index are both decided by the plan, so the planner completes the line
        // width itself. Recording finished line widths here instead would be
        // circular -- the keys do not exist yet when pass 1 runs.
        uint64_t accessionBytes;
    };

    uint64_t chunkIdx;
    // Index into the sorted input file list. Needed to reconstruct which source
    // file every key came from when the .lookup file is written.
    uint64_t fileIdx;
    uint64_t seqCount;
    // Nucleotide detection votes, aggregated across chunks by the planner so the
    // database type is decided from the whole input rather than one chunk.
    uint64_t nuclVotes;
    uint64_t sampleCount;
    std::vector<Bucket> buckets;

    ChunkHistogram() : chunkIdx(0), fileIdx(0), seqCount(0), nuclVotes(0), sampleCount(0) {}

    void write(const std::string &path) const;
    // Exits if the file is missing or malformed: a chunk histogram that cannot be
    // read means the plan would silently place sequences at wrong offsets.
    static ChunkHistogram read(const std::string &path);
};

// Where one worker must put every sequence of its chunk in pass 2.
class ChunkPlan {
public:
    struct Entry {
        uint64_t length;
        uint64_t count;
        // First key, and first data/header byte, of this chunk's run of
        // sequences of this length. The worker consumes them in input order,
        // advancing each cursor as it writes.
        uint64_t keyStart;
        uint64_t dataOffset;
        uint64_t hdrOffset;
        // Where this bucket's `.lookup` lines go, and how many bytes they must
        // occupy. The width is carried rather than re-derived so pass 2 can check
        // the text it built against what the planner reserved: an offset scheme
        // that is wrong by a byte would otherwise corrupt a neighbouring bucket
        // silently.
        uint64_t lookupOffset;
        uint64_t lookupBytes;
    };

    uint64_t chunkIdx;
    uint64_t fileIdx;
    std::vector<Entry> entries;

    ChunkPlan() : chunkIdx(0), fileIdx(0) {}

    void write(const std::string &path) const;
    static ChunkPlan read(const std::string &path);
};

// Totals for the finished database, produced together with the per-chunk plans.
struct LengthRankedTotals {
    uint64_t seqCount;
    uint64_t dataBytes;
    uint64_t headerBytes;
    uint64_t lookupBytes;
    uint64_t maxSeqLen;
    uint64_t nuclVotes;
    uint64_t sampleCount;

    LengthRankedTotals()
        : seqCount(0), dataBytes(0), headerBytes(0), lookupBytes(0), maxSeqLen(0), nuclVotes(0),
          sampleCount(0) {}
};

// Turns the per-chunk histograms into per-chunk placement plans.
//
// Input histograms may arrive in any order; they are ordered by chunkIdx here so
// the result does not depend on directory listing order. plans is resized to one
// entry per histogram, indexed the same way as the sorted chunk order.
LengthRankedTotals buildLengthRankedPlan(std::vector<ChunkHistogram> &histograms,
                                         std::vector<ChunkPlan> &plans);

#endif
