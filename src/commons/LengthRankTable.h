#ifndef MMSEQS_LENGTHRANKTABLE_H
#define MMSEQS_LENGTHRANKTABLE_H

#include <cstdint>
#include <string>
#include <vector>

// Sequence length from a database key, in space proportional to the number of
// distinct lengths rather than to the database.
//
// This exists to delete a field from the k-mer record. Stock kmermatcher keeps
// `seqkey_to_len[dbKeySize]` (kmermatcher.cpp:1236), an array sized by *key
// space* -- 2-4 TB at 1e12 -- so the distributed map carried a 2-byte `seqLen` in
// every k-mer record instead. At 21 k-mers per sequence that is 42 B/seq, or
// 42 TB at 1e12, spent re-stating something the key already determines.
//
// It already determines it because `createdbparallel` assigns keys **longest
// sequence first**: key == global length rank. So all sequences of one length
// occupy a contiguous key range, and the entire key -> length map is one
// (firstKey, length) pair per distinct length, ordered.
//
// Sparse rather than an array indexed by length: the entry count is then bounded
// by the number of distinct lengths, not by the longest sequence. For proteins
// capped at 65535 residues that is at most 1 MB either way, but a dense array
// would scale with a single long contig and this does not.
//
// The pairs come free: `createdbparallel` pass 1 already histograms lengths in
// order to compute byte offsets, and `buildLengthRankedPlan` hands them over as it
// assigns the keys they describe -- so the table cannot drift from the key order
// it encodes.
class LengthRankTable {
public:
    // One distinct length and where its keys start, longest first.
    struct Run {
        uint64_t length;
        uint64_t firstKey;
        uint64_t count;
    };

    // On disk and in memory: ordered by ascending firstKey, which is the same as
    // descending length.
    struct Entry {
        uint64_t firstKey;
        uint32_t length;
        uint32_t reserved;
    };

    LengthRankTable() : maxLength(0), entryCount(0) {}

    // <db>.lenrank
    static std::string fileName(const std::string &dbPath);

    // Writes to a private sibling and renames, so an interrupted build never
    // leaves a short file that a later stage reads as a valid (and wrong) table.
    // Runs must be sorted by descending length and must tile key space.
    static void write(const std::string &dbPath, const std::vector<Run> &runs, uint64_t entryCount);

    static bool exists(const std::string &dbPath);

    // Exits if the file is missing, truncated, or fails its magic/version check.
    // A silently wrong length would mis-rank k-mer group centres -- a different
    // clustering with no error anywhere, which is the failure this pipeline is
    // least able to detect, so it is never tolerated.
    void open(const std::string &dbPath);

    bool isOpen() const { return entries.empty() == false; }

    uint64_t getEntryCount() const { return entryCount; }
    uint64_t getMaxLength() const { return maxLength; }
    size_t getRunCount() const { return entries.size(); }

    // Length of the sequence with this key. Exits if the key is outside the
    // database.
    unsigned int lengthOf(uint64_t key) const;

    // Non-fatal form, for validating the table against a database rather than
    // trusting it.
    bool tryLengthOf(uint64_t key, unsigned int &length) const;

    // The 16-bit `pos` and the historic 16-bit `seqLen` of a k-mer record cap
    // sequences here. Nothing enforced it before: a longer sequence wrapped
    // silently and its k-mers were extracted at wrong positions.
    static const uint64_t MAX_SEQUENCE_LENGTH = 65535;

    static const uint64_t MAGIC = 0x4B4E41524E454C4DULL;
    static const uint32_t VERSION = 1;

private:
    std::vector<Entry> entries;
    uint64_t maxLength;
    uint64_t entryCount;
};

#endif
