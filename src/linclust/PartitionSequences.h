#ifndef MMSEQS_PARTITIONSEQUENCES_H
#define MMSEQS_PARTITIONSEQUENCES_H

#include <cstdint>
#include <string>
#include <vector>

// The sequences one k-mer partition needs, fetched by key into a local arena.
//
// This is the piece that lets the alignment run inside the reduce without a
// resident index. Stock align2clust opens the whole sequence DB with USE_INDEX
// (Align2clust.cpp:400), which is 24 B per sequence -- 2.4 TB at 1e11 and 24 TB
// at 1e12, on a node that has 2 TB. Here only the keys a partition actually
// touches are read, addressed directly through the dense companion index.
//
// Reads are issued in ascending key order, which matters more than it looks:
// dense keys mean ascending key == ascending file offset, so a scattered set of
// sequences is fetched as a forward scan rather than as random seeks. Adjacent
// requests are coalesced into single reads for the same reason.
//
// Sizing is the open question this class exists to answer. A partition touches
// roughly `kmersPerSequence / P` of the database, so the arena is small at test
// scale and large at target scale; getBytes() is reported so the extrapolation
// is measured rather than guessed.
class PartitionSequences {
public:
    explicit PartitionSequences(const std::string &dbName);
    ~PartitionSequences();

    // Fetches every key in `sortedKeys` (ascending, unique) into the arena,
    // replacing whatever was there before.
    void load(const std::vector<uint64_t> &sortedKeys);

    // Residues for a key, or NULL if it was not loaded. The returned pointer is
    // valid until the next load().
    const char *get(uint64_t key, unsigned int *length) const;

    uint64_t getBytes() const { return arena.size(); }
    uint64_t getCount() const { return keys.size(); }
    // Bytes actually read from disk, including coalesced gaps.
    uint64_t getBytesRead() const { return bytesRead; }

private:
    PartitionSequences(const PartitionSequences &);
    PartitionSequences &operator=(const PartitionSequences &);

    void readAt(int fd, void *dst, size_t length, size_t offset, const char *what);

    std::string dbName;
    int dataFd;
    int indexFd;
    uint64_t entryCount;
    uint64_t firstKey;

    std::vector<uint64_t> keys;     // ascending, parallel to offsets/lengths
    std::vector<uint64_t> offsets;  // into arena
    std::vector<uint32_t> lengths;  // residues, without the newline/terminator
    std::vector<char> arena;
    uint64_t bytesRead;
};

#endif
