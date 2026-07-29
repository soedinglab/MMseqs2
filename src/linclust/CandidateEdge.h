#ifndef MMSEQS_CANDIDATEEDGE_H
#define MMSEQS_CANDIDATEEDGE_H

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

// One (representative, member) pair proposed by the distributed reduce.
//
// Packed binary rather than a prefilter DB. A `pref` DB stores each hit as ASCII
// with a per-representative index entry, and at 1e11 sequences the index alone is
// per-key state no single node can hold. 15 bytes per edge is also ~4x smaller
// than the text form, which matters directly: this is the largest intermediate
// the pipeline writes.
struct __attribute__((__packed__)) CandidateEdge {
    // 48-bit keys, as in KmerRecord: 2.8e14 keys, past the 1e12 target, and two
    // bytes cheaper per side than a full 64-bit key.
    uint8_t repBytes[6];
    uint8_t memberBytes[6];
    int16_t diagonal;
    // Nucleotide strand. Stock carries it in bit 63 of the representative key,
    // which a 48-bit key has no room for.
    uint8_t reverseStrand;
    // How many k-mers put this pair on this diagonal, saturating at 255. Stock's
    // prefilter score, used downstream to rank and filter candidates.
    uint8_t score;

    uint64_t getRep() const { return get(repBytes); }
    uint64_t getMember() const { return get(memberBytes); }
    void setRep(uint64_t key) { set(repBytes, key); }
    void setMember(uint64_t key) { set(memberBytes, key); }

private:
    static uint64_t get(const uint8_t *bytes) {
        uint64_t value = 0;
        for (int i = 5; i >= 0; i--) {
            value = (value << 8) | bytes[i];
        }
        return value;
    }
    static void set(uint8_t *bytes, uint64_t value) {
        for (int i = 0; i < 6; i++) {
            bytes[i] = static_cast<uint8_t>(value & 0xFF);
            value >>= 8;
        }
    }
};

// Buffered writer for one partition's edge file.
//
// One file per partition, written whole, so a partition redone after a crash
// simply overwrites its file. That makes the reduce idempotent, unlike the map,
// whose per-worker shards can retain a dead worker's partial output.
class EdgeWriter {
public:
    EdgeWriter(const std::string &path, size_t bufferRecords = 1024 * 1024);
    ~EdgeWriter();

    void append(const CandidateEdge &edge);
    void close();

    uint64_t getEdgeCount() const { return edgeCount; }

    static std::string partitionPath(const std::string &dir, unsigned int partition);

private:
    EdgeWriter(const EdgeWriter &);
    EdgeWriter &operator=(const EdgeWriter &);

    void flush();

    std::string path;
    FILE *file;
    std::vector<CandidateEdge> buffer;
    size_t bufferRecords;
    uint64_t edgeCount;
    bool closed;
};

#endif
