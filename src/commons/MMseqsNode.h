#ifndef MMSEQS_NODE_H
#define MMSEQS_NODE_H

#include "IndexTypes.h"

#include <cstdio>
#include <string>
#include <vector>

class DBWriter;

// a member of a cluster with its representative, the record of most exchanges between nodes
struct ClusterPair {
    DBKeyType repKey;
    DBKeyType memberKey;

    static bool compareByMember(const ClusterPair &first, const ClusterPair &second) {
        return first.memberKey < second.memberKey;
    }

    // the size of the pair at record, 0 if bytes do not hold all of it
    static size_t recordSize(const char *, size_t bytes) {
        return (bytes >= sizeof(ClusterPair)) ? sizeof(ClusterPair) : 0;
    }

    // the representative itself first, so it is the first line of its entry
    static bool compareByRep(const ClusterPair &first, const ClusterPair &second) {
        if (first.repKey != second.repKey) {
            return first.repKey < second.repKey;
        }
        if ((first.memberKey == first.repKey) != (second.memberKey == second.repKey)) {
            return first.memberKey == first.repKey;
        }
        return first.memberKey < second.memberKey;
    }
};

class MMseqsNode {
public:
    static const unsigned int MASTER = 0;

    static unsigned int rank;
    static unsigned int count;

    static void init();
    static inline bool isMaster() {
        return rank == MASTER;
    }

    static void publishFile(const std::string &file, const std::string &content);
    static std::string waitForFile(const std::string &file);
    // waits for prefix + node of every node and returns their content in node order; the files are removed
    static std::string gatherFiles(const std::string &prefix);
    // publishFile and gatherFiles for files of this run, so a file that a killed run left is not gathered
    static void publishRunFile(const std::string &file, const std::string &content);
    static std::string waitForRunFile(const std::string &file);
    static std::string gatherRunFiles(const std::string &prefix);
    // a log that one node appends to while another reads it: its synced byte count and a key are published as <log>.progress, read from a new open
    static void publishLog(FILE *file, const std::string &log, size_t key, size_t bytes);
    static bool readLogProgress(const std::string &log, size_t &key, size_t &bytes);
    static void waitForLog(const std::string &log, size_t &key, size_t &bytes);
    static void readLog(const std::string &log, size_t from, size_t to, char *data);
    // every node publishes content as file + "." + node; returns all contents in node order, gathered by the master into file
    static std::string barrier(const std::string &file, const std::string &content = "");
    // a barrier that returns the smallest memoryLimit of all nodes, so they all cut the same splits; with threads all nodes must run with as many
    static size_t getSmallestMemoryLimit(const std::string &file, size_t memoryLimit, int threads = 0);

    // exchange through files: a node appends flushes sorted by destination to <prefix>.<node>, and the index line of a flush lists
    // this node and where the bytes of every node start in the file, then its end; a node reads its part of every flush as a block
    struct ExchangeBlock {
        FILE *file;
        size_t pos;
        size_t end;
        std::vector<char> buffer;
        size_t bufferPos;
        size_t bufferEnd;
        // end of the records that the buffer holds in full and start of the last of them
        size_t recordEnd;
        size_t lastRecord;
    };
    // nodeOffsets: where the bytes of every node start in the flush, then the flush size
    static void addExchangeFlush(std::string &indexData, size_t &writtenBytes, const std::vector<size_t> &nodeOffsets);
    // closes the file of this node and returns the index lines of all nodes once every node has closed its file
    static std::string closeExchange(FILE *file, const std::string &fileName, const std::string &indexFile, const std::string &indexData);
    static void openExchange(const std::string &index, const std::string &prefix, std::vector<ExchangeBlock> &blocks, std::vector<FILE *> &files);
    // moves the buffered records up to last of every block to records, recordOffsets their starts then the end; returns their count, 0 once none is left
    static size_t readExchangeRange(std::vector<ExchangeBlock> &blocks, size_t bufferMemory, size_t (*recordSize)(const char *, size_t),
                                    bool (*compare)(const ClusterPair &, const ClusterPair &), const ClusterPair &last,
                                    std::vector<char> &records, std::vector<size_t> *recordOffsets, const std::string &prefix);
    static void closeExchangeFiles(std::vector<FILE *> &files, const std::string &prefix);

    // sorts a flush of pairs, appends it to the exchange file of this node and adds its index line; nodeKeyOffsets: keys of every node
    static void writePairs(ClusterPair *pairs, size_t pairCount, bool (*compare)(const ClusterPair &, const ClusterPair &), const std::vector<size_t> &nodeKeyOffsets,
                           FILE *file, const std::string &fileName, std::string &indexData, size_t &writtenBytes);
    // the pairs that this node sends through the exchange file <exchange>.<node>, buffered and written by writePairs once pairCapacity are buffered
    struct ClusterPairWriter {
        ClusterPairWriter(const std::string &exchange, size_t pairCapacity, bool (*compare)(const ClusterPair &, const ClusterPair &), const std::vector<size_t> &nodeKeyOffsets);
        ~ClusterPairWriter();
        void add(DBKeyType repKey, DBKeyType memberKey) {
            if (pairCount == pairCapacity) {
                flush();
            }
            pairs[pairCount].repKey = repKey;
            pairs[pairCount].memberKey = memberKey;
            pairCount++;
        }
        void flush();
        // writes the buffered pairs and returns the index lines of all nodes once every node has closed its file
        std::string close();

        const std::string exchange;
        const std::string fileName;
        FILE *file;
        ClusterPair *pairs;
        const size_t pairCapacity;
        size_t pairCount;
        bool (*compare)(const ClusterPair &, const ClusterPair &);
        const std::vector<size_t> &nodeKeyOffsets;
        std::string indexData;
        size_t writtenBytes;
    };
    // the representative of entry p becomes the one received for p in the exchange receivedIndex, receivedPrefix if given
    static std::string sendClusterPairs(const std::string &clu, bool withRepresentative, bool (*compare)(const ClusterPair &, const ClusterPair &), const std::vector<size_t> &nodeKeyOffsets,
                                        const std::string &exchange, const std::string &receivedIndex = "", const std::string &receivedPrefix = "");
    // once every node has read lastExchange, a node removes its file of it; the master waits for all nodes and removes the indexes of exchanges
    static void removeExchanges(const std::string &prefix, const std::string &lastExchange, const std::vector<std::string> &exchanges);
    // ORs the bits of keys [from, to) of a file, bit k at fileOffset = key fileFirstKey + k, into bit key - firstKey of bits, a word longer than its keys
    static void readKeyBits(int fd, const std::string &fileName, size_t fileOffset, size_t fileFirstKey, size_t from, size_t to, size_t firstKey, uint64_t *bits);
    // the representative bits of keys [firstKey, firstKey + keyCount) from the files <db>.bits.<node>, each {firstKey, keyCount} then its bits
    static uint64_t *readRepresentativeBits(const std::string &db, const std::vector<size_t> &nodeKeyOffsets, size_t firstKey, size_t keyCount);
    // --createdb-mode 3 cluster databases: <clu>.index.<node> holds uint64 (key, offset) per entry in key order, then (DB_KEY_INVALID, data size)
    static void writeClusterIndexEntry(FILE *file, const std::string &fileName, DBKeyType key, size_t offset);
    // closes writer of <db>.<node> of a --createdb-mode 3 result and removes its text index; ends and closes indexFile, the cluster database index
    static void closeClusterWriter(DBWriter &writer, FILE *indexFile, const std::string &indexFileName);
    // writes the pairs of an exchange in rep order as the cluster database <clu>.<node>; singletonBits: bit k makes key firstKey + k a singleton
    static void writeClusterPairs(const std::string &index, const std::string &prefix, const std::string &clu, const uint64_t *singletonBits, size_t firstKey, size_t keyCount);
};

#endif
