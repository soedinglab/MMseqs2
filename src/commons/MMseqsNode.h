#ifndef MMSEQS_NODE_H
#define MMSEQS_NODE_H

#include "IndexTypes.h"

#include <cstdio>
#include <queue>
#include <string>
#include <vector>

// a member of a cluster with its representative, the record of most exchanges between nodes
struct ClusterPair {
    DBKeyType repKey;
    DBKeyType memberKey;

    static bool compareByMember(const ClusterPair &first, const ClusterPair &second) {
        return first.memberKey < second.memberKey;
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

struct QueuedPair {
    ClusterPair pair;
    size_t block;
};

struct CompareQueuedPair {
    bool (*compare)(const ClusterPair &, const ClusterPair &);
    bool operator()(const QueuedPair &first, const QueuedPair &second) const {
        return compare(second.pair, first.pair);
    }
};
typedef std::priority_queue<QueuedPair, std::vector<QueuedPair>, CompareQueuedPair> PairQueue;

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
    // every node publishes content as file + "." + node; returns all contents in node order, gathered by the master into file
    static std::string barrier(const std::string &file, const std::string &content = "");

    // exchange through files: a node appends flushes sorted by destination to <prefix>.<node>, and the index line of a flush lists
    // this node and where the bytes of every node start in the file, then its end; a node reads its part of every flush as a block
    struct ExchangeBlock {
        FILE *file;
        size_t pos;
        size_t end;
        std::vector<char> buffer;
        size_t bufferPos;
        size_t bufferEnd;
    };
    // nodeOffsets: where the bytes of every node start in the flush, then the flush size
    static void addExchangeFlush(std::string &indexData, size_t &writtenBytes, const std::vector<size_t> &nodeOffsets);
    // closes the file of this node and returns the index lines of all nodes once every node has closed its file
    static std::string closeExchange(FILE *file, const std::string &fileName, const std::string &indexFile, const std::string &indexData);
    static void openExchange(const std::string &index, const std::string &prefix, std::vector<ExchangeBlock> &blocks, std::vector<FILE *> &files);
    // the next bytes of a block, valid until the next read of the block; NULL at its end
    static const char *readExchange(ExchangeBlock &block, size_t bytes, const std::string &prefix);
    static void closeExchangeFiles(std::vector<FILE *> &files, const std::string &prefix);

    // sorts a flush of pairs, appends it to the exchange file of this node and adds its index line; nodeKeyOffsets: keys of every node
    static void writePairs(ClusterPair *pairs, size_t pairCount, bool (*compare)(const ClusterPair &, const ClusterPair &), const std::vector<size_t> &nodeKeyOffsets,
                           FILE *file, const std::string &fileName, std::string &indexData, size_t &writtenBytes);
    // queues the first pair of every block that the nodes sent to this node
    static void openPairs(PairQueue &queue, std::vector<ExchangeBlock> &blocks, std::vector<FILE *> &files, const std::string &index, const std::string &prefix);
    static void queueNextPair(PairQueue &queue, std::vector<ExchangeBlock> &blocks, size_t block, const std::string &prefix);
    // the representative of entry p becomes the one received for p in the exchange receivedIndex, receivedPrefix if given
    static std::string sendClusterPairs(const std::string &clu, bool withRepresentative, bool (*compare)(const ClusterPair &, const ClusterPair &), const std::vector<size_t> &nodeKeyOffsets,
                                        const std::string &exchange, const std::string &receivedIndex = "", const std::string &receivedPrefix = "");
    // once every node has read lastExchange, a node removes its file of it; the master waits for all nodes and removes the indexes of exchanges
    static void removeExchanges(const std::string &prefix, const std::string &lastExchange, const std::vector<std::string> &exchanges);
    // --createdb-mode 3 cluster databases: <clu>.index.<node> holds uint64 (key, offset) per entry in key order, then (DB_KEY_INVALID, data size)
    static void writeClusterIndexEntry(FILE *file, const std::string &fileName, DBKeyType key, size_t offset);
};

#endif
