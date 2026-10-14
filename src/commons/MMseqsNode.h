#ifndef MMSEQS_NODE_H
#define MMSEQS_NODE_H

#include <string>

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
};

#endif
