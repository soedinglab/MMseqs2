#ifndef DBCONCAT_H
#define DBCONCAT_H

#include "DBReader.h"

class DBConcat : public DBReader<unsigned int> {

public:
    DBConcat(const std::string &dataFileNameA, const std::string &indexFileNameA,
             const std::string &dataFileNameB, const std::string &indexFileNameB,
             const std::string &dataFileNameC, const std::string &indexFileNameC,
             unsigned int threads, int dataMode = USE_DATA | USE_INDEX | USE_LOOKUP, bool write = true,
             bool preserveKeysA = false, bool preserveKeysB = false, bool takeLargerEntry = false);

    ~DBConcat();

    unsigned int dbAKeyMap(unsigned int);
    unsigned int dbBKeyMap(unsigned int);

private:
    size_t indexSizeA;
    size_t indexSizeB;

    std::pair<unsigned int, unsigned int> *keysA, *keysB;

    bool sameDatabase;

    struct compareFirstEntry {
        bool operator()(const std::pair<unsigned int, unsigned int> &lhs,
                        const std::pair<unsigned int, unsigned int> &rhs) const {
            return (lhs.first < rhs.first);
        }
    };

    struct compareKeyToFirstEntry {
        bool operator()(const unsigned int &lhs, const std::pair<unsigned int, unsigned int> &rhs) const {
            return (lhs <= rhs.first);
        }
    };

};


#endif
