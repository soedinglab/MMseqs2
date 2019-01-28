#ifndef DBCONCAT_H
#define DBCONCAT_H

#include "DBReader.h"

class DBConcat : public DBReader<unsigned int> {

public:
    DBConcat(const std::string &dataFileNameA, const std::string &indexFileNameA,
             const std::string &dataFileNameB, const std::string &indexFileNameB,
             const std::string &dataFileNameC, const std::string &indexFileNameC,
             unsigned int threads, int dataMode = USE_DATA | USE_INDEX,
             bool preserveKeysA = false, bool preserveKeysB = false, bool takeLargerEntry = false);

    ~DBConcat();

    void concat(bool write = true);

    unsigned int dbAKeyMap(unsigned int);
    unsigned int dbBKeyMap(unsigned int);

private:
    std::string dataFileNameA;
    std::string indexFileNameA;
    std::string dataFileNameB;
    std::string indexFileNameB;
    std::string dataFileNameC;
    std::string indexFileNameC;

    size_t indexSizeA;
    size_t indexSizeB;

    std::pair<unsigned int, unsigned int> *keysA, *keysB;

    unsigned int threads;

    bool sameDatabase;

    bool preserveKeysA; // do not change the keys of DBA
    bool preserveKeysB; // do not change the keys of DBA
    bool takeLargerEntry; // do not write empty entries

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
