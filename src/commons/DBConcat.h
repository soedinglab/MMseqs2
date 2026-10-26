#ifndef DBCONCAT_H
#define DBCONCAT_H

#include "IndexTypes.h"

#include <string>
#include <utility>

class DBConcat {
public:
    DBConcat(const std::string &dataFileNameA, const std::string &indexFileNameA,
             const std::string &dataFileNameB, const std::string &indexFileNameB,
             const std::string &dataFileNameC, const std::string &indexFileNameC,
             unsigned int threads, bool write = true, bool preserveKeysA = false, bool preserveKeysB = false, bool takeLargerEntry = false, size_t trimRight = 0);

    ~DBConcat();

    DBKeyType dbAKeyMap(DBKeyType);
    DBKeyType dbBKeyMap(DBKeyType);

private:
    size_t indexSizeA;
    size_t indexSizeB;

    std::pair<DBKeyType, DBKeyType> *keysA, *keysB;

    bool sameDatabase;

    struct compareFirstEntry {
        bool operator()(const std::pair<DBKeyType, DBKeyType> &lhs,
                        const std::pair<DBKeyType, DBKeyType> &rhs) const {
            return (lhs.first < rhs.first);
        }
    };

    struct compareKeyToFirstEntry {
        bool operator()(const DBKeyType &lhs, const std::pair<DBKeyType, DBKeyType> &rhs) const {
            return (lhs <= rhs.first);
        }
    };
};

#endif
