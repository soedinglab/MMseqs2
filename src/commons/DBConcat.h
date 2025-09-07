
#ifndef DBCONCAT_H
#define DBCONCAT_H
#include "MMseqsTypes.h"
#include <string>
#include <utility>

class DBConcat {
public:
    DBConcat(const std::string &dataFileNameA, const std::string &indexFileNameA,
             const std::string &dataFileNameB, const std::string &indexFileNameB,
             const std::string &dataFileNameC, const std::string &indexFileNameC,
             unsigned int threads, bool write = true, bool preserveKeysA = false, bool preserveKeysB = false, bool takeLargerEntry = false, size_t trimRight = 0);

    ~DBConcat();

    KeyType dbAKeyMap(KeyType);
    KeyType dbBKeyMap(KeyType);

private:
    size_t indexSizeA;
    size_t indexSizeB;

    std::pair<KeyType, KeyType> *keysA, *keysB;

    bool sameDatabase;

    struct compareFirstEntry {
        bool operator()(const std::pair<KeyType, KeyType> &lhs,
                        const std::pair<KeyType, KeyType> &rhs) const {
            return (lhs.first < rhs.first);
        }
    };

    struct compareKeyToFirstEntry {
        bool operator()(const KeyType &lhs, const std::pair<KeyType, KeyType> &rhs) const {
            return (lhs <= rhs.first);
        }
    };
};

#endif
