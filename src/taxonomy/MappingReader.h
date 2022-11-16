#ifndef MAPPING_READER_H
#define MAPPING_READER_H

#include "Util.h"
#include "MemoryMapped.h"
#include <algorithm>

class MappingReader {
public:
    static std::pair<char *, size_t> serialize(const MappingReader &reader) {
        size_t serialized_size = reader.magicLen + reader.count * sizeof(Pair);
        char* data = (char*)malloc(serialized_size);
        memcpy(data, reader.magic, reader.magicLen);
        memcpy(data + reader.magicLen, reader.entries, reader.count * sizeof(Pair));
        return std::make_pair(data, serialized_size);
    }

    MappingReader(const std::string &db, const bool dbInput = true) {
        std::string input = dbInput ? db + "_mapping" : db;
        file = new MemoryMapped(input, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
        if (!file->isValid()) {
            delete file;
            file = NULL;
            Debug(Debug::ERROR) << db << "_mapping does not exist. Please create the taxonomy mapping!\n";
            EXIT(EXIT_FAILURE);
        }
        char *data = (char *) file->getData();
        size_t dataSize = file->size();
        if (file->size() > magicLen && memcmp(data, magic, magicLen) == 0) {
            entries = reinterpret_cast<Pair*>(data + magicLen);
            count = (dataSize - magicLen) / sizeof(Pair);
            return;
        }
        std::vector<std::pair<unsigned int, unsigned int>> mapping;
        size_t currPos = 0;
        const char *cols[3];
        size_t isSorted = true;
        unsigned int prevId = 0;
        while (currPos < dataSize) {
            Util::getWordsOfLine(data, cols, 2);
            unsigned int id = Util::fast_atoi<size_t>(cols[0]);
            isSorted *= (id >= prevId);
            unsigned int taxid = Util::fast_atoi<size_t>(cols[1]);
            data = Util::skipLine(data);
            mapping.push_back(std::make_pair(id, taxid));
            currPos = data - (char *) file->getData();
            prevId = id;
        }
        file->close();
        delete file;
        file = NULL;
        if (mapping.size() == 0) {
            Debug(Debug::ERROR) << db << "_mapping is empty. Rerun createtaxdb to recreate taxonomy mapping.\n";
            EXIT(EXIT_FAILURE);
        }

        count = mapping.size();
        entries = new Pair[count];
        for (size_t i = 0; i < count; ++i) {
            entries[i].dbkey = mapping[i].first;
            entries[i].taxon = mapping[i].second;
        }
        if (isSorted == false) {
            std::stable_sort(entries, entries + count, compareTaxa);
        }
    }

    ~MappingReader() {
        if (file != NULL) {
            file->close();
            delete file;
        } else {
            delete[] entries;
        }
    }

    unsigned int lookup(unsigned int key) {
        unsigned int taxon = 0;
        // match dbKey to its taxon based on mapping
        Pair val;
        val.dbkey = key;
        Pair* end = entries + count;
        Pair* found = std::upper_bound(entries, end, val, compareTaxa);
        if (found == end || found->dbkey != key) {
            taxon = 0;
        } else {
            taxon = found->taxon;
        }
        return taxon;
    }

private:
    MemoryMapped* file;
    struct __attribute__((__packed__)) Pair{
        unsigned int dbkey;
        unsigned int taxon;
    };
    Pair* entries;
    size_t count;
    //                    T  A   X   M  Version
    const char magic[5] = {19, 0, 23, 12, 0};
    const size_t magicLen = 5;
    static bool compareTaxa(const Pair &lhs, const Pair &rhs) {
        return (lhs.dbkey <= rhs.dbkey);
    }
};

#endif
