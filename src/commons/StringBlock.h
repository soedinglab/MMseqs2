#ifndef STRING_BLOCK_H
#define STRING_BLOCK_H

#include <cstddef>
#include <cstdlib>
#include <cstring>

class StringBlock {
public:
    StringBlock(size_t byteCapacity_ = 32 * 1024, size_t entryCapacity_ = 1024) {
        byteCapacity = byteCapacity_;
        data = (char*)malloc(byteCapacity * sizeof(char));

        entryCapacity = entryCapacity_;
        offsets = (size_t*)malloc(entryCapacity * sizeof(size_t));

        offsets[0] = 0;
        entryCount = 0;
        externalData = false;
    }

    ~StringBlock() {
        if (externalData == false) {
            free(data);
            free(offsets);
        }
    }

    const char* getString(size_t idx) const {
        if (idx >= entryCount) {
            return NULL;
        }
        return data + offsets[idx];
    }

    size_t append(const char* string, size_t length) {
        size_t nextOffset = offsets[entryCount];
        if (nextOffset + length >= byteCapacity) {
            byteCapacity = (nextOffset + length + 1) * 1.5;
            data = (char*)realloc(data, byteCapacity * sizeof(char));
        }
        memcpy(data + nextOffset, string, length);
        data[nextOffset + length] = '\0';
        entryCount++;

        if (entryCount >= entryCapacity) {
            entryCapacity = (entryCount + 1) * 1.5;
            offsets = (size_t*)realloc(offsets, entryCapacity * sizeof(size_t));
        }
        offsets[entryCount] = nextOffset + length + 1;
        return entryCount - 1;
    }

    void compact();

    static size_t memorySize(const StringBlock& block) {
        return 3 * sizeof(size_t) + block.byteCapacity * sizeof(char) + block.entryCapacity * sizeof(size_t);
    }

    static char* serialize(const StringBlock &block) {
        char* mem = (char*) malloc(memorySize(block));
        char* p = mem;
        memcpy(p, &block.byteCapacity, sizeof(size_t));
        p += sizeof(size_t);
        memcpy(p, &block.entryCapacity, sizeof(size_t));
        p += sizeof(size_t);
        memcpy(p, &block.entryCount, sizeof(size_t));
        p += sizeof(size_t);
        memcpy(p, block.data, block.byteCapacity * sizeof(char));
        p += block.byteCapacity * sizeof(char);
        memcpy(p, block.offsets, block.entryCapacity * sizeof(size_t));
        p += block.entryCapacity * sizeof(size_t);
        return mem;
    }

    static StringBlock* unserialize(const char* mem) {
        const char* p = mem;
        size_t byteCapacity = *((size_t*)p);
        p += sizeof(size_t);
        size_t entryCapacity = *((size_t*)p);
        p += sizeof(size_t);
        size_t entryCount = *((size_t*)p);
        p += sizeof(size_t);
        char* data = (char*)p;
        p += byteCapacity * sizeof(char);
        size_t* offsets = (size_t*)p;
        p += entryCapacity * sizeof(size_t);
        return new StringBlock(byteCapacity, entryCapacity, entryCount, data, offsets);
    }

private:
    StringBlock(size_t byteCapacity, size_t entryCapacity, size_t entryCount, char* data, size_t* offsets)
            : byteCapacity(byteCapacity), entryCapacity(entryCapacity), entryCount(entryCount), data(data), offsets(offsets), externalData(true) {};
    size_t byteCapacity;
    size_t entryCapacity;
    size_t entryCount;

    char* data;
    size_t* offsets;
    bool externalData;

    struct SortBlockByIndex {
        SortBlockByIndex(char* data, size_t* offsets) : data(data), offsets(offsets) {}
        bool operator() (size_t i, size_t j) const {
            return strcmp(data + offsets[i], data + offsets[j]) < 0;
        }
        char* data;
        size_t* offsets;
    };
};

#endif
