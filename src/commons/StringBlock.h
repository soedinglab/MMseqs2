#ifndef STRING_BLOCK_H
#define STRING_BLOCK_H

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <numeric>

#include "FastSort.h"

template<typename T>
class StringBlock {
public:
    StringBlock(size_t byteCapacity_ = 32 * 1024, T entryCapacity_ = 1024) {
        byteCapacity = byteCapacity_;
        data = (char*)malloc(byteCapacity * sizeof(char));

        entryCapacity = entryCapacity_;
        offsets = (T*)malloc(entryCapacity * sizeof(T));

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

    const char* getString(T idx) const {
        if (idx >= entryCount) {
            return NULL;
        }
        return data + offsets[idx];
    }

    size_t append(const char* string, size_t length) {
        T nextOffset = offsets[entryCount];
        if (nextOffset + length >= byteCapacity) {
            byteCapacity = (nextOffset + length + 1) * 1.5;
            data = (char*)realloc(data, byteCapacity * sizeof(char));
        }
        memcpy(data + nextOffset, string, length);
        data[nextOffset + length] = '\0';
        entryCount++;

        if (entryCount >= entryCapacity) {
            entryCapacity = (entryCount + 1) * 1.5;
            offsets = (T*)realloc(offsets, entryCapacity * sizeof(T));
        }
        offsets[entryCount] = nextOffset + length + 1;
        return entryCount - 1;
    }

    void compact() {
        T* indices = new T[entryCount];
        std::iota(indices, indices + entryCount, 0);
        SORT_SERIAL(indices, indices + entryCount, SortBlockByIndex(data, offsets));
        size_t unique = 1;
        size_t totalLength = strlen(getString(indices[0]));
        for (size_t i = 1; i < entryCount; ++i) {
            if (strcmp(getString(indices[i]), getString(indices[i - 1])) == 0) {
                offsets[indices[i]] = offsets[indices[i - 1]];
            } else {
                unique++;
                totalLength += strlen(getString(indices[i]));
            }
        }
        char* newData = (char*)malloc((totalLength + unique) * sizeof(char));
        T* newOffsets = (T*)malloc(entryCount * sizeof(T));
        T offset = 0;
        for (T i = 0; i < entryCount; ++i) {
            if (i != 0 && strcmp(getString(indices[i]), getString(indices[i - 1])) == 0) {
                newOffsets[indices[i]] = newOffsets[indices[i - 1]];
            } else {
                newOffsets[indices[i]] = offset;
                size_t length = strlen(getString(indices[i]));
                memcpy(newData + offset, getString(indices[i]), length);
                newData[offset + length] = '\0';
                offset += length + 1;
            }
        }
        free(data);
        data = newData;
        free(offsets);
        offsets = newOffsets;
        entryCapacity = entryCount;
        byteCapacity = (totalLength + unique) * sizeof(char);
        delete[] indices;
    }

    static size_t memorySize(const StringBlock& block) {
        return  sizeof(size_t) + 2 * sizeof(T) + block.byteCapacity * sizeof(char) + block.entryCapacity * sizeof(T);
    }

    static char* serialize(const StringBlock &block) {
        char* mem = (char*) malloc(memorySize(block));
        char* p = mem;
        memcpy(p, &block.byteCapacity, sizeof(size_t));
        p += sizeof(size_t);
        memcpy(p, &block.entryCapacity, sizeof(T));
        p += sizeof(T);
        memcpy(p, &block.entryCount, sizeof(T));
        p += sizeof(T);
        memcpy(p, block.data, block.byteCapacity * sizeof(char));
        p += block.byteCapacity * sizeof(char);
        memcpy(p, block.offsets, block.entryCapacity * sizeof(T));
        p += block.entryCapacity * sizeof(T);
        return mem;
    }

    static StringBlock* unserialize(const char* mem) {
        const char* p = mem;
        size_t byteCapacity = *((size_t*)p);
        p += sizeof(size_t);
        size_t entryCapacity = *((T*)p);
        p += sizeof(T);
        size_t entryCount = *((T*)p);
        p += sizeof(T);
        char* data = (char*)p;
        p += byteCapacity * sizeof(char);
        T* offsets = (T*)p;
        p += entryCapacity * sizeof(T);
        return new StringBlock<T>(byteCapacity, entryCapacity, entryCount, data, offsets);
    }

private:
    StringBlock(size_t byteCapacity, T entryCapacity, T entryCount, char* data, T* offsets)
            : byteCapacity(byteCapacity), entryCapacity(entryCapacity), entryCount(entryCount), data(data), offsets(offsets), externalData(true) {};
    size_t byteCapacity;
    T entryCapacity;
    T entryCount;

    char* data;
    T* offsets;
    bool externalData;

    struct SortBlockByIndex {
        SortBlockByIndex(char* data, T* offsets) : data(data), offsets(offsets) {}
        bool operator() (T i, T j) const {
            return strcmp(data + offsets[i], data + offsets[j]) < 0;
        }
        char* data;
        T* offsets;
    };
};

#endif
