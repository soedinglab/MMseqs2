#include "StringBlock.h"

#include <algorithm>
#include <numeric>


void StringBlock::compact() {
    size_t* indices = new size_t[entryCount];
    std::iota(indices, indices + entryCount, 0);
    std::sort(indices, indices + entryCount, SortBlockByIndex(data, offsets));
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
    size_t* newOffsets = (size_t*)malloc(entryCount * sizeof(size_t));
    size_t offset = 0;
    for (size_t i = 0; i < entryCount; ++i) {
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
