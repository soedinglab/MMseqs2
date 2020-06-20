//
// Created by Martin Steinegger on 6/21/20.
//

#ifndef MMSEQS_MEMORYTRACKER_H
#define MMSEQS_MEMORYTRACKER_H
#include "stddef.h"
class MemoryTracker{
public:
    static size_t getSize() { return totalMemorySizeInst;};
protected:
    static size_t totalMemorySizeInst;
    static void incrementMemory(size_t memorySize) { totalMemorySizeInst+=memorySize; }
    static void decrementMemory(size_t memorySize) { totalMemorySizeInst-=memorySize; }
};
#endif //MMSEQS_MEMORYTRACKER_H
