#include "DynamicArray.h"

DynamicArray::DynamicArray(){
    size = 0;
    capacity = 0;
}

DynamicArray::DynamicArray(int initialCapacity){
    size = 0;
    capacity = initialCapacity;
    entries = new int[initialCapacity];
}

DynamicArray::~DynamicArray(){
    delete[] entries;
}

void DynamicArray::pushBack(int entry){
    if (size == capacity)
        expand();
    entries[size++] = entry;
}

void DynamicArray::expand(){
    capacity = capacity*2 + 1;
    int* entriesNew = new int[capacity];
    memcpy(entriesNew, entries, sizeof(int)*size);
    if (size > 0)
        delete[] entries;
    entries = entriesNew;

}

void DynamicArray::removeDuplicates(){
    if (size == 0)
        return;
    std::sort(entries, entries+size);
    // remove duplicates in-place
    int boundary = 1;
    for (int i = 1; i < size; i++){
        if (entries[i] != entries[i-1])
            entries[boundary++] = entries[i];
    }
    size = boundary;
}

void DynamicArray::shrinkToFit(){
    if (size == capacity)
        return;
    int* entriesNew = new int[size];
    memcpy(entriesNew, entries, sizeof(int)*size);
    delete[] entries;
    entries = entriesNew;
    capacity = size;
}

int* DynamicArray::getEntries(){
    return entries;
}

int DynamicArray::getSize(){
    return size;
}

int DynamicArray::getCapacity(){
    return capacity;
}
