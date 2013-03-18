#ifndef DYNAMIC_ARRAY_H
#define DYNAMIC_ARRAY_H

#include <cstring>
#include <cstdlib>
#include <algorithm>

class DynamicArray{

    public:
        DynamicArray();

        DynamicArray(int initialCapacity);

        ~DynamicArray();

        // append an entry
        void pushBack(int entry);

        // append an array of entries
        void addEntries(int* entry, int size);

        // expand the entries array
        void expand();

        // remove duplicate entries
        void removeDuplicates();

        // shrink the capacity of the array to its actual size
        void shrinkToFit();

        // get the pointer to the underlying data structure
        int* getEntries();

        // get the number of entries in the array
        int getSize();

        // get the maximum size before the array has to be enlarged
        int getCapacity();

        // reset the size: (but not the capacity!)
        void clear();

    private:

        // number of entries in the array
        int size;

        // maximum number of entries before the array has to be enlarged
        int capacity;

        int* entries;

};

#endif
