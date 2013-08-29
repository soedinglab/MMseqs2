#include "IndexTable.h"

IndexTable::IndexTable (int alphabetSize, int kmerSize)
{
    this->alphabetSize = alphabetSize;
    this->kmerSize = kmerSize;

    tableSize = ipow(alphabetSize, kmerSize);

    sizes = new int[tableSize];
    currPos = new int[tableSize];
    for (int i = 0; i < tableSize; i++){
        sizes[i] = 0;
        currPos[i] = 0;
    }

    table = new int*[tableSize];

    idxer = new Indexer(alphabetSize, kmerSize);
}

IndexTable::~IndexTable(){
    for (int i = 0; i < tableSize; i++){
        if (sizes[i] > 0)
            delete[] table[i];
    }
    delete[] table;
    delete[] sizes;
    delete idxer;
}

void IndexTable::addKmerCount (Sequence* s){
    unsigned int kmerIdx;
    s->resetCurrPos();
    idxer->reset();

    while(s->hasNextKmer(kmerSize)){
        kmerIdx = idxer->getNextKmerIndex(s->nextKmer(kmerSize), kmerSize);
        sizes[kmerIdx]++;
    }

}

void IndexTable::init(){
    for (int i = 0; i < tableSize; i++){
        if (sizes[i] > 0)
            table[i] = new int[sizes[i]];
    }
}

void IndexTable::addSequence (Sequence* s){
    // iterate over all k-mers of the sequence and add the id of s to the sequence list of the k-mer (tableDummy)
    unsigned int kmerIdx;
    s->resetCurrPos();
    this->s = s;
    idxer->reset();

    while(s->hasNextKmer(kmerSize)){
        kmerIdx = idxer->getNextKmerIndex(s->nextKmer(kmerSize), kmerSize);
        table[kmerIdx][currPos[kmerIdx]++] = s->id;
    }
}

void IndexTable::removeDuplicateEntries(){

    delete[] currPos;
    for (int e = 0; e < tableSize; e++){
        if (sizes[e] == 0)
            continue;
        int* entries = table[e];
        int size = sizes[e];

        std::sort(entries, entries+size);
        // remove duplicates in-place
        int boundary = 1;
        for (int i = 1; i < size; i++){
            if (entries[i] != entries[i-1])
                entries[boundary++] = entries[i];
        }

        size = boundary;
        // no duplicates found
        if (size == sizes[e])
            continue;
        // copy the remaining entries to a smaller array
        int* entriesNew = new int[size];
        memcpy(entriesNew, entries, sizeof(int)*size);
        delete[] entries;
        table[e] = entriesNew;
        sizes[e] = size;
    }

}

void IndexTable::print(){
    int* testKmer = new int[kmerSize];
    for (int i = 0; i < tableSize; i++){
        if (sizes[i] > 0){
            idxer->printKmer(testKmer, i, kmerSize, s->int2aa);
            std::cout << "\n";
            for (int j = 0; j < sizes[i]; j++){
                std::cout << "\t" << table[i][j] << "\n";
            }
        }
    }
    delete[] testKmer;
}

int* IndexTable::getDBSeqList (int kmer, int* matchedListSize){
    *matchedListSize = sizes[kmer];
    return table[kmer];
}

int IndexTable::ipow (int base, int exponent){
    int res = 1;
    for (int i = 0; i < exponent; i++)
        res = res*base;
    return res;
}
