#include "IndexTable.h"

IndexTable::IndexTable (int alphabetSize, int kmerSize)
{
    this->alphabetSize = alphabetSize;
    this->kmerSize = kmerSize;

    tableSize = ipow(alphabetSize, kmerSize);

    table = new DynamicArray* [tableSize];
    for (int i = 0; i < tableSize; i++){
        table[i] = new DynamicArray();
    }

    idxer = new Indexer(alphabetSize, kmerSize);
}

void IndexTable::addSequence (Sequence* s){
    // iterate over all k-mers of the sequence and add the id of s to the sequence list of the k-mer (tableDummy)
    unsigned int kmerIdx;
    s->resetCurrPos();
    this->s = s;
    idxer->reset();

    while(s->hasNextKmer(kmerSize)){
        const int* kmer = s->nextKmer(kmerSize);
        kmerIdx = idxer->getNextKmerIndex(kmer, kmerSize);
        table[kmerIdx]->pushBack(s->id);
    }
}

void IndexTable::checkSizeAndCapacity(){

    long sizes = 0;
    long capacities = 0;
    for (int i = 0; i < tableSize; i++){
        sizes += table[i]->getSize();
        capacities += table[i]->getCapacity();
    }
    std::cout << sizes << " " << capacities;
}

void IndexTable::reduceMemoryUsage(){
    for (int i = 0; i < tableSize; i++){
        table[i]->shrinkToFit();
    }
}

void IndexTable::init(){
    for (int i = 0; i < tableSize; i++){
        // remove duplicate sequence entries
        table[i]->removeDuplicates();
        table[i]->shrinkToFit();
    }
}

int* IndexTable::getDBSeqList (int kmer, int* matchedListSize){
    *matchedListSize = table[kmer]->getSize();
    return table[kmer]->getEntries();
}

int IndexTable::ipow (int base, int exponent){
    int res = 1;
    for (int i = 0; i < exponent; i++)
        res = res*base;
    return res;
}
