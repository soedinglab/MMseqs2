#include "IndexTableGlobal.h"

IndexTableGlobal::IndexTableGlobal (int alphabetSize, int kmerSize, int skip)
                    : IndexTable(alphabetSize, kmerSize, skip) { }

IndexTableGlobal::~IndexTableGlobal(){ }

void IndexTableGlobal::initMemory(){
    // allocate memory for the sequence id lists
    // tablesSizes is added to put the Size of the entry infront fo the memory
    entries = new unsigned int[tableEntriesNum+tableSize];
}

void IndexTableGlobal::init(){
    unsigned int* it = entries;
    // set the pointers in the index table to the start of the list for a certain k-mer
    for (size_t i = 0; i < tableSize; i++){
        table[i] = it;
        table[i][0] = sizes[i];
        it += sizes[i] + 1; // +1 for sizes element
    }
}

void IndexTableGlobal::addSequence (Sequence* s){
    // iterate over all k-mers of the sequence and add the id of s to the sequence list of the k-mer (tableDummy)
    unsigned int kmerIdx;
    this->size++; // amount of sequences added
    s->resetCurrPos();
    this->s = s;
    idxer->reset();
    
    while(s->hasNextKmer()){
        kmerIdx = idxer->int2index(s->nextKmer(), 0, kmerSize);
        table[kmerIdx][currPos[kmerIdx]++] = s->getId();
        for (int i = 0; i < skip && s->hasNextKmer(); i++){
            idxer->getNextKmerIndex(s->nextKmer(), kmerSize);
        }
    }
}

void IndexTableGlobal::removeDuplicateEntries(){
    
    delete[] currPos;
    
    for (size_t e = 0; e < tableSize; e++){
        if (sizes[e] == 0)
            continue;
        unsigned int* entries = table[e] + 1; // because of size at position 1
        int size = sizes[e];
        
        std::sort(entries, entries + size);
        // remove duplicates in-place
        int boundary = 1;
        for (int i = 1; i < size; i++){
            if (entries[i] != entries[i-1])
                entries[boundary++] = entries[i];
        }
        size = boundary;
        table[e][0]  = size;
    }
    delete[] sizes;

}

void IndexTableGlobal::print(char * int2aa){
    for (size_t i = 0; i < tableSize; i++){
        if (table[i][0] > 0){
            idxer->printKmer(i, kmerSize, int2aa);
            std::cout << "\n";
            for (unsigned int j = 0; j < table[i][0]; j++){
                std::cout << "\t" << table[i][j+1] << "\n";
            }
        }
    }
}

