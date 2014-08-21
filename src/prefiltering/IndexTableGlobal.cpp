#include "IndexTableGlobal.h"

IndexTableGlobal::IndexTableGlobal (int alphabetSize, int kmerSize, int skip)
                    : IndexTable(alphabetSize, kmerSize, skip) { }

IndexTableGlobal::~IndexTableGlobal(){ delete[] entries; }

void IndexTableGlobal::initMemory(){
    // allocate memory for the sequence id lists
    // tablesSizes is added to put the Size of the entry infront fo the memory
    entries = new char[(tableEntriesNum+tableSize)*sizeof(int)];
    
    std::fill_n(currPos, tableSize, 1); // needed because of size at beginning

}

void IndexTableGlobal::init(){
    unsigned int* it = (unsigned int *) entries;
    // set the pointers in the index table to the start of the list for a certain k-mer
    for (size_t i = 0; i < tableSize; i++){
        table[i] = (char *) it;
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
        unsigned int* row = (unsigned int *) table[kmerIdx];
        row[currPos[kmerIdx]++] = s->getId();
        for (int i = 0; i < skip && s->hasNextKmer(); i++){
            idxer->getNextKmerIndex(s->nextKmer(), kmerSize);
        }
    }
}

void IndexTableGlobal::removeDuplicateEntries(){
    delete[] currPos;
    currPos = NULL;
    
    for (size_t e = 0; e < tableSize; e++){
        if (sizes[e] == 0)
            continue;
        unsigned int* row = (unsigned int *) table[e];
        unsigned int* entries = row + 1; // because of size at position 1
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
    sizes = NULL;
}

void IndexTableGlobal::print(char * int2aa){
    for (size_t i = 0; i < tableSize; i++){
        unsigned int* row = (unsigned int *) table[i];
        if (row[0] > 0){
            idxer->printKmer(i, kmerSize, int2aa);
            std::cout << "\n";
            for (unsigned int j = 0; j < row[0]; j++){
                std::cout << "\t" << row[j+1] << "\n";
            }
        }
    }
}

void IndexTableGlobal::initTableByExternalData(uint64_t tableEntriesNum, unsigned short * sizes,
                                               unsigned int * pentries, unsigned int sequenzeCount){
    this->tableEntriesNum = tableEntriesNum;
    this->size = sequenzeCount;
    initMemory();
    memcpy ( this->entries , pentries, sizeof(unsigned int) * (this->tableEntriesNum  + this->tableSize));
    char* it = this->entries;
    // set the pointers in the index table to the start of the list for a certain k-mer
    for (size_t i = 0; i < tableSize; i++){
        table[i] = it;
        it += (sizes[i] + 1) * sizeof(int); // +1 for sizes element
    }
    delete [] this->sizes;
    this->sizes = NULL;
}