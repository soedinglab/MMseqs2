#include "IndexTableLocal.h"

IndexTableLocal::IndexTableLocal (int alphabetSize, int kmerSize, int skip)
: IndexTable(alphabetSize, kmerSize, skip) { }

IndexTableLocal::~IndexTableLocal(){ }

void IndexTableLocal::initMemory(){
    // allocate memory for the sequence id lists
    // tablesSizes is added to put the Size of the entry infront fo the memory
    entries = new char [tableEntriesNum * sizeof(IndexEntryLocal) +
                        tableSize * sizeof(int)];
    std::fill_n(currPos, tableSize, 0); // needed because of size at beginning

}

void IndexTableLocal::init(){
    char * it = entries;
    // set the pointers in the index table to the start of the list for a certain k-mer
    for (size_t i = 0; i < tableSize; i++){
        table[i] = (char *) it;
        table[i][0] = sizes[i];
        it += (sizes[i] * sizeof(IndexEntryLocal) + sizeof(unsigned int)); // jump to next kmer entrie; +1 for sizes element
    }
    delete[] sizes;
    sizes = NULL;
}

void IndexTableLocal::addSequence (Sequence* s){
    // iterate over all k-mers of the sequence and add the id of s to the sequence list of the k-mer (tableDummy)
    unsigned int kmerIdx;
    this->size++; // amount of sequences added
    s->resetCurrPos();
    this->s = s;
    idxer->reset();
    
    while(s->hasNextKmer()){
        kmerIdx = idxer->int2index(s->nextKmer(), 0, kmerSize);

        IndexEntryLocal * entry = (IndexEntryLocal *) (table[kmerIdx]
                                                       + (currPos[kmerIdx]++ * sizeof(IndexEntryLocal))
                                                       + sizeof(int));
        entry->seqId      = s->getId();
        entry->position_j = s->getCurrentPosition();
        
        //unsigned int* row = (unsigned int *) table[kmerIdx];
        //row[currPos[kmerIdx]++] = s->getId();
        for (int i = 0; i < skip && s->hasNextKmer(); i++){
            idxer->getNextKmerIndex(s->nextKmer(), kmerSize);
        }
    }
}

void IndexTableLocal::removeDuplicateEntries(){
    delete[] currPos;
    currPos = NULL;
    //  we dont remove entries for local search
}

void IndexTableLocal::print(char * int2aa){
    for (size_t i = 0; i < tableSize; i++){
        const unsigned int entrieSize = ((unsigned int*)table[i])[0];
        if (entrieSize > 0){
            idxer->printKmer(i, kmerSize, int2aa);
            std::cout << "\n";
            IndexEntryLocal * entries = (IndexEntryLocal *) (table[i] + sizeof(unsigned int));
            for (unsigned int j = 0; j < entrieSize; j++){
                std::cout << "\t(" << entries[i].seqId << ", " << entries[i].position_j << ")\n";
            }
        }
    }
}

void IndexTableLocal::initTableByExternalData(uint64_t tableEntriesNum, unsigned short * sizes,
                                              unsigned int * pentries, unsigned int sequenzeCount){
    this->tableEntriesNum = tableEntriesNum;
    this->size = sequenzeCount;
    initMemory();
    memcpy ( this->entries, pentries, sizeof(IndexEntryLocal) * this->tableEntriesNum + this->tableSize * sizeof(unsigned int));
    char* it = this->entries;
    // set the pointers in the index table to the start of the list for a certain k-mer
    for (size_t i = 0; i < tableSize; i++){
        table[i] = it;
        it += (sizes[i] * sizeof(IndexEntryLocal) + sizeof(unsigned int)); // jump to next kmer entrie; +1 for sizes element
    }
    delete [] this->sizes;
    this->sizes = NULL;
}