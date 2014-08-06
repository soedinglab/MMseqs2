#include "IndexTable.h"
/*
IndexTable::IndexTable (int alphabetSize, int kmerSize, int skip)
{
    this->alphabetSize = alphabetSize;
    this->kmerSize = kmerSize;
    this->size = 0;
    this->skip = skip;
    
    tableSize = ipow(alphabetSize, kmerSize);
    
    sizes = new unsigned short[tableSize];
    memset(sizes, 0, sizeof(short) * tableSize);
    
    currPos = new int[tableSize];
    std::fill_n(currPos, tableSize, 0); // needed because of size at beginning
    
    table = new unsigned int*[tableSize];
    
    idxer = new Indexer(alphabetSize, kmerSize);
    
    this->tableEntriesNum = 0;
}
IndexTable::~IndexTable(){
    delete[] entries;
    delete[] table;
    delete idxer;
}

unsigned short * IndexTable::getSizes(){
    if(sizes == NULL)
        std::cerr << "AAAAAH" << std::endl;
    return sizes;
}

void IndexTable::addKmerCount (Sequence* s){
    unsigned int kmerIdx;
    s->resetCurrPos();
    idxer->reset();
    
    while(s->hasNextKmer()){
        kmerIdx = idxer->int2index(s->nextKmer(), 0, kmerSize);
        sizes[kmerIdx]++;
        tableEntriesNum++;
        for (int i = 0; i < skip && s->hasNextKmer(); i++){
            idxer->getNextKmerIndex(s->nextKmer(), kmerSize);
        }
    }
    this->s = s;
}

void IndexTable::initMemory(){
    // allocate memory for the sequence id lists
    // tablesSizes is added to put the Size of the entry infront fo the memory
    entries = new char [tableEntriesNum * sizeof(IndexEntry) +
                                         tableSize * sizeof(int)];
}


void IndexTable::init(){
    char * it = entries;
    // set the pointers in the index table to the start of the list for a certain k-mer
    for (size_t i = 0; i < tableSize; i++){
        table[i] = it;
        int * row = table[i];
        row[0] = sizes[i];
        it += (sizes[i] * sizeof(IndexEntry)) + sizeof(int); // jump to next kmer entrie; +1 for sizes element
    }
    
    delete[] sizes;
    sizes = NULL;

}

KmerEntry * IndexTable::getEntries(){
    return entries;
}

void IndexTable::addSequence (Sequence* s){
    // iterate over all k-mers of the sequence and add the id of s to the sequence list of the k-mer (tableDummy)
    unsigned int kmerIdx;
    this->size++; // amount of sequences added
    s->resetCurrPos();
    this->s = s;
    idxer->reset();
    
    while(s->hasNextKmer()){
        kmerIdx = idxer->int2index(s->nextKmer(), 0, kmerSize);
        IndexEntryLocal entry = *(table[kmerIdx]+(currPos[kmerIdx]++ * sizeof(IndexEntry)) + sizeof(int));
        entry.seqId      = s->getId();
        entry.position_j = s->getCurrentPosition();
        // if skip is activated skip next positions
        for (int i = 0; i < skip && s->hasNextKmer(); i++){
            idxer->getNextKmerIndex(s->nextKmer(), kmerSize);
        }
    }
    // current position not needed
    delete[] currPos;

}

void IndexTable::removeDuplicateEntries(){
    
   // we dont remove in case of local search

}

void IndexTable::print(char * int2aa){
    for (size_t i = 0; i < tableSize; i++){
        const unsigned int tableSize = ((unsigned int*)table[i])[0];
        if ( > 0){
            idxer->printKmer(i, kmerSize, int2aa);
            std::cout << "\n";
            IndexEntry * entries = table[i] + sizeof(unsigned int);
            for (unsigned int j = 0; j < tableSize; j++){
                std::cout << "\t(" << entries[i].seqId << ", " << entries[i].position_j << ")\n";
            }
        }
    }
}


void IndexTable::initTableByExternalData(uint64_t tableEntriesNum,
                                         unsigned short * sizes, unsigned int * pentries, unsigned int sequenzeCount){
    this->tableEntriesNum = tableEntriesNum;
    this->size = sequenzeCount;
    initMemory();
    memcpy ( this->entries , pentries, sizeof(unsigned int) * (this->tableEntriesNum  + this->tableSize));
    unsigned int* it = this->entries;
    // set the pointers in the index table to the start of the list for a certain k-mer
    for (size_t i = 0; i < tableSize; i++){
        table[i] = it;
        it += sizes[i] + 1; // +1 for sizes element
    }
    delete [] this->sizes;
    this->sizes = NULL;
}

int IndexTable::ipow (int base, int exponent){
    int res = 1;
    for (int i = 0; i < exponent; i++)
        res = res * base;
    return res;
}*/
