	#include "IndexTableGlobal.h"

IndexTableGlobal::IndexTableGlobal (int alphabetSize, int kmerSize, int skip)
                    : IndexTable(alphabetSize, kmerSize, skip, sizeof(int)) { }

IndexTableGlobal::~IndexTableGlobal(){ delete[] entries; }


void IndexTableGlobal::addSequence (Sequence* s){
    // iterate over all k-mers of the sequence and add the id of s to the sequence list of the k-mer (tableDummy)
    unsigned int kmerIdx;
    this->size++; // amount of sequences added
    s->resetCurrPos();
    this->s = s;
    idxer->reset();
    
    while(s->hasNextKmer()){
        kmerIdx = idxer->int2index(s->nextKmer(), 0, kmerSize);
        table[kmerIdx][0] = s->getId();
        table[kmerIdx] += sizeof(unsigned int); // move one element further
        for (int i = 0; i < skip && s->hasNextKmer(); i++){
            idxer->getNextKmerIndex(s->nextKmer(), kmerSize);	
        }
    }
}


// removes duplicate entries
void IndexTableGlobal::removeDuplicateEntries(){
    unsigned int* entriesWrite = (unsigned int *) table[0];
    unsigned int writePostion = 0; // must be global to pack memory sequential
    
    // Example:
    // 1 2 2 3 | 2 3 3
    // ^         ^
    // will be transformed into
    // 1 2 3 | 2 3 | x x
    // ^       ^
    for (size_t e = 0; e < tableSize; e++){
        const ptrdiff_t size = (table[e + 1] - table[e]) / sizeof(unsigned int); // size of an entry
        unsigned int* entriesRead = (unsigned int *) table[e];
        // update table pointer position
        table[e] = (char *) (entriesWrite + writePostion);

        if (size == 0)
            continue;
        
        std::sort(entriesRead, entriesRead + size);
        // remove duplicates in-place
        entriesWrite[writePostion++] = entriesRead[0]; // copy first element
        for (int i = 1; i < size; i++){
            if (entriesRead[i] != entriesRead[i-1])
                entriesWrite[writePostion++] = entriesRead[i];
        }
    }
    // set pointer to last element
    table[tableSize] = (char *) (entriesWrite + writePostion);

}

void IndexTableGlobal::print(char * int2aa){
    for (size_t i = 0; i < tableSize; i++){
        unsigned int * row = (unsigned int *) table[i];
        ptrdiff_t size = (table[i+1] - table[i])/sizeof(unsigned int);
        
        if (size > 0){
            idxer->printKmer(i, kmerSize, int2aa);
//            std::cout << " " << size << " " << nextRow << " ";
            std::cout << " size(= " << size << ")\n";
            for (unsigned int j = 0; j < size; j++){
                std::cout << "\t" << row[j] << "\n";
            }
        }
    }
}

