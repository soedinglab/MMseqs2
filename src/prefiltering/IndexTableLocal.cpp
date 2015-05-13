#include "IndexTableLocal.h"

IndexTableLocal::IndexTableLocal (int alphabetSize, int kmerSize, int skip)
: IndexTable(alphabetSize, kmerSize, skip, sizeof(IndexEntryLocal)) { }

IndexTableLocal::~IndexTableLocal(){ }


void IndexTableLocal::addSequence (Sequence* s){
    // iterate over all k-mers of the sequence and add the id of s to the sequence list of the k-mer (tableDummy)
    unsigned int kmerIdx;
    this->size++; // amount of sequences added
    s->resetCurrPos();
    this->s = s;
    idxer->reset();
    
    while(s->hasNextKmer()){
        kmerIdx = idxer->int2index(s->nextKmer(), 0, kmerSize);

        IndexEntryLocal * entry = (IndexEntryLocal *) (table[kmerIdx]);
        table[kmerIdx] += sizeof(IndexEntryLocal);
        entry->seqId      = s->getId();
        entry->position_j = s->getCurrentPosition() % 256;
        
        //unsigned int* row = (unsigned int *) table[kmerIdx];
        //row[currPos[kmerIdx]++] = s->getId();
        for (int i = 0; i < skip && s->hasNextKmer(); i++){
            idxer->getNextKmerIndex(s->nextKmer(), kmerSize);
        }
    }
}

void IndexTableLocal::removeDuplicateEntries(){
    //  we dont remove entries for local search
}

void IndexTableLocal::print(char * int2aa){
    for (size_t i = 0; i < tableSize; i++){
        ptrdiff_t entrieSize = (table[i+1] - table[i]) / sizeof(IndexEntryLocal);

        if (entrieSize > 0){
            idxer->printKmer(i, kmerSize, int2aa);
            Debug(Debug::INFO) << "\n";
            IndexEntryLocal * entries = (IndexEntryLocal *) table[i];
            for (unsigned int j = 0; j < entrieSize; j++){
                Debug(Debug::INFO) << "\t(" << entries[i].seqId << ", " << entries[i].position_j << ")\n";
            }
        }
    }
}
