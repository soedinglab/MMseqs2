#include "IndexTable.h"

IndexTable::IndexTable (int alphabetSize, int kmerSize, int skip)
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

    this->skip = skip;

    this->tableEntriesNum = 0;
}
IndexTable::~IndexTable(){
    delete[] entries;
    delete[] currPos;
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
        tableEntriesNum++;
        for (int i = 0; i < skip && s->hasNextKmer(kmerSize); i++){
            idxer->getNextKmerIndex(s->nextKmer(kmerSize), kmerSize);
        }
    }
    this->s = s;
}

void IndexTable::init(){
    // allocate memory for the sequence id lists
    entries = new int[tableEntriesNum];
    int* it = entries;
    // set the pointers in the index table to the start of the list for a certain k-mer
    for (int i = 0; i < tableSize; i++){
        if (sizes[i] > 0){
            table[i] = it;
            it += sizes[i];
        }
    }
}

void IndexTable::addSequence (Sequence* s){
    // iterate over all k-mers of the sequence and add the id of s to the sequence list of the k-mer (tableDummy)
    unsigned int kmerIdx;
    s->resetCurrPos();
    this->s = s;
    idxer->reset();

    int* workspace = new int(kmerSize);
    while(s->hasNextKmer(kmerSize)){
        kmerIdx = idxer->getNextKmerIndex(s->nextKmer(kmerSize), kmerSize);
        if(currPos[kmerIdx] >= sizes[kmerIdx]){
#pragma omp critical
            {
                std::cout << "ERROR\n";
                std::cout << "CurrPos: " << currPos[kmerIdx] << ", sizes: " << sizes[kmerIdx] << "\n";
                std::cout << "kmerIdx: " << kmerIdx << ", kmer: ";
                idxer->printKmer(workspace, kmerIdx, kmerSize, s->int2aa);
                std::cout << "\nSequence:\n";
                std::cout << "id: " << s->getId() << ", DB id: " << s->getDbKey() << "\n";
            }
        }
        table[kmerIdx][currPos[kmerIdx]++] = s->getId();
        for (int i = 0; i < skip && s->hasNextKmer(kmerSize); i++){
            idxer->getNextKmerIndex(s->nextKmer(kmerSize), kmerSize);
        }
    }
    delete[] workspace;
}

/*void IndexTable::removeDuplicateEntries(){

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

}*/

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
