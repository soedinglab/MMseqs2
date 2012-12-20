#include "IndexTable.h"

IndexTable::IndexTable (int alphabetSize_, int kmerSize_):
    alphabetSize(alphabetSize_),
    kmerSize(kmerSize_)
{
    tableSize = ipow(alphabetSize, kmerSize);

    std::cout << "table size: " << tableSize << "\n";
    std::cout << "int size: " << sizeof(int) << "\n";
    std::cout << "size of list: " << sizeof(std::list<int>) << "\n";
    
    table = new int* [tableSize];
    
    tableDummy = new std::list<int>* [tableSize];
    for (int i = 0; i < tableSize; i++){
        tableDummy[i] = new std::list<int>();
    }

    idxer = new Indexer(alphabetSize, kmerSize);
    
    listSizes = new int[tableSize];
}

void IndexTable::addSequence (Sequence* s){
    // iterate over all k-mers of the sequence and add the id of s to the sequence list of the k-mer (tableDummy)
    const int* kmer;
    int kmerIdx;

    s->reset();
    while(s->hasNextKmer(kmerSize)){
        kmer = s->nextKmer(kmerSize);
        kmerIdx = idxer->int2index(kmer, 0, kmerSize);
        tableDummy[kmerIdx]->push_back(s->id);
    }
}

void IndexTable::init(){
    for (int i = 0; i < tableSize; i++){
        std::list<int>* l = tableDummy[i];
        int* seqList = new int [l->size()];

        std::list<int>::iterator it;
        int j = 0;
        for (it = l->begin(); it != l->end(); it++){
            seqList[j] = *it;
            j++;
        }

        listSizes[i] = j;
        table[i] = seqList;

        delete l;
    }
    delete[] tableDummy;
}

int* IndexTable::getDBSeqList (int kmer, int* matchedListSize){
    *matchedListSize = listSizes[kmer];
    return table[kmer];
}

int IndexTable::ipow (int base, int exponent){
    int res = 1;
    for (int i = 0; i < exponent; i++)
        res = res*base;
    return res;
}
