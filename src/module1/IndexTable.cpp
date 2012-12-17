#include "IndexTable.h"

IndexTable::IndexTable (size_t alphabetSize_, size_t kmerSize_):
    alphabetSize(alphabetSize_),
    kmerSize(kmerSize_)
{
    tableSize = ipow(alphabetSize, kmerSize);

    std::cout << "table size: " << tableSize << "\n";
    std::cout << "size_t size: " << sizeof(size_t) << "\n";
    std::cout << "size of list: " << sizeof(std::list<size_t>) << "\n";
    
    table = new size_t* [tableSize];
    
    tableDummy = new std::list<size_t>* [tableSize];
    for (size_t i = 0; i < tableSize; i++){
        tableDummy[i] = new std::list<size_t>();
    }

    idxer = new Indexer(alphabetSize, kmerSize);
    
    listSizes = new size_t[tableSize];
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
    for (size_t i = 0; i < tableSize; i++){
        std::list<size_t>* l = tableDummy[i];
        size_t* seqList = new size_t [l->size()];

        std::list<size_t>::iterator it;
        size_t j = 0;
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

size_t* IndexTable::getDBSeqList (size_t kmer, size_t* matchedListSize){
    *matchedListSize = listSizes[kmer];
    return table[kmer];
}

size_t IndexTable::ipow (size_t base, size_t exponent){
    size_t res = 1;
    for (size_t i = 0; i < exponent; i++)
        res = res*base;
    return res;
}
