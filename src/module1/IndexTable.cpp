#include "IndexTable.h"

IndexTable::IndexTable (int alphabetSize_, int kmerSize_):
    alphabetSize(alphabetSize_),
    kmerSize(kmerSize_)
{
    tableSize = ipow(alphabetSize, kmerSize);

    std::cout << "table size: " << tableSize << "\n";
    std::cout << "int size: " << sizeof(int) << "\n";
    std::cout << "vector size: " << sizeof(std::vector<int>) << "\n";
    
    table = new int* [tableSize];
    
    tableDummy = new std::vector<int>* [tableSize];
    tableDummy[0] = new std::vector<int>();
    std::cout << "Initial capacity: " << tableDummy[0]->capacity() << "\n";
    std::cout << "Initial sizeof(): " << sizeof(*tableDummy[0]) << "\n";
    for (int i = 0; i < 10; i++){
        tableDummy[0]->push_back(i);
    }
    std::cout << "Capacity: " << tableDummy[0]->capacity() << "\n";
    std::cout << "sizeof(): " << sizeof(*tableDummy[0]) << "\n";


    for (int i = 0; i < tableSize; i++){
        tableDummy[i] = new std::vector<int>();
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
        std::vector<int>* v = tableDummy[i];
        int* seqList = new int [v->size()];

        std::vector<int>::iterator it;
        int j = 0;
        for (it = v->begin(); it != v->end(); it++){
            seqList[j] = *it;
            j++;
        }

        listSizes[i] = j;
        table[i] = seqList;

        delete v;
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
