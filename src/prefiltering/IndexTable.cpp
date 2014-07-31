#include "IndexTable.h"

IndexTable::IndexTable (int alphabetSize, int kmerSize, int skip)
{
    this->alphabetSize = alphabetSize;
    this->kmerSize = kmerSize;
    this->size = 0;
    this->skip = skip;

    tableSize = ipow(alphabetSize, kmerSize);

    sizes = new int[tableSize];
    memset(sizes, 0, sizeof(int) * tableSize);
    
    currPos = new int[tableSize];
    memset(currPos, 0, sizeof(int) * tableSize);

    table = new int*[tableSize];

    idxer = new Indexer(alphabetSize, kmerSize);

    this->tableEntriesNum = 0;
}
IndexTable::~IndexTable(){
    delete[] entries;
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
//    std::list<int>* sizesList = new std::list<int>();
    // allocate memory for the sequence id lists
    entries = new int[tableEntriesNum];
    int* it = entries;
    // set the pointers in the index table to the start of the list for a certain k-mer
    for (int i = 0; i < tableSize; i++){
//        if (i < tableSize - 1)
//            sizesList->push_back(sizes[i]);
        if (sizes[i] > 0){
            table[i] = it;
            it += sizes[i];
        }
    }
/*    sizesList->sort();
    std::cout << sizesList->front() << " " << sizesList->back() << "\n";
    std::list<int>::iterator itt = sizesList->end();
    itt--;
    itt--;
    int maxVal = *itt; 
    int maxPos = maxVal / 10 + 1;
    int* ar = new int[maxPos];
    memset(ar, 0, sizeof(int) * (maxPos));

    for (std::list<int>::iterator it = sizesList->begin(); it != sizesList->end(); it++){
        ar[*it/10]++;
    }
    std::ofstream outf;
    outf.open("/cluster/user/maria/test/kmer_cnt.out");
    for (int i = 0; i < maxPos; i++)
        outf << (i*10) << " " << ar[i] << "\n";
    outf.close();
    delete sizesList;
    delete[] ar;*/
}

void IndexTable::addSequence (Sequence* s){
    // iterate over all k-mers of the sequence and add the id of s to the sequence list of the k-mer (tableDummy)
    unsigned int kmerIdx;
    this->size++; // amount of sequences added
    s->resetCurrPos();
    this->s = s;
    idxer->reset();

    while(s->hasNextKmer(kmerSize)){
        kmerIdx = idxer->getNextKmerIndex(s->nextKmer(kmerSize), kmerSize);
        table[kmerIdx][currPos[kmerIdx]++] = s->getId();
        for (int i = 0; i < skip && s->hasNextKmer(kmerSize); i++){
            idxer->getNextKmerIndex(s->nextKmer(kmerSize), kmerSize);
        }
    }
}

void IndexTable::removeDuplicateEntries(){

    delete[] currPos;
//    std::list<int>* sizesList = new std::list<int>();
    
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
        sizes[e] = size;
//        if (e < tableSize - 1)
//            sizesList->push_back(sizes[e]);

    }
/*    sizesList->sort();
    std::cout << sizesList->front() << " " << sizesList->back() << "\n";
    std::list<int>::iterator itt = sizesList->end();
    itt--;
    itt--;
    int maxVal = *itt;
    int maxPos = maxVal / 10 + 1;
    int* ar = new int[maxPos];
    memset(ar, 0, sizeof(int) * (maxPos));

    for (std::list<int>::iterator it = sizesList->begin(); it != sizesList->end(); it++){
        ar[*it/10]++;
    }
    std::ofstream outf;
    outf.open("/cluster/user/maria/test/kmer_cnt_with_rm.out");
    for (int i = 0; i < maxPos; i++)
        outf << (i*10) << " " << ar[i] << "\n";
    outf.close();
    delete sizesList;
    delete[] ar;
*/

}

void IndexTable::print(){
    for (int i = 0; i < tableSize; i++){
        if (sizes[i] > 0){
            idxer->printKmer(i, kmerSize, s->int2aa);
            std::cout << "\n";
            for (int j = 0; j < sizes[i]; j++){
                std::cout << "\t" << table[i][j] << "\n";
            }
        }
    }
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
