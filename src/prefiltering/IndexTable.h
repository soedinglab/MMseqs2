#ifndef INDEX_TABLE_H
#define INDEX_TABLE_H

//
// Written by Martin Steinegger martin.steinegger@mpibpc.mpg.de and Maria Hauser mhauser@genzentrum.lmu.de
//
// Abstract: Index table stores the list of DB sequences containing a certain k-mer, for each k-mer.
//

#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>
#include <sys/mman.h>
#include <new>
#include <DBReader.h>
#ifdef OPENMP
#include <omp.h>
#endif
#include "omptl/omptl_algorithm"
#include "Sequence.h"
#include "Indexer.h"
#include "Debug.h"
#include "Util.h"
#include "SequenceLookup.h"
#include "MathUtil.h"

// IndexEntryLocal is an entry with position and seqId for a kmer
// structure needs to be packed or it will need 8 bytes instead of 6
struct __attribute__((__packed__)) IndexEntryLocal {
    unsigned int seqId;
    unsigned short position_j;

    static bool comapreByIdAndPos(IndexEntryLocal first, IndexEntryLocal second){
        if(first.seqId < second.seqId )
            return true;
        if(second.seqId < first.seqId )
            return false;
        if(first.position_j < second.position_j )
            return true;
        if(second.position_j < first.position_j )
            return false;
        return false;
    }
};

class IndexTable {

public:
    IndexTable (int alphabetSize, int kmerSize) {
        this->alphabetSize = alphabetSize;
        this->kmerSize = kmerSize;
        this->size = 0;
        this->sizeOfEntry = sizeof(IndexEntryLocal);
        tableSize = MathUtil::ipow(alphabetSize, kmerSize);
        table = new(std::nothrow) char*[tableSize + 1]; // 1 + needed for the last pointer to calculate the size
        Util::checkAllocation(table, "Could not allocate table memory in IndexTable");
        memset(table, 0, sizeof(char * ) * (tableSize + 1)); // set all pointers to 0
        idxer = new Indexer(alphabetSize, kmerSize);
        this->tableEntriesNum = 0;
        entries = NULL;
        sequenceLookup = NULL;
    }

    virtual ~IndexTable(){
        deleteEntries();
        delete[] table;
        delete idxer;
        if(sequenceLookup != NULL){
            delete sequenceLookup;
        }
    }

    void deleteEntries(){
        if(entries != NULL && externalData == false){
            delete [] entries;
            entries = NULL;
        }else{
            munlock(entries, tableEntriesNum);
        }
    }

    // count k-mers in the sequence, so enough memory for the sequence lists can be allocated in the end
    size_t addKmerCount (Sequence* s, Indexer * idxer){
        unsigned int kmerIdx;
        s->resetCurrPos();
        //idxer->reset();
        size_t countKmer = 0;
        while(s->hasNextKmer()){
            kmerIdx = idxer->int2index(s->nextKmer(), 0, kmerSize);
            //table[kmerIdx] += 1;
            // size increases by one
            __sync_fetch_and_add( (int *) &table[kmerIdx], 1 );
            countKmer++;
        }
        return countKmer;
    }

    inline  char * getTable(unsigned int kmer){
        return table[kmer];
    }

    // get list of DB sequences containing this k-mer
    template<typename T> inline T* getDBSeqList (int kmer, size_t* matchedListSize){
        const ptrdiff_t diff =  (table[kmer + 1] - table[kmer]) / sizeof( T );
        *matchedListSize = diff;
        return (T *) table[kmer];
    }

    // get pointer to entries array
    char * getEntries(){
        return entries;
    }

    // init the arrays for the sequence lists
    void initMemory(size_t sequenzeCount, size_t tableEntriesNum, size_t aaDbSize, SequenceLookup * seqLookup, size_t dbSize) {
        this->tableEntriesNum = tableEntriesNum;
        this->size = dbSize; // amount of sequences added

        if(seqLookup != NULL){
            sequenceLookup = seqLookup;
        }
        // allocate memory for the sequence id lists
        // tablesSizes is added to put the Size of the entry infront fo the memory
        entries = new(std::nothrow) char [(tableEntriesNum + 1) * this->sizeOfEntry]; // +1 for table[tableSize] pointer address
        externalData = false;
        Util::checkAllocation(entries, "Could not allocate entries memory in IndexTable::initMemory");
    }

    // allocates memory for index tables
    void init(){
        char * it = entries;
        // set the pointers in the index table to the start of the list for a certain k-mer
        for (size_t i = 0; i < tableSize; i++){
            const size_t entriesCount = (size_t) table[i];
            table[i] = (char *) it;
            it += (entriesCount * this->sizeOfEntry);
        }
        table[tableSize] = it;
    }

    // init index table with external data (needed for index readin)
    void initTableByExternalData(size_t sequenzeCount, size_t tableEntriesNum,
                                 char * entries, size_t * entriesSize, SequenceLookup * lookup) {
        this->tableEntriesNum = tableEntriesNum;
        this->size = sequenzeCount;
        //        initMemory(sequenzeCount, tableEntriesNum, seqDataSize);
        if(lookup != NULL){
            sequenceLookup = lookup;
        }
        this->entries = entries;
        Debug(Debug::WARNING) << "Cache database  \n";
        char* it = this->entries;
        // set the pointers in the index table to the start of the list for a certain k-mer
        magicByte = 0; // read each entry to keep them in memory
        for (size_t i = 0; i < tableSize; i++){
            size_t entrySize = entriesSize[i] * this->sizeOfEntry;
            table[i] = it;
            magicByte += *table[i];
            it += entrySize;
        }
        table[tableSize] = it;
        externalData = true;
        Debug(Debug::WARNING) << "Read IndexTable ... Done\n";
    }

    void revertPointer(){
        for(size_t i = tableSize - 1; i > 0; i--){
            table[i] = table[i-1];
        }
        table[0] = entries;
    }

    void printStatisitic(char * int2aa){
        size_t minKmer = 0;

        size_t entries = 0;
        double avgKmer = 0;

        size_t emptyKmer = 0;
        const size_t top_N = 10;
        std::pair<size_t, size_t> topElements[top_N];
        for(size_t j =0; j < top_N; j++)
            topElements[j].first = 0;

        for(size_t i = 0; i < tableSize - 1; i++){
            const ptrdiff_t size =  (table[i + 1] - table[i]) / this->sizeOfEntry;
            minKmer = std::min(minKmer, (size_t)size);
            entries += size;
            if(size == 0){
                emptyKmer++;
            }
            if(((size_t)size) < topElements[top_N-1].first)
                continue;
            for(size_t j =0; j < top_N; j++){
                if (topElements[j].first < ((size_t)size)) {
                    topElements[j].first = size;
                    topElements[j].second = i;
                    break;
                }
            }
        }
        avgKmer = ((double)entries) / ((double)tableSize);
        Debug(Debug::WARNING) << "DB statistic\n";
        Debug(Debug::WARNING) << "Entries:         " << entries << "\n";
        Debug(Debug::WARNING) << "DB Size:         " << entries * sizeOfEntry + tableSize * sizeof(int *) << " (byte)\n";
        Debug(Debug::WARNING) << "Avg Kmer Size:   " << avgKmer << "\n";
        Debug(Debug::WARNING) << "Top " << top_N << " Kmers\n   ";
        for(size_t j =0; j < top_N; j++){
            Debug(Debug::WARNING) << "\t";

            this->idxer->printKmer(topElements[j].second, kmerSize, int2aa);
            Debug(Debug::WARNING) << "\t\t" << topElements[j].first << "\n";
        }
        Debug(Debug::WARNING) << "Min Kmer Size:   " << minKmer << "\n";
        Debug(Debug::WARNING) << "Empty list: " << emptyKmer << "\n";
        Debug(Debug::WARNING) << "\n";

    }

    // FUNCTIONS TO OVERWRITE
    // add k-mers of the sequence to the index table
    void addSequence (Sequence* s, Indexer * idxer, size_t aaFrom, size_t aaSize){
        // iterate over all k-mers of the sequence and add the id of s to the sequence list of the k-mer (tableDummy)
        s->resetCurrPos();
        idxer->reset();
        while(s->hasNextKmer()){
            const int * kmer = s->nextKmer();
            unsigned int kmerIdx = idxer->int2index(kmer, 0, kmerSize);
            if(kmerIdx >= aaFrom  && kmerIdx < aaFrom + aaSize){
                // if region got masked do not add kmer
                if((table[kmerIdx+1] - table[kmerIdx]) == 0)
                    continue;
                IndexEntryLocal * entry = (IndexEntryLocal *) (table[kmerIdx]);
                entry->seqId      = s->getId();
                entry->position_j = s->getCurrentPosition();
                table[kmerIdx] += sizeof(IndexEntryLocal);
            }
        }
    }

    // prints the IndexTable
    void print(char * int2aa) {
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
    };

    // get amount of sequences in Index
    size_t getSize() {  return size; };

    // returns the size of  table entries
    int64_t getTableEntriesNum(){ return tableEntriesNum; };

    // returns table size
    size_t getTableSize(){ return tableSize; };

    // returns table
    char ** getTable(){ return table; };

    // returns the size of the entry (int for global) (IndexEntryLocal for local)
    size_t getSizeOfEntry() { return sizeOfEntry; }

    SequenceLookup *getSequenceLookup(){ return sequenceLookup; }

    int getKmerSize() {
        return kmerSize;
    }

    int getAlphabetSize() {
        return alphabetSize;
    }

    static int computeKmerSize(size_t aaSize) {
        return aaSize < getUpperBoundAACountForKmerSize(6) ? 6 : 7;
    }


    static size_t getUpperBoundAACountForKmerSize(int kmerSize) {
        switch(kmerSize){
            case 6:
                return 3350000000;
            case 7:
                return (UINT_MAX - 1); // UINT_MAX is often reserved as safe flag
        }
        return 0;
    }


protected:
    // number of entries in all sequence lists
    int64_t tableEntriesNum; // must be 64bit

    // alphabetSize**kmerSize
    size_t tableSize;

    // Index table: contains pointers to the k-mer start position in entries array
    // Needed for fast access
    char** __restrict table;

    // Index table entries: ids of sequences containing a certain k-mer, stored sequentially in the memory
    char* entries;

    Indexer* idxer;

    int alphabetSize;

    int kmerSize;

    // amount of sequences in Index
    size_t size;

    // entry size
    size_t sizeOfEntry;

    // sequence lookup
    SequenceLookup *sequenceLookup;

    // external data from mmap
    bool externalData;

    // magic byte to avoid compiler optimisation
    size_t magicByte;
};




#endif
