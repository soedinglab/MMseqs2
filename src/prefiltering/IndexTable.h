#ifndef INDEX_TABLE_H
#define INDEX_TABLE_H

//
// Written by Martin Steinegger martin.steinegger@campus.lmu.de and Maria Hauser mhauser@genzentrum.lmu.de
// 
// Abstract: Index table stores the list of DB sequences containing a certain k-mer, for each k-mer. 
//

#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>

#include "Sequence.h"
#include "Indexer.h"
#include "Util.h"


class IndexTable {

    public:

        IndexTable (int alphabetSize, int kmerSize, int skip) {
            this->alphabetSize = alphabetSize;
            this->kmerSize = kmerSize;
            this->size = 0;
            this->skip = skip;
        
            tableSize = Util::ipow(alphabetSize, kmerSize);
        
            sizes = new unsigned short[tableSize];
            memset(sizes, 0, sizeof(short) * tableSize);
        
            currPos = new int[tableSize];
            std::fill_n(currPos, tableSize, 1); // needed because of size at beginning
        
            table = new char*[tableSize];
        
            idxer = new Indexer(alphabetSize, kmerSize);
        
            this->tableEntriesNum = 0;
        }

        virtual ~IndexTable(){
            delete[] entries;
            delete[] table;
            delete idxer; }

        // count k-mers in the sequence, so enough memory for the sequence lists can be allocated in the end
        void addKmerCount (Sequence* s){
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
    
        // get list of DB sequences containing this k-mer
        inline unsigned int* getDBSeqList (int kmer, int* matchedListSize){
            unsigned int * __restrict tabPosition = (unsigned int *) table[kmer];
            *matchedListSize = tabPosition[0];
            return tabPosition + 1;
        }
    
        // get pointer to sizes array
        unsigned short * getSizes(){
            if(sizes == NULL)
                std::cerr << "AAAAAH" << std::endl;
            return sizes;
        }
    
        // get pointer to entries array
        char * getEntries(){
            return entries;
        }
    
        // FUNCTIONS TO OVERWRITE
        // add k-mers of the sequence to the index table
        virtual void addSequence (Sequence* s) = 0;

        // removes dublicates
        virtual void removeDuplicateEntries() = 0;

        // init the arrays for the sequence lists 
        virtual void init() = 0;
    
        // allocates memory for index tables
        virtual void initMemory() = 0;
    
        // prints the IndexTable
        virtual void print(char * int2aa) = 0;
    
        // init index table with external data (needed for index readin)
        virtual void initTableByExternalData(uint64_t tableEntriesNum, unsigned short * sizes,
                                         unsigned int * pentries, unsigned int sequenzeCount) = 0;
    
        // get amount of sequences in Index
        unsigned int getSize() {  return size; };
    
        // returns the size of  table entries
        int64_t getTableEntriesNum(){ return tableEntriesNum; };
    
        // returns table size
        unsigned int getTableSize(){ return tableSize; };

    protected:
        // number of entries in all sequence lists
        int64_t tableEntriesNum; // must be 64bit
    
        // alphabetSize**kmerSize
        unsigned int tableSize;
    
        // Index table: contains pointers to the point in the entries array where starts the list of sequence ids for a certain k-mer
        char** __restrict table;

        // Index table entries: ids of sequences containing a certain k-mer, stored sequentially in the memory
        char* entries;

        unsigned short* sizes;

        // only for init: current position in the DB id array of the index table where the next sequence id can be written
        int* currPos;

        Indexer* idxer;
    
        int alphabetSize;

        int kmerSize;

        Sequence* s;

        // number of skipped k-mers
        int skip;
    
        // amount of sequences in Index
        unsigned int size;
    


};

#endif
