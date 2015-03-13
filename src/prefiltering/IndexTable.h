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
#include "Debug.h"
#include "Util.h"

// IndexEntryLocal is an entry with possition and seqId for a kmer
// strucutre needs to be packed or it will need 8 bytes instead of 6
struct __attribute__((__packed__)) IndexEntryLocal {
    unsigned int seqId;
    unsigned short position_j;
};


class IndexTable {

    public:

        IndexTable (int alphabetSize, int kmerSize, int skip, size_t sizeOfEntry) {
            this->alphabetSize = alphabetSize;
            this->kmerSize = kmerSize;
            this->size = 0;
            this->skip = skip;
            this->sizeOfEntry = sizeOfEntry;
            tableSize = Util::ipow(alphabetSize, kmerSize);

            table = new char*[tableSize + 1]; // 1 + needed for the last pointer to calculate the size
            memset(table, 0, sizeof(char * ) * (tableSize + 1));

            idxer = new Indexer(alphabetSize, kmerSize);
        
            this->tableEntriesNum = 0;
            
            entries = NULL;
        }

        virtual ~IndexTable(){
            deleteEntries();
            delete[] table;
            delete idxer;
        }
    
        void deleteEntries(){
            if(entries != NULL){
                delete [] entries;
                entries = NULL;
            }
        }

        // count k-mers in the sequence, so enough memory for the sequence lists can be allocated in the end
        void addKmerCount (Sequence* s){
            unsigned int kmerIdx;
            s->resetCurrPos();
            idxer->reset();
        
            while(s->hasNextKmer()){
                kmerIdx = idxer->int2index(s->nextKmer(), 0, kmerSize);
                table[kmerIdx] += 1; // size increases by one
                tableEntriesNum++;
                for (int i = 0; i < skip && s->hasNextKmer(); i++){
                    idxer->getNextKmerIndex(s->nextKmer(), kmerSize);
                }
            }
            this->s = s;
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
        void initMemory(){
            // allocate memory for the sequence id lists
            // tablesSizes is added to put the Size of the entry infront fo the memory
            entries = new char [(tableEntriesNum + 1) * this->sizeOfEntry]; // +1 for table[tableSize] pointer address
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
        void initTableByExternalData(uint64_t tableEntriesNum, size_t * sizes,
                                     char * pentries, size_t sequenzeCount){
            this->tableEntriesNum = tableEntriesNum;
            this->size = sequenzeCount;
            initMemory();
            Debug(Debug::WARNING) << "Copy " << this->tableEntriesNum
                                  << " Entries (" <<  this->tableEntriesNum*this->sizeOfEntry  << " byte)\n";
            memcpy ( this->entries , pentries, this->tableEntriesNum * this->sizeOfEntry );
            Debug(Debug::WARNING) << "Setup Sizes  \n";
            char* it = this->entries;
            // set the pointers in the index table to the start of the list for a certain k-mer
            for (size_t i = 0; i < tableSize; i++){
                table[i] = it;
                it += sizes[i] * this->sizeOfEntry;
            }
            table[tableSize] = it;
            Debug(Debug::WARNING) << "Read IndexTable ... Done\n";


        }
    
        void revertPointer(){
            //TODO maybe not - 1
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
            Debug(Debug::WARNING) << "Empty Kmer list: " << emptyKmer << "\n";
            Debug(Debug::WARNING) << "\n";

        }
    
        // FUNCTIONS TO OVERWRITE
        // add k-mers of the sequence to the index table
        virtual void addSequence (Sequence* s) = 0;

        // removes dublicates
        virtual void removeDuplicateEntries() = 0;
    
        // prints the IndexTable
        virtual void print(char * int2aa) = 0;
    
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
    protected:
        // number of entries in all sequence lists
        int64_t tableEntriesNum; // must be 64bit
    
        // alphabetSize**kmerSize
        size_t tableSize;
    
        // Index table: contains pointers to the point in the entries array where starts the list of sequence ids for a certain k-mer
        char** __restrict table;

        // Index table entries: ids of sequences containing a certain k-mer, stored sequentially in the memory
        char* entries;

        Indexer* idxer;
    
        int alphabetSize;

        int kmerSize;

        Sequence* s;

        // number of skipped k-mers
        int skip;
    
        // amount of sequences in Index
        size_t size;
    
        // entry size
        size_t sizeOfEntry;
    
};

#endif
