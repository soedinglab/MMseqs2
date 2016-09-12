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
#include <Log.h>
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
    
    IndexTable (int alphabetSize, int kmerSize, bool hasSequenceLookup) {
        this->alphabetSize = alphabetSize;
        this->kmerSize = kmerSize;
        this->size = 0;
        this->hasSequenceLookup = hasSequenceLookup;
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
    void addKmerCount (Sequence* s){
        unsigned int kmerIdx;
        s->resetCurrPos();
        idxer->reset();
        
        while(s->hasNextKmer()){
            kmerIdx = idxer->int2index(s->nextKmer(), 0, kmerSize);
            table[kmerIdx] += 1; // size increases by one
            tableEntriesNum++;
        }
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
    void initMemory(size_t sequenzeCount, size_t tableEntriesNum, size_t aaDbSize) {
        this->tableEntriesNum = tableEntriesNum;
        if(hasSequenceLookup == true){
            sequenceLookup = new SequenceLookup(sequenzeCount, aaDbSize);
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
        if(hasSequenceLookup == true && lookup != NULL){
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
    void addSequence (Sequence* s){
        // iterate over all k-mers of the sequence and add the id of s to the sequence list of the k-mer (tableDummy)
        unsigned int kmerIdx;
        this->size++; // amount of sequences added
        s->resetCurrPos();
        idxer->reset();
        while(s->hasNextKmer()){
            kmerIdx = idxer->int2index(s->nextKmer(), 0, kmerSize);
            // if region got masked do not add kmer
            if((table[kmerIdx+1] - table[kmerIdx]) == 0)
                continue;
            IndexEntryLocal * entry = (IndexEntryLocal *) (table[kmerIdx]);
            entry->seqId      = s->getId();
            entry->position_j = s->getCurrentPosition();
            table[kmerIdx] += sizeof(IndexEntryLocal);
        }
        if(hasSequenceLookup == true){
            sequenceLookup->addSequence(s);
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
    
    
    struct __attribute__((__packed__)) TableEntry{
        unsigned int kmer;
        unsigned int seqId;
        unsigned short position_j;
    };
    
    static void fillEntries(DBReader<unsigned int> *dbr, size_t aaSize, size_t dbFrom, size_t dbTo,
                            Sequence *seq, SequenceLookup * sequenceLookup,
                            IndexTable::TableEntry *entries, size_t tableEntrieNum, size_t alphabetSize, size_t threads) {
#pragma omp parallel
        {
            Sequence s(seq->getMaxLen(), seq->aa2int, seq->int2aa,
                       seq->getSeqType(), seq->getKmerSize(), seq->isSpaced(), false);
            Indexer idxer(alphabetSize, seq->getKmerSize());
            
            size_t threadFrom, threadSize;
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            Util::decomposeDomainByAminoAcid(aaSize, dbr->getSeqLens() + dbFrom, dbTo - dbFrom,
                                             thread_idx, threads, &threadFrom, &threadSize);
            size_t kmerStartPos = 0;
            for(size_t id = 0; id < threadFrom; id++ ){
                kmerStartPos +=  Util::overlappingKmers(static_cast<int>(dbr->getSeqLens(id)) - 2, s.getEffectiveKmerSize());
            }
            std::cout << thread_idx << "\t" << threadFrom << "\t" << threadSize << "\t" << kmerStartPos << std::endl;
            TableEntry * threadEntry = entries + kmerStartPos;
            for (unsigned int id = threadFrom; id < (threadFrom + threadSize); id++){
                s.resetCurrPos();
                Log::printProgress(id);
                unsigned int qKey = dbr->getDbKey(id);
                s.mapSequence(id, qKey, sequenceLookup->getSequence(dbFrom + id));
                idxer.reset();
                while(s.hasNextKmer()){
                    size_t kmerIdx = idxer.int2index(s.nextKmer(), 0, s.getKmerSize());
                    threadEntry->kmer       = kmerIdx;
                    threadEntry->seqId      = s.getId();
                    threadEntry->position_j = s.getCurrentPosition();
                    threadEntry++;
                }
            }
        }
        
    }
    
    
    static bool compareTableEntryByKmer(TableEntry first, TableEntry second){
        return (first.kmer < second.kmer) ? true : false;
    }

    static void parallelFillDatabase(DBReader<unsigned int>* dbr, Sequence* seq, IndexTable * indexTable,
                                     BaseMatrix *subMat, size_t dbFrom, size_t dbTo, size_t threads)
    {
        size_t dbSize = dbTo - dbFrom;
        size_t * sequenceOffSet = new size_t[dbSize];
        size_t aaDbSize = 0;
        size_t tableEntriesNum = 0;
        sequenceOffSet[0] = 0;
        for (unsigned int id = dbFrom; id < dbTo; id++){
            int seqLen = std::max(static_cast<int>(dbr->getSeqLens(id)) - 2, 0);
            aaDbSize += seqLen; // remove /n and /0
            size_t idFromNull = (id - dbFrom);
            if(id < dbTo - 1){
                sequenceOffSet[idFromNull + 1] =  sequenceOffSet[idFromNull] + seqLen;
            }
            tableEntriesNum += Util::overlappingKmers(seqLen, seq->getEffectiveKmerSize());
        }
        
        SequenceLookup * sequenceLookup = new SequenceLookup(dbSize, aaDbSize);
        size_t maskedResidues = 0;
        size_t aaCount = 0;
#pragma omp parallel
        {
            Sequence s(seq->getMaxLen(), seq->aa2int, seq->int2aa,
                       seq->getSeqType(), seq->getKmerSize(), seq->isSpaced(), false);
#pragma omp for schedule(dynamic, 100) reduction(+: aaCount, maskedResidues)
            for (unsigned int id = dbFrom; id < dbTo; id++){
                s.resetCurrPos();
                Log::printProgress(id - dbFrom);
                char* seqData = dbr->getData(id);
                unsigned int qKey = dbr->getDbKey(id);
                s.mapSequence(id - dbFrom, qKey, seqData);
                maskedResidues += Util::maskLowComplexity(subMat, &s, s.L, 12, 3,
                                                          indexTable->getAlphabetSize(), s.aa2int['X']);
                sequenceLookup->addSequence(&s, sequenceOffSet[id-dbFrom]);
                aaCount += s.L;
            }
        }
        delete [] sequenceOffSet;
        dbr->remapData();
        if(aaCount != aaDbSize){
            Debug(Debug::ERROR) << "Index table: aaCount (" << aaCount << ") != aaDbSize (" << aaDbSize << ")\n";
        }
        Debug(Debug::INFO) << "\nIndex table: Masked residues: " << maskedResidues << "\n";
        // if memory is enough to accommodate 10 bytes pro entry
        TableEntry * sortableEntries = (TableEntry *) malloc(sizeof(TableEntry) * tableEntriesNum);
        fillEntries(dbr, aaDbSize, dbFrom, dbTo, seq, sequenceLookup, sortableEntries, tableEntriesNum, subMat->alphabetSize, threads);
        Debug(Debug::INFO) << "\nIndex table: sort residues ... ";
        omptl::sort(sortableEntries, sortableEntries + tableEntriesNum, compareTableEntryByKmer);
        Debug(Debug::INFO) << "Done";
        
        IndexEntryLocal * indexEntries = (IndexEntryLocal *) sortableEntries;
        Debug(Debug::INFO) << "\nIndex table: Count and rewrite entries  ... ";
        for(size_t pos = 0; pos < tableEntriesNum; pos++){
            indexTable->table[sortableEntries[pos].kmer] += 1;
            indexEntries[pos].seqId = sortableEntries[pos].seqId;
            indexEntries[pos].position_j = sortableEntries[pos].position_j;
        }
        Debug(Debug::INFO) << "Done";
        char * entries =  (char *) realloc(indexEntries, sizeof(IndexEntryLocal) * tableEntriesNum);
        indexTable->initTableByExternalData(dbSize, tableEntriesNum, entries, (size_t *)indexTable->table, sequenceLookup);
//        Indexer idx(subMat->alphabetSize, seq->getKmerSize());
#pragma omp for schedule(static)
        for(size_t pos = 0; pos < indexTable->tableSize; pos++){
            size_t size;
            IndexEntryLocal * indexEntry = indexTable->getDBSeqList<IndexEntryLocal>(pos, &size);
            if(size > 0){
                std::sort(indexEntry, indexEntry + size, IndexEntryLocal::comapreByIdAndPos);
//                for(size_t i = 0; i < size; i++){
//                    idx.printKmer(pos, seq->getKmerSize(), subMat->int2aa);
//                    std::cout << "\t" << indexEntry[i].seqId << "\t";
//                    std::cout << indexEntry[i].position_j << std::endl;
//                }
            }
        }
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
    
    // is sequence lookup needed
    bool hasSequenceLookup;
    
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
