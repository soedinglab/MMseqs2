//
// Written by Maria Hauser mhauser@genzentrum.lmu.de
//
// Test class for k-mer generation and index table testing.
//


#include <cstdio>
#include <iostream>

#include "SubstitutionMatrix.h"
#include "IndexTable.h"
#include "QueryScore.h"

int main (int argc, const char * argv[])
{
    int kmerSize = 6;

    SubstitutionMatrix* sm = new SubstitutionMatrix("/cluster/user/maria/kClust2/data/blosum62.out");

    std::cout << "Sequence (id 0):\n";
    char* sequence = "MIPAEAGRPSLADS";
    std::cout << sequence << "\n\n";
    Sequence* s = new Sequence (10000, sm->aa2int, sm->int2aa);
    s->mapSequence(sequence);
    std::cout << "Int sequence:\n";
    for (int i = 0; i < s->L; i++)
        std::cout << s->int_sequence[i] << " ";
    std::cout << "\n";

    ///////////////////////////////////////////
    // k-mer generation test
    ///////////////////////////////////////////

    Indexer* idxer = new Indexer(20, kmerSize);

    int* kmer;
    unsigned int kmerIdx;

    int* testKmer = new int[kmerSize];
    std::cout << "Pos:\tkmer idx:\tint k-mer:\tchar k-mer:\n";
    
    for (int pos = 0; pos < (s->L-kmerSize); pos++){
        kmer = s->int_sequence + pos;
        kmerIdx = idxer->getNextKmerIndex(kmer, kmerSize);
        
        std::cout << pos << "\t" << kmerIdx << "\t";
        // test: reverse conversion index -> int k-mer -> char k-mer
        idxer->index2int(testKmer, kmerIdx, kmerSize);
        std::cout << "\t";
        for (int i = 0; i < kmerSize; i++)
            std::cout << testKmer[i] << " ";
        std::cout << "\t";
        for (int i = 0; i < kmerSize; i++)
            std::cout << sm->int2aa[testKmer[i]];
        std::cout << "\n";
    }

    /////////////////////////////////////////////
    // Index table test
    ////////////////////////////////////////////

    std::cout << "\nTesting index table!\n";
    std::cout << "Initial allocation...\n";
    IndexTable* it = new IndexTable(20, kmerSize);
    std::cout << " done.\n";
    it->addSequence(s);
    // adding next sequence
/*    std::cout << "Adding another sequence (id 1):\n";
    sequence = "MSSAEAGRPSLADS";
    std::cout << sequence << "\n";
    s->setId(1);
    s->mapSequence(sequence);
    it->addSequence(s);
*/
    std::cout << "INIT Index Table...\n";
    it->init();
    std::cout << " done.";
    exit(0);

    std::cout << "\nSequence list for k-mer index " << kmerIdx << ":\n";
    int listSize = 0;
    int* seqList = it->getDBSeqList(kmerIdx, &listSize);
    if (listSize > 0){
        for (int i = 0; i < listSize-1; i++)
            std::cout << seqList[i] << ",";
        std::cout << seqList[listSize-1] << "\n";
    }

    std::cout << "\nSequence list for k-mer index 0:\n";
    seqList = it->getDBSeqList(0, &listSize);
    if (listSize > 0){
        for (int i = 0; i < listSize-1; i++)
            std::cout << seqList[i] << ",";
        std::cout << seqList[listSize-1] << "\n";
    }
    std::cout << "\n";

    //////////////////////////////////////////////////////
    // Query Score test
    /////////////////////////////////////////////////////

    std::cout << "Testing QueryScore! (each exact k-mer match has the score 1)\n";
    QueryScore* qs = new QueryScore(5, 3);
    
    // Simulation of the index table match step
    for (int pos = 0; pos < (s->L-kmerSize); pos++){
        kmer = s->int_sequence + pos;
        kmerIdx = idxer->getNextKmerIndex(kmer, kmerSize);
        
        listSize = 0;
        seqList = it->getDBSeqList(kmerIdx, &listSize);

        qs->addScores(seqList, listSize, 1);
        qs->addScores(seqList, listSize, 2);
    }

    // get the result from the QueryScore
    std::cout << "Prefiltering results for the sequence with id 1:\n";
    std::cout << "seqID\tscore\n";
    std::list<hit_t>* res = qs->getResult(12);
    std::list<hit_t>::iterator iter;
    for (iter = res->begin(); iter != res->end(); iter++){
        std::cout << iter->seqId << "\t" << iter->prefScore << "\n";
    }

    return 0;
}

