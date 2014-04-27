//
// Written by Maria Hauser mhauser@genzentrum.lmu.de
//
// Test class for k-mer generation and index table testing.
//


#include <cstdio>
#include <iostream>

#include "../commons/SubstitutionMatrix.h"
#include "../prefiltering/IndexTable.h"
#include "../prefiltering/QueryScore.h"
#include "../prefiltering/QueryScoreGlobal.h"

int main (int argc, const char * argv[])
{
    int kmerSize = 7;
    int alphabetSize = 21;

    SubstitutionMatrix* sm = new SubstitutionMatrix("/cluster/user/maria/mmseqs/data/blosum62.out", 8.0);

    std::cout << "Sequence (id 0):\n";
    char* sequence = "MIPAEAGRPSLADS";
    std::cout << sequence << "\n\n";
    Sequence* s = new Sequence (10000, sm->aa2int, sm->int2aa, Sequence::AMINO_ACIDS);
    s->mapSequence(0, "TEST", sequence);
    std::cout << "Int sequence:\n";
    for (int i = 0; i < s->L; i++)
        std::cout << s->int_sequence[i] << " ";
    std::cout << "\n\n";

    Sequence* s1 = new Sequence (10000, sm->aa2int, sm->int2aa, Sequence::AMINO_ACIDS);
    sequence = "MSSAEAGRPSLADS";
    s1->mapSequence(1, "TEST1", sequence);
    std::cout << "Sequence (id 1):\n";
    std::cout << sequence << "\n\n";

    ///////////////////////////////////////////
    // k-mer generation test
    ///////////////////////////////////////////

    Indexer* idxer = new Indexer(21, kmerSize);

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

    std::cout << "Pos:\tkmer idx:\tint k-mer:\tchar k-mer:\n";
    idxer->reset();
    for (int pos = 0; pos < (s1->L-kmerSize); pos++){
        kmer = s1->int_sequence + pos;
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
    IndexTable* it = new IndexTable(alphabetSize, kmerSize, 0);
    it->addKmerCount(s);
    it->addKmerCount(s1);
    it->init();
    it->addSequence(s);
    it->addSequence(s1);
    std::cout << " done.\n";

    for (int kmerIdx = 0; kmerIdx < pow(alphabetSize, kmerSize); kmerIdx++){
        int listSize = 0;
        int* seqList = it->getDBSeqList(kmerIdx, &listSize);
        if (listSize > 0){
            std::cout << "\nSequence list for k-mer index " << kmerIdx << " (";
            idxer->printKmer(kmerIdx, kmerSize, sm->int2aa);
            std::cout << ")\n";
            std::cout << "size: " << listSize << "\n";
            for (int i = 0; i < listSize-1; i++)
                std::cout << seqList[i] << ",";
            std::cout << seqList[listSize-1] << "\n";

        }
    }


    //////////////////////////////////////////////////////
    // Query Score test
    /////////////////////////////////////////////////////

    short unsigned int seqLens[2];
    seqLens[0] = 14;
    seqLens[1] = 14;
    std::cout << "Testing QueryScore! (each exact k-mer match has the score 1)\n";
    QueryScore* qs = new QueryScoreGlobal(2, &seqLens[0], kmerSize, 125, 8.1433e-07, 5.0);
    
    // Simulation of the index table match step
    for (int pos = 0; pos < (s->L-kmerSize); pos++){
        kmer = s->int_sequence + pos;
        kmerIdx = idxer->getNextKmerIndex(kmer, kmerSize);
        
        int listSize = 0;
        int* seqList = it->getDBSeqList(kmerIdx, &listSize);

        qs->addScores(seqList, listSize, 1);
    }

    // get the result from the QueryScore
    std::cout << "Prefiltering results for the sequence with id 1:\n";
    std::cout << "seqID\tscore\n";
    std::pair<hit_t *,size_t> res = qs->getResult(12,100000000);
    for (int i = 0; i < res.second; i++ ){
        std::cout << res.first[i].seqId << "\t" << res.first[i].prefScore << "\n";
    }

    return 0;
}

