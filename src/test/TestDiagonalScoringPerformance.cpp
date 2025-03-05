#include <iostream>
#include <list>
#include <algorithm>
#include <cmath>

#include "SequenceLookup.h"
#include "SubstitutionMatrix.h"
#include "UngappedAlignment.h"
#include "ExtendedSubstitutionMatrix.h"
#include "FileUtil.h"

#include "kseq.h"
#include <unistd.h> // read

KSEQ_INIT(int, read)

#include "Clustering.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Parameters.h"

const char* binary_name = "test_diagonalscoringperformance";
DEFAULT_PARAMETER_SINGLETON_INIT

int main (int, const char**) {
    size_t kmer_size = 6;
    Parameters& par = Parameters::getInstance();
    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 8.0, 0.0);
    SubstitutionMatrix::print(subMat.subMatrix,subMat.num2aa,subMat.alphabetSize);

    std::string S1 = "AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL";
    const char* S1char = S1.c_str();
//    std::cout << S1char << "\n\n";
    Sequence s1(10000,  0, &subMat, kmer_size, true, false);
    s1.mapSequence(0,0,S1char, S1.size());

    std::string S2 = "MLKIRYSSAFKKDLKPFQHDKSAISVINTVLKLLATGKPLPREYKEHSLKGDYIGYLECHGKPDLLLIYKRTEQEVFLYRVGSHAKLF";
    const char* S2char = S2.c_str();
//    std::cout << S2char << "\n\n";
    Sequence s2(10000,  0, &subMat, kmer_size, true, false);
    s2.mapSequence(0,0,S2char, S2.size());

    FILE *fasta_file = FileUtil::openFileOrDie("/Users/mad/Documents/databases/mmseqs_benchmark/benchmarks/clustering_benchmark/db/db_full.fas", "r", true);
    kseq_t *seq = kseq_init(fileno(fasta_file));
    size_t dbEntrySize = 0;
    size_t dbCnt = 0;
    while (kseq_read(seq) >= 0) {
        dbEntrySize += seq->seq.l;
        dbCnt += 1;
    }
    SequenceLookup lookup(dbCnt*10, dbEntrySize*10);
    //kseq_destroy(seq);
    Sequence dbSeq(40000,  0, &subMat, kmer_size, true, false);
    size_t id = 0;
    size_t maxLen = 0;
    for(size_t i = 0; i < 10; i++){
        fclose(fasta_file);
        fasta_file = FileUtil::openFileOrDie("/Users/mad/Documents/databases/mmseqs_benchmark/benchmarks/clustering_benchmark/db/db_full.fas", "r", true);
        kseq_rewind(seq);
        while (kseq_read(seq) >= 0) {
            dbSeq.mapSequence(id,id,seq->seq.s, seq->seq.l);
            maxLen = std::max(seq->seq.l, maxLen);
//        if(id == 202423){
//            std::cout << seq->seq.s << std::endl;
//        }
            lookup.addSequence(&dbSeq);

            id += 1;
        }
    }
    kseq_destroy(seq);
    std::cout << maxLen << std::endl;
    UngappedAlignment matcher(maxLen, &subMat, &lookup);
    CounterResult hits[16000];
    hits[0].id =142424;
    hits[0].diagonal = 50;
    hits[1].id = 191382;
    hits[1].diagonal = 4;
    hits[2].id = 135950;
    hits[2].diagonal = 4;
    hits[3].id = 63969;
    hits[3].diagonal = 4;
    hits[4].id = 244188;
    hits[4].diagonal = 4;

    for(size_t i = 5; i < 16; i++) {
        hits[i].id = 159147;
        hits[i].diagonal = 31;
    }

    float * compositionBias = new float[s1.L];
    SubstitutionMatrix::calcLocalAaBiasCorrection(&subMat, s1.numSequence, s1.L, compositionBias, 1.0);


    matcher.createProfile(&s1, compositionBias);
    matcher.align(hits, 16);
    std::cout << (int)hits[0].count << " ";
    std::cout << (int)hits[1].count << " ";
    std::cout << (int)hits[2].count << " ";
    std::cout << (int)hits[3].count << std::endl;

    matcher.align(hits, 1);
    matcher.align(hits + 1, 1);
    matcher.align(hits + 2, 1);
    matcher.align(hits + 3, 1);

    std::cout << (int)hits[0].count<< " ";
    std::cout << (int)hits[1].count<< " ";
    std::cout << (int)hits[2].count<< " ";
    std::cout << (int)hits[3].count<< std::endl;
    for(size_t i = 0; i < 10000; i++){
        for(int j = 1; j < 16000; j++){
            hits[j].id = rand()%dbCnt;
            hits[j].diagonal =  rand()%s1.L;
        }
        //   std::reverse(hits, hits+1000);
        matcher.align(hits, 16000);
    }
//    std::cout << ExtendedSubstitutionMatrix::calcScore(s1.sequence, s1.sequence,s1.L, subMat.subMatrix) << " " << (int)hits[0].diagonalScore <<  std::endl;
//    std::cout << (int)hits[0].diagonalScore <<  std::endl;
    for(int i = 0; i < 1000; i++){
        std::cout << hits[i].id << "\t" << (int) hits[i].diagonal  << "\t" << (int)hits[i].count <<  std::endl;
    }
    delete [] compositionBias;
}
