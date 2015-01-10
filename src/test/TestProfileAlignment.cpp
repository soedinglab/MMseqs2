//  main.cpp
//  forautocompl
//
//  Created by Martin Steinegger on 26.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//
#include <iostream>
#include <smith_waterman_sse2.h>
#include "Sequence.h"
#include "Indexer.h"
#include "ExtendedSubstitutionMatrix.h"
#include "SubstitutionMatrix.h"
#include "ReducedMatrix.h"
#include "KmerGenerator.h"
#include "BaseMatrix.h"
#include "../alignment/smith_waterman_sse2.h"
int main (int argc, const char * argv[])
{

    const size_t kmer_size=6;

    SubstitutionMatrix subMat("/Users/mad/Documents/workspace/mmseqs/data/blosum62.out",2.0);
    std::cout << "Subustitution matrix:\n";
    SubstitutionMatrix::print(subMat.subMatrix,subMat.int2aa,subMat.alphabetSize);
    //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);
    std::cout << "\n";

    std::cout << "subMatrix:\n";
    //    ReducedMatrix subMat(subMat.probMatrix, 20);
    //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);
    std::cout << "\n";
    //static const char ref_seq[40] = {'C', 'A', 'G', 'C', 'C', 'T', 'T', 'T', 'C', 'T', 'G', 'A', 'C', 'C', 'C', 'G', 'G', 'A', 'A', 'A', 'T',
    //						'C', 'A', 'A', 'A', 'A', 'T', 'A', 'G', 'G', 'C', 'A', 'C', 'A', 'A', 'C', 'A', 'A', 'A', '\0'};
    //static const char read_seq[16] = {'C', 'T', 'G', 'A', 'G', 'C', 'C', 'G', 'G', 'T', 'A', 'A', 'A', 'T', 'C', '\0'};	// read sequence
    std::string profile = "#\n"
            "NULL   3762     5295    4223    4226    4397    3752    5261    3875    4122    3320    5317    4501    4723    4910    4303    4114    4286    3783    6139    4901\n"
            "HMM    A        C       D       E       F       G       H       I       K       L       M       N       P       Q       R       S       T       V       W       Y\n"
            "M->M     M->I    M->D    I->M    I->I    D->M    D->D    Neff    Neff_I  Neff_D\n"
            "       0        *       *       0       *       0       *       *       *       *\n"
            "M 1    4271     5972    5972    5387    4387    5164    6387    3328    4972    2271    2480    5650    5972    5164    5164    4972    4650    3448    6972    5387    1\n"
            "       0        *       *       *       *       *       *       1000    0       0\n"
            "\n"
            "S 2    3153     5715    4420    4317    5590    3965    5715    5087    4268    4590    6175    4268    5175    5005    4715    2077    3590    4651    7590    5853    2\n"
            "       0        *       *       *       *       *       *       1000    0       0\n"
            "\n"
            "E 3    4205     7478    3449    1671    6063    4893    5363    5478    3741    4815    6478    4671    5256    3976    4420    4205    4671    5063    8063    6063    3\n"
            "       0        *       *       *       *       *       *       1000    0       0\n"
            "\n"
            "I 4    4459     5829    5954    5829    4507    5713    6829    1836    5507    2568    4770    6244    6244    6414    5829    5326    4659    2507    7414    5507    4\n"
            "       0        *       *       *       *       *       *       1000    0       0\n"
            "\n"
            "L 5    4542     5969    6161    5721    4161    5647    6799    3123    5384    1369    4268    6268    6268    5969    5445    5384    4924    3414    7161    5445    5\n"
            "       0        *       *       *       *       *       *       1000    0       0\n"
            "\n"
            "I 6    4459     5829    5954    5829    4507    5713    6829    1836    5507    2568    4770    6244    6244    6414    5829    5326    4659    2507    7414    5507    6\n"
            "       0        *       *       *       *       *       *       1000    0       0\n"
            "\n"
            "V 7    3891     5698    5921    5506    4862    5418    6921    2599    5258    2951    4982    6046    5921    6046    5599    4982    4377    1819    7506    5599    7\n"
            "       0        *       *       *       *       *       *       1000    0       0\n"
            "\n"
//            "P 8    4174     6981    5107    4759    6244    4866    6244    5396    4659    4866    6566    5566    959     5566    5396    4566    4866    4981    8566    6566    8\n"
//            "       0        *       *       *       *       *       *       1000    0       0\n"
//            "\n"
            "//\n";

    std::cout << "Sequence (id 0):\n";
    //const char* sequence = read_seq;
    const char* sequence = profile.c_str();
    std::cout << sequence << "\n\n";
    Sequence* s = new Sequence(10000, &subMat, Sequence::HMM_PROFILE, kmer_size, true);
    s->mapSequence(0,"lala",sequence);
    s->printProfile();
    Sequence* dbSeq = new Sequence(10000, &subMat, Sequence::AMINO_ACIDS, kmer_size, true);
    //dbSeq->mapSequence(1,"lala2",ref_seq);
    const char* sequence2 = "MSEILLLIVP";
    dbSeq->mapSequence(1,"lala2",sequence2);
    SmithWaterman aligner(15000, subMat.alphabetSize);
    int8_t * tinySubMat = new int8_t[subMat.alphabetSize*subMat.alphabetSize];

    aligner.ssw_init(s, s->getAlignmentProfile(), subMat.alphabetSize, 2);
    int32_t maskLen = s->L / 2;
    int gap_open = 10;
    int gap_extend = 1;
    s_align * alignment = aligner.ssw_align(dbSeq->int_sequence, dbSeq->L, gap_open, gap_extend, 2, 0, 0, maskLen);
    if(alignment->cigar){
        std::cout << "Cigar" << std::endl;

        int32_t targetPos = alignment->dbStartPos1, queryPos = alignment->qStartPos1;
        for (int32_t c = 0; c < alignment->cigarLen; ++c) {
            char letter = SmithWaterman::cigar_int_to_op(alignment->cigar[c]);
            uint32_t length = SmithWaterman::cigar_int_to_len(alignment->cigar[c]);
            for (uint32_t i = 0; i < length; ++i){
                if (letter == 'M') {
                    fprintf(stdout,"%c",subMat.int2aa[dbSeq->int_sequence[targetPos]]);
                    if (dbSeq->int_sequence[targetPos] == s->int_sequence[queryPos]){
                        fprintf(stdout, "|");
                    }
                    else fprintf(stdout, "*");
                    fprintf(stdout,"%c",subMat.int2aa[s->int_sequence[queryPos]]);
                    ++queryPos;
                    ++targetPos;
                } else {
                    if (letter == 'I'){
                        fprintf(stdout,"%c",subMat.int2aa[s->int_sequence[queryPos]]);
                        ++queryPos;
                    } else{
                        fprintf(stdout,"  %c",subMat.int2aa[dbSeq->int_sequence[targetPos]]);
                        ++targetPos;
                    }
                }
                std::cout << std::endl;
            }
        }
    }
    std::cout << alignment->score1  << " "<< alignment->qEndPos1 << " " << alignment->qStartPos1  << " "<< alignment->qEndPos1 << " "
            << alignment->dbStartPos1 << " "<< alignment->dbEndPos1 << std::endl;
    delete [] tinySubMat;
    delete [] alignment->cigar;
    delete alignment;
    delete s;
    delete dbSeq;
    return 0;
}

