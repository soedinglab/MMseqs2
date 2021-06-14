//  main.cpp
//  forautocompl
//
//  Created by Martin Steinegger on 26.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//
#include <iostream>
#include <cmath>

#include "StripedSmithWaterman.h"
#include "Sequence.h"
#include "Indexer.h"
#include "ExtendedSubstitutionMatrix.h"
#include "SubstitutionMatrix.h"
#include "ReducedMatrix.h"
#include "KmerGenerator.h"
#include "BaseMatrix.h"
#include "StripedSmithWaterman.h"
#include "Parameters.h"
#include "Matcher.h"

const char* binary_name = "test_alignment";

int main (int, const char**) {
    const size_t kmer_size=6;

    char buffer[1024];
    const Matcher::result_t result(1351, 232, 1.0, 1.0, 0.99, 0.000000001, 20,
                                   3, 15, 22, 4, 18, 354, "MMMMMIIMMMMDDMMMMMM");
    size_t len = Matcher::resultToBuffer(buffer, result, true, false);
    std::cout << std::string(buffer, len) << std::endl;

    Parameters& par = Parameters::getInstance();
    par.initMatrices();
    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 2.0, -0.0f);
    std::cout << "Substitution matrix:\n";
    SubstitutionMatrix::print(subMat.subMatrix,subMat.num2aa,subMat.alphabetSize);
    SubstitutionMatrix::print(subMat.subMatrix,subMat.num2aa,subMat.alphabetSize);
//    for(int i = 0; i < 255; i++){
//        std::cout << i << "\t" << MathUtil::convertCharToFloat(i) << std::endl;
//    }
    //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);
    std::cout << "\n";

    std::cout << "subMatrix:\n";
    //    ReducedMatrix subMat(subMat.probMatrix, 20);
    //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);
    std::cout << "\n";
    //static const char ref_seq[40] = {'C', 'A', 'G', 'C', 'C', 'T', 'T', 'T', 'C', 'T', 'G', 'A', 'C', 'C', 'C', 'G', 'G', 'A', 'A', 'A', 'T',
    //						'C', 'A', 'A', 'A', 'A', 'T', 'A', 'G', 'G', 'C', 'A', 'C', 'A', 'A', 'C', 'A', 'A', 'A', '\0'};
    //static const char read_seq[16] = {'C', 'T', 'G', 'A', 'G', 'C', 'C', 'G', 'G', 'T', 'A', 'A', 'A', 'T', 'C', '\0'};	// read sequence
//	std::string tim = "APRKFFVGGNWKMNGKRKSLGELIHTLDGAKLSADTEVVCGAPSIYLDFARQKLDAKIGVAAQNCYKVPKGAFTGEISPAMIKDIGAAWVILGH"
//                      "SERRHVFGESDELIGQKVAHALAEGLGVIACIGEKLDEREAGITEKVVFQETKAIADNVKDWSKVVLAYEPVWAIGTGKTATPQQAQEVHEKLR"
//			          "GWLKTHVSDAVAVQSRIIYGGSVTGGNCKELASQHDVDGFLVGGASLKPEFVDIINAKH";
//    std::string tim1 = "LAEVGDARSLLEDDLVDLPDARFFKAMGREFVKLMLQGEASEAIKAPRAAAAVLPKQYTRDEDGDGVNLVLLVERVLEVPDECRLYIIGVAARVAGATVVYATGSRKKDAALPIANDETHLTAVLAKGESLPPPPENPMSADRVRWEHIQRIYEMCDRNVSETARRLNMHRRTLQRILAKRSPR";
//    std::string tim2 = "EMDLAFVELGADRSLLLVDDDEPFLKRLAKAMEKRGFVLETAQSVAEGKAIAQARPPAYAVVDLRLEDGNGLDVVEVLRERRPDCRIVVLTGYGAIATAVAAVKIGATDYLSKPADANEVTHALLAKGESLPPPPENPMSADRVRWEHIQRIYEMCDRNVSETARRLNMHRRTLQRILAKRSPR";
    std::string tim1 = "'AAAGGTGACCGGGCACGGTGGCCCATGCCTATAATCCCAGCACTTTGGGAGGCCCAGGCAGGTGGATCACTTGAGGTCAGGAGTTCGAGACCAGCCTGGC";
    std::string tim2 = "'GATTGAAAAACTCCCAGGCTGGACACGGTGGCCCATGCCTGTAATCCCAGCACTCTGGGAGGCTGAGGTGGGCTGATCCCTTGAGGTCAGGAGTTCGAGACCATCCTGGAAAATGTGGCA";

    std::cout << "Sequence (id 0):\n";
    //const char* sequence = read_seq;
    std::cout << tim1 << "\n\n";

    Sequence* s = new Sequence(10000, 0, &subMat, kmer_size, true, false);
    s->mapSequence(0,0,tim1.c_str(), tim1.size());
    Sequence* dbSeq = new Sequence(10000, 0, &subMat, kmer_size, true, false);
    //dbSeq->mapSequence(1,"lala2",ref_seq);
    dbSeq->mapSequence(1,1,tim2.c_str(), tim2.size());
    SmithWaterman aligner(15000, subMat.alphabetSize, true);
    int8_t * tinySubMat = new int8_t[subMat.alphabetSize*subMat.alphabetSize];
    for (int i = 0; i < subMat.alphabetSize; i++) {
        for (int j = 0; j < subMat.alphabetSize; j++) {
            std::cout << ( i*subMat.alphabetSize + j) << " " << subMat.subMatrix[i][j] << " ";

            tinySubMat[i*subMat.alphabetSize + j] = (int8_t)subMat.subMatrix[i][j];
        }
        std::cout << std::endl;
    }
    float sum=0.0;
    for(int i = 0; i < subMat.alphabetSize; i++){
        sum += subMat.subMatrix[i][i];
    }
    std::cout << "Test: " << sum/ subMat.alphabetSize << std::endl;
    aligner.ssw_init(s, tinySubMat, &subMat, 2);
    int32_t maskLen = s->L / 2;
    int gap_open = 11;
    int gap_extend = 1;
    float seqId = 1.0;
    int aaIds = 0;
    EvalueComputation evalueComputation(100000, &subMat, gap_open, gap_extend);
    s_align alignment = aligner.ssw_align(dbSeq->numSequence, dbSeq->L, gap_open, gap_extend, 2, 10000, &evalueComputation, 0, 0.0, maskLen);
    if(alignment.cigar){
        std::cout << "Cigar" << std::endl;

        int32_t targetPos = alignment.dbStartPos1, queryPos = alignment.qStartPos1;
        for (int32_t c = 0; c < alignment.cigarLen; ++c) {
            char letter = SmithWaterman::cigar_int_to_op(alignment.cigar[c]);
            uint32_t length = SmithWaterman::cigar_int_to_len(alignment.cigar[c]);
            for (uint32_t i = 0; i < length; ++i){
                if (letter == 'M') {
                    fprintf(stdout,"%c",subMat.num2aa[s->numSequence[queryPos]]);

                    if (dbSeq->numSequence[targetPos] == s->numSequence[queryPos]){
                        fprintf(stdout, "|");
                        aaIds++;
                    }
                    else fprintf(stdout, "*");
                    fprintf(stdout,"%c",subMat.num2aa[dbSeq->numSequence[targetPos]]);

                    ++queryPos;
                    ++targetPos;
                } else {
                    if (letter == 'I'){
                        fprintf(stdout,"%c",subMat.num2aa[s->numSequence[queryPos]]);
                        fprintf(stdout, " ");
                        fprintf(stdout, "|");
                        ++queryPos;
                    }else{
                        fprintf(stdout, "|");
                        fprintf(stdout, " ");
                        fprintf(stdout,"%c",subMat.num2aa[dbSeq->numSequence[targetPos]]);
                        ++targetPos;
                    };
                }
                std::cout << std::endl;
            }
        }
    }
    std::cout <<  alignment.score1 << " " << alignment.qStartPos1  << "-"<< alignment.qEndPos1 << " "
    << alignment.dbStartPos1 << "-"<< alignment.dbEndPos1 << std::endl;
    seqId = (float)aaIds/(float)(std::min(s->L, dbSeq->L)); //TODO
    int32_t qAlnLen = std::max(alignment.qEndPos1 - alignment.qStartPos1, static_cast<int32_t>(1));
    int32_t dbAlnLen = std::max(alignment.dbEndPos1 - alignment.dbStartPos1, static_cast<int32_t>(1));
    if(alignment.cigar){
        seqId =  static_cast<float>(aaIds) / static_cast<float>(std::max(qAlnLen, dbAlnLen));
    }else{
        seqId = Matcher::estimateSeqIdByScorePerCol(alignment.score1, qAlnLen, dbAlnLen);
    }

    std::cout << "Seqid: "<< seqId << " aaIds " << aaIds <<std::endl;
    double evalue =  pow (2,-(double)alignment.score1/2);


    //score* 1/lambda
    //
    std::cout << evalue << std::endl;
    std::cout << (s->L) << std::endl;
    std::cout << (dbSeq->L) << std::endl;
    // calcuate stop score
    const double qL = static_cast<double>(s->L);
    const double dbL = static_cast<double>(dbSeq->L);
    size_t dbSize = (qL * 11000000 * dbL);
    dbSize = 1.74e+10;
    std::cout << dbSize << std::endl;
    std::cout << evalue * dbSize << std::endl;

    double lambda= 0.267;
//    double K= 0.041;
//    double Kmn=(qL * seqDbSize * dbSeq->L);
    double Kmn=1.74e+12;
    std::cout << dbSize/Kmn<< " " <<  Kmn * exp(-(alignment.score1 * lambda)) << std::endl;
    delete [] tinySubMat;
    delete [] alignment.cigar;
    delete s;
    delete dbSeq;
    return 0;
}

