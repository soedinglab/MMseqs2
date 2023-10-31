//  main.cpp
//  forautocompl
//
//  Created by Martin Steinegger on 26.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//
#include <iostream>
#include <cmath>
#include <ksw2.h>

//#include "gaba_wrap.h"			/* multiple target configuration: gcc example.c libgaba.a */

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <ksw2/ksw2.h>
#include <DistanceCalculator.h>
#include <SubstitutionMatrix.h>
#include <NucleotideMatrix.h>
#include <Sequence.h>
#include <BandedNucleotideAligner.h>

const char* binary_name = "test_ksw2";


/**
 * @fn random_base
 *
 * @brief return a character randomly from {'A', 'C', 'G', 'T'}.
 */
char random_base(void)
{
     char const table[4] = {'A', 'C', 'G', 'T'};
//    char const table[4] = { 0x01, 0x02, 0x04, 0x08 };
    return(table[rand() % 4]);
}

/**
 * @fn generate_random_sequence
 *
 * @brief malloc memory of size <len>, then fill it randomely with {'A', 'C', 'G', 'T'}.
 *
 * @param[in] len : the length of sequence to generate.
 * @return a pointer to the generated sequence.
 */
char *generate_random_sequence(int len)
{
    int i;
    char *seq;		/** a pointer to sequence */
    seq = (char *)malloc(sizeof(char) * (len + 32 + 1));
    if(seq == NULL) { return NULL; }
    for(i = 0; i < len; i++) {
        seq[i] = random_base();
    }
    seq[len] = '\0';
    return seq;
}

/**
 * @fn generate_mutated_sequence
 *
 * @brief take a DNA sequence represented in ASCII, mutate the sequence at the given probability.
 *
 * @param[in] seq : a pointer to the string of {'A', 'C', 'G', 'T'}.
 * @param[in] x : mismatch rate
 * @param[in] d : insertion / deletion rate
 * @param[in] bw : bandwidth (the alignment path of the input sequence and the result does not go out of the band)
 *
 * @return a pointer to mutated sequence.
 */
char *generate_mutated_sequence(char *seq, int len, double x, double d, int bw)
{
    int i, j, wave = 0;			/** wave is q-coordinate of the alignment path */
    char *mutated_seq;

    if(seq == NULL) { return NULL; }
    mutated_seq = (char *)malloc(sizeof(char) * (len + 32 + 1));
    if(mutated_seq == NULL) { return NULL; }
    for(i = 0, j = 0; i < len; i++) {
        if(((double)rand() / (double)RAND_MAX) < x) {
            mutated_seq[i] = random_base();	j++;	/** mismatch */
        } else if(((double)rand() / (double)RAND_MAX) < d) {
            if(rand() & 0x01 && wave > -bw+1) {
                mutated_seq[i] = (j < len) ? seq[j++] : random_base();
                j++; wave--;						/** deletion */
            } else if(wave < bw-2) {
                mutated_seq[i] = random_base();
                wave++;								/** insertion */
            } else {
                mutated_seq[i] = (j < len) ? seq[j++] : random_base();
            }
        } else {
            mutated_seq[i] = (j < len) ? seq[j++] : random_base();
        }
    }
    mutated_seq[len] = '\0';
    return mutated_seq;
}

/**
 * @struct params
 */
struct params {
    int64_t len;
    int64_t cnt;
    double x;
    double d;
    char **pa;
    char **pb;
};


void printSeq(char * seq, int len){
    for(int i = 0; i < len; i++){
        // char const table[4] = {'A', 'C', 'G', 'T'};
        switch(seq[i]){
            case 0x01:
                std::cout << "A";
                break;
            case 0x02:
                std::cout << "C";
                break;
            case 0x04:
                std::cout << "G";
                break;
            case 0x08:
                std::cout << "T";
                break;
        }
    }
    printf("\n");
}

struct alignment_t {
    char * cigar;
    int len;
    alignment_t(int bufferSize){
        len = 0;
        cigar = new char[bufferSize];
    }
};


//static int64_t convert_to_cigar(void *_b, int64_t len, char c) {
//    alignment_t *aln = (alignment_t*)_b;
//    for(int64_t i = 0; i < len; i++){
//        aln->cigar[aln->len + i] = c;
//    }
//    aln->len += len;
//    return len+1;
//}


void align(const char * qSeqAscii, uint8_t *qseq, uint8_t *qseqrev, const int qLen,
           const char * tSeqAscii, uint8_t *tseq, uint8_t *tseqrev, const int tLen,
           short diagonal, SubstitutionMatrix::FastMatrix &fastMatrix,
           int sc_mch, int sc_mis, int gapo, int gape)
{
    int8_t i, a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
    int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };

    unsigned short distanceToDiagonal = abs(diagonal);
    DistanceCalculator::LocalAlignment alignment;
    int qUngappedStartPos, qUngappedEndPos, dbUngappedStartPos, dbUngappedEndPos;
    if (diagonal >= 0){
        int diagonalLen = std::min(tLen, qLen - distanceToDiagonal);
        alignment = DistanceCalculator::computeSubstitutionStartEndDistance(
                qSeqAscii + distanceToDiagonal,
                tSeqAscii,
                diagonalLen, fastMatrix.matrix);
        qUngappedStartPos = alignment.startPos + distanceToDiagonal;
        qUngappedEndPos = alignment.endPos + distanceToDiagonal;
        dbUngappedStartPos = alignment.startPos;
        dbUngappedEndPos = alignment.endPos;
    }else{
        int diagonalLen = std::min(tLen - distanceToDiagonal, qLen);
        alignment = DistanceCalculator::computeSubstitutionStartEndDistance(qSeqAscii,
                                                                            tSeqAscii +
                                                                            distanceToDiagonal,
                                                                            diagonalLen,
                                                                            fastMatrix.matrix);
        qUngappedStartPos = alignment.startPos;
        qUngappedEndPos = alignment.endPos;
        dbUngappedStartPos = alignment.startPos + distanceToDiagonal;
        dbUngappedEndPos = alignment.endPos + distanceToDiagonal;
    }

//    printf("%d\t%d\t%d\t%d\t%d\n", alignment.score, qUngappedStartPos, qUngappedEndPos, dbUngappedStartPos,dbUngappedEndPos);

    // get middle position of ungapped alignment
    int qStartRev = (qLen - (qUngappedStartPos + (qUngappedEndPos - qUngappedStartPos)/2)) - 1;
    int tStartRev = (tLen - (dbUngappedStartPos + (dbUngappedEndPos - dbUngappedStartPos)/2)) - 1;


    std::string queryRev(qSeqAscii);
    std::reverse(queryRev.begin(), queryRev.end());
    std::string targetRev(tSeqAscii);
    std::reverse(targetRev.begin(), targetRev.end());

    for(int pos = qStartRev; pos < qLen; pos++){
        std::cout << queryRev[pos];
    }
    std::cout << std::endl;
    for(int pos = tStartRev; pos < tLen; pos++){
        std::cout << targetRev[pos];
    }
    std::cout << std::endl;
    ksw_extz_t ez;
    int flag = 0;
    flag |= KSW_EZ_SCORE_ONLY;
    flag |= KSW_EZ_EXTZ_ONLY;

    ksw_extz2_sse(0, qLen - qStartRev, qseqrev + qStartRev, tLen - tStartRev, tseqrev + tStartRev, 5, mat, gapo, gape, -1, -1, flag, &ez);
    printf("%s\t%s\t%d", "query", "target", ez.score);
    printf("\t%d\t%d\t%d\n", ez.max, ez.max_t, ez.max_q);
//    printf("\t%d\t%d\t%d\t%d", ez.mqe, ez.mte, ez.mqe_t, ez.mte_q);

    for(int pos = 0; pos <= ez.max_t; pos++){
        std::cout << targetRev[tStartRev+pos];
    }
    std::cout << std::endl;
    for(int pos = 0; pos <= ez.max_q; pos++){
        std::cout << queryRev[qStartRev+pos];
    }
    std::cout << std::endl;
    int qStartPos = qLen - ( qStartRev + ez.max_q ) - 1;
    int tStartPos = tLen - ( tStartRev + ez.max_t ) - 1;

    int alignFlag = 0;
    alignFlag |= KSW_EZ_EXTZ_ONLY;

    ksw_extz_t ezAlign;
    ezAlign.cigar = new uint32_t[qLen+tLen];
    ksw_extz2_sse(0, qLen-qStartPos, qseq+qStartPos, tLen-tStartPos, tseq+tStartPos, 5,
                  mat, gapo, gape, -1, -1, alignFlag, &ezAlign);
    printf("%s\t%s\t%d", "query", "target", ezAlign.score);
    printf("\t%d\t%d\t%d\t%d\t%d\n", ezAlign.max, tStartPos, ezAlign.max_t, qStartPos, ezAlign.max_q);

    for(int pos = tStartPos; pos <= tStartPos+ ezAlign.max_t; pos++){
        std::cout << tSeqAscii[pos];
    }
    std::cout << std::endl;
    for(int pos = qStartPos; pos <= qStartPos+ ezAlign.max_q; pos++){
        std::cout << qSeqAscii[pos];
    }
    std::cout << std::endl;

    for (i = 0; i < ezAlign.n_cigar; ++i) // print CIGAR
        printf("%d%c", ezAlign.cigar[i]>>4, "MID"[ezAlign.cigar[i]&0xf]);
    std::cout << std::endl;


    std::string queryAln;
    std::string targetAln;
    std::string middleAln;
    int queryPos = qStartPos;
    int targetPos = tStartPos;
    // int aaIds = 0;
    for(int i = 0; i < ezAlign.n_cigar; i++){
        int len = ezAlign.cigar[i]>>4;
        switch("MID"[ezAlign.cigar[i]&0xf]){
            case 'M': {
                for(int pos = 0; pos < len; pos++){
                queryAln.push_back(qSeqAscii[queryPos]);
                if (qSeqAscii[queryPos] == tSeqAscii[targetPos]) {
                    middleAln.push_back('|');
                    // aaIds++;
                } else {
                    middleAln.push_back('*');
                }
                targetAln.push_back(tSeqAscii[targetPos]);
                ++queryPos;
                ++targetPos;
                }
            }
                break;
            case 'I': {
                for(int pos = 0; pos < len; pos++) {
                    targetAln.push_back('-');
                    queryAln.push_back(qSeqAscii[queryPos]);
                    middleAln.push_back(' ');
                    ++queryPos;
                }
                break;
            }
            case 'D':{
                for(int pos = 0; pos < len; pos++) {
                    queryAln.push_back('-');
                    middleAln.push_back(' ');
                    targetAln.push_back(tSeqAscii[targetPos]);
                    ++targetPos;
                }
                break;
            }
        }
    }
            std::cout << queryAln << std::endl;
        std::cout << middleAln << std::endl;
        std::cout << targetAln << std::endl;
//        std::cout << static_cast<float>(aaIds)/ static_cast<float>(alignment.len) << std::endl;

    if(ezAlign.cigar) {
        delete [] ezAlign.cigar;
    }
}


int main (int, const char**) {
    int64_t i;
    struct params p;
    /** set defaults */
    p.len = 1000;
    p.cnt = 1000;
    p.x = 0.1;
    p.d = 0.1;
    p.pa = p.pb = NULL;


//    fprintf(stderr, "len\t%" PRId64 "\ncnt\t%" PRId64 "\nx\t%f\nd\t%f\n", p.len, p.cnt, p.x, p.d);

    /** init sequence */
//    queryStr = generate_random_sequence(p.len);
//    targetStr = generate_mutated_sequence(queryStr, p.len, p.x, p.d, 8);
//    std::string query  = "GCTCCGGAAGTCACAGTTTCAATCCCAAAACTGATCGATGCTCTCTCCATGC";
//    std::string target = "AAAAATCCGGAACAGTTTCAATCCCACTGATCGATGCTCTCTACACCATGC";
        std::string query = generate_random_sequence(p.len);
//        std::string target = generate_mutated_sequence((char*)query.c_str(), (int) query.size(), p.x, p.d, 8);
    short diagonal = 0;

//    std::string query  =    "GCTCCGGAAGTCACAGTTTCAATCCCAAAACTGATCGATGCTCTCTCCATGC";
//    std::string target =     "AAAAATCCGGAACAGTTTCAATCCCACTGATCGATGCTCTCTACACCATGCAAAAAA";
//    short diagonal = 15-14;

    NucleotideMatrix subMat("blosum62.out", 2.0, -0.0f);
    BandedNucleotideAligner aligner((BaseMatrix*)&subMat, 10000,  5, 1, 40);
    EvalueComputation evalueComputation(100000, &subMat, 7, 1);
    
    
//    SubstitutionMatrix::FastMatrix fastMatrix = SubstitutionMatrix::createAsciiSubMat(subMat);
    const size_t kmer_size=6;
    Sequence* queryObj = new Sequence(10000, 0, &subMat, kmer_size, true, false);
    queryObj->mapSequence(0, 0, query.c_str(), query.size());
    aligner.initQuery(queryObj);

    Sequence* targetObj = new Sequence(10000, 0, &subMat, kmer_size, true, false);
    for(i = 0; i < p.cnt; i++) {
        std::string target = generate_mutated_sequence((char*)query.c_str(), (int) query.size(), p.x, p.d, 8);
        targetObj->mapSequence(1, 1, target.c_str(), target.size());
        std::string backtrace;
        s_align alignment = aligner.align(targetObj,diagonal, false, backtrace, &evalueComputation);
        std::string queryAln;
        std::string targetAln;
        std::string middleAln;
//        int aaIds = 0;
//        if(alignment.cigar){
//            int32_t targetPos = alignment.dbStartPos1, queryPos = alignment.qStartPos1;
//            for (int32_t c = 0; c < alignment.cigarLen; ++c) {
//                char letter = SmithWaterman::cigar_int_to_op(alignment.cigar[c]);
//                uint32_t length = SmithWaterman::cigar_int_to_len(alignment.cigar[c]);
//                for (uint32_t i = 0; i < length; ++i){
//                    if (letter == 'M') {
//                        queryAln.push_back(subMat.num2aa[queryObj->sequence[queryPos]]);
//                        if (targetObj->sequence[targetPos] == queryObj->sequence[queryPos]){
//                            middleAln.push_back('|');
//                            aaIds++;
//                        }
//                        else {
//                            middleAln.push_back('*');
//                        }
//                        targetAln.push_back(subMat.num2aa[targetObj->sequence[targetPos]]);
//                        ++queryPos;
//                        ++targetPos;
//                    } else {
//                        if (letter == 'I'){
//                            queryAln.push_back(subMat.num2aa[queryObj->sequence[queryPos]]);
//                            middleAln.push_back(' ');
//                            targetAln.push_back('-');
//                            ++queryPos;
//                        }else{
//                            queryAln.push_back('-');
//                            middleAln.push_back(' ');
//                            targetAln.push_back(subMat.num2aa[targetObj->sequence[targetPos]]);
//                            ++targetPos;
//                        };
//                    }
//                }
//            }
//        }
//        std::cout << queryAln << std::endl;
//        std::cout << middleAln << std::endl;
//        std::cout << targetAln << std::endl;
        std::cout <<  alignment.score1 << " " << alignment.qStartPos1  << "-"<< alignment.qEndPos1 << " "
        << alignment.dbStartPos1 << "-"<< alignment.dbEndPos1 << std::endl;
    }
    

//    std::string query  =    "CCGCTCCGGAAGTCACAGTTTCAATCCCAAAACTGATCGATGCTCTCTCCATGC";
//    std::string target = "AAAAAAAGCTCCGGAAGTCACAGTTTCAATCCCAAAACTGATCGATGCTCTCTCCATGCAAAAAAAA";


//    std::string queryRev(query);
//    std::reverse(queryRev.begin(), queryRev.end());
//    std::string targetRev(target);
//    std::reverse(targetRev.begin(), targetRev.end());
//
//
//
//    for(i = 0; i < p.cnt; i++) {
////        flag |= KSW_EZ_APPROX_DROP;
////        flag |= KSW_EZ_APPROX_MAX;
////        std::cout << queryRev.c_str()+7 << std::endl;
////        std::cout << targetRev.c_str()+16 << std::endl;
////        align( targetRev.c_str()+16, queryRev.c_str()+7, flag, 2, -3, 5, 1);
//
//        short diagonal = 15-14;
//        int tLen = target.size();
//        int qLen = query.size();
//
//        uint8_t c[256];
//        memset(c, 4, 256);
//        c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
//        c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
//        uint8_t * ts = (uint8_t*)malloc(tLen);
//        uint8_t * qs = (uint8_t*)malloc(qLen);
//        uint8_t * tsr = (uint8_t*)malloc(tLen);
//        uint8_t * qsr = (uint8_t*)malloc(qLen);
//        for (i = 0; i < tLen; ++i) {
//            ts[i] = c[(uint8_t)target[i]]; // encode to 0/1/2/3
//            tsr[i] = c[(uint8_t)targetRev[i]]; // encode to 0/1/2/3
//        }
//        for (i = 0; i < qLen; ++i) {
//            qs[i] = c[(uint8_t)query[i]];
//            qsr[i] = c[(uint8_t)queryRev[i]];
//        }
//
//        align(query.c_str(), qs, qsr, qLen, target.c_str(), ts, tsr, tLen,
//              diagonal, fastMatrix, 2, -3, 5, 1);
//        free(ts);
//        free(qs);
//        free(tsr);
//        free(qsr);
//
////        align( targetRev.c_str()+14, queryRev.c_str()+7, flag, 2, -3, 5, 1);
//
//
////        std::string queryAln;
////        std::string targetAln;
////        std::string middleAln;
////        int queryPos = queryStart;
////        int targetPos = targetStart;
////        int aaIds = 0;
////        for(int i = 0; i < alignment.len; i++){
////            if(alignment.cigar[i]=='M'){
////                queryAln.push_back(query[queryPos]);
////                if (query[queryPos] == target[targetPos]){
////                    middleAln.push_back('|');
////                    aaIds++;
////                }
////                else {
////                    middleAln.push_back('*');
////                }
////                targetAln.push_back(target[targetPos]);
////                ++queryPos;
////                ++targetPos;
////            }else if(alignment.cigar[i]=='I'){
////                targetAln.push_back('-');
////                queryAln.push_back(query[queryPos]);
////                middleAln.push_back(' ');
////                ++queryPos;
////            }else if(alignment.cigar[i]=='D'){
////                queryAln.push_back('-');
////                middleAln.push_back(' ');
////                targetAln.push_back(target[targetPos]);
////                ++targetPos;
////            }
////        }
////        std::cout << queryAln << std::endl;
////        std::cout << middleAln << std::endl;
////        std::cout << targetAln << std::endl;
////        std::cout << static_cast<float>(aaIds)/ static_cast<float>(alignment.len) << std::endl;
//
//
//        //        std::cout << r->rapos << std::endl;
//        //        gaba_dp_dump_cigar_reverse(c, p.len, r->path->array, 0, r->path->len);
//
//        //        std::cout << c << std::endl;
//        //11M1I11M1D58M1I3M
//        //3M1I58M1D11M1I11M
//    }

    delete queryObj;
    delete targetObj;
    return 0;
}
