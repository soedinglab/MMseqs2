#include <iostream>
#include <list>
#include <algorithm>
#include <math.h>
#include <SequenceLookup.h>
#include <SubstitutionMatrix.h>
#include <DiagonalMatcher.h>
#include <QueryScore.h>
#include <ExtendedSubstitutionMatrix.h>


#include "Clustering.h"
#include "SetElement.h"

#include "DBReader.h"
#include "DBWriter.h"


int main(int argc, char **argv)
{

    size_t kmer_size = 6;
    SubstitutionMatrix subMat("/Users/mad/Documents/workspace/mmseqs/data/blosum62.out",
                              8.0, -0.2);
    SubstitutionMatrix::print(subMat.subMatrix,subMat.int2aa,subMat.alphabetSize);

    std::string S1 = "PQITLWQG";
    const char* S1char = S1.c_str();
    std::cout << S1char << "\n\n";
    Sequence s1(10000, subMat.aa2int, subMat.int2aa, 0, kmer_size, true);
    s1.mapSequence(0,0,S1char);
    //                0123456789
    std::string S2 = "XXXXXXXXXPQITLWQG";
    const char* S2char = S2.c_str();
    std::cout << S2char << "\n\n";
    Sequence s2(10000, subMat.aa2int, subMat.int2aa, 0, kmer_size, true);
    s2.mapSequence(1,1,S2char);
    std::string S3 = "PQITLWQGXXXXXXXXX";
    const char* S3char = S3.c_str();
    std::cout << S3char << "\n\n";
    Sequence s3(10000, subMat.aa2int, subMat.int2aa, 0, kmer_size, true);
    s3.mapSequence(2,2, S3char);
    std::string S4 = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
            "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
            "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
            "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXPQITLWQG";
    const char* S4char = S4.c_str();
    std::cout << S4char << "\n\n";
    Sequence s4(10000, subMat.aa2int, subMat.int2aa, 0, kmer_size, true);
    s4.mapSequence(3,3, S4char);



    std::string S5 = "QDELTAGPCATVHVITVQMAKSGELQAIAPEVAQSLAEFFAVLADPNRLRLLSLLARSELCVGDLAQAIGVSESAVSHQLRSLRNLRLVSYRKQGRHVYYQLQDHHIVALYQNALDHLQECR";
    const char* S5char = S5.c_str();
    std::cout << S5char << "\n\n";
    Sequence s5(10000, subMat.aa2int, subMat.int2aa, 0, kmer_size, true);
    s5.mapSequence(4,4, S5char);
//    s5.reverse();

    std::string S6 = "DARAAADLAPVAVETTAVLSTPSLLDMEATLKRAETAGFHAGLDITRLGKDFEADLTRGFAEPNHEVALKWTPVPAFVPAVRSFPAAAAAAAAAATAVATTATGTATEPAADTSTTNAASATSAAGADSAEGKADAADPVAKGRYEDCRKNGFGEIGDFGPSQHKSCISEGCGHLPCMVCGGDGKMHESNKRCNKCAKPSGCQVREYPPGRMCRVFVYRDALLPSDGGRLGNQPACRHGCAYPQAVHRPNSCLSPEQGSLRFSCYKCGQEGHIAKDCCPKTCHVGSNAHPYFASPDGDKDCRYRPSRDRVRDDARSSAGRYPGDSTGGGGYGYRGRYDYYDGDDGRRGSDSYPYPSRYGPGGGSYGYGSGASSSGGGGRSGGYSASSYGASGSYGGGYTYVSAGSGSSYSASGPASPYGSSGSGQSGSGGSYGYSGGGGYYQASSSGSSQASYSGSSPASYGPGSSYGGAGSQGGGYAAAGTAYAYASGASGSPPQGGNYTGRYTSRPQS";
    const char* S6char = S6.c_str();
    std::cout << S6char << "\n\n";
    Sequence s6(10000, subMat.aa2int, subMat.int2aa, 0, kmer_size, true);
    s6.mapSequence(5,5, S6char);


    std::string S7 = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
            "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
            "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
            "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
            "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
            "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
            "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
            "XXXXXXXXXXXXXXXXXXXXXXPQITLWQG";
    const char* S7char = S7.c_str();
    std::cout << S7char << "\n\n";
    Sequence s7(10000, subMat.aa2int, subMat.int2aa, 0, kmer_size, true);
    s7.mapSequence(6,6, S7char);




    SequenceLookup lookup(6, s1.L + s2.L + s3.L + s4.L + s5.L + s6.L + s7.L );
    lookup.addSequence(&s1);
    lookup.addSequence(&s2);
    lookup.addSequence(&s3);
    lookup.addSequence(&s4);
    lookup.addSequence(&s5);
    lookup.addSequence(&s6);
    lookup.addSequence(&s7);



    float * compositionBias = new float[1000];
    hit_t hits[32];
    DiagonalMatcher matcher(1000, &subMat, &lookup);

    SubstitutionMatrix::calcLocalAaBiasCorrection(&subMat, s5.int_sequence, s5.L, compositionBias);
    memset(compositionBias, 0.0, sizeof(float)*s5.L);
//    std::cout << compositionBias[74] << std::endl;
//    std::cout << compositionBias[79] << std::endl;
//    std::cout << compositionBias[80] << std::endl;
//    std::cout << compositionBias[84] << std::endl;
//    std::cout << compositionBias[90] << std::endl;

//    for(size_t i = 0; i < s5.L; i++){
//        std::cout << compositionBias[i] << " ";
//    }
    std::cout << std::endl;
    for(size_t i = 0; i < 16; i++){
        hits[i].seqId = s6.getId();
        hits[i].diagonal = 65296;
        hits[i].diagonalScore = 0;
    }
    matcher.processQuery(&s5, compositionBias, std::make_pair(hits, 16));
    std::cout << hits[0].diagonal << " " <<  (int)hits[0].diagonalScore << std::endl;



//    SubstitutionMatrix::calcLocalAaBiasCorrection(&subMat, s1.int_sequence, s1.L, compositionBias);
//
//    hits[0].seqId = s1.getId();
//    hits[0].diagonal = 0;
//    hits[0].diagonalScore = 0;
//    matcher.processQuery(&s1, compositionBias, std::make_pair(hits, 1));
//    std::cout << ExtendedSubstitutionMatrix::calcScore(s1.int_sequence, s1.int_sequence,s1.L, subMat.subMatrix) << " " << (int)hits[0].diagonalScore <<  std::endl;
//
//    for(int i = 0; i < 16; i++){
//        hits[i].seqId = s1.getId();
//        hits[i].diagonal = 0;
//    }
//    matcher.processQuery(&s1, compositionBias, std::make_pair(hits, 16));
//    std::cout << ExtendedSubstitutionMatrix::calcScore(s1.int_sequence, s1.int_sequence,s1.L, subMat.subMatrix) << " " << (int)hits[0].diagonalScore <<  std::endl;
//
//
//    hits[0].seqId = s1.getId();
//    hits[0].diagonal = 9;
//    hits[0].diagonalScore = 0;
//    matcher.processQuery(&s2, compositionBias, std::make_pair(hits, 1));
//    std::cout << ExtendedSubstitutionMatrix::calcScore(s1.int_sequence, s1.int_sequence,s1.L, subMat.subMatrix) << " " << (int)hits[0].diagonalScore <<  std::endl;
//
//    for(int i = 0; i < 16; i++){
//        hits[i].seqId = s1.getId();
//        hits[i].diagonal = 9;
//    }
//    hits[0].diagonalScore = 0;
//    matcher.processQuery(&s2, compositionBias, std::make_pair(hits, 16));
//    std::cout << ExtendedSubstitutionMatrix::calcScore(s1.int_sequence, s1.int_sequence,s1.L, subMat.subMatrix) << " " << (int)hits[0].diagonalScore <<  std::endl;
//
//    for(int i = 0; i < 16; i++){
//        hits[i].seqId = s2.getId();
//        hits[i].diagonal = 9;
//    }
//    hits[0].diagonalScore = 0;
//    matcher.processQuery(&s1, compositionBias, std::make_pair(hits, 16));
//    std::cout << ExtendedSubstitutionMatrix::calcScore(s1.int_sequence, s1.int_sequence,s1.L, subMat.subMatrix) << " " << (int)hits[0].diagonalScore <<  std::endl;
//
//    hits[0].diagonalScore = 0;
//    matcher.processQuery(&s1, compositionBias, std::make_pair(hits, 1));
//    std::cout << ExtendedSubstitutionMatrix::calcScore(s1.int_sequence, s1.int_sequence,s1.L, subMat.subMatrix) << " " << (int)hits[0].diagonalScore <<  std::endl;
//
//    for(int i = 0; i < 16; i++){
//        hits[i].seqId = s2.getId();
//        hits[i].diagonal = 9;
//    }
//    hits[0].diagonalScore = 0;
//    matcher.processQuery(&s3, compositionBias, std::make_pair(hits, 16));
//    std::cout << ExtendedSubstitutionMatrix::calcScore(s1.int_sequence, s1.int_sequence,s1.L, subMat.subMatrix) << " " << (int)hits[0].diagonalScore <<  std::endl;
//
//    hits[0].diagonalScore = 0;
//    matcher.processQuery(&s3, compositionBias, std::make_pair(hits, 1));
//    std::cout << ExtendedSubstitutionMatrix::calcScore(s1.int_sequence, s1.int_sequence,s1.L, subMat.subMatrix) << " " << (int)hits[0].diagonalScore <<  std::endl;
//
//
//    hits[0].seqId = s4.getId();
//    hits[0].diagonalScore = 0;
//    hits[0].diagonal = 256;
//    matcher.processQuery(&s1, compositionBias, std::make_pair(hits, 1));
//    std::cout << ExtendedSubstitutionMatrix::calcScore(s1.int_sequence, s1.int_sequence,s1.L, subMat.subMatrix) << " " << (int)hits[0].diagonalScore <<  std::endl;
//
//
//    hits[0].seqId = s1.getId();
//    hits[0].diagonalScore = 0;
//    hits[0].diagonal = 256;
//    matcher.processQuery(&s4,compositionBias, std::make_pair(hits, 1));
//    std::cout << ExtendedSubstitutionMatrix::calcScore(s1.int_sequence, s1.int_sequence,s1.L, subMat.subMatrix) << " " << (int)hits[0].diagonalScore <<  std::endl;
//
//    hits[0].seqId = s7.getId();
//    hits[0].diagonalScore = 0;
//    hits[0].diagonal = 512;
//    matcher.processQuery(&s1,compositionBias, std::make_pair(hits, 16));
//    std::cout << ExtendedSubstitutionMatrix::calcScore(s1.int_sequence, s1.int_sequence,s1.L, subMat.subMatrix) << " " << (int)hits[0].diagonalScore <<  std::endl;
//
//    hits[0].seqId = s1.getId();
//    hits[0].diagonalScore = 0;
//    hits[0].diagonal = 512;
//    matcher.processQuery(&s7,compositionBias, std::make_pair(hits, 16));
//    std::cout << ExtendedSubstitutionMatrix::calcScore(s1.int_sequence, s1.int_sequence,s1.L, subMat.subMatrix) << " " << (int)hits[0].diagonalScore <<  std::endl;
//
//
//    hits[0].seqId = s7.getId();
//    hits[0].diagonalScore = 0;
//    hits[0].diagonal = 0;
//    matcher.processQuery(&s7, compositionBias, std::make_pair(hits, 16));
//    std::cout << ExtendedSubstitutionMatrix::calcScore(s1.int_sequence, s1.int_sequence,s1.L, subMat.subMatrix) << " " << (int)hits[0].diagonalScore <<  std::endl;

    delete [] compositionBias;
}
