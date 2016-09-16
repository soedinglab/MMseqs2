//  main.cpp
//  forautocompl
//
//  Created by Martin Steinegger on 26.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//
#include <iostream>
#include <smith_waterman_sse2.h>
#include <MsaFilter.h>
#include "PSSMCalculator.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "MultipleAlignment.h"
#include "Parameters.h"
int main (int argc, const char * argv[])
{
    Parameters& par = Parameters::getInstance();

    const size_t kmer_size=6;

    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0, 0);
    std::cout << "Subustitution matrix:\n";
    SubstitutionMatrix::print(subMat.subMatrix,subMat.int2aa,subMat.alphabetSize);
    //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);
    std::cout << "\n";

    std::cout << "subMatrix:\n";
    //    ReducedMatrix subMat(subMat.probMatrix, 20);
    //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);
    std::cout << "\n";
    std::cout << "Sequence (id 0):\n";
    std::string S1 = "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF";
    const char* S1char = S1.c_str();
    std::cout << S1char << "\n\n";
    Sequence s1(10000, subMat.aa2int, subMat.int2aa, 0, kmer_size, true, true);
    s1.mapSequence(0,0,S1char);
    std::string S2 = "PQFSLWKRPVVTAYIEGQPVEVLLDTGADDSIVAGIELGNNIVGGIGGFINTLEYKNVEIEVLNKKVRATIMTGDTPINIFGRNILTALGMSLNL";
    const char* S2char = S2.c_str();
    std::cout << S2char << "\n\n";
    Sequence s2(10000, subMat.aa2int, subMat.int2aa, 0, kmer_size, true, true);
    s2.mapSequence(1,1,S2char);
    std::string S3 = "PQFHLWKRPVVTAGQPVEVLLDTGADDSIVTGIELGPHYTPKIVGGIGGFINTKEYKNVEVEVLGKRIKGTIMTGDTPINIFGRNLLTALGMSLNF";
    const char* S3char = S3.c_str();
    std::cout << S3char << "\n\n";
    Sequence s3(10000, subMat.aa2int, subMat.int2aa, 0, kmer_size, true, true);
    s3.mapSequence(2,2, S3char);
    std::string S4 = "LAMTMEHKDRPLVRVILTNTGSHPVKQRSVYITALLDTGADDTVISEEDWPTDWPVMEAANPQIHGIGGGIPVRKSRDMIELGVINRDGSLERPLLLFPLVAMTPVNILGRDCLQGLGLRLTNL";
    const char* S4char = S4.c_str();
    std::cout << S4char << "\n\n";
    Sequence s4(10000, subMat.aa2int, subMat.int2aa, 0, kmer_size, true, true);
    s4.mapSequence(3,3, S4char);
    std::string S5 = "PQFSLWKRPVVTAYIEGQPVEVLLDTGADDSIVAGIELGNNYSPKIVGGIGGFINTLEYKNVEIEVLNKKVRATIMTGDTPINIFGRNILTALGMSLNL";
    const char* S5char = S5.c_str();
    std::cout << S5char << "\n\n";
    Sequence s5(10000, subMat.aa2int, subMat.int2aa, 0, kmer_size, true, true);
    s5.mapSequence(4,4, S5char);

    std::vector<Sequence *> seqSet;
    seqSet.push_back(&s2);
    seqSet.push_back(&s3);
    seqSet.push_back(&s4);
    //seqSet.push_back(s5);
    Matcher * aligner = new Matcher(10000, &subMat, 100000 ,seqSet.size(), false);
    MultipleAlignment msaAligner(1000, 10, &subMat, aligner);
    MultipleAlignment::MSAResult res = msaAligner.computeMSA(&s1, seqSet, true);
    MsaFilter filter(1000, 10000, &subMat);
    MsaFilter::MsaFilterResult msafilter = filter.filter((const char**)res.msaSequence, res.setSize, res.centerLength, 0, 0, -20.0f, 50, 100);
    std::cout << msafilter.setSize << std::endl;
    for(size_t i = 0; i < res.setSize; i++){
        std::cout << "Included sequence=" << (int) msafilter.keep[i] << std::endl;
    }
    MultipleAlignment::print(res, &subMat);
    PSSMCalculator pssm(&subMat, 1000, 1.0, 1.5);
    pssm.computePSSMFromMSA(res.setSize, res.centerLength, (const char**)res.msaSequence, false);
    pssm.printProfile(res.centerLength);
    pssm.printPSSM(res.centerLength);
    MultipleAlignment::deleteMSA(&res);
    delete aligner;
    return 0;
}

//PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF
//                     ALLDTGADDTVISEEDWPTDWPVMEAANPQIHGIGGGIPVRKSRDMIELGVINRDGSLERPLLLFPLVAMTPVNILGRDCLQGLGLRLTNL
