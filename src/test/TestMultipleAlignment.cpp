//  main.cpp
//  forautocompl
//
//  Created by Martin Steinegger on 26.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//
#include <iostream>
#include <StripedSmithWaterman.h>
#include <MsaFilter.h>
#include "PSSMCalculator.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "MultipleAlignment.h"
#include "Parameters.h"

const char* binary_name = "test_multiplealignment";

int main(int, const char**) {
    Parameters& par = Parameters::getInstance();
    par.initMatrices();

    const size_t kmer_size=6;

    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 2.0, 0);
    std::cout << "Subustitution matrix:\n";
    SubstitutionMatrix::print(subMat.subMatrix,subMat.num2aa,subMat.alphabetSize);
    //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);
    std::cout << "\n";

    std::cout << "subMatrix:\n";
    //    ReducedMatrix subMat(subMat.probMatrix, 20);
    //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);
    std::cout << "\n";
    EvalueComputation evaluer(100000, &subMat, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
    Matcher * aligner = new Matcher(Parameters::DBTYPE_AMINO_ACIDS, 10000, &subMat, &evaluer, false, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
    std::vector<Matcher::result_t> alnResults;
    std::vector<std::vector<unsigned char>> seqSet;
    std::cout << "Sequence (id 0):\n";
    std::string S1 = "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF";
    const char* S1char = S1.c_str();
    std::cout << S1char << "\n\n";
    Sequence q(10000, 0, &subMat, kmer_size, true, true);
    q.mapSequence(0,0,S1char, S1.size());
    aligner->initQuery(&q);
    std::string S2 = "PQFSLWKRPVVTAYIEGQPVEVLLDTGADDSIVAGIELGNNIVGGIGGFINTLEYKNVEIEVLNKKVRATIMTGDTPINIFGRNILTALGMSLNL";
    const char* S2char = S2.c_str();
    std::cout << S2char << "\n\n";
    Sequence s(10000, 0, &subMat, kmer_size, true, true);
    s.mapSequence(1,1,S2char, S2.size());
    seqSet.emplace_back(s.numSequence, s.numSequence + s.L);
    alnResults.emplace_back(aligner->getSWResult(&s, INT_MAX, false, 0, 0.0, FLT_MAX, Matcher::SCORE_COV_SEQID, 0, false));
    std::string S3 = "PQFHLWKRPVVTAGQPVEVLLDTGADDSIVTGIELGPHYTPKIVGGIGGFINTKEYKNVEVEVLGKRIKGTIMTGDTPINIFGRNLLTALGMSLNF";
    const char* S3char = S3.c_str();
    std::cout << S3char << "\n\n";
    s.mapSequence(2,2, S3char, S3.size());
    seqSet.emplace_back(s.numSequence, s.numSequence + s.L);
    alnResults.emplace_back(aligner->getSWResult(&s, INT_MAX, false, 0, 0.0, FLT_MAX, Matcher::SCORE_COV_SEQID, 0, false));
    std::string S4 = "LAMTMEHKDRPLVRVILTNTGSHPVKQRSVYITALLDTGADDTVISEEDWPTDWPVMEAANPQIHGIGGGIPVRKSRDMIELGVINRDGSLERPLLLFPLVAMTPVNILGRDCLQGLGLRLTNL";
    const char* S4char = S4.c_str();
    std::cout << S4char << "\n\n";
    s.mapSequence(3,3, S4char, S4.size());
    seqSet.emplace_back(s.numSequence, s.numSequence + s.L);
    alnResults.emplace_back(aligner->getSWResult(&s, INT_MAX, false, 0, 0.0, FLT_MAX, Matcher::SCORE_COV_SEQID, 0, false));
    std::string S5 = "PQFSLWKRPVVTAYIEGQPVEVLLDTGADDSIVAGIELGNNYSPKIVGGIGGFINTLEYKNVEIEVLNKKVRATIMTGDTPINIFGRNILTALGMSLNL";
    const char* S5char = S5.c_str();
    std::cout << S5char << "\n\n";
    s.mapSequence(4,4, S5char, S5.size());
    seqSet.emplace_back(s.numSequence, s.numSequence + s.L);
    alnResults.emplace_back(aligner->getSWResult(&s, INT_MAX, false, 0, 0.0, FLT_MAX, Matcher::SCORE_COV_SEQID, 0, false));

    MultipleAlignment msaAligner(1000, &subMat);
    MultipleAlignment::MSAResult res = msaAligner.computeMSA(&s, seqSet, alnResults, true);
    MsaFilter filter(1000, 10000, &subMat, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
    size_t filterSetSize = filter.filter(res, alnResults, 0, 0, -20.0, 50, 100);
    std::cout << "Filtered:" << filterSetSize << std::endl;
    MultipleAlignment::print(res, &subMat);
    PSSMCalculator pssm(&subMat, 1000, 5, 1.0, 1.5);
    pssm.computePSSMFromMSA(filterSetSize, res.centerLength, (const char**)res.msaSequence, false);
    pssm.printProfile(res.centerLength);
    pssm.printPSSM(res.centerLength);
    MultipleAlignment::deleteMSA(&res);
    delete aligner;
    return 0;
}

//PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF
//                     ALLDTGADDTVISEEDWPTDWPVMEAANPQIHGIGGGIPVRKSRDMIELGVINRDGSLERPLLLFPLVAMTPVNILGRDCLQGLGLRLTNL
