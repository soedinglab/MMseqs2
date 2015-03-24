//  main.cpp
//  forautocompl
//
//  Created by Martin Steinegger on 26.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//
#include <iostream>
#include <smith_waterman_sse2.h>
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "MultipleAlignment.h"
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
    std::cout << "Sequence (id 0):\n";
    std::string S1 = "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF";
    const char* S1char = S1.c_str();
    std::cout << S1char << "\n\n";
    Sequence* s1 = new Sequence(10000, subMat.aa2int, subMat.int2aa, 0, kmer_size, true);
    s1->mapSequence(0,"s1",S1char);
    std::string S2 = "PQFSLWKRPVVTAYIEGQPVEVLLDTGADDSIVAGIELGNNIVGGIGGFINTLEYKNVEIEVLNKKVRATIMTGDTPINIFGRNILTALGMSLNL";
    const char* S2char = S2.c_str();
    std::cout << S2char << "\n\n";
    Sequence* s2 = new Sequence(10000, subMat.aa2int, subMat.int2aa, 0, kmer_size, true);
    s2->mapSequence(1,"s2",S2char);
    std::string S3 = "PQFHLWKRPVVTAGQPVEVLLDTGADDSIVTGIELGPHYTPKIVGGIGGFINTKEYKNVEVEVLGKRIKGTIMTGDTPINIFGRNLLTALGMSLNF";
    const char* S3char = S3.c_str();
    std::cout << S3char << "\n\n";
    Sequence* s3 = new Sequence(10000, subMat.aa2int, subMat.int2aa, 0, kmer_size, true);
    s3->mapSequence(2,"s3", S3char);
    std::string S4 = "LAMTMEHKDRPLVRVILTNTGSHPVKQRSVYITALLDTGADDTVISEEDWPTDWPVMEAANPQIHGIGGGIPVRKSRDMIELGVINRDGSLERPLLLFPLVAMTPVNILGRDCLQGLGLRLTNL";
    const char* S4char = S4.c_str();
    std::cout << S4char << "\n\n";
    Sequence* s4 = new Sequence(10000, subMat.aa2int, subMat.int2aa, 0, kmer_size, true);
    s4->mapSequence(3,"s4", S4char);
    std::string S5 = "PQFSLWKRPVVTAYIEGQPVEVLLDTGADDSIVAGIELGNNYSPKIVGGIGGFINTLEYKNVEIEVLNKKVRATIMTGDTPINIFGRNILTALGMSLNL";
    const char* S5char = S5.c_str();
    std::cout << S5char << "\n\n";
    Sequence* s5 = new Sequence(10000, subMat.aa2int, subMat.int2aa, 0, kmer_size, true);
    s5->mapSequence(4,"s5", S5char);

    std::vector<Sequence *> seqSet;
    seqSet.push_back(s2);
    seqSet.push_back(s3);
    seqSet.push_back(s4);
    //seqSet.push_back(s5);
    MultipleAlignment msaAligner(1000,10,&subMat);
    MultipleAlignment::MSAResult res = msaAligner.computeMSA(s1,seqSet);
    MultipleAlignment::print(res);
    msaAligner.computePSSMFromMSA(res);
    return 0;
}

//PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF
//                     ALLDTGADDTVISEEDWPTDWPVMEAANPQIHGIGGGIPVRKSRDMIELGVINRDGSLERPLLLFPLVAMTPVNILGRDCLQGLGLRLTNL