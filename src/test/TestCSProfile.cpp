//  main.cpp
//  forautocompl
//
//  Created by Martin Steinegger on 26.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <ProfileStates.h>
#include <CSProfile.h>
#include "Parameters.h"
#include "StripedSmithWaterman.h"
#include "MsaFilter.h"
#include "PSSMCalculator.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "MultipleAlignment.h"

const char* binary_name = "test_csprofile";


int main (int argc, const char * argv[])
{
    Parameters& par = Parameters::getInstance();
    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0, 0.0);
    

    
    std::cout << "Subustitution matrix:";
    SubstitutionMatrix::print(subMat.subMatrix,subMat.int2aa,subMat.alphabetSize);
    //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);

    std::string tim1 = "MAAWCSPRWLRVAVGTPRLPAAAGRGVQQPQGGVVAASLCRKLCVSAFGLSMGAHGPRAL"
            "LTLRPGVRLTGTKSFPFVCTASFHTSASLAKDDYYQILGVPRNASQKDIKKAYYQLAKKY"
            "HPDTNKDDPKAKEKFSQLAEAYEVLSDEVKRKQYDAYGSAGFDPGASSSGQGYWRGGPSV"
            "DPEELFRKIFGEFSSSPFGDFQNVFDQPQEYIMELTFNQAAKGVNKEFTVNIMDTCERCD"
            "GKGNEPGTKVQHCHYCSGSGMETINTGPFVMRSTCRRCGGRGSIITNPCVVCRGAGQAKQ"
            "KKRVTVPVPAGVEDGQTVRMPVGKREIFVTFRVQKSPVFRRDGADIHSDLFISIAQAILG"
            "GTAKAQGLYETINVTIPAGIQTDQKIRLTGKGIPRINSYGYGDHYIHIKIRVPKRLSSRQ"
            "QNLILSYAEDETDVEGTVNGVTHTSTGKRSTGN";
    Sequence* s = new Sequence(10000, 0, &subMat, 6, true, false);
    s->mapSequence(0,0,tim1.c_str());
     //std::string libraryString((const char *)Library1_lib, Library1_lib_len);
    CSProfile * ps = new CSProfile(10000);

    
    //std::ifstream libFile;
    //libFile.open ("/Users/mad/Documents/workspace/csblast/data/K4000.crf");
    
//    std::string libData((std::istreambuf_iterator<char>(libFile)),
//                 std::istreambuf_iterator<char>());
    //ps->read(libData);

    float * profile = ps->computeProfile(s, 1.0, 0.9);
    printf("Query profile of sequence\n");
    printf("Pos ");
    for(size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
        printf("%3c ", subMat.int2aa[aa]);
    }
    printf("\n");
    for(int i = 0; i < s->L; i++){
        printf("%3c ", subMat.int2aa[s->int_sequence[i]]);
        for(size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++){
            const float aaProb = profile[i * Sequence::PROFILE_AA_SIZE + aa];

            float logProb = MathUtil::flog2(aaProb / subMat.pBack[aa]);
            const float pssmVal = 2.0 * logProb ;
            char finalPssmVal = static_cast<char>((pssmVal < 0.0) ? pssmVal - 0.5 : pssmVal + 0.5);

            printf("%3d ", (int)finalPssmVal);
//            printf("%3d ", profile_score[i * profile_row_size + aa] );
        }
        printf("\n");
    }
    delete s;
    delete ps;

}

//PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF
//                     ALLDTGADDTVISEEDWPTDWPVMEAANPQIHGIGGGIPVRKSRDMIELGVINRDGSLERPLLLFPLVAMTPVNILGRDCLQGLGLRLTNL
