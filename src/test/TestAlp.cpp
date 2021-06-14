/* $Id: $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's offical duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================*/

/*****************************************************************************

File name: gumbelparams_app.t.cpp

Author: Sergey Sheetlin

Contents: pairwise alignment algorithms

******************************************************************************/


#include <vector>
#include<fstream>
#include<string>
#include <SubstitutionMatrix.h>
#include <EvalueComputation.h>

#include "sls_alignment_evaluer.hpp"

const char* binary_name = "test_alp";

using namespace Sls;
using namespace std;

int main(int, const char**) {
    error ee_error("",0);


    long alphabetSize=25;//a number of letters in the amino acid alphabet

    //allocating two dimensional array alphabetSize X alphabetSize to store the scoring matrix
    long int i;
    long int **S=new long int *[alphabetSize];

    for(i=0;i<alphabetSize;i++)
    {
        S[i]=new long int [alphabetSize];
    };

    //assigning BLOSUM62 matrix
    S[0][0]=4; S[0][1]=-1; S[0][2]=-2; S[0][3]=-2; S[0][4]=0; S[0][5]=-1; S[0][6]=-1; S[0][7]=0; S[0][8]=-2; S[0][9]=-1; S[0][10]=-1; S[0][11]=-1; S[0][12]=-1; S[0][13]=-2; S[0][14]=-1; S[0][15]=1; S[0][16]=0; S[0][17]=-3; S[0][18]=-2; S[0][19]=0; S[0][20]=-2; S[0][21]=-1; S[0][22]=-1; S[0][23]=-1; S[0][24]=-4	;
    S[1][0]=-1; S[1][1]=5; S[1][2]=0; S[1][3]=-2; S[1][4]=-3; S[1][5]=1; S[1][6]=0; S[1][7]=-2; S[1][8]=0; S[1][9]=-3; S[1][10]=-2; S[1][11]=2; S[1][12]=-1; S[1][13]=-3; S[1][14]=-2; S[1][15]=-1; S[1][16]=-1; S[1][17]=-3; S[1][18]=-2; S[1][19]=-3; S[1][20]=-1; S[1][21]=-2; S[1][22]=0; S[1][23]=-1; S[1][24]=-4	;
    S[2][0]=-2; S[2][1]=0; S[2][2]=6; S[2][3]=1; S[2][4]=-3; S[2][5]=0; S[2][6]=0; S[2][7]=0; S[2][8]=1; S[2][9]=-3; S[2][10]=-3; S[2][11]=0; S[2][12]=-2; S[2][13]=-3; S[2][14]=-2; S[2][15]=1; S[2][16]=0; S[2][17]=-4; S[2][18]=-2; S[2][19]=-3; S[2][20]=4; S[2][21]=-3; S[2][22]=0; S[2][23]=-1; S[2][24]=-4	;
    S[3][0]=-2; S[3][1]=-2; S[3][2]=1; S[3][3]=6; S[3][4]=-3; S[3][5]=0; S[3][6]=2; S[3][7]=-1; S[3][8]=-1; S[3][9]=-3; S[3][10]=-4; S[3][11]=-1; S[3][12]=-3; S[3][13]=-3; S[3][14]=-1; S[3][15]=0; S[3][16]=-1; S[3][17]=-4; S[3][18]=-3; S[3][19]=-3; S[3][20]=4; S[3][21]=-3; S[3][22]=1; S[3][23]=-1; S[3][24]=-4	;
    S[4][0]=0; S[4][1]=-3; S[4][2]=-3; S[4][3]=-3; S[4][4]=9; S[4][5]=-3; S[4][6]=-4; S[4][7]=-3; S[4][8]=-3; S[4][9]=-1; S[4][10]=-1; S[4][11]=-3; S[4][12]=-1; S[4][13]=-2; S[4][14]=-3; S[4][15]=-1; S[4][16]=-1; S[4][17]=-2; S[4][18]=-2; S[4][19]=-1; S[4][20]=-3; S[4][21]=-1; S[4][22]=-3; S[4][23]=-1; S[4][24]=-4	;
    S[5][0]=-1; S[5][1]=1; S[5][2]=0; S[5][3]=0; S[5][4]=-3; S[5][5]=5; S[5][6]=2; S[5][7]=-2; S[5][8]=0; S[5][9]=-3; S[5][10]=-2; S[5][11]=1; S[5][12]=0; S[5][13]=-3; S[5][14]=-1; S[5][15]=0; S[5][16]=-1; S[5][17]=-2; S[5][18]=-1; S[5][19]=-2; S[5][20]=0; S[5][21]=-2; S[5][22]=4; S[5][23]=-1; S[5][24]=-4	;
    S[6][0]=-1; S[6][1]=0; S[6][2]=0; S[6][3]=2; S[6][4]=-4; S[6][5]=2; S[6][6]=5; S[6][7]=-2; S[6][8]=0; S[6][9]=-3; S[6][10]=-3; S[6][11]=1; S[6][12]=-2; S[6][13]=-3; S[6][14]=-1; S[6][15]=0; S[6][16]=-1; S[6][17]=-3; S[6][18]=-2; S[6][19]=-2; S[6][20]=1; S[6][21]=-3; S[6][22]=4; S[6][23]=-1; S[6][24]=-4	;
    S[7][0]=0; S[7][1]=-2; S[7][2]=0; S[7][3]=-1; S[7][4]=-3; S[7][5]=-2; S[7][6]=-2; S[7][7]=6; S[7][8]=-2; S[7][9]=-4; S[7][10]=-4; S[7][11]=-2; S[7][12]=-3; S[7][13]=-3; S[7][14]=-2; S[7][15]=0; S[7][16]=-2; S[7][17]=-2; S[7][18]=-3; S[7][19]=-3; S[7][20]=-1; S[7][21]=-4; S[7][22]=-2; S[7][23]=-1; S[7][24]=-4	;
    S[8][0]=-2; S[8][1]=0; S[8][2]=1; S[8][3]=-1; S[8][4]=-3; S[8][5]=0; S[8][6]=0; S[8][7]=-2; S[8][8]=8; S[8][9]=-3; S[8][10]=-3; S[8][11]=-1; S[8][12]=-2; S[8][13]=-1; S[8][14]=-2; S[8][15]=-1; S[8][16]=-2; S[8][17]=-2; S[8][18]=2; S[8][19]=-3; S[8][20]=0; S[8][21]=-3; S[8][22]=0; S[8][23]=-1; S[8][24]=-4	;
    S[9][0]=-1; S[9][1]=-3; S[9][2]=-3; S[9][3]=-3; S[9][4]=-1; S[9][5]=-3; S[9][6]=-3; S[9][7]=-4; S[9][8]=-3; S[9][9]=4; S[9][10]=2; S[9][11]=-3; S[9][12]=1; S[9][13]=0; S[9][14]=-3; S[9][15]=-2; S[9][16]=-1; S[9][17]=-3; S[9][18]=-1; S[9][19]=3; S[9][20]=-3; S[9][21]=3; S[9][22]=-3; S[9][23]=-1; S[9][24]=-4	;
    S[10][0]=-1; S[10][1]=-2; S[10][2]=-3; S[10][3]=-4; S[10][4]=-1; S[10][5]=-2; S[10][6]=-3; S[10][7]=-4; S[10][8]=-3; S[10][9]=2; S[10][10]=4; S[10][11]=-2; S[10][12]=2; S[10][13]=0; S[10][14]=-3; S[10][15]=-2; S[10][16]=-1; S[10][17]=-2; S[10][18]=-1; S[10][19]=1; S[10][20]=-4; S[10][21]=3; S[10][22]=-3; S[10][23]=-1; S[10][24]=-4	;
    S[11][0]=-1; S[11][1]=2; S[11][2]=0; S[11][3]=-1; S[11][4]=-3; S[11][5]=1; S[11][6]=1; S[11][7]=-2; S[11][8]=-1; S[11][9]=-3; S[11][10]=-2; S[11][11]=5; S[11][12]=-1; S[11][13]=-3; S[11][14]=-1; S[11][15]=0; S[11][16]=-1; S[11][17]=-3; S[11][18]=-2; S[11][19]=-2; S[11][20]=0; S[11][21]=-3; S[11][22]=1; S[11][23]=-1; S[11][24]=-4	;
    S[12][0]=-1; S[12][1]=-1; S[12][2]=-2; S[12][3]=-3; S[12][4]=-1; S[12][5]=0; S[12][6]=-2; S[12][7]=-3; S[12][8]=-2; S[12][9]=1; S[12][10]=2; S[12][11]=-1; S[12][12]=5; S[12][13]=0; S[12][14]=-2; S[12][15]=-1; S[12][16]=-1; S[12][17]=-1; S[12][18]=-1; S[12][19]=1; S[12][20]=-3; S[12][21]=2; S[12][22]=-1; S[12][23]=-1; S[12][24]=-4	;
    S[13][0]=-2; S[13][1]=-3; S[13][2]=-3; S[13][3]=-3; S[13][4]=-2; S[13][5]=-3; S[13][6]=-3; S[13][7]=-3; S[13][8]=-1; S[13][9]=0; S[13][10]=0; S[13][11]=-3; S[13][12]=0; S[13][13]=6; S[13][14]=-4; S[13][15]=-2; S[13][16]=-2; S[13][17]=1; S[13][18]=3; S[13][19]=-1; S[13][20]=-3; S[13][21]=0; S[13][22]=-3; S[13][23]=-1; S[13][24]=-4	;
    S[14][0]=-1; S[14][1]=-2; S[14][2]=-2; S[14][3]=-1; S[14][4]=-3; S[14][5]=-1; S[14][6]=-1; S[14][7]=-2; S[14][8]=-2; S[14][9]=-3; S[14][10]=-3; S[14][11]=-1; S[14][12]=-2; S[14][13]=-4; S[14][14]=7; S[14][15]=-1; S[14][16]=-1; S[14][17]=-4; S[14][18]=-3; S[14][19]=-2; S[14][20]=-2; S[14][21]=-3; S[14][22]=-1; S[14][23]=-1; S[14][24]=-4	;
    S[15][0]=1; S[15][1]=-1; S[15][2]=1; S[15][3]=0; S[15][4]=-1; S[15][5]=0; S[15][6]=0; S[15][7]=0; S[15][8]=-1; S[15][9]=-2; S[15][10]=-2; S[15][11]=0; S[15][12]=-1; S[15][13]=-2; S[15][14]=-1; S[15][15]=4; S[15][16]=1; S[15][17]=-3; S[15][18]=-2; S[15][19]=-2; S[15][20]=0; S[15][21]=-2; S[15][22]=0; S[15][23]=-1; S[15][24]=-4	;
    S[16][0]=0; S[16][1]=-1; S[16][2]=0; S[16][3]=-1; S[16][4]=-1; S[16][5]=-1; S[16][6]=-1; S[16][7]=-2; S[16][8]=-2; S[16][9]=-1; S[16][10]=-1; S[16][11]=-1; S[16][12]=-1; S[16][13]=-2; S[16][14]=-1; S[16][15]=1; S[16][16]=5; S[16][17]=-2; S[16][18]=-2; S[16][19]=0; S[16][20]=-1; S[16][21]=-1; S[16][22]=-1; S[16][23]=-1; S[16][24]=-4	;
    S[17][0]=-3; S[17][1]=-3; S[17][2]=-4; S[17][3]=-4; S[17][4]=-2; S[17][5]=-2; S[17][6]=-3; S[17][7]=-2; S[17][8]=-2; S[17][9]=-3; S[17][10]=-2; S[17][11]=-3; S[17][12]=-1; S[17][13]=1; S[17][14]=-4; S[17][15]=-3; S[17][16]=-2; S[17][17]=11; S[17][18]=2; S[17][19]=-3; S[17][20]=-4; S[17][21]=-2; S[17][22]=-2; S[17][23]=-1; S[17][24]=-4	;
    S[18][0]=-2; S[18][1]=-2; S[18][2]=-2; S[18][3]=-3; S[18][4]=-2; S[18][5]=-1; S[18][6]=-2; S[18][7]=-3; S[18][8]=2; S[18][9]=-1; S[18][10]=-1; S[18][11]=-2; S[18][12]=-1; S[18][13]=3; S[18][14]=-3; S[18][15]=-2; S[18][16]=-2; S[18][17]=2; S[18][18]=7; S[18][19]=-1; S[18][20]=-3; S[18][21]=-1; S[18][22]=-2; S[18][23]=-1; S[18][24]=-4	;
    S[19][0]=0; S[19][1]=-3; S[19][2]=-3; S[19][3]=-3; S[19][4]=-1; S[19][5]=-2; S[19][6]=-2; S[19][7]=-3; S[19][8]=-3; S[19][9]=3; S[19][10]=1; S[19][11]=-2; S[19][12]=1; S[19][13]=-1; S[19][14]=-2; S[19][15]=-2; S[19][16]=0; S[19][17]=-3; S[19][18]=-1; S[19][19]=4; S[19][20]=-3; S[19][21]=2; S[19][22]=-2; S[19][23]=-1; S[19][24]=-4	;
    S[20][0]=-2; S[20][1]=-1; S[20][2]=4; S[20][3]=4; S[20][4]=-3; S[20][5]=0; S[20][6]=1; S[20][7]=-1; S[20][8]=0; S[20][9]=-3; S[20][10]=-4; S[20][11]=0; S[20][12]=-3; S[20][13]=-3; S[20][14]=-2; S[20][15]=0; S[20][16]=-1; S[20][17]=-4; S[20][18]=-3; S[20][19]=-3; S[20][20]=4; S[20][21]=-3; S[20][22]=0; S[20][23]=-1; S[20][24]=-4	;
    S[21][0]=-1; S[21][1]=-2; S[21][2]=-3; S[21][3]=-3; S[21][4]=-1; S[21][5]=-2; S[21][6]=-3; S[21][7]=-4; S[21][8]=-3; S[21][9]=3; S[21][10]=3; S[21][11]=-3; S[21][12]=2; S[21][13]=0; S[21][14]=-3; S[21][15]=-2; S[21][16]=-1; S[21][17]=-2; S[21][18]=-1; S[21][19]=2; S[21][20]=-3; S[21][21]=3; S[21][22]=-3; S[21][23]=-1; S[21][24]=-4	;
    S[22][0]=-1; S[22][1]=0; S[22][2]=0; S[22][3]=1; S[22][4]=-3; S[22][5]=4; S[22][6]=4; S[22][7]=-2; S[22][8]=0; S[22][9]=-3; S[22][10]=-3; S[22][11]=1; S[22][12]=-1; S[22][13]=-3; S[22][14]=-1; S[22][15]=0; S[22][16]=-1; S[22][17]=-2; S[22][18]=-2; S[22][19]=-2; S[22][20]=0; S[22][21]=-3; S[22][22]=4; S[22][23]=-1; S[22][24]=-4	;
    S[23][0]=-1; S[23][1]=-1; S[23][2]=-1; S[23][3]=-1; S[23][4]=-1; S[23][5]=-1; S[23][6]=-1; S[23][7]=-1; S[23][8]=-1; S[23][9]=-1; S[23][10]=-1; S[23][11]=-1; S[23][12]=-1; S[23][13]=-1; S[23][14]=-1; S[23][15]=-1; S[23][16]=-1; S[23][17]=-1; S[23][18]=-1; S[23][19]=-1; S[23][20]=-1; S[23][21]=-1; S[23][22]=-1; S[23][23]=-1; S[23][24]=-4	;
    S[24][0]=-4; S[24][1]=-4; S[24][2]=-4; S[24][3]=-4; S[24][4]=-4; S[24][5]=-4; S[24][6]=-4; S[24][7]=-4; S[24][8]=-4; S[24][9]=-4; S[24][10]=-4; S[24][11]=-4; S[24][12]=-4; S[24][13]=-4; S[24][14]=-4; S[24][15]=-4; S[24][16]=-4; S[24][17]=-4; S[24][18]=-4; S[24][19]=-4; S[24][20]=-4; S[24][21]=-4; S[24][22]=-4; S[24][23]=-4; S[24][24]=1	;

    //allocation of arrays that store the background probabilities: "letterFreqs1" for sequence #1 and "letterFreqs2" for sequence #2
    double *letterFreqs1=new double[alphabetSize];
    double *letterFreqs2=new double[alphabetSize];

    //assigning background probabilities (Robinson and Robinson frequencies for the alphabet ARNDCQEGHILKMFPSTWYVBJZX*)
    letterFreqs1[0]=0.07805; letterFreqs1[1]=0.05129; letterFreqs1[2]=0.04487; letterFreqs1[3]=0.05364; letterFreqs1[4]=0.01925; letterFreqs1[5]=0.04264; letterFreqs1[6]=0.06295; letterFreqs1[7]=0.07377; letterFreqs1[8]=0.02199; letterFreqs1[9]=0.05142; letterFreqs1[10]=0.09019; letterFreqs1[11]=0.05744; letterFreqs1[12]=0.02243; letterFreqs1[13]=0.03856; letterFreqs1[14]=0.05203; letterFreqs1[15]=0.0712; letterFreqs1[16]=0.05841; letterFreqs1[17]=0.0133; letterFreqs1[18]=0.03216; letterFreqs1[19]=0.06441; letterFreqs1[20]=0; letterFreqs1[21]=0; letterFreqs1[22]=0; letterFreqs1[23]=0; letterFreqs1[24]=0;

    //the probabilities are symmetric
    for(i=0;i<alphabetSize;i++)
    {
        letterFreqs2[i]=letterFreqs1[i];
    };

    //maximum allowed calculation time in seconds
    double max_seconds=60.0;

    //setting better mnemonics for the scoring matrix array
//    const long ** substitutionScoreMatrix=(const long **)S;

    //creating an object to test the functions (and store the Gumbel parameters)
    AlignmentEvaluer evaluer;

    evaluer.set_gapped_computation_parameters_simplified(max_seconds);
    const double lambdaTolerance = 0.01;
    const double kTolerance = 0.05;
    const double maxMegabytes = 500;
    const long randomSeed = 42;

    Parameters& par = Parameters::getInstance();
    par.initMatrices();
    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 2.0, 0.0);
    long ** tmpMat = new long *[subMat.alphabetSize];
    long * tmpMatData = new long[subMat.alphabetSize*subMat.alphabetSize];

    for(int i = 0; i < subMat.alphabetSize; i++) {
        tmpMat[i] = &tmpMatData[i * subMat.alphabetSize];
        for (int j = 0; j < subMat.alphabetSize; j++) {
            tmpMat[i][j] = subMat.subMatrix[i][j];
        }
        letterFreqs1[i] = subMat.pBack[i];
        letterFreqs2[i] = subMat.pBack[i];
        std::cout << subMat.pBack[i] << std::endl;
    }
    //the function "initGapless" is called for BLOSUM62 matrix and Robinson and Robinson frequencies; the "AlignmentEvaluer" object is initialized
    evaluer.initGapped(
            subMat.alphabetSize-1,
            tmpMat,
            letterFreqs1,
            letterFreqs2,
            11,1,11,1, false,lambdaTolerance,kTolerance,max_seconds,maxMegabytes,randomSeed);


std::cout << std::setprecision(20) <<
    evaluer.parameters().lambda <<"\t" <<
    evaluer.parameters().K <<"\t" <<
    evaluer.parameters().a_J<<"\t" <<
    evaluer.parameters().b_J<<"\t" <<
    evaluer.parameters().a_I<<"\t" <<
    evaluer.parameters().b_I<<"\t" <<
    evaluer.parameters().alpha_J<<"\t" <<
    evaluer.parameters().beta_J<<"\t" <<
    evaluer.parameters().alpha_I<<"\t" <<
    evaluer.parameters().beta_I<<"\t" <<
    evaluer.parameters().sigma<<"\t" <<
    evaluer.parameters().tau<<"\t" << std::endl;


    evaluer.initGapless(
            subMat.alphabetSize-1,
            tmpMat,
            letterFreqs1,
            letterFreqs2,max_seconds);


    std::cout << std::setprecision(20) <<
              evaluer.parameters().lambda <<"\t" <<
              evaluer.parameters().K <<"\t" <<
              evaluer.parameters().a_J<<"\t" <<
              evaluer.parameters().b_J<<"\t" <<
              evaluer.parameters().a_I<<"\t" <<
              evaluer.parameters().b_I<<"\t" <<
              evaluer.parameters().alpha_J<<"\t" <<
              evaluer.parameters().beta_J<<"\t" <<
              evaluer.parameters().alpha_I<<"\t" <<
              evaluer.parameters().beta_I<<"\t" <<
              evaluer.parameters().sigma<<"\t" <<
              evaluer.parameters().tau<<"\t" << std::endl;



//    BlastScoreUtils::BlastStat stats = BlastScoreUtils::getAltschulStatsForMatrix("blosum62", 11, 1, true);
//    double kmnByLen[350];
//    for(int len = 0; len < 350; len++){
//        kmnByLen[len] = BlastScoreUtils::computeKmn(len, stats.K, stats.lambda, stats.alpha, stats.beta, 3500000000, 10000000);
//    }
//    double d = 0.1;
//    for(size_t i = 10; i < 300; i+=20) {
//        for (size_t j = 10; j < 350; j ++) {
//            double area =  evaluer.area(i, j, static_cast<double>(3500000000));
//
//            double evalPerArea  =  evaluer.evaluePerArea(i);
//            double evalue_alp = area * evalPerArea;
//            double evalue = BlastScoreUtils::computeEvalue(i, kmnByLen[j], stats.lambda);
//            std::cout << "s: " << i << "\tl: " << j << "\t" << evalue_alp << "\t" << evalue << std::endl;
//
//        }
//    }
//    printf("%f", d);
//    for(size_t i = 100; i < 200; i+=50){
//        printf("%d\t%d\t%f\n", i , 10, static_cast<float>(evaluer.area(i, 10, 10000000)));
//        printf("%d\t%d\t%f\n", i , 50, static_cast<float>(evaluer.area(i, 50, 10000000)));
//        printf("%d\t%d\t%f\n", i , 100, static_cast<float>(evaluer.area(i, 100, 10000000)));
//        printf("%d\t%d\t%f\n", i , 200, static_cast<float>(evaluer.area(i, 200, 10000000)));
//        printf("%d\t%d\t%f\n", i , 350, static_cast<float>(evaluer.area(i, 350, 10000000)));
//    }
    //outputting the gapless parameters into a file
//    string st_gapless="gapless_lib.out";
//    ofstream f_gapless(st_gapless.data());
//    if(!f_gapless)
//    {
//        throw error("Error - file "+st_gapless+" is not found\n",3);
//    };

    //the "operator <<" is called for the object what results outputting the gapless Gumbel parameters into the file "gapless_lib.out"


    //memory release
    delete[]letterFreqs2;
    delete[]letterFreqs1;

    for(i=0;i<alphabetSize;i++)
    {
        delete[]S[i];
    };
    delete[]S;


}


