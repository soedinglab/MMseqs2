//
// Created by mad on 10/26/15.
//

#include <iostream>
#include <cstring>
#include <math.h>
#include "tantan.h"
#include "SubstitutionMatrix.h"
#include "Sequence.h"
#include "Parameters.h"

const char* binary_name = "test_tantan";

int main (int, const char**) {
    const size_t kmer_size = 6;

    Parameters& par = Parameters::getInstance();
    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0);
    std::cout << "Substitution matrix:";
    SubstitutionMatrix::print(subMat.subMatrix, subMat.num2aa, subMat.alphabetSize);

    const char *ref = "MTLHSNSTTSSLFPNISSSWIHSPSDAGLPPGTVTHFGSYNVSRAAGNFSSPDGTTDDPLGGHTVWQVVFIAFLTGILALVTIIGNILVIVSFKVNKQLKTVNNYFLLSLACADLIIGVISMNLFTTYIIMNRWALGNLACDLWLAIDYVASNASVMNLLVISFDRYFSITRPLTYRAKRTTKRAGVMIGLAWVISFVLWAPAILFWQYFVGKRTVPPGECFIQFLSEPTITFGTAIAAFYMPVTIMTILYWRIYKETEKRTKELAGLQASGTEAETENFVHPTGSSRSCSSYELQQQSMKRSNRRKYGRCHFWFTTKSWKPSSEQMDQDHSSSDSWNNNDAAASLENSASSDEEDIGSETRAIYSIVLKLPGHSTILNSTKLPSSDNLQVPEEELGMVDLERKADKLQAQKSVDDGGSFPKSFSKLPIQLESAVDTAKTSDVNSSVGKSTATLPLSFKEATLAKRFALKTRSQITKRKRMSLVKEKKAAQTLSAILLAFIITWTPYNIMVLVNTFCDSCIPKTFWNLGYWLCYINSTVNPVCYALCNKTFRTTFKMLLLCQCDKKKRRKQQYQRQSVIFHKRAPEQAL";
    const size_t len = strlen(ref);
    Sequence refSeq(10000, 0, &subMat, kmer_size, false, true);
    refSeq.mapSequence(0, 0, ref, strlen(ref));

    char hardMaskTable[256];
    std::fill_n(hardMaskTable, 256, subMat.aa2num[(int) 'X']);
    double probMatrix[21][21];

    const double *probMatrixPointers[64];

    for (int i = 0; i < 21; ++i){
        probMatrixPointers[i] = probMatrix[i];
        for(int j = 0; j < 21; ++j){
            probMatrix[i][j]  = exp(0.324032 * subMat.subMatrix[i][j]);
            //std::cout << probMatrix[i][j] << "\t";
        }
        //std::cout << std::endl;
    }
    char  refInt[100000];


    for(size_t i = 0; i < 100000; i++){
        for(int i = 0; i < refSeq.L; i++){
            refInt[i] = (char) refSeq.numSequence[i];
        }
        tantan::maskSequences(refInt, refInt+len, 50 /*options.maxCycleLength*/,
                              probMatrixPointers,
                              0.005 /*options.repeatProb*/, 0.05 /*options.repeatEndProb*/,
                              0.9 /*options.repeatOffsetProbDecay*/,
                              0, 0,
                              0.5 /*options.minMaskProb*/, hardMaskTable);
    }
    for(int i = 0; i < refSeq.L; i++){
//        refInt[i] = (char) refSeq.sequence[i];
        std::cout << subMat.num2aa[(int)refInt[i]];

    }
    std::cout << std::endl;
    return 0;
}
