//
// Created by mad on 10/26/15.
//

#include <iostream>

#include "SubstitutionMatrix.h"
#include "Sequence.h"
#include "Parameters.h"

const char* binary_name = "test_compositionbias";

void calcLocalAaBiasCorrection(Sequence* seq, SubstitutionMatrix * m){
    const int windowSize = 40;
    short * composition = new short[seq->L];
    // calculate local amino acid bias
    for (int i = 0; i < seq->L; i++){
        const int minPos = std::max(0, (i - windowSize/2));
        const int maxPos = std::min(seq->L, (i + windowSize/2));
        const int _2d = maxPos - minPos;
        // negative score for the amino acids in the neighborhood of i
        int sumSubScores = 0;
        short * subMat = m->subMatrix[seq->numSequence[i]];
        for (int j = minPos; j < maxPos; j++){
            sumSubScores += (j != i) ? subMat[seq->numSequence[j]] : 0;
        }
        float deltaS_i = (float) sumSubScores;
        deltaS_i /= -1.0 * _2d;
        // positive score for the background score distribution for i
        for (int a = 0; a < m->alphabetSize; a++){
            deltaS_i += m->pBack[a] * subMat[a];
        }
        composition[i] = deltaS_i;
        std::cout << i << ": "  << _2d << " "  << sumSubScores << " " << deltaS_i << std::endl;
    }

    for (int i = 0; i < seq->L-6; i++) {
        short score = 0;
        for (int j = 0; j < 6; j++) {
            score += composition[i+j];
        }
        std::cout << i << " : " << score << std::endl;
    }

        delete [] composition;
}



int main (int, const char**) {
    const size_t kmer_size = 6;

    Parameters& par = Parameters::getInstance();
    par.initMatrices();
    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 8.0, 0);
    std::cout << "Substitution matrix:";
    SubstitutionMatrix::print(subMat.subMatrix, subMat.num2aa, subMat.alphabetSize);

    const char *ref = "MDDVKIERLKRLNEDVLEDLIEVYMRGYEGLEEYGGEGRDYARDYIKWCWKKAPDGFFVAKVGDRIVGFIVCDRDWYSRYEGKIVGAIHEFVVDKGWQGKGIGKKLLTKCLEFLGKYNDTIELWVGEKNFGAMRLYEKFGFKKVGKSGIWIRMVRRQLS";
    Sequence refSeq(10000, 0, &subMat, kmer_size, false, true);
    refSeq.mapSequence(0, 0, ref, strlen(ref));


    calcLocalAaBiasCorrection(&refSeq, &subMat);

    return 0;
}
