#ifndef NUCLEOTIDE_MATRIX_H
#define NUCLEOTIDE_MATRIX_H

#include "SubstitutionMatrix.h"

class NucleotideMatrix : public SubstitutionMatrix {

    public:
        NucleotideMatrix(const char *scoringMatrixFileName_, float bitFactor, float scoreBias);

        virtual ~NucleotideMatrix();

        using BaseMatrix::getBitFactor;

    void setupLetterMapping();

};

#endif
