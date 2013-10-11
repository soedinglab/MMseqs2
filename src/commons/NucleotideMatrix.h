#ifndef NUCLEOTIDE_MATRIX_H
#define NUCLEOTIDE_MATRIX_H

#include "BaseMatrix.h"

class NucleotideMatrix : public BaseMatrix {

    public:
        NucleotideMatrix();

        virtual ~NucleotideMatrix();

        using BaseMatrix::getBitFactor;
};

#endif
