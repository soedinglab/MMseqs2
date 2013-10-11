#ifndef SUBSTITUTION_MATRIX_H
#define SUBSTITUTION_MATRIX_H

// Written by Maria Hauser mhauser@genzentrum.lmu.de
//  
// Represents a simple amino acid substitution matrix. 
// Can read a straight-forward substitution matrix file (.mat ending of the file) and probabilities file (.out ending of the file).
// If a probabilities file is given, it calculates a biased matrix (produces shorter = more precise alignments).
//


#include <iostream>
#include <cmath>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <algorithm>

#include "BaseMatrix.h"

class SubstitutionMatrix: public BaseMatrix {

    public:
        SubstitutionMatrix(const char* scoringMatrixFileName_, float bitFactor);

        virtual ~SubstitutionMatrix();

        virtual float getBitFactor() {return bitFactor; }

    private:

        const char* scoringMatrixFileName;

        void readProbMatrix();

        float bitFactor;

};

#endif
