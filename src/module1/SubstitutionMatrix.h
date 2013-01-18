#ifndef SUBSTITUTION_MATRIX_H
#define SUBSTITUTION_MATRIX_H

// Wrote by Maria Hauser mhauser@genzentrum.lmu.de
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

class SubstitutionMatrix {

    public:
        SubstitutionMatrix(char* scoringMatrixFileName, size_t alphabetSize);

        ~SubstitutionMatrix();

        void print();

        int* aa2int;

        char* int2aa;

        short** scMatrix;

        double* pBackground;
    
        double ** probMatrix;

        const size_t ALPHABET_SIZE;

    private:

        void readScoringMatrix();

        void readBiasedScoringMatrix(double bitFactor, double scoringBias);

        inline double _log2 (double x) { return log10(x)/0.301029996; }

        char* scoringMatrixFileName;

        void ltrim(std::string &s) {
            s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        }

        void rtrim(std::string &s) {
            s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        }

        void trim(std::string &s) {
            rtrim(s);
            ltrim(s);
        }

        std::vector<std::string>& split(const std::string &s, char delim, std::vector<std::string> &elems);

        std::vector<std::string> split(const std::string &s, char delim);
};

#endif
