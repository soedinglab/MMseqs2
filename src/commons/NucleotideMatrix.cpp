#include "NucleotideMatrix.h"
#include <climits>
NucleotideMatrix::NucleotideMatrix(const char *scoringMatrixFileName_, float bitFactor, float scoreBias)
        : SubstitutionMatrix("nucleotide.out", bitFactor, scoreBias) {
//    this->alphabetSize = 5;
    setupLetterMapping();
    matrixName = "nucleotide.out";
}


void NucleotideMatrix::setupLetterMapping(){
    for(int letter = 0; letter < UCHAR_MAX; letter++){
        char upperLetter = toupper(static_cast<char>(letter));
        switch(upperLetter){
            case 'A':
            case 'T':
            case 'G':
            case 'C':
                this->aa2int[static_cast<int>(letter)] = this->aa2int[static_cast<int>(upperLetter)];
                break;
            default:
                this->aa2int[static_cast<int>(letter)] = this->aa2int[(int)'X'];
                break;
        }
    }
}


NucleotideMatrix::~NucleotideMatrix(){
}

