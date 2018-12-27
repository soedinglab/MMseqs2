#include "NucleotideMatrix.h"
#include <climits>

// TODO think about making the matrix dynamic
NucleotideMatrix::NucleotideMatrix(const char* /*scoringMatrixFileName_*/, float bitFactor, float scoreBias)
        : SubstitutionMatrix("nucleotide.out", bitFactor, scoreBias) {
//    this->alphabetSize = 5;
    matrixName = "nucleotide.out";
    setupLetterMapping();
    reverseLookup = new int[alphabetSize];
    reverseLookup[aa2int[static_cast<int>('A')]] = aa2int[static_cast<int>('T')];
    reverseLookup[aa2int[static_cast<int>('G')]] = aa2int[static_cast<int>('C')];
    reverseLookup[aa2int[static_cast<int>('C')]] = aa2int[static_cast<int>('G')];
    reverseLookup[aa2int[static_cast<int>('T')]] = aa2int[static_cast<int>('A')];
    reverseLookup[aa2int[static_cast<int>('X')]] = aa2int[static_cast<int>('X')];
}


void NucleotideMatrix::setupLetterMapping(){
    for(int letter = 0; letter < UCHAR_MAX; letter++){
        char upperLetter = toupper(static_cast<char>(letter));
        /*
         * R.................A or G
         * Y.................C or T
         * S.................G or C
         * W.................A or T
         * K.................G or T
         * M.................A or C
         * B.................C or G or T
         * D.................A or G or T
         * H.................A or C or T
         * V.................A or C or G
         */
        switch(upperLetter){
            case 'A':
            case 'T':
            case 'G':
            case 'C':
                this->aa2int[static_cast<int>(letter)] = this->aa2int[static_cast<int>(upperLetter)];
                break;
            case 'U':
            case 'W':
                this->aa2int[static_cast<int>(letter)] = this->aa2int[static_cast<int>('T')];
                break;
            case 'K':
            case 'B':
            case 'D':
            case 'V':
            case 'R':
            case 'S':
                this->aa2int[static_cast<int>(letter)] = this->aa2int[static_cast<int>('G')];
                break;
            case 'M':
            case 'Y':
            case 'H':
                this->aa2int[static_cast<int>(letter)] = this->aa2int[static_cast<int>('C')];
                break;
            default:
                this->aa2int[static_cast<int>(letter)] = this->aa2int[(int)'X'];
                break;
        }
    }
}


NucleotideMatrix::~NucleotideMatrix(){
    delete [] reverseLookup;
}

