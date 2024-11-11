#ifndef CONVERT_CUH
#define CONVERT_CUH

#include <cctype>
#include <cstring>

namespace cudasw4{

struct ConvertAA_20{
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    char operator()(const char& AA) {
        // ORDER of AminoAcids following MMseqs2
        // lower-case has upper-case encoding
        switch(AA & ~32) { // bit twiddling to turn all AA to uppercase
            case 'A': return 0;
            case 'C': return 1;
            case 'D':
            case 'B': return 2;
            case 'Z':
            case 'E': return 3;
            case 'F': return 4;
            case 'G': return 5;
            case 'H': return 6;
            case 'I': return 7;
            case 'K': return 8;
            case 'J':
            case 'L': return 9;
            case 'M': return 10;
            case 'N': return 11;
            case 'P': return 12;
            case 'Q': return 13;
            case 'R': return 14;
            case 'S': return 15;
            case 'T': return 16;
            case 'V': return 17;
            case 'W': return 18;
            case 'Y': return 19;
            default:  return 20;
        }
    }
};

struct InverseConvertAA_20{
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    char operator()(const char& AA) {
        switch(AA) {
            case 0:  return 'A';
            case 1:  return 'C';
            case 2:  return 'D';
            case 3:  return 'E';
            case 4:  return 'F';
            case 5:  return 'G';
            case 6:  return 'H';
            case 7:  return 'I';
            case 8:  return 'K';
            case 9:  return 'L';
            case 10: return 'M';
            case 11: return 'N';
            case 12: return 'P';
            case 13: return 'Q';
            case 14: return 'R';
            case 15: return 'S';
            case 16: return 'T';
            case 17: return 'V';
            case 18: return 'W';
            case 19: return 'Y';
            default: return 'X';
        }
    }
};


struct ConvertAA_20_CaseSensitive{
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    char operator()(const char& AA) {
        switch (AA) {
            case 'A': return 0;
            case 'C': return 1;
            case 'B':
            case 'D': return 2;
            case 'Z':
            case 'E': return 3;
            case 'F': return 4;
            case 'G': return 5;
            case 'H': return 6;
            case 'I': return 7;
            case 'K': return 8;
            case 'J':
            case 'L': return 9;
            case 'M': return 10;
            case 'N': return 11;
            case 'P': return 12;
            case 'Q': return 13;
            case 'R': return 14;
            case 'S': return 15;
            case 'T': return 16;
            case 'V': return 17;
            case 'W': return 18;
            case 'Y': return 19;

            case 'a': return 32;
            case 'c': return 33;
            case 'b':
            case 'd': return 34;
            case 'z':
            case 'e': return 35;
            case 'f': return 36;
            case 'g': return 37;
            case 'h': return 38;
            case 'i': return 39;
            case 'k': return 40;
            case 'j':
            case 'l': return 41;
            case 'm': return 42;
            case 'n': return 43;
            case 'p': return 44;
            case 'q': return 45;
            case 'r': return 46;
            case 's': return 47;
            case 't': return 48;
            case 'v': return 49;
            case 'w': return 50;
            case 'y': return 51;
            default: return 20;
        }
    }
};

struct InverseConvertAA_20_CaseSensitive{
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    char operator()(const char& AA) {
        switch (AA) {
            case 0:  return 'A';
            case 1:  return 'C';
            case 2:  return 'D';
            case 3:  return 'E';
            case 4:  return 'F';
            case 5:  return 'G';
            case 6:  return 'H';
            case 7:  return 'I';
            case 8:  return 'K';
            case 9:  return 'L';
            case 10: return 'M';
            case 11: return 'N';
            case 12: return 'P';
            case 13: return 'Q';
            case 14: return 'R';
            case 15: return 'S';
            case 16: return 'T';
            case 17: return 'V';
            case 18: return 'W';
            case 19: return 'Y';
            case 32: return 'a';
            case 33: return 'c';
            case 34: return 'd';
            case 35: return 'e';
            case 36: return 'f';
            case 37: return 'g';
            case 38: return 'h';
            case 39: return 'i';
            case 40: return 'k';
            case 41: return 'l';
            case 42: return 'm';
            case 43: return 'n';
            case 44: return 'p';
            case 45: return 'q';
            case 46: return 'r';
            case 47: return 's';
            case 48: return 't';
            case 49: return 'v';
            case 50: return 'w';
            case 51: return 'y';
            default: return '-';
        }
    }
};

struct CaseSensitive_to_CaseInsensitive{
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    char operator()(const char& AA) {
        return AA % 32;
    }

    //vectorized for 4 values packed in a single int
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    unsigned int operator()(const unsigned int& packed4) {
        constexpr unsigned int mod32mask = 0x1F1F1F1F;
        return packed4 & mod32mask;
    }
};

/*
    Map lower-case encoded letters to "invalid letter"
*/
struct ClampToInvalid{
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    char operator()(const char& AA) {
        return AA < 20 ? AA : 20;
    }

    //vectorized for 4 values packed in a single int
    #ifdef __CUDACC__
    __device__
    #endif
    unsigned int operator()(const unsigned int& packed4) {
        #ifdef __CUDA_ARCH__

        constexpr unsigned int mask20 = 0x14141414; // decimal 20 per byte
        return __vminu4(packed4, mask20);

        #else

        char asChar[4];
        std::memcpy(&asChar[0], &packed4, sizeof(unsigned int));
        asChar[0] = operator()(asChar[0]);
        asChar[1] = operator()(asChar[1]);
        asChar[2] = operator()(asChar[2]);
        asChar[3] = operator()(asChar[3]);

        unsigned int result;
        std::memcpy(&result, &asChar[0], sizeof(unsigned int));
        return result;

        #endif
    }

};




struct ConvertAA_20_mmseqs_to_ncbi{
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    char operator()(const char& encodedAA) {
        //             A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
        //   (NCBI)    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
        //   (mmseqs)  0 14 11  2  1 13  3  5  6  7  9  8 10  4 12 15 16 18 19 17        
        constexpr std::array<int, 21> mmseqsToNcbi{
            0, //A
            4, //C
            3, //D
            6, //E
            13, //F
            7, //G
            8, //H
            9, //I
            11, //K
            10, //L
            12, //M
            2, //N
            14, //P
            5, //Q
            1, //R
            15, //S
            16, //T
            19, //V
            17, //W
            18, //Y
            20, //else
        };

        return mmseqsToNcbi[encodedAA];
    }
};

struct ConvertAA_20_ncbi_to_mmseqs{
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    char operator()(const char& encodedAA) {
        //             A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
        //   (NCBI)    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
        //   (mmseqs)  0 14 11  2  1 13  3  5  6  7  9  8 10  4 12 15 16 18 19 17        
        constexpr std::array<int, 21> ncbiToMMseqs{
            0,
            14,
            11,
            2,
            1,
            13,
            3,
            5,
            6,
            7,
            9,
            8,
            10,
            4,
            12,
            15,
            16,
            18,
            19,
            17,
            20, //else
        };

        return ncbiToMMseqs[encodedAA];
    }
};



} //namespace cudasw4

#endif