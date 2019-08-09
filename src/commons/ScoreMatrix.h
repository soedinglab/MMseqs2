#ifndef SCOREMATRIX_H
#define SCOREMATRIX_H

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <simd.h>
#include <cmath>

// <match score, k-mer index>
struct ScoreMatrix {
    size_t elementSize;
    size_t rowSize;
    short *score;
    unsigned int *index;

    ScoreMatrix() : elementSize(0), rowSize(0), score(NULL), index(NULL) {}

    ScoreMatrix(short *scoreMatrix, unsigned int *indexMatrix, size_t elementSize, size_t rowSize) :
            elementSize(elementSize), rowSize(rowSize), score(scoreMatrix), index(indexMatrix) {}

    bool isValid() {
        return score != NULL && index != NULL;
    }

    static size_t size(const ScoreMatrix &mat) {
        size_t memSize = (mat.elementSize * mat.rowSize) * (sizeof(short) + sizeof(unsigned int));
        return memSize;
    }

    static char *serialize(const ScoreMatrix &mat) {
        char *data = (char *) malloc(size(mat));
        char *p = data;
        size_t scoreSize = mat.elementSize * mat.rowSize * sizeof(short);
        memcpy(p, mat.score, scoreSize);
        p += scoreSize;
        size_t indexSize = mat.elementSize * mat.rowSize * sizeof(unsigned int);
        memcpy(p, mat.index, indexSize);

        return data;
    }

    static ScoreMatrix unserialize(const char *data, size_t alphabetSize, size_t kmerSize) {
        size_t size = pow(alphabetSize, kmerSize);
        size_t row_size = size / MAX_ALIGN_INT;
        row_size = (row_size + 1) * MAX_ALIGN_INT; // for SIMD memory alignment

        size_t scoreSize = size * row_size * sizeof(short);
        return ScoreMatrix((short *) data, (unsigned int *) (data + scoreSize), size, row_size);
    }

    static ScoreMatrix unserializeCopy(const char *data, size_t alphabetSize, size_t kmerSize) {
        size_t size = pow(alphabetSize, kmerSize);
        size_t row_size = size / MAX_ALIGN_INT;
        row_size = (row_size + 1) * MAX_ALIGN_INT; // for SIMD memory alignment

        size_t scoreSize = size * row_size * sizeof(short);
        short * score = (short *) mem_align(MAX_ALIGN_INT, (size * (row_size)) * sizeof(short));
        unsigned int * index = (unsigned int *)mem_align(MAX_ALIGN_INT, (size * (row_size)) * sizeof(unsigned int));
        memcpy(score, data, scoreSize);
        memcpy(index, data + scoreSize, size * row_size * sizeof(unsigned int));
        return ScoreMatrix(score, index, size, row_size);
    }
};

#endif
