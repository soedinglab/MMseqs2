#ifndef SCOREMATRIX_H
#define SCOREMATRIX_H

#include <cstddef>
#include <cstdlib>
#include <cstring>

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

    static ScoreMatrix unserialize(const char *data, size_t alphabetSize, size_t kmerSize);

    static ScoreMatrix unserializeCopy(const char *data, size_t alphabetSize, size_t kmerSize);
};

#endif
