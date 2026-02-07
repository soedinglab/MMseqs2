#include "ScoreMatrix.h"

#include <cmath>
#include <cstring>

#include "simd.h"

ScoreMatrix ScoreMatrix::unserialize(const char *data, size_t alphabetSize, size_t kmerSize) {
    size_t size = pow(alphabetSize, kmerSize);
    size_t row_size = size / MAX_ALIGN_INT;
    row_size = (row_size + 1) * MAX_ALIGN_INT; // for SIMD memory alignment

    size_t scoreSize = size * row_size * sizeof(short);
    return ScoreMatrix((short *) data, (unsigned int *) (data + scoreSize), size, row_size);
}

ScoreMatrix ScoreMatrix::unserializeCopy(const char *data, size_t alphabetSize, size_t kmerSize) {
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
