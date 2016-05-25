#ifndef MMSEQS_COMPRESSEDA3M_H
#define MMSEQS_COMPRESSEDA3M_H

#include "Matcher.h"

class DBConcat;

class CompressedA3M {
public:
    static std::string fromAlignmentResult(const std::vector<Matcher::result_t>& alignment, DBConcat* referenceDBr);
};


#endif //MMSEQS_COMPRESSEDA3M_H
