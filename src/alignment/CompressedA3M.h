#ifndef MMSEQS_COMPRESSEDA3M_H
#define MMSEQS_COMPRESSEDA3M_H

#include "Matcher.h"
#include "DBReader.h"

class DBConcat;

class CompressedA3M {
public:
    static std::string fromAlignmentResult(const std::vector<Matcher::result_t>& alignment, DBConcat& referenceDBr);

    static std::string extractA3M(const char *data, size_t data_size,
                                  DBReader<unsigned int>& sequenceReader,
                                  DBReader<unsigned int>& headerReader);
};


#endif //MMSEQS_COMPRESSEDA3M_H
