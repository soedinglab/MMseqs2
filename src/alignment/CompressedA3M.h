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
                                  DBReader<unsigned int>& headerReader, int thread_idx);

    static void extractMatcherResults(unsigned int &key, std::vector<Matcher::result_t> &results,
                                      const char *data, size_t dataSize, DBReader<unsigned int>& sequenceReader, bool skipFirst);
};

#endif
