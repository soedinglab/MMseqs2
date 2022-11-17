#ifndef RESULT2PROFILE_H
#define RESULT2PROFILE_H

#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"

#include <unordered_map>

class StatsComputer {
public:
    StatsComputer(const Parameters &par);
    ~StatsComputer();

    int run();
private:
    int stat;

    std::string queryDb;
    std::string queryDbIndex;

    std::string targetDb;
    std::string targetDbIndex;

    const bool tsvOut;

    DBReader<unsigned int> *resultReader;
    DBWriter *statWriter;

    int threads;

    template<typename T>
    struct PerSequence {
        typedef T(*type)(const char *);
    };

    template<typename T>
    int sequenceWise(typename PerSequence<T>::type call, bool onlyResultDb = false);

    int countNumberOfLines();
    int meanValue();
    int minValue();
    int maxValue();
    int sumValue();
};


#endif
