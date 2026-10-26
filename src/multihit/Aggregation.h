#ifndef MMSEQS_AGGREGATION_H
#define MMSEQS_AGGREGATION_H

#include "DBReader.h"
#include "DBWriter.h"

#include <vector>
#include <map>

class Aggregation {
public:
    Aggregation(const std::string &targetDbName, const std::string &resultDbName, const std::string &outputDbName,
                unsigned int threads, unsigned int compressed);

    virtual ~Aggregation();

    int run();
    virtual void prepareInput(DBKeyType querySetKey, unsigned int thread_idx) = 0;
    virtual std::string aggregateEntry(std::vector<std::vector<std::string>> &dataToAggregate, DBKeyType querySetKey, DBKeyType targetSetKey, unsigned int thread_idx) = 0;

protected:
    std::string resultDbName;
    std::string outputDbName;
    DBReader<DBKeyType> *targetSetReader;
    unsigned int threads;
    unsigned int compressed;

    void buildMap(char *data, int thread_idx, std::map<DBKeyType, std::vector<std::vector<std::string>>> &dataToAggregate);
};

#endif
