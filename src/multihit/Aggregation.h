#ifndef MMSEQS_AGGREGATION_H
#define MMSEQS_AGGREGATION_H

#include "DBReader.h"
#include "DBWriter.h"

#include <vector>
#include <map>

class Aggregation {
public:
    Aggregation(const std::string &targetDbName, const std::string &resultDbName, const std::string &outputDbName,
                unsigned int threads);

    virtual ~Aggregation();

    int run();

    virtual std::string aggregateEntry(std::vector<std::vector<std::string>> &dataToAggregate, unsigned int querySetKey, unsigned int targetSetKey) = 0;
protected:
    std::string resultDbName;
    std::string outputDbName;
    DBReader<unsigned int> *targetSetReader;
    unsigned int threads;

    void buildMap(char *data, std::map<unsigned int, std::vector<std::vector<std::string>>> &dataToAggregate);
};

#endif
