#ifndef MMSEQS_AGGREGATIONFUNCTIONS_H
#define MMSEQS_AGGREGATIONFUNCTIONS_H

#include "DBReader.h"
#include "DBWriter.h"

#include <vector>
#include <map>

class Aggregation {
public:
    Aggregation(const std::string &targetDbName, const std::string &resultDbName, const std::string &outputDbName, unsigned int threads);
    int run();

    virtual std::string aggregateEntry(std::vector<std::vector<std::string>> &dataToAggregate, unsigned int querySetKey, unsigned int targetSetKey) = 0;

protected:
    std::string resultDbName;
    std::string outputDbName;
    DBReader<unsigned int> *targetSetReader;
    unsigned int threads;

    bool buildMap(std::stringstream &data, std::map<unsigned int, std::vector<std::vector<std::string>>> &dataToAggregate);

};

class BestHitBySetFilter : public Aggregation {
public :
    BestHitBySetFilter(const std::string &targetDbName, const std::string &resultDbName, const std::string &outputDbName, bool simpleBestHitMode, unsigned int threads);
    std::string aggregateEntry(std::vector<std::vector<std::string>> &dataToAggregate, unsigned int querySetKey, unsigned int targetSetKey);

private:
    DBReader<unsigned int> *targetSizeReader;
    bool simpleBestHitMode;
};


class PvalueAggregator : public Aggregation {
public:
    PvalueAggregator(std::string queryDbName, std::string targetDbName, const std::string &resultDbName,
                     const std::string &outputDbName, float alpha, unsigned int threads);

    std::string aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate, unsigned int querySetKey, unsigned int targetSetKey);

private:
    double alpha;
    DBReader<unsigned int> *querySizeReader;
    DBReader<unsigned int> *targetSizeReader;
};


class SetSummaryAggregator : public Aggregation {
public:
    SetSummaryAggregator(const std::string &queryDbName, const std::string &targetDbName,
                              const std::string &resultDbName, const std::string &outputDbName, bool shortOutput,
                              float alpha, unsigned int threads);
    ~SetSummaryAggregator();

    std::string aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate, unsigned int querySetKey, unsigned int targetSetKey);

private:
    DBReader<unsigned int> *querySizeReader;
    DBReader<unsigned int> *targetSourceReader;
    DBReader<unsigned int> *targetSizeReader;
    float alpha;
    bool shortOutput;
};

#endif //MMSEQS_AGGREGATIONFUNCTIONS_H
