#ifndef MMSEQS_AGGREGATIONFUNCTIONS_H
#define MMSEQS_AGGREGATIONFUNCTIONS_H

#include <vector>
#include <map>
#include "DBReader.h"
#include "DBWriter.h"

typedef struct {
    unsigned int querySetKey;
    unsigned int targetSetKey;
} aggregParams;


class Aggregation {
protected:
    size_t targetSetColumn;
    std::string inputDBname;
    std::string outputDBname;
    unsigned int nbrThread;
public:
    Aggregation(std::string argInputDBname, std::string argOutputDBname, unsigned int argNbrThread, size_t argTargetSetColumn);
    bool buildMap(std::stringstream& dataStream, std::map<unsigned int, std::vector<std::vector<std::string> > > &dataToAggregate);


    void runAggregate();
    virtual std::string aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate, aggregParams* params) = 0;
};

class BestHitAggregator : public Aggregation{
private:
    size_t scoreColumn; // Field where to retrieve score values
    std::string targetSetSizeName;
    DBReader<unsigned int> *targetSetSizeDB;
    bool simpleBestHitMode;
public :
    BestHitAggregator(std::string argInputDBname, std::string argOutputDBname, std::string argTargetSetSizeName,
                      size_t argTargetColumn, unsigned int argNbrThread,  size_t argScoreColumn=1, bool argSimpleBestHitMode = false); // argScoreColumn=3 for evalue

    std::string aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate, aggregParams* params) override;
};


class PvalAggregator : public Aggregation{
private:
    double alpha;
    size_t scoreColumn;
    std::string querySetSizeDBname;
    std::string targetSetSizeDBname;
    DBReader<unsigned int>* querySetSizeDB ;
    DBReader<unsigned int>* targetSetSizeDB ;
public:
    PvalAggregator(std::string argInputDBname, std::string argOutputDBname, unsigned int arg_nbrThread,
                   std::string argQuerySetSizeDBname, std::string argTargetSetSizeDBname, size_t argTargetSetColumn, float alpha, size_t argScoreColumn=1);

    std::string aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate, aggregParams* params) override;

};


class GetHitDistance : public Aggregation {
private:

    std::string querySetSizeDBname;
    std::string targetSetSizeDBname;
    DBReader<unsigned int>* querySetSizeDB;
    DBReader<unsigned int>* targetSetGenomes;
    DBReader<unsigned int>* targetSetSizeDB;
    float alpha;
public:
    GetHitDistance(std::string arg_inputDBname, std::string arg_outputDBname, std::string argTargetSetSizeDB, std::string argQuerySetSizeDB,std::string argTargetGenomeDB,
                             unsigned int arg_nbrThread, float argAlpha, size_t arg_targetColumn=12);
    ~GetHitDistance();
    std::string aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate, aggregParams* params) override;
};

#endif //MMSEQS_AGGREGATIONFUNCTIONS_H
