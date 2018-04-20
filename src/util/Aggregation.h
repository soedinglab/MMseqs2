//
// Created by clovisjr on 23/02/18.
//

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
    size_t targetColumn ;
    std::string inputDBname ;
    std::string outputDBname ;
    unsigned int nbrThread ;
public:
    Aggregation(std::string arg_inputDBname, std::string arg_outputDBname, unsigned int arg_nbrThread, size_t arg_targetColumn=0) ;
    bool buildMap(std::stringstream& dataStream, std::map<unsigned int, std::vector<std::vector<std::string> > > &dataToAggregate);
    void runAggregate();
    virtual std::string aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate, aggregParams* params) = 0;
};

class bestHitAggregator : public Aggregation{
private:
    size_t pValColumn ; // Field where to retrieve score values
    std::string querySetSizeName ;
    DBReader<unsigned int> *querySetSize ;
public :
    bestHitAggregator(std::string arg_inputDBname, std::string arg_outputDBname, std::string arg_querySetSizeName="",
                      unsigned int arg_nbrThread=1, size_t arg_targetColumn=10, size_t arg_scoreColumn=3,
                      DBReader<unsigned int> * querySetSize=nullptr);

    std::string aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate, aggregParams* params) override;

} ;


class pvalAggregator : public Aggregation{
private:
    size_t pValColumn ;
    std::string querySetSizeDBname ;
    std::string targetSetSizeDBname ;
    DBReader<unsigned int>* querySetSizeDB  ;
    DBReader<unsigned int>* targetSetSizeDB ;
    size_t targetGlobalSize;

public:
    pvalAggregator(std::string arg_inputDBname, std::string arg_outputDBname, unsigned int arg_nbrThread,
                   std::string arg_querySetSizeDBname, std::string arg_targetSetSizeDBname,
                   DBReader<unsigned int> * querySetSizeDB=nullptr, DBReader<unsigned int> *targetSetSizeDB=nullptr,
                   size_t arg_targetColumn=10, size_t arg_scoreColumn=3) ;

    std::string aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate, aggregParams* params) override;

};



#endif //MMSEQS_AGGREGATIONFUNCTIONS_H
