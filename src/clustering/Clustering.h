#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <list>
#include <string>
#include <map>

#include "DBReader.h"
#include "DBWriter.h"
#include "SetElement.h"

class Clustering {
public:
    Clustering(const std::string &seqDB, const std::string &seqDBIndex,
               const std::string &alnResultsDB, const std::string &alnResultsDBIndex,
               const std::string &outDB, const std::string &outDBIndex,
               int validateClustering, unsigned int maxIteration,
               int similarityScoreType, int threads);

    void run(int mode);


    ~Clustering();

private:
    // check if every element is member in only one cluster
    bool validate_result(std::list<set *> *ret, unsigned int uniqu_element_count);

    void writeData(DBWriter *dbw, const std::map<unsigned int, std::vector<unsigned int>> &ret);

    DBReader<unsigned int> *seqDbr;
    DBReader<unsigned int> *alnDbr;

    int validate;

    //values for affinity clustering
    unsigned int maxIteration;
    int similarityScoreType;

    int threads;
    std::string outDB;
    std::string outDBIndex;
};

#endif
