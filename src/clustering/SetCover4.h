//
// Created by lars on 28.07.15.
//

//
// Created by lars on 08.06.15.
//

#ifndef MMSEQS_SETCOVER4_H
#define MMSEQS_SETCOVER4_H


#include <DBReader.h>
#include <DBWriter.h>
#include <set>
#include <list>
#include "SetElement.h"

class SetCover4 {
public:
    SetCover4(DBReader * seqDbr, DBReader * alnDbr, float seqIdThr, float coverage);

    std::list<set *>  execute();
private:
    DBReader* seqDbr;

    DBReader* alnDbr;

    float seqIdThr;

    float coverage;
//datastructures
    int * clustersizes;
    float *clusterweights;
    std::set<int>*orderedClustersizes;
    std::set<int>*orderedClusterweights;
    unsigned int maxClustersize;
    float maxClusterweight;
    unsigned int dbSize;
    int spreadingfactor=10;
    int weightslength;
    float offset=0.8;

//methods
    void insertCluster(int clusterid, int clustersize);

    void insertClusterWeight(int clusterid, float clusterweight);
    void eraseClusterWeight(int clusterid, float clusterweight);
};



#endif //MMSEQS_SETCOVER4_H

