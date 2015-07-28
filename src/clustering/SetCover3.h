//
// Created by lars on 08.06.15.
//

#ifndef MMSEQS_SETCOVER3_H
#define MMSEQS_SETCOVER3_H


#include <DBReader.h>
#include <DBWriter.h>
#include <set>
#include <list>
#include "SetElement.h"

class SetCover3 {
public:
    SetCover3(DBReader * seqDbr, DBReader * alnDbr, float seqIdThr, float coverage);

    std::list<set *>  execute();
private:
    DBReader* seqDbr;

    DBReader* alnDbr;

    float seqIdThr;

    float coverage;
//datastructures
    int * clustersizes;
    std::set<int>*orderedClustersizes;
    unsigned int maxClustersize;
    unsigned int dbSize;


//methods
    void insertCluster(int clusterid, int clustersize);

};



#endif //MMSEQS_SETCOVER3_H
