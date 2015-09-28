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
    SetCover3(DBReader * seqDbr, DBReader * alnDbr, float seqIdThr, float coverage,int threads,int scoretype, int maxiterations);

    std::list<set *>  execute(int mode);
private:
    DBReader* seqDbr;

    DBReader* alnDbr;

    float seqIdThr;

    float coverage;

    int threads;
    int scoretype;
//datastructures
    int * clustersizes;
    unsigned int maxClustersize;
    unsigned int dbSize;

    //
    int* sorted_clustersizes;
    int* clusterid_to_arrayposition;
    int* borders_of_set;

//methods

    void initClustersizes();

    void removeClustersize(int clusterid);

    void decreaseClustersize(int clusterid);
//for connected component
    int maxiterations;
};



#endif //MMSEQS_SETCOVER3_H
