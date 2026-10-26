//
// Created by lars on 08.06.15.
//

#ifndef MMSEQS_CLUSTERINGALGORITHMS_H
#define MMSEQS_CLUSTERINGALGORITHMS_H

#include <set>
#include <list>
#include <vector>
#include <unordered_map>

#include "DBReader.h"

class ClusteringAlgorithms {
public:
    ClusteringAlgorithms(DBReader<DBKeyType>* seqDbr, DBReader<DBKeyType>* alnDbr, int threads,int scoretype, int maxiterations, DBKeyType *keyToSet, size_t *sourceOffsets, DBKeyType **sourceLookupTable, DBKeyType *sourceList, size_t sourceLen, bool needSET);
    ~ClusteringAlgorithms();
    std::pair<DBKeyType, DBKeyType> * execute(int mode);
private:
    DBReader<DBKeyType>* seqDbr;

    DBReader<DBKeyType>* alnDbr;

    bool needSET;
    int threads;
    int scoretype;
//datastructures
    size_t maxClustersize;
    size_t dbSize;
    // signed (uses a -1 sentinel); 4 bytes by default, 8 under MMSEQS_INT64_IDS
#ifdef MMSEQS_INT64_IDS
    int64_t * clustersizes;
#else
    int32_t * clustersizes;
#endif
    DBLocalId* sorted_clustersizes;
    DBLocalId* clusterid_to_arrayposition;   // position in sorted_clustersizes (<= dbSize)
    size_t* borders_of_set;                  // small (maxClustersize+1); kept size_t
    DBKeyType* keyToSet;
    size_t* sourceOffsets;
    DBKeyType** sourceLookupTable;
    DBKeyType* sourceList;
    size_t sourceLen;

//methods

    void initClustersizes();

    void removeClustersize(DBLocalId clusterid);

    void decreaseClustersize(DBLocalId clusterid);
//for connected component
    int maxiterations;


    void setCover(DBLocalId **elementLookup, unsigned short ** elementScoreLookupTable,
                  DBLocalId *assignedcluster, short *bestscore, size_t *offsets);

    void greedyIncremental(DBLocalId **elementLookupTable, size_t *elementOffsets,
                           size_t n, DBLocalId *assignedcluster) ;


    void greedyIncrementalLowMem(DBLocalId *assignedcluster) ;


    void readInClusterData(DBLocalId **elementLookupTable, DBLocalId *&elements,
                           unsigned short **scoreLookupTable, unsigned short *&scores,
                           size_t *elementOffsets, size_t totalElementCount)  ;

};



#endif //MMSEQS_CLUSTERINGALGORITHMS_H
