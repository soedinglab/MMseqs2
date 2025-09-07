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
    ClusteringAlgorithms(DBReader<KeyType>* seqDbr, DBReader<KeyType>* alnDbr, int threads, int scoretype, int maxiterations,
                         KeyType *keyToSet, size_t *sourceOffsets, KeyType **sourceLookupTable,
                         KeyType *sourceList, KeyType sourceLen, bool needSET);
    ~ClusteringAlgorithms();
    std::pair<KeyType, KeyType> * execute(int mode);
private:
    DBReader<KeyType>* seqDbr;

    DBReader<KeyType>* alnDbr;

    bool needSET;
    int threads;
    int scoretype;
//datastructures
    KeyType maxClustersize;
    KeyType dbSize;
    int * clustersizes;
    KeyType* sorted_clustersizes;
    KeyType* clusterid_to_arrayposition;
    KeyType* borders_of_set;
    KeyType* keyToSet;
    size_t* sourceOffsets;
    KeyType** sourceLookupTable;
    KeyType* sourceList;
    KeyType sourceLen;

//methods

    void initClustersizes();

    void removeClustersize(KeyType clusterid);

    void decreaseClustersize(KeyType clusterid);
//for connected component
    int maxiterations;


    void setCover(KeyType **elementLookup, unsigned short ** elementScoreLookupTable,
                  KeyType *assignedcluster, short *bestscore, size_t *offsets);

    void greedyIncremental(KeyType **elementLookupTable, size_t *elementOffsets,
                           size_t n, KeyType *assignedcluster) ;


    void greedyIncrementalLowMem(KeyType *assignedcluster) ;


    void readInClusterData(KeyType **elementLookupTable, KeyType *&elements,
                           unsigned short **scoreLookupTable, unsigned short *&scores,
                           size_t *elementOffsets, size_t totalElementCount)  ;

};



#endif //MMSEQS_CLUSTERINGALGORITHMS_H
