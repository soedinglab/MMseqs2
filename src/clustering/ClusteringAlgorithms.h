//
// Created by lars on 08.06.15.
//

#ifndef MMSEQS_CLUSTERINGALGORITHMS_H
#define MMSEQS_CLUSTERINGALGORITHMS_H

#include <set>
#include <list>
#include <vector>
#include <map>

#include "DBReader.h"
#include "SetElement.h"

class ClusteringAlgorithms {
public:
    ClusteringAlgorithms(DBReader<unsigned int>* seqDbr, DBReader<unsigned int>* alnDbr, int threads,int scoretype, int maxiterations);
    ~ClusteringAlgorithms();
    std::map<unsigned int, std::vector<unsigned int>> execute(int mode);
private:
    DBReader<unsigned int>* seqDbr;

    DBReader<unsigned int>* alnDbr;


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
    // all results sets
    set *sets;

    void setCover(unsigned int **elementLookup, unsigned short ** elementScoreLookupTable,  int *assignedcluster, short *bestscore, size_t *offsets);

    void greedyIncremental(unsigned int **elementLookupTable, size_t *elementOffsets,
                           unsigned short **elementScoreLookupTable, size_t elementCount,
                           size_t n, int *assignedcluster, unsigned short *scoreelements) ;

    void readInClusterData(unsigned int **elementLookupTable, unsigned int *&elements,
                           unsigned short ** elementScoreLookupTable, unsigned short *&scoreelements,
                           size_t *elementOffsets, size_t elementCount)  ;

};



#endif //MMSEQS_CLUSTERINGALGORITHMS_H
