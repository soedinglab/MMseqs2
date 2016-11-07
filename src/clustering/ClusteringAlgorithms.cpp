//
// Created by lars on 08.06.15.
//
#include <sys/time.h>
#include "ClusteringAlgorithms.h"
#include "Util.h"
#include "Debug.h"
#include "AlignmentSymmetry.h"
#include <queue>
#include <algorithm>
#include <new>
#include <climits>

ClusteringAlgorithms::ClusteringAlgorithms(DBReader<unsigned int>* seqDbr, DBReader<unsigned int>* alnDbr,
                                           int threads, int scoretype, int maxiterations){
    this->seqDbr=seqDbr;
    if(seqDbr->getSize() != alnDbr->getSize()){
        Debug(Debug::ERROR) << "Sequence db size != result db size\n";
        EXIT(EXIT_FAILURE);
    }
    this->alnDbr=alnDbr;
    this->dbSize=alnDbr->getSize();
    this->threads=threads;
    this->scoretype=scoretype;
    this->maxiterations=maxiterations;
    ///time
    this->clustersizes=new int[dbSize];
    this->sets = new set[dbSize];
    std::fill_n(clustersizes, dbSize, 0);
}

ClusteringAlgorithms::~ClusteringAlgorithms(){
    delete [] clustersizes;
    delete [] sets;
}

std::map<unsigned int, std::vector<unsigned int>>  ClusteringAlgorithms::execute(int mode) {
    memset(sets, 0, sizeof(set)*(dbSize));
    const char * data = alnDbr->getData();
    const size_t dataSize = alnDbr->getDataSize();
    // init data
    const size_t elementCount = Util::countLines(data, dataSize);
    unsigned int * elements = new(std::nothrow) unsigned int[elementCount];
    Util::checkAllocation(elements, "Could not allocate elements memory in ClusteringAlgorithms::execute");
    unsigned int ** elementLookupTable = new(std::nothrow) unsigned int*[dbSize];
    Util::checkAllocation(elementLookupTable, "Could not allocate elementLookupTable memory in ClusteringAlgorithms::execute");
    unsigned short **scoreLookupTable = new(std::nothrow) unsigned short *[dbSize];
    Util::checkAllocation(scoreLookupTable, "Could not allocate scoreLookupTable memory in ClusteringAlgorithms::execute");
    unsigned short *score = NULL;
    size_t *elementOffsets = new(std::nothrow) size_t[dbSize + 1];
    Util::checkAllocation(elementOffsets, "Could not allocate elementOffsets memory in ClusteringAlgorithms::execute");
    elementOffsets[dbSize] = 0;

    unsigned int *assignedcluster = new(std::nothrow) unsigned int[dbSize];
    Util::checkAllocation(assignedcluster, "Could not allocate assignedcluster memory in ClusteringAlgorithms::execute");
    std::fill_n(assignedcluster, dbSize, UINT_MAX);
    short *bestscore = new(std::nothrow) short[dbSize];
    Util::checkAllocation(bestscore, "Could not allocate bestscore memory in ClusteringAlgorithms::execute");
    std::fill_n(bestscore, dbSize, SHRT_MIN);
    readInClusterData(elementLookupTable, elements, scoreLookupTable, score, elementOffsets, elementCount);
    //time
    if (mode==2){
        greedyIncremental(elementLookupTable, elementOffsets,
                          dbSize, assignedcluster);
    }else {
        ClusteringAlgorithms::initClustersizes();
        if (mode == 1) {
            setCover(elementLookupTable, scoreLookupTable, assignedcluster, bestscore, elementOffsets);
        } else if (mode == 3) {
            Debug(Debug::INFO) << "connected component mode" << "\n";
            for (int cl_size = dbSize - 1; cl_size >= 0; cl_size--) {
                unsigned int representative = sorted_clustersizes[cl_size];
                if (assignedcluster[representative] == UINT_MAX) {
                    assignedcluster[representative] = representative;
                    std::queue<int> myqueue;
                    myqueue.push(representative);
                    std::queue<int> iterationcutoffs;
                    iterationcutoffs.push(0);
                    //delete clusters of members;
                    while (!myqueue.empty()) {
                        int currentid = myqueue.front();
                        int iterationcutoff = iterationcutoffs.front();
                        assignedcluster[currentid] = representative;
                        myqueue.pop();
                        iterationcutoffs.pop();
                        size_t elementSize = (elementOffsets[currentid + 1] - elementOffsets[currentid]);
                        for (size_t elementId = 0; elementId < elementSize; elementId++) {
                            unsigned int elementtodelete = elementLookupTable[currentid][elementId];
                            if (assignedcluster[elementtodelete] == UINT_MAX && iterationcutoff < maxiterations) {
                                myqueue.push(elementtodelete);
                                iterationcutoffs.push((iterationcutoff + 1));
                            }
                            assignedcluster[elementtodelete] = representative;
                        }
                    }

                }
            }
        }
        //delete unnecessary datastructures
        delete [] sorted_clustersizes;
        delete [] clusterid_to_arrayposition;
        delete [] borders_of_set;
    }

    delete [] elementLookupTable;
    delete [] elements;
    delete [] elementOffsets;
    delete [] scoreLookupTable;
    delete [] bestscore;
    delete [] score;

    std::map<unsigned int, std::vector<unsigned int>> retMap;
    for(size_t i = 0; i < dbSize; i++) {
        if(assignedcluster[i] == UINT_MAX){
            Debug(Debug::ERROR) << "there must be an error: " << seqDbr->getDbKey(i) <<
            " is not assigned to a cluster\n";
            continue;
        }

        // make sure the representative is always the first entry
        if(retMap.find(assignedcluster[i]) == retMap.end()) {
            retMap[assignedcluster[i]].push_back(assignedcluster[i]);
        }

        // and don't add it a second time
        if(i != assignedcluster[i]) {
            retMap[assignedcluster[i]].push_back(i);
        }
    }
    delete [] assignedcluster;
    return retMap;
}

void ClusteringAlgorithms::initClustersizes(){
    int * setsize_abundance=new int[maxClustersize+1];

    std::fill_n(setsize_abundance,maxClustersize+1,0);
    //count how often a set size occurs
    for (unsigned int i = 0; i < dbSize; ++i) {
        setsize_abundance[clustersizes[i]]++;
    }
    //compute offsets
    borders_of_set= new unsigned int[maxClustersize+1];
    borders_of_set[0] = 0;
    for (unsigned int i = 1; i < maxClustersize+1; ++i) {
        borders_of_set[i] = borders_of_set[i-1] + setsize_abundance[i-1];
    }
    //fill array
    sorted_clustersizes = new(std::nothrow)  unsigned int[dbSize + 1];
    Util::checkAllocation(sorted_clustersizes, "Could not allocate sorted_clustersizes memory in ClusteringAlgorithms::initClustersizes");

    std::fill_n(sorted_clustersizes, dbSize+1, 0);
    clusterid_to_arrayposition = new(std::nothrow)  unsigned int[dbSize + 1];
    Util::checkAllocation(clusterid_to_arrayposition, "Could not allocate sorted_clustersizes memory in ClusteringAlgorithms::initClustersizes");

    std::fill_n(clusterid_to_arrayposition, dbSize + 1, 0);
    //reuse setsize_abundance as offset counter
    std::fill_n(setsize_abundance, maxClustersize + 1, 0);
    for (unsigned int i = 0; i < dbSize; ++i) {
        int position=borders_of_set[clustersizes[i]] + setsize_abundance[clustersizes[i]];
        sorted_clustersizes[position] = i;
        clusterid_to_arrayposition[i] = position;
        setsize_abundance[clustersizes[i]]++;
    }
    delete [] setsize_abundance;
}


void ClusteringAlgorithms::removeClustersize(int clusterid){
    clustersizes[clusterid]=0;
    sorted_clustersizes[clusterid_to_arrayposition[clusterid]] = UINT_MAX;
    clusterid_to_arrayposition[clusterid]=UINT_MAX;
}

void ClusteringAlgorithms::decreaseClustersize(int clusterid){
    const unsigned int oldposition=clusterid_to_arrayposition[clusterid];
    const unsigned int newposition=borders_of_set[clustersizes[clusterid]];
    const unsigned int swapid=sorted_clustersizes[newposition];
    if(swapid != UINT_MAX){
        clusterid_to_arrayposition[swapid]=oldposition;
    }
    sorted_clustersizes[oldposition]=swapid;

    sorted_clustersizes[newposition]=clusterid;
    clusterid_to_arrayposition[clusterid]=newposition;
    borders_of_set[clustersizes[clusterid]]++;
    clustersizes[clusterid]--;
}

void ClusteringAlgorithms::setCover(unsigned int **elementLookupTable, unsigned short ** elementScoreLookupTable,
                                    unsigned int *assignedcluster, short *bestscore, size_t *newElementOffsets) {
    for (int cl_size = dbSize - 1; cl_size >= 0; cl_size--) {
        const unsigned int representative = sorted_clustersizes[cl_size];
        if (representative == UINT_MAX) {
            continue;
        }
//          Debug(Debug::INFO)<<alnDbr->getDbKey(representative)<<"\n";
        removeClustersize(representative);
        assignedcluster[representative] = representative;

        //delete clusters of members;
        size_t elementSize = (newElementOffsets[representative + 1] - newElementOffsets[representative]);
        for (size_t elementId = 0; elementId < elementSize; elementId++) {
            const unsigned int elementtodelete = elementLookupTable[representative][elementId];
            // float seqId = elementScoreTable[representative][elementId];
            const short seqId = elementScoreLookupTable[representative][elementId];
            //  Debug(Debug::INFO)<<seqId<<"\t"<<bestscore[elementtodelete]<<"\n";
            // becareful of this criteria
            if (seqId > bestscore[elementtodelete]) {
                assignedcluster[elementtodelete] = representative;
                bestscore[elementtodelete] = seqId;
            }
            //Debug(Debug::INFO)<<bestscore[elementtodelete]<<"\n";
            if (elementtodelete == representative) {
                continue;
            }
            if (clustersizes[elementtodelete] < 1) {
                continue;
            }
            removeClustersize(elementtodelete);
        }

        for (size_t elementId = 0; elementId < elementSize; elementId++) {
            bool representativefound = false;
            const unsigned int elementtodelete = elementLookupTable[representative][elementId];
            const unsigned int currElementSize = (newElementOffsets[elementtodelete + 1] -
                                                  newElementOffsets[elementtodelete]);
            if (elementtodelete == representative) {
                clustersizes[elementtodelete] = -1;
                continue;
            }
            if (clustersizes[elementtodelete] < 0) {
                continue;
            }
            clustersizes[elementtodelete] = -1;
            //decrease clustersize of sets that contain the element
            for (size_t elementId2 = 0; elementId2 < currElementSize; elementId2++) {
                const unsigned int elementtodecrease = elementLookupTable[elementtodelete][elementId2];
                if (representative == elementtodecrease) {
                    representativefound = true;
                }
                if (clustersizes[elementtodecrease] == 1) {
                    Debug(Debug::ERROR) << "there must be an error: " << seqDbr->getDbKey(elementtodelete) <<
                    " deleted from " << seqDbr->getDbKey(elementtodecrease) <<
                    " that now is empty, but not assigned to a cluster\n";
                } else if (clustersizes[elementtodecrease] > 0) {
                    decreaseClustersize(elementtodecrease);
                }
            }
            if (!representativefound) {
                Debug(Debug::ERROR) << "error with cluster:\t" << seqDbr->getDbKey(representative) <<
                "\tis not contained in set:\t" << seqDbr->getDbKey(elementtodelete) << ".\n";
            }
        }
    }
}

void ClusteringAlgorithms::greedyIncremental(unsigned int **elementLookupTable, size_t *elementOffsets,
                                             size_t n, unsigned int *assignedcluster) {
    for(size_t i = 0; i < n; i++) {
        // seqDbr is descending sorted by length
        // the assumption is that clustering is B -> B (not A -> B)
        Debug::printProgress(i);
        if(assignedcluster[i] == UINT_MAX){
            size_t elementSize = (elementOffsets[i + 1] - elementOffsets[i]);
            for (size_t elementId = 0; elementId < elementSize; elementId++) {
                const unsigned int currElm = elementLookupTable[i][elementId];
                if(assignedcluster[currElm] == currElm){
                    assignedcluster[i] = currElm;
                    break;
                }
            }
            if(assignedcluster[i] == UINT_MAX) {
                assignedcluster[i] = i;
            }
        }
    }
}

void ClusteringAlgorithms::readInClusterData(unsigned int **elementLookupTable, unsigned int *&elements,
                                             unsigned short **scoreLookupTable, unsigned short *&scores,
                                             size_t *elementOffsets, size_t totalElementCount) {
    //time
    struct timeval start, end;
    gettimeofday(&start, NULL);
#pragma omp parallel for schedule(dynamic, 1000)
    for(size_t i = 0; i < dbSize; i++) {
        const unsigned int clusterId = seqDbr->getDbKey(i);
        const size_t alnId = alnDbr->getId(clusterId);
        const char *data = alnDbr->getData(alnId);
        const size_t dataSize = alnDbr->getSeqLens(alnId);
        elementOffsets[i] = Util::countLines(data, dataSize);
    }

    // make offset table
    AlignmentSymmetry::computeOffsetFromCounts(elementOffsets, dbSize);
    // set element edge pointers by using the offset table
    AlignmentSymmetry::setupPointers<unsigned int>(elements, elementLookupTable, elementOffsets, dbSize,
                                                   totalElementCount);
    // fill elements
    AlignmentSymmetry::readInData(alnDbr, seqDbr, elementLookupTable, NULL, 0, elementOffsets);
    Debug(Debug::WARNING) << "\nSort entries.\n";
    AlignmentSymmetry::sortElements(elementLookupTable, elementOffsets, dbSize);
    Debug(Debug::WARNING) << "\nFind missing connections.\n";

    size_t *newElementOffsets = new size_t[dbSize + 1];
    memcpy(newElementOffsets, elementOffsets, sizeof(size_t) * (dbSize + 1));

    // findMissingLinks detects new possible connections and updates the elementOffsets with new sizes
    const size_t symmetricElementCount = AlignmentSymmetry::findMissingLinks(elementLookupTable,
                                                                       newElementOffsets, dbSize,
                                                                       threads);
    // resize elements
    delete[] elements;
    elements = new(std::nothrow) unsigned int[symmetricElementCount];
    Util::checkAllocation(elements, "Could not allocate elements memory in readInClusterData");
    std::fill_n(elements, symmetricElementCount, UINT_MAX);
    // init score vector
    scores = new(std::nothrow) unsigned short[symmetricElementCount];
    Util::checkAllocation(scores, "Could not allocate scores memory in readInClusterData");
    std::fill_n(scores, symmetricElementCount, 0);
    Debug(Debug::WARNING) << "\nFound " << symmetricElementCount - totalElementCount << " new connections.\n";
    AlignmentSymmetry::setupPointers<unsigned int>  (elements, elementLookupTable, newElementOffsets, dbSize, symmetricElementCount);
    AlignmentSymmetry::setupPointers<unsigned short>(scores, scoreLookupTable, newElementOffsets, dbSize, symmetricElementCount);
    //time
    Debug(Debug::WARNING) << "\nReconstruct initial order.\n";
    alnDbr->remapData(); // need to free memory
    AlignmentSymmetry::readInData(alnDbr, seqDbr, elementLookupTable, scoreLookupTable, scoretype, elementOffsets);
    alnDbr->remapData(); // need to free memory
    Debug(Debug::WARNING) << "\nAdd missing connections.\n";
    AlignmentSymmetry::addMissingLinks(elementLookupTable, elementOffsets, newElementOffsets, dbSize, scoreLookupTable);
    maxClustersize = 0;
    for (size_t i = 0; i < dbSize; i++) {
        size_t elementCount = newElementOffsets[i + 1] - newElementOffsets[i];
        maxClustersize = std::max((unsigned int) elementCount, maxClustersize);
        clustersizes[i] = elementCount;
    }

    memcpy(elementOffsets, newElementOffsets, sizeof(size_t) * (dbSize + 1));
    delete [] newElementOffsets;
    //time
    gettimeofday(&end, NULL);
    size_t sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for Read in: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
}