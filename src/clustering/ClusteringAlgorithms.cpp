//
// Created by lars on 08.06.15.
//

#include <sys/time.h>
#include "Log.h"
#include "ClusteringAlgorithms.h"
#include "Util.h"
#include "Debug.h"
#include "AlignmentSymmetry.h"
#include <queue>
#include <algorithm>


ClusteringAlgorithms::ClusteringAlgorithms(DBReader<unsigned int>* seqDbr, DBReader<unsigned int>* alnDbr,int threads, int scoretype, int maxiterations){
    this->seqDbr=seqDbr;
    this->alnDbr=alnDbr;
    this->dbSize=alnDbr->getSize();
    this->threads=threads;
    this->scoretype=scoretype;
    this->maxiterations=maxiterations;
    ///time
    this->clustersizes=new int[dbSize];
    this->sets = new set[dbSize];

    std::fill_n(clustersizes,dbSize,0);
}

ClusteringAlgorithms::~ClusteringAlgorithms(){
    delete [] clustersizes;
    delete [] sets;
}

std::map<unsigned int, std::vector<unsigned int>>  ClusteringAlgorithms::execute(int mode) {

    memset(sets, 0, sizeof(set)*(dbSize));

    const char * data = alnDbr->getData();
    size_t dataSize = alnDbr->getDataSize();
    // init data
    size_t elementCount = Util::countLines(data, dataSize);
    unsigned int * elements = new unsigned int[elementCount];
    unsigned int ** elementLookupTable = new unsigned int*[dbSize];
    unsigned short **elementScoreLookupTable = new unsigned short *[dbSize];
    unsigned short *scoreelements = new unsigned short[elementCount];;

    size_t *elementOffsets = new size_t[dbSize + 1];
    elementOffsets[dbSize]=0;


    std::list<set *> result;
    size_t n = seqDbr->getSize();
    int *assignedcluster = new int[n];
    std::fill_n(assignedcluster, n, -1);
    short *bestscore = new short[n];
    std::fill_n(bestscore, n, -10);
    readInClusterData(elementLookupTable, elements, elementScoreLookupTable, scoreelements, elementOffsets, elementCount);
    //time
    if (mode==2){
        greedyIncremental(elementLookupTable, elementOffsets, elementScoreLookupTable,
                          elementCount, n, assignedcluster, scoreelements);
    }else {
        ClusteringAlgorithms::initClustersizes();

        //time
        //delete from beginning
        if (mode == 1) {
            setCover(elementLookupTable, elementScoreLookupTable, assignedcluster, bestscore, elementOffsets);
        } else if (mode == 3) {
            Debug(Debug::INFO) << "connected component mode" << "\n";
            for (int cl_size = dbSize - 1; cl_size >= 0; cl_size--) {
                int representative = sorted_clustersizes[cl_size];
                if (assignedcluster[representative] == -1) {
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
                            if (assignedcluster[elementtodelete] == -1 && iterationcutoff < maxiterations) {
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
    delete [] bestscore;

    std::map<unsigned int, std::vector<unsigned int>> retMap;
    for(size_t i = 0; i < n; i++) {
        if(assignedcluster[i] == -1){
            Debug(Debug::ERROR) << "there must be an error: " << seqDbr->getDbKey(i) <<
            " is not assigned to a cluster\n";
            continue;
        }
        retMap[assignedcluster[i]].push_back(i);
    }

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
    borders_of_set= new int [maxClustersize+1];
    borders_of_set[0]=0;
    for (unsigned int i = 1; i < maxClustersize+1; ++i) {
        borders_of_set[i]=borders_of_set[i-1]+setsize_abundance[i-1];
    }
    //fill array
    sorted_clustersizes =new int [dbSize+1];
    std::fill_n(sorted_clustersizes,dbSize+1,0);
    clusterid_to_arrayposition=new int[dbSize+1];
    std::fill_n(clusterid_to_arrayposition,dbSize+1,0);
    //reuse setsize_abundance as offset counter
    std::fill_n(setsize_abundance,maxClustersize+1,0);
    for (unsigned int i = 0; i < dbSize; ++i) {
        int position=borders_of_set[clustersizes[i]]+setsize_abundance[clustersizes[i]];
        sorted_clustersizes[position]=i;
        clusterid_to_arrayposition[i]=position;
        setsize_abundance[clustersizes[i]]++;
    }
    delete [] setsize_abundance;
}


void ClusteringAlgorithms::removeClustersize(int clusterid){
    clustersizes[clusterid]=0;
    sorted_clustersizes[clusterid_to_arrayposition[clusterid]]=-1;
    clusterid_to_arrayposition[clusterid]=-1;
}

void ClusteringAlgorithms::decreaseClustersize(int clusterid){
    int oldposition=clusterid_to_arrayposition[clusterid];
    int newposition=borders_of_set[clustersizes[clusterid]];
    int swapid=sorted_clustersizes[newposition];
    if(swapid!=-1){
        clusterid_to_arrayposition[swapid]=oldposition;
    }
    sorted_clustersizes[oldposition]=swapid;

    sorted_clustersizes[newposition]=clusterid;
    clusterid_to_arrayposition[clusterid]=newposition;
    borders_of_set[clustersizes[clusterid]]++;
    clustersizes[clusterid]--;
}

void ClusteringAlgorithms::setCover(unsigned int **elementLookupTable, unsigned short ** elementScoreLookupTable,
                                    int *assignedcluster, short *bestscore, size_t *newElementOffsets) {
    for (int cl_size = dbSize - 1; cl_size >= 0; cl_size--) {
        int representative = sorted_clustersizes[cl_size];
        if (representative < 0) {
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
            short seqId = elementScoreLookupTable[representative][elementId];
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
                                             unsigned short **elementScoreLookupTable, size_t elementCount,
                                             size_t n, int *assignedcluster, unsigned short *scoreelements) {
    scoreelements=new unsigned short[elementCount];
    std::fill_n(scoreelements, elementCount, 0);
    elementScoreLookupTable= new unsigned short*[dbSize];
    AlignmentSymmetry::setupElementLookupPointerShort(scoreelements, elementScoreLookupTable, elementOffsets, dbSize);
    for(size_t i = 0; i < n; i++) {
        // seqDbr is descending sorted by length
        // the assumption is that clustering is B -> B (not A -> B)
        Log::printProgress(i);
        if(assignedcluster[i]==-1){
            size_t elementSize = (elementOffsets[i + 1] - elementOffsets[i]);
            for (size_t elementId = 0; elementId < elementSize; elementId++) {
                unsigned int id = elementLookupTable[i][elementId];
                if(assignedcluster[id]==id){
                    assignedcluster[i]=id;
                    break;
                }
            }
            if(assignedcluster[i]==-1) {
                assignedcluster[i]=i;
            }
        }
    }
    delete [] scoreelements;
    delete [] elementScoreLookupTable;
}

void ClusteringAlgorithms::readInClusterData(unsigned int **elementLookupTable, unsigned int *&elements,
                                             unsigned short ** elementScoreLookupTable, unsigned short *&scoreelements,
                                             size_t *elementOffsets, size_t elementCount) {

    //time
    struct timeval start, end;
    gettimeofday(&start, NULL);
#pragma omp parallel
    {
#pragma omp for schedule(dynamic, 1000)
        for(size_t i = 0; i < dbSize; i++) {
            unsigned int clusterId = seqDbr->getDbKey(i);
            const size_t alnId = alnDbr->getId(clusterId);
            char *data = alnDbr->getData(alnId);
            size_t dataSize = alnDbr->getSeqLens()[alnId];
            size_t elementCount = Util::countLines(data, dataSize);
            elementOffsets[i] = elementCount;
        }

    }
    // make offset table
    AlignmentSymmetry::computeOffsetTable(elementOffsets, dbSize);

    // set element edge pointers by using the offset table
    AlignmentSymmetry::setupElementLookupPointer(elements, elementLookupTable, elementOffsets, dbSize);

    // fill elements
    AlignmentSymmetry::readInData(alnDbr, seqDbr, elementLookupTable);
    // set element edge pointers by using the offset table
    AlignmentSymmetry::setupElementLookupPointer(elements, elementLookupTable, elementOffsets, dbSize);

    Debug(Debug::WARNING) << "\nSort entries.\n";
    // sort each element vector for bsearch
#pragma omp parallel for schedule(dynamic, 1000)
    for (size_t i = 0; i < dbSize; i++) {
        Log::printProgress(i);
        std::sort(elementLookupTable[i], elementLookupTable[i] + (elementOffsets[i + 1] - elementOffsets[i]));
    }

    //time
    Debug(Debug::WARNING) << "\nFind missing connections.\n";
    size_t *newElementOffsets = new size_t[dbSize + 1];
    memcpy(newElementOffsets, elementOffsets, sizeof(size_t) * (dbSize + 1));
    size_t newElementCount = AlignmentSymmetry::findMissingLinks(elementLookupTable, newElementOffsets, dbSize,
                                                                 threads);
    //time
    Debug(Debug::WARNING) << "\nFind missing connections.\n";

    // resize elements
    delete[] elements;
    elements = new unsigned int[newElementCount];
    std::fill_n(elements, newElementCount, 0);
    delete[] scoreelements;
    scoreelements = new unsigned short[newElementCount];
    std::fill_n(scoreelements, newElementCount, 0);

    Debug(Debug::WARNING) << "\nFound " << newElementCount - elementCount << " new connections.\n";
    AlignmentSymmetry::setupElementLookupPointer(elements, elementLookupTable, newElementOffsets, dbSize);
    AlignmentSymmetry::setupElementLookupPointerShort(scoreelements, elementScoreLookupTable, newElementOffsets,
                                                      dbSize);
    //time
    Debug(Debug::WARNING) << "\nReconstruct initial order.\n";
    alnDbr->remapData();
    seqDbr->remapData();
    AlignmentSymmetry::readInData(alnDbr, seqDbr, elementLookupTable, elementScoreLookupTable, scoretype);
    // set element edge pointers by using the offset table
    AlignmentSymmetry::setupElementLookupPointer(elements, elementLookupTable, newElementOffsets, dbSize);
    AlignmentSymmetry::setupElementLookupPointerShort(scoreelements, elementScoreLookupTable, newElementOffsets,
                                                      dbSize);

    Debug(Debug::WARNING) << "\nAdd missing connections.\n";
    AlignmentSymmetry::addMissingLinks(elementLookupTable, elementOffsets, dbSize, elementScoreLookupTable);
    AlignmentSymmetry::setupElementLookupPointer(elements, elementLookupTable, newElementOffsets, dbSize);
    AlignmentSymmetry::setupElementLookupPointerShort(scoreelements, elementScoreLookupTable, newElementOffsets,
                                                      dbSize);
    maxClustersize = 0;
    for (size_t i = 0; i < dbSize; i++) {
        elementCount = newElementOffsets[i + 1] - newElementOffsets[i];
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