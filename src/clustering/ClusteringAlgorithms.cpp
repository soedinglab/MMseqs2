//
// Created by lars on 08.06.15.
//

#include <sys/time.h>
#include <Log.h>
#include "ClusteringAlgorithms.h"
#include "Util.h"
#include "Debug.h"
#include "AffinityClustering.h"
#include "AlignmentSymmetry.h"
#include <queue>

ClusteringAlgorithms::ClusteringAlgorithms(DBReader<unsigned int>* seqDbr, DBReader<unsigned int>* alnDbr,int threads, int scoretype, int maxiterations){
    this->seqDbr=seqDbr;
    this->alnDbr=alnDbr;
    this->dbSize=alnDbr->getSize();
    this->threads=threads;
    this->scoretype=scoretype;
    this->maxiterations=maxiterations;
}

std::list<set *>  ClusteringAlgorithms::execute(int mode) {

    set* sets = new set[dbSize];
    memset(sets, 0, sizeof(set)*(dbSize));
    //time
    struct timeval start, end;
    gettimeofday(&start, NULL);
    ///time
    clustersizes=new int[dbSize];
    std::fill_n(clustersizes,dbSize,0);


    const char * data = alnDbr->getData();
    size_t dataSize = alnDbr->getDataSize();
    size_t elementCount = Util::count_lines(data, dataSize);
    unsigned int * elements = new unsigned int[elementCount];
    unsigned int ** elementLookupTable = new unsigned int*[dbSize];
    size_t *elementOffsets = new size_t[dbSize + 1];
    elementOffsets[dbSize]=0;
    unsigned short *scoreelements;
    unsigned short **elementScoreLookupTable;
    unsigned int * seqDbrIdToalnDBrId= new unsigned int[dbSize];


    std::list<set *> result;
    size_t n = seqDbr->getSize();
    int *assignedcluster = new int[n];
    std::fill_n(assignedcluster, n, -1);
    short *bestscore = new short[n];
    std::fill_n(bestscore, n, -10);
#pragma omp parallel
    {

#pragma omp for schedule(dynamic, 1000)
    for(size_t i = 0; i < dbSize; i++) {
        unsigned int clusterId = seqDbr->getDbKey(i);
        const size_t alnId = alnDbr->getId(clusterId);
        seqDbrIdToalnDBrId[i]=alnId;
        char *data = alnDbr->getData(alnId);
        size_t dataSize = alnDbr->getSeqLens()[alnId];
        size_t elementCount = Util::count_lines(data, dataSize);
        elementOffsets[i] = elementCount;
    }

    }
    //time
    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for Parallel read in: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);

    //time
    // make offset table
    AlignmentSymmetry::computeOffsetTable(elementOffsets, dbSize);
    //time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for offset table: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time
    // set element edge pointers by using the offset table
    AlignmentSymmetry::setupElementLookupPointer(elements, elementLookupTable, elementOffsets, dbSize);
    //time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for lookuppointer: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time
    // fill elements
    AlignmentSymmetry::readInData(alnDbr, seqDbr, elementLookupTable);
    // set element edge pointers by using the offset table
    AlignmentSymmetry::setupElementLookupPointer(elements, elementLookupTable, elementOffsets, dbSize);
    if (mode==2){
        scoreelements=new unsigned short[elementCount];
        std::fill_n(scoreelements, elementCount, 0);
        elementScoreLookupTable= new unsigned short*[dbSize];
        AlignmentSymmetry::setupElementLookupPointerShort(scoreelements, elementScoreLookupTable, elementOffsets, dbSize);
    }
    //time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for Read in: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time
    if (mode==2){

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





    }else {
        //time
        Debug(Debug::WARNING) << "\nSort entries.\n";
        // sort each element vector for bsearch
#pragma omp parallel for schedule(dynamic, 1000)
        for (size_t i = 0; i < dbSize; i++) {
            Log::printProgress(i);
            std::sort(elementLookupTable[i], elementLookupTable[i] + (elementOffsets[i + 1] - elementOffsets[i]));
        }
        //time
        gettimeofday(&end, NULL);
        sec = end.tv_sec - start.tv_sec;
        Debug(Debug::INFO) << "\nTime for Sorting entries: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
        gettimeofday(&start, NULL);
        //time
        Debug(Debug::WARNING) << "\nFind missing connections.\n";
        size_t *newElementOffsets = new size_t[dbSize + 1];
        memcpy(newElementOffsets, elementOffsets, sizeof(size_t) * (dbSize + 1));
        size_t newElementCount = AlignmentSymmetry::findMissingLinks(elementLookupTable, newElementOffsets, dbSize,
                                                                     threads);
        delete[] elements;
        //time
        gettimeofday(&end, NULL);
        sec = end.tv_sec - start.tv_sec;
        Debug(Debug::INFO) << "\nTime for finding missing connections: " << (sec / 60) << " m " << (sec % 60) <<
        "s\n\n";
        gettimeofday(&start, NULL);
        //time
        Debug(Debug::WARNING) << "\nFind missing connections.\n";

        // resize elements
        elements = new unsigned int[newElementCount];
        std::fill_n(elements, newElementCount, 0);
        unsigned short *scoreelements = new unsigned short[newElementCount];
        std::fill_n(scoreelements, newElementCount, 0);
        unsigned short **elementScoreLookupTable = new unsigned short *[dbSize];
        Debug(Debug::WARNING) << "\nFound " << newElementCount - elementCount << " new connections.\n";
        AlignmentSymmetry::setupElementLookupPointer(elements, elementLookupTable, newElementOffsets, dbSize);
        AlignmentSymmetry::setupElementLookupPointerShort(scoreelements, elementScoreLookupTable, newElementOffsets,
                                                          dbSize);
        //time
        gettimeofday(&end, NULL);
        sec = end.tv_sec - start.tv_sec;
        Debug(Debug::INFO) << "\nTime: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
        gettimeofday(&start, NULL);
        //time
        Debug(Debug::WARNING) << "\nReconstruct initial order.\n";
        alnDbr->remapData();
        seqDbr->remapData();
        AlignmentSymmetry::readInData(alnDbr, seqDbr, elementLookupTable, elementScoreLookupTable, scoretype);
        // set element edge pointers by using the offset table
        AlignmentSymmetry::setupElementLookupPointer(elements, elementLookupTable, newElementOffsets, dbSize);
        AlignmentSymmetry::setupElementLookupPointerShort(scoreelements, elementScoreLookupTable, newElementOffsets,
                                                          dbSize);
//time
        gettimeofday(&end, NULL);
        sec = end.tv_sec - start.tv_sec;
        Debug(Debug::INFO) << "\nTime: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
        gettimeofday(&start, NULL);
        //time
        Debug(Debug::WARNING) << "\nAdd missing connections.\n";
        AlignmentSymmetry::addMissingLinks(elementLookupTable, elementOffsets, dbSize, elementScoreLookupTable);
        AlignmentSymmetry::setupElementLookupPointer(elements, elementLookupTable, newElementOffsets, dbSize);
        AlignmentSymmetry::setupElementLookupPointerShort(scoreelements, elementScoreLookupTable, newElementOffsets,
                                                          dbSize);
        //time
        gettimeofday(&end, NULL);
        sec = end.tv_sec - start.tv_sec;
        Debug(Debug::INFO) << "\nTime: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
        gettimeofday(&start, NULL);
        //time

        maxClustersize = 0;
        for (size_t i = 0; i < dbSize; i++) {
            elementCount = newElementOffsets[i + 1] - newElementOffsets[i];
            maxClustersize = std::max((unsigned int) elementCount, maxClustersize);
            clustersizes[i] = elementCount;


        }
        ClusteringAlgorithms::initClustersizes();
        //time
        gettimeofday(&end, NULL);
        sec = end.tv_sec - start.tv_sec;
        Debug(Debug::INFO) << "\nTime for Maximum determination: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
        gettimeofday(&start, NULL);




        //time
        gettimeofday(&end, NULL);
        sec = end.tv_sec - start.tv_sec;
        Debug(Debug::INFO) << "\nTime for Clustersize insertion " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
        gettimeofday(&start, NULL);
        //time
        //delete from beginning
        if (mode == 1) {
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
                        size_t elementSize = (newElementOffsets[currentid + 1] - newElementOffsets[currentid]);
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
        delete [] elementLookupTable;
    delete [] elements;
    delete [] elementOffsets;
    delete [] seqDbrIdToalnDBrId;
    delete [] bestscore;
    delete [] sorted_clustersizes;
    delete [] clusterid_to_arrayposition;
    delete [] borders_of_set;

    }
//time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for Cluster computation: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time

    for(size_t i = 0; i < n; i++) {
        AffinityClustering::add_to_set(i,&sets[assignedcluster[i]],assignedcluster[i]);
    }

    for(size_t i = 0; i < n; i++) {
        set * max_set = &sets[i];
        if (max_set->elements == NULL)
            continue;
        result.push_back(max_set); // O(1)
    }
    //time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for preparing cluster sets: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time
    return result;
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
