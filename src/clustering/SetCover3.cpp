//
// Created by lars on 08.06.15.
//

#include <sys/time.h>
#include "SetCover3.h"
#include "Util.h"
#include "Debug.h"
#include "AffinityClustering.h"
#include "AlignmentSymmetry.h"

SetCover3::SetCover3(DBReader * seqDbr, DBReader * alnDbr, float seqIdThr, float coverage, int threads){
    this->seqDbr=seqDbr;
    this->alnDbr=alnDbr;
    this->seqIdThr=seqIdThr;
    this->coverage=coverage;
    this->dbSize=alnDbr->getSize();
    this->threads=threads;
}

std::list<set *>  SetCover3::execute() {
    //time
    struct timeval start, end;
    gettimeofday(&start, NULL);
    ///time
    clustersizes=new int[dbSize];
    memset(clustersizes,0,dbSize);

    const char * data = alnDbr->getData();
    size_t dataSize = alnDbr->getDataSize();
    size_t elementCount = Util::count_lines(data, dataSize);
    unsigned int * elements = new unsigned int[elementCount];
    unsigned int ** elementLookupTable = new unsigned int*[dbSize];
    size_t *elementOffsets = new size_t[dbSize + 1];
    int* maxClustersizes=new int [threads];
    unsigned int * seqDbrIdToalnDBrId= new unsigned int[dbSize];
#pragma omp parallel
    {
    int thread_idx = 0;
#ifdef OPENMP
    thread_idx = omp_get_thread_num();
#endif
        maxClustersizes[thread_idx]=0;
#pragma omp for schedule(dynamic, 1000)
    for(size_t i = 0; i < dbSize; i++) {
        char *clusterId = seqDbr->getDbKey(i);
        const size_t alnId = alnDbr->getId(clusterId);
        seqDbrIdToalnDBrId[i]=alnId;
        char *data = alnDbr->getData(alnId);
        size_t dataSize = alnDbr->getSeqLens()[alnId];
        size_t elementCount = Util::count_lines(data, dataSize);
        elementOffsets[i] = elementCount;
        maxClustersizes[thread_idx]=std::max((int)elementCount,maxClustersizes[thread_idx]);
        clustersizes[i]=elementCount;
    }

    }
    //time
    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for Parallel read in: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time
    int maxClustersize=0;
    for (int j = 0; j < threads; ++j) {
        maxClustersize=std::max(maxClustersize,maxClustersizes[j]);
    }
    //time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for Maximum determination: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
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
    //time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for Read in: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time

    std::list<set *> result;
    size_t n = seqDbr->getSize();
    int* assignedcluster=new int[n];
    memset(assignedcluster, -1, sizeof(int)*(n));
    float* bestscore=new float[n];
    memset(bestscore, -10, sizeof(float)*(n));
    char *similarity = new char[255+1];

    set* sets = new set[n];
    memset(sets, 0, sizeof(set *)*(n));


    //save ordered clustersizes
    orderedClustersizes=new std::set<int>[maxClustersize+1];
    for (int i = 0; i < dbSize; i++) {
        SetCover3::insertCluster(i,clustersizes[i]);
    }
    //time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for Clustersize insertion " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time
    //delete from beginning
    for (int cl_size = maxClustersize; cl_size > 0; cl_size--) {
        while(orderedClustersizes[cl_size].size()>0) {
            int representative =*orderedClustersizes[cl_size].begin();
//          Debug(Debug::INFO)<<alnDbr->getDbKey(representative)<<"\n";
            orderedClustersizes[cl_size].erase(representative);
            clustersizes[representative]=0;
            assignedcluster[representative]=representative;
            char *data = alnDbr->getData(seqDbrIdToalnDBrId[representative]);
            //delete clusters of members;
            size_t elementSize = (elementOffsets[representative +1] - elementOffsets[representative]);
            for(size_t elementId = 0; elementId < elementSize; elementId++) {

                const unsigned int elementtodelete = elementLookupTable[representative][elementId];
                const unsigned int currElementSize = (elementOffsets[elementtodelete +1] - elementOffsets[elementtodelete]);

                bool representativefound=false;
                //


                Util::parseByColumnNumber(data, similarity, 4); //column 4 = sequence identity
                data = Util::skipLine(data);
                float seqId = atof(similarity);
                if(seqId>bestscore[elementtodelete]) {
                    assignedcluster[elementtodelete] = representative;
                    bestscore[elementtodelete] = seqId;
                }
                if(elementtodelete== representative){
                    continue;
                }
                if(clustersizes[elementtodelete]<1){
                    continue;
                }

                orderedClustersizes[clustersizes[elementtodelete]].erase(elementtodelete);
                clustersizes[elementtodelete]=0;
                //decrease clustersize of sets that contain the element
                for(size_t elementId2 = 0; elementId2 < currElementSize; elementId2++) {
                    const unsigned int elementtodecrease = elementLookupTable[elementtodelete][elementId2];
                    if(representative == elementtodecrease){
                        representativefound=true;
                    }
                    if(clustersizes[elementtodecrease]==1) {
                        Debug(Debug::ERROR)<<"there must be an error: "<<alnDbr->getDbKey(elementtodelete)<<" deleted from "<<alnDbr->getDbKey(elementtodecrease)<<" that now is empty, but not assigned to a cluster\n";
                    }else if (clustersizes[elementtodecrease]>0) {
                        orderedClustersizes[clustersizes[elementtodecrease]].erase(elementtodecrease);
                        clustersizes[elementtodecrease]--;
                        orderedClustersizes[clustersizes[elementtodecrease]].insert(elementtodecrease);
                    }
                }
                if(!representativefound){
                    Debug(Debug::ERROR)<<"error with cluster:\t"<<alnDbr->getDbKey(representative)<<"\tis not contained in set:\t"<<alnDbr->getDbKey(elementtodelete)<<".\n";
                }

            }


        }

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


void SetCover3::insertCluster(int clusterid, int clustersize){
    //if(orderedClustersizes[clustersize] ==NULL){
    orderedClustersizes[clustersize].insert(clusterid);
}
