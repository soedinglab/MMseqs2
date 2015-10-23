//
// Created by lars on 28.07.15.
//
//
// Created by lars on 08.06.15.
//

#include <sys/time.h>
#include <Log.h>
#include "SetCover4.h"
#include "Util.h"
#include "Debug.h"
#include "AffinityClustering.h"
#include "AlignmentSymmetry.h"

SetCover4::SetCover4(DBReader * seqDbr, DBReader * alnDbr, float seqIdThr, float coverage,int threads,int scoretype){

    this->seqDbr=seqDbr;
    this->alnDbr=alnDbr;
    this->seqIdThr=seqIdThr;
    this->coverage=coverage;
    this->dbSize=alnDbr->getSize();
    this->threads=threads;
    this->scoretype=scoretype;

}

std::list<set *>  SetCover4::execute() {
    set* sets = new set[dbSize];
    memset(sets, 0, sizeof(set *)*(dbSize));
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

    unsigned int * seqDbrIdToalnDBrId= new unsigned int[dbSize];
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
#pragma omp for schedule(dynamic, 1000)
        for(size_t i = 0; i < dbSize; i++) {
            const char *clusterId = seqDbr->getDbKey(i).c_str();
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
    //time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for Read in: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time

    //time
    Debug(Debug::WARNING) << "\nSort entries.\n";
    // sort each element vector for bsearch
#pragma omp parallel for schedule(dynamic, 1000)
    for(size_t i = 0; i < dbSize; i++) {
        Log::printProgress(i);
        std::sort(elementLookupTable[i], elementLookupTable[i] + (elementOffsets[i+1] - elementOffsets[i]));
    }
    //time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for Sorting entries: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time
    Debug(Debug::WARNING) << "\nFind missing connections.\n";
    size_t * newElementOffsets = new size_t[dbSize + 1];
    memcpy(newElementOffsets, elementOffsets, sizeof(size_t) * (dbSize + 1));
    size_t newElementCount = AlignmentSymmetry::findMissingLinks(elementLookupTable, newElementOffsets, dbSize,threads);
    delete [] elements;
    //time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for finding missing connections: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time
    Debug(Debug::WARNING) << "\nFind missing connections.\n";

    // resize elements
    elements = new unsigned int[newElementCount];
    std::fill_n(elements, newElementCount, 0);
    unsigned short *scoreelements=new unsigned short[newElementCount];
    std::fill_n(scoreelements, newElementCount, 0);
    unsigned short **elementScoreLookupTable= new unsigned short*[dbSize];
    Debug(Debug::WARNING) << "\nFound "<< newElementCount - elementCount << " new connections.\n";
    AlignmentSymmetry::setupElementLookupPointer(elements, elementLookupTable, newElementOffsets, dbSize);
    AlignmentSymmetry::setupElementLookupPointerShort(scoreelements, elementScoreLookupTable, newElementOffsets, dbSize);
    //time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time
    Debug(Debug::WARNING) << "\nReconstruct initial order.\n";
    alnDbr->remapData();
    seqDbr->remapData();
    AlignmentSymmetry::readInData(alnDbr, seqDbr, elementLookupTable,elementScoreLookupTable,scoretype);
    // set element edge pointers by using the offset table
    AlignmentSymmetry::setupElementLookupPointer(elements, elementLookupTable, newElementOffsets, dbSize);
    AlignmentSymmetry::setupElementLookupPointerShort(scoreelements, elementScoreLookupTable, newElementOffsets, dbSize);
//time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time
    Debug(Debug::WARNING) << "\nAdd missing connections.\n";
    AlignmentSymmetry::addMissingLinks(elementLookupTable, elementOffsets, dbSize,elementScoreLookupTable);
    AlignmentSymmetry::setupElementLookupPointer(elements, elementLookupTable, newElementOffsets, dbSize);
    AlignmentSymmetry::setupElementLookupPointerShort(scoreelements, elementScoreLookupTable, newElementOffsets, dbSize);
    //time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time





    ///////////////////////////////////////////////////////
    std::list<set *> result;
    size_t n = seqDbr->getSize();
    int* assignedcluster=new int[n];
    memset(assignedcluster, -1, sizeof(int)*(n));
    float* bestscore=new float[n];
    memset(bestscore, -10, sizeof(float)*(n));


    //initialise
    clustersizes=new int[dbSize];
    memset(clustersizes,0,dbSize);
    clusterweights=new float[dbSize];
    memset(clusterweights,0,dbSize);
    maxClustersize=0;
    maxClusterweight=0;

//read in data and determine biggest cluster

    for(size_t i = 0; i < dbSize; i++) {
        elementCount=newElementOffsets[i +1] - newElementOffsets[i];
        maxClustersize= std::max((unsigned int) elementCount, maxClustersize);
        clustersizes[i] = elementCount;


    }
    for (int i = 0; i < dbSize; i++) {
        unsigned int currentclustersize=0;
        unsigned int currentclusterweight=0;
        size_t elementSize = (newElementOffsets[i + 1] - newElementOffsets[i]);
        for (size_t elementId = 0; elementId < elementSize; elementId++) {
            short seqId = elementScoreLookupTable[i][elementId]+offset;
            currentclusterweight+=seqId;
            currentclustersize++;
        }
        maxClustersize=std::max(currentclustersize,maxClustersize);
        maxClusterweight=std::max(currentclusterweight,maxClusterweight);
        clustersizes[i]=currentclustersize;
        clusterweights[i]=currentclusterweight;
    }
    //save ordered clustersizes
    orderedClustersizes=new std::set<int>[maxClustersize+1];
    weightslength=(maxClusterweight+1)*spreadingfactor;
    orderedClusterweights=new std::set<int>[weightslength];
    for (int i = 0; i < dbSize; i++) {
        SetCover4::insertCluster(i,clustersizes[i]);
        SetCover4::insertClusterWeight(i,clusterweights[i]);
    }
    //delete from beginning
    for (int cl_size = weightslength-1; cl_size >= 0; cl_size--) {
        while(orderedClusterweights[cl_size].size()>0) {
            int representative =*orderedClusterweights[cl_size].begin();
//          Debug(Debug::INFO)<<alnDbr->getDbKey(representative)<<"\n";
            SetCover4::eraseClusterWeight(representative,clusterweights[representative]);
            orderedClustersizes[clustersizes[representative]].erase(representative);
            clustersizes[representative]=0;
            assignedcluster[representative]=representative;

            //delete clusters of members;
            size_t elementSize = (newElementOffsets[representative + 1] - newElementOffsets[representative]);
            for (size_t elementId = 0; elementId < elementSize; elementId++) {

                const unsigned int elementtodelete = elementLookupTable[representative][elementId];
                bool representativefound=false;
                //
                unsigned int seqId = elementScoreLookupTable[representative][elementId]+offset;
                unsigned int oldseqId=bestscore[elementtodelete];
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
                SetCover4::eraseClusterWeight(elementtodelete,clusterweights[elementtodelete]);
                orderedClustersizes[clustersizes[elementtodelete]].erase(elementtodelete);

                clustersizes[elementtodelete]=0;
                clusterweights[elementtodelete]=0;
                //decrease clustersize of sets that contain the element
                const unsigned int currElementSize = (newElementOffsets[elementtodelete + 1] -
                                                      newElementOffsets[elementtodelete]);
                for (size_t elementId2 = 0; elementId2 < currElementSize; elementId2++) {
                    const unsigned int elementtodecrease = elementLookupTable[elementtodelete][elementId2];
                    unsigned int seqIdOfElement =elementScoreLookupTable[elementtodelete][elementId2]+offset;
                    if(representative == elementtodecrease){
                        representativefound=true;
                    }
                    if(clustersizes[elementtodecrease]==1) {
                        Debug(Debug::ERROR)<<"there must be an error: "<<alnDbr->getDbKey(elementtodelete)<<" deleted from "<<alnDbr->getDbKey(elementtodecrease)<<" that now is empty, but not assigned to a cluster\n";
                    }else if (clustersizes[elementtodecrease]>0) {
                        orderedClustersizes[clustersizes[elementtodecrease]].erase(elementtodecrease);
                        SetCover4::eraseClusterWeight(elementtodecrease,clusterweights[elementtodecrease]);
                        clustersizes[elementtodecrease]--;
                        clusterweights[elementtodecrease]+=std::min(oldseqId,seqIdOfElement)-std::max(std::min(oldseqId,seqIdOfElement),std::min(seqId,seqIdOfElement));
                        orderedClustersizes[clustersizes[elementtodecrease]].insert(elementtodecrease);
                        SetCover4::insertClusterWeight(elementtodecrease,clusterweights[elementtodecrease]);
                    }


                }
                if(!representativefound){
                    Debug(Debug::ERROR)<<"error with cluster:\t"<<alnDbr->getDbKey(representative)<<"\tis not contained in set:\t"<<alnDbr->getDbKey(elementtodelete)<<".\n";
                }

                            }


        }

    }


    for(size_t i = 0; i < n; i++) {
        AffinityClustering::add_to_set(i,&sets[assignedcluster[i]],assignedcluster[i]);
    }
    for(size_t i = 0; i < n; i++) {
        set * max_set = &sets[i];
        if (max_set->elements == NULL)
            continue;
        result.push_back(max_set); // O(1)
    }
    return result;

}


void SetCover4::insertCluster(int clusterid, int clustersize){
    //if(orderedClustersizes[clustersize] ==NULL){
    orderedClustersizes[clustersize].insert(clusterid);
}

void SetCover4::insertClusterWeight(int clusterid,float clusterweight){
    orderedClusterweights[(int)(clusterweight*1)*spreadingfactor].insert(clusterid);
}
void SetCover4::eraseClusterWeight(int clusterid,float clusterweight){
    orderedClusterweights[(int)(clusterweight*1)*spreadingfactor].erase(clusterid);
}

