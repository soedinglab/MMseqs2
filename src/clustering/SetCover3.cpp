//
// Created by lars on 08.06.15.
//

#include <sys/time.h>
#include <Log.h>
#include "SetCover3.h"
#include "Util.h"
#include "Debug.h"
#include "AffinityClustering.h"
#include "AlignmentSymmetry.h"
#include <queue>

SetCover3::SetCover3(DBReader * seqDbr, DBReader * alnDbr, float seqIdThr, float coverage, int threads){
    this->seqDbr=seqDbr;
    this->alnDbr=alnDbr;
    this->seqIdThr=seqIdThr;
    this->coverage=coverage;
    this->dbSize=alnDbr->getSize();
    this->threads=threads;
}

std::list<set *>  SetCover3::execute(int mode) {

    set* sets = new set[dbSize];
    memset(sets, 0, sizeof(set *)*(dbSize));
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
        char *clusterId = seqDbr->getDbKey(i);
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
    memset(elements, 0, sizeof( unsigned int) * (newElementCount ));
    unsigned short *scoreelements=new unsigned short[newElementCount];
    memset(scoreelements, 0, sizeof( unsigned short) * (newElementCount ));
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
    AlignmentSymmetry::readInData(alnDbr, seqDbr, elementLookupTable,elementScoreLookupTable);
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

    maxClustersize=0;
    for(size_t i = 0; i < dbSize; i++) {
        elementCount=newElementOffsets[i +1] - newElementOffsets[i];
        maxClustersize= std::max((unsigned int) elementCount, maxClustersize);
        clustersizes[i] = elementCount;


    }
    SetCover3::initClustersizes();
    //time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for Maximum determination: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);


    std::list<set *> result;
    size_t n = seqDbr->getSize();
    int* assignedcluster=new int[n];
    memset(assignedcluster, -1, sizeof(int)*(n));
    short* bestscore=new short[n];
    memset(bestscore, -10, sizeof(short)*(n));



    //time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for Clustersize insertion " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time
    //delete from beginning
    if(mode==1) {
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
    }else if (mode==3){
        Debug(Debug::INFO)<<"connected component mode"<<"\n";
        for (int cl_size = dbSize - 1; cl_size >= 0; cl_size--) {
            int representative = sorted_clustersizes[cl_size];
            if(assignedcluster[representative]==-1){
                assignedcluster[representative] = representative;
                std::queue<int> myqueue;
                    myqueue.push(representative);
                //delete clusters of members;
                while(!myqueue.empty()){
                    int currentid=myqueue.front();
                    assignedcluster[currentid]=representative;
                    myqueue.pop();
                    size_t elementSize = (newElementOffsets[currentid + 1] - newElementOffsets[currentid]);
                    for (size_t elementId = 0; elementId < elementSize; elementId++) {
                        unsigned int elementtodelete = elementLookupTable[currentid][elementId];
                        if(assignedcluster[elementtodelete]==-1){
                            myqueue.push(elementtodelete);
                        }
                        assignedcluster[elementtodelete]=representative;
                    }
                }

            }
        }

    }else if (mode==2){
        Debug(Debug::INFO)<<"connected component mode"<<"\n";
        int* ranks=new int[n];
        memset(ranks, 0, sizeof(int)*(n));
        int* incommingconnections=new int[n];
        memset(incommingconnections, 0, sizeof(int)*(n));
        int* connactionswithsamerank=new int[n];
        memset(connactionswithsamerank, 0, sizeof(int)*(n));
        int connectioncutoff=2;
        for (int cl_size = dbSize - 1; cl_size >= 0; cl_size--) {
            //  for (int cl_size =0 ;cl_size<dbSize ; cl_size++) {
            int representative = sorted_clustersizes[cl_size];
            if(assignedcluster[representative]==-1){
                assignedcluster[representative] = representative;
                connactionswithsamerank[representative]=connectioncutoff+1;
                incommingconnections[representative]=connectioncutoff+1;
                std::queue<int> myqueue;
                std::queue<int> myqueue2;
                myqueue.push(representative);
                int elementSize = (newElementOffsets[representative + 1] - newElementOffsets[representative]);
                elementSize=std::min(elementSize,std::max(elementSize/5,10));
                connectioncutoff=elementSize/2;
                for (size_t elementId = 0; elementId < elementSize; elementId++) {
                    unsigned int elementtodelete = elementLookupTable[representative][elementId];
                    if (assignedcluster[elementtodelete] == -1) {
                        incommingconnections[elementtodelete]=connectioncutoff+1;
                        assignedcluster[elementtodelete] = representative;
                        ranks[elementtodelete] = 0;
                    }
                }

                //delete clusters of members;
                while(!myqueue.empty()) {
                    while (!myqueue.empty()) {
                        int currentid = myqueue.front();

                        myqueue.pop();
                        if(incommingconnections[currentid]<=connectioncutoff){
                            assignedcluster[currentid]=-1;
                             incommingconnections[currentid]=0;
                            continue;
                        }
                        assignedcluster[currentid] = representative;
                        size_t elementSize = (newElementOffsets[currentid + 1] - newElementOffsets[currentid]);
                        for (size_t elementId = 0; elementId < elementSize; elementId++) {
                            unsigned int elementtodelete = elementLookupTable[currentid][elementId];
                            if (assignedcluster[elementtodelete] == -1) {
                                incommingconnections[elementtodelete]++;
                            } else if (ranks[elementtodelete] == ranks[currentid]) {
                                connactionswithsamerank[currentid]++;
                            }

                        }
                        if (connactionswithsamerank[currentid] > connectioncutoff&& incommingconnections[currentid]>connectioncutoff) {
                            myqueue2.push(currentid);
                        }else{
                          //  assignedcluster[currentid]=-1;
                           // incommingconnections[currentid]=0;
                        }


                    }
                    while (!myqueue2.empty()) {
                        int currentid = myqueue2.front();
                        assignedcluster[currentid] = representative;
                        myqueue2.pop();
                        size_t elementSize = (newElementOffsets[currentid + 1] - newElementOffsets[currentid]);
                        for (size_t elementId = 0; elementId < elementSize; elementId++) {
                            unsigned int elementtodelete = elementLookupTable[currentid][elementId];
                            if (assignedcluster[elementtodelete] == -1) {
                                myqueue.push(elementtodelete);
                                ranks[elementtodelete] = ranks[currentid] + 1;
                            }
                            assignedcluster[elementtodelete] = representative;
                        }
                    }
                }

            }
        }

    }
    //delete unnecessary datastructures
   /* delete [] elementLookupTable;
    delete [] elements;
    delete [] elementOffsets;
    delete [] maxClustersizes;
    delete [] seqDbrIdToalnDBrId;
    delete [] bestscore;
    delete [] sorted_clustersizes;
    delete [] clusterid_to_arrayposition;
    delete [] borders_of_set;
*/
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




void SetCover3::initClustersizes(){
    int * setsize_abundance=new int[maxClustersize+1];
    memset(setsize_abundance,0,sizeof(int)*(maxClustersize+1));
    //count how often a set size occurs
    for (int i = 0; i < dbSize; ++i) {
        setsize_abundance[clustersizes[i]]++;
    }
    //compute offsets
    borders_of_set= new int [maxClustersize+1];
    borders_of_set[0]=0;
    for (int i = 1; i < maxClustersize+1; ++i) {
        borders_of_set[i]=borders_of_set[i-1]+setsize_abundance[i-1];
    }
    //fill array
    sorted_clustersizes =new int [dbSize];
    memset(sorted_clustersizes,0,sizeof(int)*dbSize);
    clusterid_to_arrayposition=new int[dbSize];
    memset(clusterid_to_arrayposition,0,sizeof(int)*dbSize);
    //reuse setsize_abundance as offset counter
    memset(setsize_abundance,0,sizeof(int)*(maxClustersize+1));
    for (int i = 0; i < dbSize; ++i) {
        int position=borders_of_set[clustersizes[i]]+setsize_abundance[clustersizes[i]];
        sorted_clustersizes[position]=i;
        clusterid_to_arrayposition[i]=position;
        setsize_abundance[clustersizes[i]]++;
    }
}


void SetCover3::removeClustersize(int clusterid){
    clustersizes[clusterid]=0;
    sorted_clustersizes[clusterid_to_arrayposition[clusterid]]=-1;
    clusterid_to_arrayposition[clusterid]=-1;
}

void SetCover3::decreaseClustersize(int clusterid){
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