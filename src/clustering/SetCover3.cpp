//
// Created by lars on 08.06.15.
//

#include "SetCover3.h"
#include "Util.h"
#include "Debug.h"
#include "AffinityClustering.h"

SetCover3::SetCover3(DBReader * seqDbr, DBReader * alnDbr, float seqIdThr, float coverage){
    this->seqDbr=seqDbr;
    this->alnDbr=alnDbr;
    this->seqIdThr=seqIdThr;
    this->coverage=coverage;
    this->dbSize=alnDbr->getSize();
}

std::list<set *>  SetCover3::execute() {
    std::list<set *> result;
    char *idbuffer1 = new char[255 + 1];
    char *idbuffer2 = new char[255 + 1];
    size_t n = seqDbr->getSize();
    int* assignedcluster=new int[n];
    memset(assignedcluster, -1, sizeof(int)*(n));
    float* bestscore=new float[n];
    memset(bestscore, -10, sizeof(float)*(n));
    char *similarity = new char[255+1];

    set* sets = new set[n];
    memset(sets, 0, sizeof(set *)*(n));

    //initialise
    clustersizes=new int[dbSize];
    memset(clustersizes,0,dbSize);
    maxClustersize=0;

//read in data and determine biggest cluster
    for (int i = 0; i < dbSize; i++) {
        unsigned int currentclustersize=0;
        char *data = alnDbr->getData(i);
        while (*data != '\0') {

            currentclustersize++;
            data = Util::skipLine(data);
        }
        maxClustersize=std::max(currentclustersize,maxClustersize);
        clustersizes[i]=currentclustersize;
    }
    //save ordered clustersizes
    orderedClustersizes=new std::set<int>[maxClustersize+1];
    for (int i = 0; i < dbSize; i++) {
        SetCover3::insertCluster(i,clustersizes[i]);
    }
    //delete from beginning
    for (int cl_size = maxClustersize; cl_size > 0; cl_size--) {
        while(orderedClustersizes[cl_size].size()>0) {
            int representative =*orderedClustersizes[cl_size].begin();
//          Debug(Debug::INFO)<<alnDbr->getDbKey(representative)<<"\n";
            orderedClustersizes[cl_size].erase(representative);
            clustersizes[representative]=0;
            assignedcluster[representative]=representative;

            //delete clusters of members;
            char *data = alnDbr->getData(representative);
            while (*data != '\0') {
                Util::parseKey(data, idbuffer1);
                int elementtodelete=alnDbr->getId(idbuffer1);
                bool representativefound=false;
                //
                Util::parseByColumnNumber(data, similarity, 4); //column 4 = sequence identity
                float seqId = atof(similarity);
                if(seqId>bestscore[elementtodelete]) {
                    assignedcluster[elementtodelete] = representative;
                    bestscore[elementtodelete] = seqId;
                }
                if(elementtodelete== representative){
                    data = Util::skipLine(data);
                    continue;
                }
                if(clustersizes[elementtodelete]<1){
                    data = Util::skipLine(data);
                    continue;
                }
                char *data2 = alnDbr->getData(elementtodelete);
                orderedClustersizes[clustersizes[elementtodelete]].erase(elementtodelete);
                clustersizes[elementtodelete]=0;
                //decrease clustersize of sets that contain the element
                while (*data2 != '\0') {
                    Util::parseKey(data2, idbuffer2);
                    int elementtodecrease=alnDbr->getId(idbuffer2);
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

                    data2 = Util::skipLine(data2);
                }
                if(!representativefound){
                    Debug(Debug::ERROR)<<"error with cluster:\t"<<alnDbr->getDbKey(representative)<<"\tis not contained in set:\t"<<alnDbr->getDbKey(elementtodelete)<<".\n";
                }
                data = Util::skipLine(data);
            }Util::parseKey(data, idbuffer1);


        }

    }


    for(size_t i = 0; i < n; i++) {
        AffinityClustering::add_to_set(seqDbr->getId(alnDbr->getDbKey(i)),&sets[seqDbr->getId(alnDbr->getDbKey(assignedcluster[i]))],seqDbr->getId(alnDbr->getDbKey(assignedcluster[i])));
    }
    for(size_t i = 0; i < n; i++) {
        set * max_set = &sets[i];
        if (max_set->elements == NULL)
            continue;
        result.push_back(max_set); // O(1)
    }
    return result;

}


void SetCover3::insertCluster(int clusterid, int clustersize){
    //if(orderedClustersizes[clustersize] ==NULL){
    orderedClustersizes[clustersize].insert(clusterid);
}
