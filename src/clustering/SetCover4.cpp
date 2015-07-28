//
// Created by lars on 28.07.15.
//
//
// Created by lars on 08.06.15.
//

#include "SetCover4.h"
#include "Util.h"
#include "Debug.h"
#include "AffinityClustering.h"

SetCover4::SetCover4(DBReader * seqDbr, DBReader * alnDbr, float seqIdThr, float coverage){

    this->seqDbr=seqDbr;
    this->alnDbr=alnDbr;
    this->seqIdThr=seqIdThr;
    this->coverage=coverage;
    this->dbSize=alnDbr->getSize();

}

std::list<set *>  SetCover4::execute() {
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
    clusterweights=new float[dbSize];
    memset(clusterweights,0,dbSize);
    maxClustersize=0;
    maxClusterweight=0;

//read in data and determine biggest cluster
    for (int i = 0; i < dbSize; i++) {
        unsigned int currentclustersize=0;
        float currentclusterweight=0;
        char *data = alnDbr->getData(i);
        while (*data != '\0') {
            Util::parseByColumnNumber(data, similarity, 4); //column 4 = sequence identity
            float seqId = atof(similarity)+offset;
            currentclusterweight+=seqId;
            currentclustersize++;
            data = Util::skipLine(data);
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
            char *data = alnDbr->getData(representative);
            while (*data != '\0') {
                Util::parseKey(data, idbuffer1);
                int elementtodelete=alnDbr->getId(idbuffer1);
                bool representativefound=false;
                //
                Util::parseByColumnNumber(data, similarity, 4); //column 4 = sequence identity
                float seqId = atof(similarity)+offset;
                float oldseqId=bestscore[elementtodelete];
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
                SetCover4::eraseClusterWeight(elementtodelete,clusterweights[elementtodelete]);
                orderedClustersizes[clustersizes[elementtodelete]].erase(elementtodelete);

                clustersizes[elementtodelete]=0;
                clusterweights[elementtodelete]=0;
                //decrease clustersize of sets that contain the element
                while (*data2 != '\0') {
                    Util::parseKey(data2, idbuffer2);
                    int elementtodecrease=alnDbr->getId(idbuffer2);
                    Util::parseByColumnNumber(data2, similarity, 4); //column 4 = sequence identity
                    float seqIdOfElement = atof(similarity)+offset;
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

