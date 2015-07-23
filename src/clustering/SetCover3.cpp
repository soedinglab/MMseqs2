//
// Created by lars on 08.06.15.
//

#include <Util.h>
#include <Debug.h>
#include "SetCover3.h"
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
    std::map<int,int>missed_ids; //save representative if symmetry condittion is not fulfilled
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
    for (int cl_size = maxClustersize-1; cl_size > 0; cl_size--) {
        while(orderedClustersizes[cl_size].size()>0) {
            int representative =*orderedClustersizes[cl_size].begin();
//            Debug(Debug::INFO)<<alnDbr->getDbKey(representative)<<"\n";
            orderedClustersizes[cl_size].erase(representative);
            representatives.insert(representative);
            clustersizes[representative]=0;
            //delete clusters of members;
            char *data = alnDbr->getData(representative);
            while (*data != '\0') {
                Util::parseKey(data, idbuffer1);
                int elementtodelete=alnDbr->getId(idbuffer1);
                bool representativefound=false;
                if(elementtodelete== representative|| clustersizes[elementtodelete]<1){
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
                    missed_ids.insert(std::make_pair(elementtodelete,representative));
              //      Debug(Debug::ERROR)<<"error with cluster:\t"<<alnDbr->getDbKey(representative)<<"\tis not contained in set:\t"<<alnDbr->getDbKey(elementtodelete)<<".\n";
                }
                data = Util::skipLine(data);
            }Util::parseKey(data, idbuffer1);


        }

    }
    //assign to cluster representatives
    //compute result
    set* sets = new set[representatives.size()];
    memset(sets, 0, sizeof(set *)*(representatives.size()+1));
    std::map<int,set*>sets_mapping;
    int position=0;
    for(int i :representatives){
        sets_mapping.insert(std::make_pair(i,&sets[position]));
        position++;
    }
    char *idbuffer3 = new char[255 + 1];
    for (int i = 0; i < dbSize; i++) {
        char *data = alnDbr->getData(i);
        int result=-1;
        while (*data != '\0') {
            Util::parseKey(data, idbuffer3);
            if(sets_mapping.find(alnDbr->getId(idbuffer3)) !=sets_mapping.end()){
                result=alnDbr->getId(idbuffer3);
                break;
            }
            data = Util::skipLine(data);
        }
        if(result==-1){
            if(missed_ids.find(i) !=missed_ids.end()){
                //translate to seqDBr id
                AffinityClustering::add_to_set(seqDbr->getId(alnDbr->getDbKey(i)),sets_mapping.find(missed_ids.find(i)->second)->second,seqDbr->getId(alnDbr->getDbKey(missed_ids.find(i)->second)));
            }else{
                Debug(Debug::INFO)<<"error with id:\t"<<alnDbr->getDbKey(i)<<"\t it is not assigned to a cluster.\n";
            }

        }else{
            AffinityClustering::add_to_set(seqDbr->getId(alnDbr->getDbKey(i)),sets_mapping.find(result)->second,seqDbr->getId(alnDbr->getDbKey(result)));
        }

    }
    for(int i=0;i<representatives.size();i++) {
        set* resultset=&sets[i];
        result.push_back(resultset); // O(1)
    }
    return result;

}


void SetCover3::insertCluster(int clusterid, int clustersize){
        //if(orderedClustersizes[clustersize] ==NULL){
            orderedClustersizes[clustersize].insert(clusterid);

    
}
