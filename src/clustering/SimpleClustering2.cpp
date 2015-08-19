//
// Created by lars on 28.07.15.
//

#include <Log.h>
#include "SimpleClustering2.h"
#include "Util.h"
#include "Debug.h"
#include "AffinityClustering.h"
SimpleClustering2::SimpleClustering2(DBReader * seqDbr, DBReader * alnDbr, float seqIdThr, float coverage){

    this->seqDbr=seqDbr;
    this->alnDbr=alnDbr;
    this->seqIdThr=seqIdThr;
    this->coverage=coverage;

}


std::list<set *> SimpleClustering2::execute(){
    std::list<set *> result;
    size_t n = seqDbr->getSize();
    int* assignedcluster=new int[n];
    memset(assignedcluster, -1, sizeof(int)*(n));
    float* bestscore=new float[n];
    memset(bestscore, -10, sizeof(float)*(n));
    char *similarity = new char[255+1];
    char *idbuffer1 = new char[255 + 1];


    for(size_t i = 0; i < n; i++) {
        // seqDbr is descending sorted by length
        // the assumption is that clustering is B -> B (not A -> B)
        Log::printProgress(i);
        if(assignedcluster[i]==-1){
            char *clusterId = seqDbr->getDbKey(i);
            char *data = alnDbr->getDataByDBKey(clusterId);
            while (*data != '\0') {
                Util::parseKey(data, idbuffer1);
                //filter by alignment thresholds
                Util::parseByColumnNumber(data, similarity, 4); //column 4 = sequence identity
                float seqId = atof(similarity);
                    if (seqId < this->seqIdThr) {
                    data = Util::skipLine(data);
                    continue;
                }
                int id=seqDbr->getId(idbuffer1);
                if(seqId>bestscore[id]){
                    assignedcluster[id]=i;
                    bestscore[id]=seqId;
                }


                data = Util::skipLine(data);
            }

        }




    }

    set* sets = new set[n];
    memset(sets, 0, sizeof(set *)*(n));
    for(size_t i = 0; i < n; i++) {
        AffinityClustering::add_to_set(i,&sets[assignedcluster[i]],assignedcluster[i]);
    }
    for(size_t i = 0; i < n; i++) {
        set * max_set = &sets[i];
        if (max_set->elements == NULL)
            continue;
        result.push_back(max_set); // O(1)
    }
    Debug(Debug::ERROR)<<result.size()<<"\n";
   // delete[] sets;
    delete[] idbuffer1;
    delete[] bestscore;
    delete[] assignedcluster;
    delete[] similarity;
    return result;
}

