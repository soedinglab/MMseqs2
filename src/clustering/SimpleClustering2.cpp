//
// Created by lars on 28.07.15.
//

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
    memset(assignedcluster, -1, sizeof(int)*(n+1));
    float* bestscore=new float[n];
    memset(bestscore, -10, sizeof(float)*(n+1));
    char *similarity = new char[255+1];
    char *idbuffer1 = new char[255 + 1];
    set* sets = new set[n];
    memset(sets, 0, sizeof(set *)*(n+1));

    for(size_t i = 0; i < n; i++) {
        // seqDbr is descending sorted by length
        // the assumption is that clustering is B -> B (not A -> B)
        char *clusterId = seqDbr->getDbKey(i);
        if(assignedcluster[i]==-1){
            char *data = alnDbr->getDataByDBKey(clusterId);
            while (*data != '\0') {
                Util::parseKey(data, idbuffer1);
                //filter by alignment thresholds
                Util::parseByColumnNumber(data, similarity, 4); //column 4 = sequence identity
                float seqId = atof(similarity);
                Util::parseByColumnNumber(data, similarity, 1); //column 1 = alignmentscore
                seqId = atof(std::string(similarity).c_str());
                int queryLength = strlen(seqDbr->getData(i));
                int dbSeqLength = strlen(seqDbr->getDataByDBKey(idbuffer1));
                float maxSeqLength = std::max(queryLength, dbSeqLength);

                //
                seqId = seqId / (maxSeqLength);
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

