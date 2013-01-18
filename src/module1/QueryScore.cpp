#include "QueryScore.h"

QueryScore::QueryScore (int dbSize, short prefThreshold){

    this->dbSize = dbSize;
    this->prefThreshold = prefThreshold;
    
    this->scores = new short[dbSize];
    memset (scores, 0, sizeof(short) * dbSize);

    this->lastMatchPos = new short[dbSize];
    memset (lastMatchPos, 0, sizeof(short) * dbSize);

    this->hitList = new std::list<int>();
    this->resList = new std::list<hit_t>();

    this->currQueryPos = 0;
}

QueryScore::~QueryScore (){
    delete[] scores;
    delete[] lastMatchPos;
    delete hitList;
    delete resList;
}

void QueryScore::moveToNextQueryPos(){
    this->currQueryPos++;
}

void QueryScore::addScores (int* hitList, int hitListSize, short score){
    int dbSeqId;
    for (int i = 0; i < hitListSize; i++){
        dbSeqId = hitList[i];
        // this position in the query sequence already matched this db sequence
        if (this->currQueryPos > this->lastMatchPos[dbSeqId]){
            scores[dbSeqId] += score;
            this->lastMatchPos[dbSeqId] = this->currQueryPos;
            if (scores[dbSeqId] >= this->prefThreshold)
                addElementToResults(dbSeqId);
        }
    }
}

void QueryScore::addElementToResults (int seqId){
    std::list<int>::iterator it;
    it = lower_bound(hitList->begin(), hitList->end(), seqId);
    // until now, sequence is not in the hitList
    if (*it != seqId)
        this->hitList->insert(it, seqId);
}

std::list<hit_t>* QueryScore::getResult (){
    std::list<int>::iterator it;
    for (it = hitList->begin(); it != hitList->end(); it++){
        hit_t hit = {*it, scores[*it]};
        resList->push_back(hit); 
    }
    return resList;
}

void QueryScore::reset(){
    memset (scores, 0, sizeof(short) * dbSize);
    memset (lastMatchPos, 0, sizeof(short) * dbSize);
    resList->clear();
    currQueryPos = 0;
}
