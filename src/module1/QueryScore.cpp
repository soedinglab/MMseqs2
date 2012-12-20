#include "QueryScore.h"

QueryScore::QueryScore (size_t dbSize, short prefThreshold){

    this->dbSize = dbSize;
    this->prefThreshold = prefThreshold;
    
    this->scores = new short[dbSize];
    memset (scores, 0, sizeof(short) * dbSize);

    this->lastMatchPos = new short[dbSize];
    memset (lastMatchPos, 0, sizeof(short) * dbSize);

    this->hitList = new std::list<size_t>();
    this->resList = new std::list<hit_t>();
}

QueryScore::~QueryScore (){
    delete[] scores;
    delete[] lastMatchPos;
    delete hitList;
    delete resList;
}

void QueryScore::addScores (size_t* hitList, size_t hitListSize, short score){
    size_t seqId;
    for (int i = 0; i < hitListSize; i++){
        seqId = hitList[i];
        scores[seqId] += score;
        if (scores[seqId] >= this->prefThreshold)
            addElementToResults(seqId);
    }
}

void QueryScore::addElementToResults (size_t seqId){
    std::list<size_t>::iterator it;
    it = lower_bound(hitList->begin(), hitList->end(), seqId);
    if (*it != seqId)
        this->hitList->insert(it, seqId);
}

std::list<hit_t>* QueryScore::getResult (){
    std::list<size_t>::iterator it;
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
}
