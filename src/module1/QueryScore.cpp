#include "QueryScore.h"

QueryScore::QueryScore (int dbSize, float prefThreshold){

    this->dbSize = dbSize;
    this->prefThreshold = prefThreshold;

    this->scores = new int[dbSize];
    memset (scores, 0, sizeof(int) * dbSize);

    //this->lastMatchPos = new short[dbSize];
    //memset (lastMatchPos, 0, sizeof(short) * dbSize);

    this->hitList = new DynamicArray();
    this->resList = new std::list<hit_t>();

    dbFractCnt = 0.0;
    qSeqCnt = 0;

}

QueryScore::~QueryScore (){
    delete[] scores;
    delete hitList;
    delete resList;
}

inline unsigned short sadd16(unsigned short a, unsigned short b)
{
	unsigned int s = (unsigned int)a+b;
	return -(s>>16) | (unsigned short)s;
}

void QueryScore::addScores (int* seqList, int seqListSize, short score){
    for (int i = 0; i < seqListSize; i++){
        scores[seqList[i]] = sadd16(scores[seqList[i]], score);
    }
}

std::list<hit_t>* QueryScore::getResult (int querySeqLen){
    // look for elements in scores above the prefiltering threshold
    int minScore = (int) (prefThreshold * (float)querySeqLen);
    //    int cnt = 0;
    for (int s = 0; s < dbSize; s++){
        //        if (scores[s] > 0)
        //            cnt++;
        if (scores[s] >= minScore){
            hitList->pushBack(s);
        }
    }
    //    dbFractCnt += (float)cnt/(float)dbSize;
    // remove duplicate entries
    hitList->removeDuplicates();

    // write the sequence ids and the corresponding prefiltering scores into the result list
    int* seqList = hitList->getEntries();
    for (int  i = 0; i < hitList->getSize(); i++){
        hit_t hit = {seqList[i], (float)scores[seqList[i]]/(float)querySeqLen};
        resList->push_back(hit); 
    }
    return resList;
}

void QueryScore::reset(){
    memset (scores, 0, sizeof(int) * dbSize);
    resList->clear();
    hitList->clear();
}

void QueryScore::printStats(){
    std::cout << "Average occupancy of the DB scores array: " << dbFractCnt/(double)qSeqCnt << "\n";
}

