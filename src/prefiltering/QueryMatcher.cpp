//
// Created by mad on 5/26/15.
//
#include <new>
#include <SubstitutionMatrix.h>
#include "QueryMatcher.h"
#include "Util.h"


hit_t parsePrefilterHit(char* data)
{
    hit_t result;
    char *wordCnt[255];
    size_t cols = Util::getWordsOfLine(data, wordCnt, 254);
    if (cols>=3)
    {
        result.seqId = std::stoul(wordCnt[0],NULL,10);
        result.prefScore = std::stol(wordCnt[1],NULL,10);
        result.diagonal = std::stol(wordCnt[2],NULL,10);
    } else { //error
        result.seqId = -1;
    }
    return result;
}


std::string prefilterHitToString(hit_t h)
{
    std::ostringstream resStream;
    resStream << h.seqId << '\t' << (int)h.prefScore << '\t' << h.diagonal << '\n';
    return resStream.str();
}



#define FE_1(WHAT, X) WHAT(X)
#define FE_2(WHAT, X, ...) WHAT(X)FE_1(WHAT, __VA_ARGS__)
#define FE_3(WHAT, X, ...) WHAT(X)FE_2(WHAT, __VA_ARGS__)
#define FE_4(WHAT, X, ...) WHAT(X)FE_3(WHAT, __VA_ARGS__)
#define FE_5(WHAT, X, ...) WHAT(X)FE_4(WHAT, __VA_ARGS__)
#define FE_6(WHAT, X, ...) WHAT(X)FE_5(WHAT, __VA_ARGS__)
#define FE_7(WHAT, X, ...) WHAT(X)FE_6(WHAT, __VA_ARGS__)
#define FE_8(WHAT, X, ...) WHAT(X)FE_7(WHAT, __VA_ARGS__)
#define FE_9(WHAT, X, ...) WHAT(X)FE_8(WHAT, __VA_ARGS__)
#define FE_10(WHAT, X, ...) WHAT(X)FE_9(WHAT, __VA_ARGS__)
#define FE_11(WHAT, X, ...) WHAT(X)FE_10(WHAT, __VA_ARGS__)

#define GET_MACRO(_1,_2,_3,_4,_5,_6,_7,_8,_9,_10,_11,NAME,...) NAME
#define FOR_EACH(action,...) \
  GET_MACRO(__VA_ARGS__,FE_11,FE_10,FE_9,FE_8,FE_7,FE_6,FE_5,FE_4,FE_3,FE_2,FE_1)(action,__VA_ARGS__)

QueryMatcher::QueryMatcher(BaseMatrix *m, IndexTable *indexTable,
                                               unsigned int *seqLens, short kmerThr,
                                               double kmerMatchProb, int kmerSize, size_t dbSize,
                                               unsigned int maxSeqLen, unsigned int effectiveKmerSize,
                                               size_t maxHitsPerQuery, bool aaBiasCorrection,
                                               bool diagonalScoring, unsigned int minDiagScoreThr)
{
    this->m = m;
    this->indexTable = indexTable;
    this->kmerSize = kmerSize;
    this->kmerThr = kmerThr;
    this->kmerGenerator = new KmerGenerator(kmerSize, m->alphabetSize, kmerThr);
    this->aaBiasCorrection = aaBiasCorrection;
    this->stats = new statistics_t();
    // assure that the whole database can be matched (extreme case)
    // this array will need 500 MB for 50 Mio. sequences ( dbSize * 2 * 5byte)
    this->dbSize = dbSize;
    this->counterResultSize = std::max((size_t)1000000, dbSize);
    this->maxDbMatches = dbSize * 2;
    this->resList = (hit_t *) mem_align(ALIGN_INT, MAX_RES_LIST_LEN * sizeof(hit_t) );
    this->databaseHits = new(std::nothrow) IndexEntryLocal[maxDbMatches];
    memset(databaseHits, 0, sizeof(IndexEntryLocal) * maxDbMatches);
    Util::checkAllocation(databaseHits, "Could not allocate databaseHits memory in QueryMatcher");
    this->foundDiagonals = new(std::nothrow) CounterResult[counterResultSize];
    memset(foundDiagonals, 0, sizeof(CounterResult) * counterResultSize);
    Util::checkAllocation(foundDiagonals, "Could not allocate foundDiagonals memory in QueryMatcher");
    this->lastSequenceHit = this->databaseHits + maxDbMatches;
    this->indexPointer = new(std::nothrow) IndexEntryLocal*[maxSeqLen + 1];
    Util::checkAllocation(indexPointer, "Could not allocate indexPointer memory in QueryMatcher");
    this->diagonalScoring = diagonalScoring;
    this->minDiagScoreThr = minDiagScoreThr;
    // data for histogram of score distribution
    this->scoreSizes = new unsigned int[SCORE_RANGE];
    memset(scoreSizes, 0, SCORE_RANGE * sizeof(unsigned int));
    this->maxHitsPerQuery = maxHitsPerQuery;
    // this array will need 128 * (maxDbMatches / 128) * 5byte ~ 500MB for 50 Mio. Sequences
    initDiagonalMatcher(dbSize, maxDbMatches);
//    this->diagonalMatcher = new CacheFriendlyOperations(dbSize, maxDbMatches / 128 );
    // needed for p-value calc.
    this->mu = kmerMatchProb;
    this->logMatchProb = log(kmerMatchProb);
    this->logScoreFactorial = new double[SCORE_RANGE];
    MathUtil::computeFactorial(logScoreFactorial, SCORE_RANGE);

    // initialize sequence lenghts with each seqLens[i] = L_i - k + 1
    this->seqLens = new float[dbSize];
    memset (this->seqLens, 0, dbSize * sizeof(float));
    for (size_t i = 0; i < dbSize; i++){
        if (seqLens[i] > (effectiveKmerSize - 1))
            this->seqLens[i] = static_cast<float>(seqLens[i] - effectiveKmerSize + 1);
        else
            this->seqLens[i] = 1.0f;
    }
    compositionBias = new float[maxSeqLen];
    ungappedAlignment = NULL;
    if(this->diagonalScoring == true) {
        ungappedAlignment = new UngappedAlignment(maxSeqLen, m, indexTable->getSequenceLookup());
    }
}

QueryMatcher::~QueryMatcher(){
    deleteDiagonalMatcher(activeCounter);
    free(resList);
    delete [] scoreSizes;
    delete [] databaseHits;
    delete [] indexPointer;
    delete [] foundDiagonals;
    delete [] logScoreFactorial;
    delete [] seqLens;
    delete [] compositionBias;
    if(ungappedAlignment != NULL){
        delete ungappedAlignment;
    }
    delete stats;
    delete kmerGenerator;
}

size_t QueryMatcher::evaluateBins(IndexEntryLocal **hitsByIndex,
                                            CounterResult *output,
                                            size_t outputSize,
                                            unsigned short indexFrom,
                                            unsigned short indexTo,
                                            bool computeTotalScore) {
    size_t localResultSize = 0;
#define COUNT_CASE(x) case x: localResultSize += cachedOperation##x->countElements(hitsByIndex, output, outputSize, indexFrom, indexTo, computeTotalScore); break;
    switch (activeCounter){
        FOR_EACH(COUNT_CASE,2,4,8,16,32,64,128,256,512,1024,2048)
    }
#undef COUNT_CASE
    return localResultSize;
}

std::pair<hit_t *, size_t> QueryMatcher::matchQuery (Sequence * seq, unsigned int identityId){
    seq->resetCurrPos();
//    std::cout << "Id: " << seq->getId() << std::endl;
    memset(scoreSizes, 0, SCORE_RANGE * sizeof(unsigned int));

    // bias correction
    if(aaBiasCorrection == true){
        if(seq->getSeqType() == Sequence::AMINO_ACIDS) {
            SubstitutionMatrix::calcLocalAaBiasCorrection(m, seq->int_sequence, seq->L, compositionBias);
        }else{
            memset(compositionBias, 0, sizeof(float) * seq->L);
        }
    } else {
        memset(compositionBias, 0, sizeof(float) * seq->L);
    }

    size_t resultSize = match(seq, compositionBias);
    std::pair<hit_t *, size_t > queryResult;
    if(diagonalScoring == true) {
        // write diagonal scores in count value
        ungappedAlignment->processQuery(seq, compositionBias, foundDiagonals, resultSize, 0);
        memset(scoreSizes, 0, SCORE_RANGE * sizeof(unsigned int));
        resultSize = keepMaxScoreElementOnly(foundDiagonals, resultSize);
        updateScoreBins(foundDiagonals, resultSize);
        unsigned int diagonalThr = computeScoreThreshold(scoreSizes, this->maxHitsPerQuery);
        diagonalThr = std::max(minDiagScoreThr, diagonalThr);
        queryResult = getResult(foundDiagonals, resultSize, seq->L, identityId, diagonalThr, true);
    }else{
        unsigned int thr = computeScoreThreshold(scoreSizes, this->maxHitsPerQuery);
        queryResult = getResult(foundDiagonals, resultSize, seq->L, identityId, thr, false);
    }
    if(queryResult.second > 1){
        if (identityId != UINT_MAX){
            if(diagonalScoring == true) {
                std::sort(resList + 1, resList + queryResult.second, hit_t::compareHitsByDiagonalScore);
            }else {
                std::sort(resList + 1, resList + queryResult.second, hit_t::compareHitsByPValue);
            }
        } else{
            if(diagonalScoring == true) {
                std::sort(resList, resList + queryResult.second, hit_t::compareHitsByDiagonalScore);
            } else {
                std::sort(resList, resList + queryResult.second, hit_t::compareHitsByPValue);
            }
        }
    }
    return queryResult;
}

size_t QueryMatcher::match(Sequence *seq, float *compositionBias) {
    // go through the query sequence
    size_t kmerListLen = 0;
    size_t numMatches = 0;
    size_t overflowNumMatches = 0;
    size_t overflowHitCount = 0;
    //size_t pos = 0;
    stats->diagonalOverflow = false;
    IndexEntryLocal* sequenceHits = databaseHits;
    size_t seqListSize;
    unsigned short indexStart = 0;
    unsigned short indexTo = 0;
    while(seq->hasNextKmer()){
        const int* kmer = seq->nextKmer();
        const int * pos = seq->getKmerPositons();
        float biasCorrection = 0;
        for (int i = 0; i < kmerSize; i++){
            biasCorrection += compositionBias[pos[i]];
        }
        // round bias to next higher or lower value
        short bias = static_cast<short>((biasCorrection < 0.0) ? biasCorrection - 0.5: biasCorrection + 0.5);
        short kmerMatchScore = std::max(kmerThr - bias, 0);
        // adjust kmer threshold based on composition bias
        kmerGenerator->setThreshold(kmerMatchScore);
        const ScoreMatrix kmerList = kmerGenerator->generateKmerList(kmer);
        kmerListLen += kmerList.elementSize;
        const unsigned short current_i = seq->getCurrentPosition();
        indexPointer[current_i] = sequenceHits;
        // match the index table
        for (unsigned int kmerPos = 0; kmerPos < kmerList.elementSize; kmerPos++) {
            // generate k-mer list
            const IndexEntryLocal *entries = indexTable->getDBSeqList<IndexEntryLocal>(kmerList.index[kmerPos],
                                                                                       &seqListSize);
            // detected overflow while matching
            if ((sequenceHits + seqListSize) >= lastSequenceHit) {
                stats->diagonalOverflow = true;
                // last pointer
                indexPointer[current_i + 1] = sequenceHits;
//                std::cout << "Overflow in i=" << indexStart << " IndexTo=" << i << std::endl;
                const size_t hitCount = evaluateBins(indexPointer,
                                                     foundDiagonals + overflowHitCount,
                                                     counterResultSize - overflowHitCount,
                                                     indexStart, current_i, (diagonalScoring == false));
                if(overflowHitCount != 0){ //merge lists
                    // hitCount is max. dbSize so there can be no overflow in mergeElemens
                    overflowHitCount = mergeElements(diagonalScoring, foundDiagonals, overflowHitCount +  hitCount);

                } else {
                    overflowHitCount = hitCount;
                }
                // reset pointer position
                sequenceHits = databaseHits;
                indexPointer[current_i] = databaseHits;
                indexStart = current_i;
                overflowNumMatches += numMatches;
                numMatches = 0;
                if((sequenceHits + seqListSize) >= lastSequenceHit){
                    goto outer;
                }
            };
            memcpy(sequenceHits, entries, sizeof(IndexEntryLocal) * seqListSize);
            sequenceHits += seqListSize;
            numMatches += seqListSize;
        }
        indexTo = current_i;
    }
    outer:
    indexPointer[indexTo + 1] = databaseHits + numMatches;
    size_t hitCount = evaluateBins(indexPointer, foundDiagonals + overflowHitCount,
                                   counterResultSize - overflowHitCount, indexStart, indexTo,  (diagonalScoring == false));
    //fill the output
    if(overflowHitCount != 0){ // overflow occurred
        hitCount = mergeElements(diagonalScoring, foundDiagonals, overflowHitCount + hitCount);
    }
    stats->doubleMatches = 0;
    if(diagonalScoring == false) {
        updateScoreBins(foundDiagonals, hitCount);
        stats->doubleMatches = getDoubleDiagonalMatches();
    }
    stats->kmersPerPos   = ((double)kmerListLen/(double)seq->L);
    stats->querySeqLen   = seq->L;
    stats->dbMatches     = overflowNumMatches + numMatches;
    return hitCount;
}

size_t QueryMatcher::getDoubleDiagonalMatches(){
    size_t retValue = 0;
    for(size_t i = 1; i < SCORE_RANGE; i++){
        retValue += scoreSizes[i] * i;
        //std::cout << scoreSizes[i] * i << std::endl;
    }
    return retValue;
}

void QueryMatcher::updateScoreBins(CounterResult *result, size_t elementCount) {
    for(size_t i = 0; i < elementCount; i++){
        scoreSizes[result[i].count]++;
    }
}

std::pair<hit_t *, size_t>  QueryMatcher::getResult(CounterResult * results,
                                                              size_t resultSize, const int l,
                                                              const unsigned int id,
                                                              const unsigned short thr,
                                                              const bool diagonalScoring) {
    size_t elementCounter = 0;
    if (id != UINT_MAX){
        hit_t * result = (resList + 0);
        const unsigned short rawScore  = SCORE_RANGE-1;
        result->seqId = id;
        result->prefScore = rawScore;
        result->diagonal = 0;
        //result->pScore = (((float)rawScore) - mu)/ sqrtMu;
        result->pScore = (diagonalScoring) ? 0.0 : -computeLogProbability(rawScore, seqLens[id],
                                                                          mu, logMatchProb, logScoreFactorial[rawScore]);
        elementCounter++;
    }

    for (size_t i = 0; i < resultSize; i++) {
        const unsigned int seqIdCurr = results[i].id;
        const unsigned int scoreCurr = results[i].count;
        const unsigned int diagCurr  = results[i].diagonal;
        // write result to list
        if(scoreCurr >= thr && id != seqIdCurr){
            hit_t *result = (resList + elementCounter);
            result->seqId = seqIdCurr;
            result->prefScore = scoreCurr;
            result->diagonal = diagCurr;
            //printf("%d\t%d\t%f\t%f\t%f\t%f\t%f\n", result->seqId, scoreCurr, seqLens[seqIdCurr], mu, logMatchProb, logScoreFactorial[scoreCurr]);
            result->pScore =  (diagonalScoring) ? 0.0 :  -computeLogProbability(scoreCurr, seqLens[seqIdCurr],
                                                                                mu, logMatchProb, logScoreFactorial[scoreCurr]);
            elementCounter++;
            if (elementCounter >= MAX_RES_LIST_LEN)
                break;
        }
    }
    std::pair<hit_t *, size_t>  pair = std::make_pair(this->resList, elementCounter);
    return pair;
}

void QueryMatcher::initDiagonalMatcher(size_t dbsize, unsigned int maxDbMatches) {
#define INIT(x)   cachedOperation##x = new CacheFriendlyOperations<x>(dbsize, maxDbMatches/x); \
                  activeCounter = x;
    if(dbsize/2 < L2_CACH_SIZE){
        INIT(2)
    }else if(dbsize/4 < L2_CACH_SIZE){
        INIT(4)
    }else if(dbsize/8 < L2_CACH_SIZE){
        INIT(8)
    }else if(dbsize/16 < L2_CACH_SIZE){
        INIT(16)
    }else if(dbsize/32 < L2_CACH_SIZE){
        INIT(32)
    }else if(dbsize/64 < L2_CACH_SIZE){
        INIT(64)
    }else if(dbsize/128 < L2_CACH_SIZE){
        INIT(128)
    }else if(dbsize/256 < L2_CACH_SIZE){
        INIT(256)
    }else if(dbsize/512 < L2_CACH_SIZE){
        INIT(512)
    }else if(dbsize/1024 < L2_CACH_SIZE){
        INIT(1024)
    }else {
        INIT(2048)
    }
#undef INIT
}

void QueryMatcher::deleteDiagonalMatcher(unsigned int activeCounter){
#define DELETE_CASE(x) case x: delete cachedOperation##x; break;
    switch (activeCounter){
        FOR_EACH(DELETE_CASE,2,4,8,16,32,64,128,256,512,1024,2048)
    }
#undef DELETE_CASE
}

size_t QueryMatcher::mergeElements(bool diagonalScoring, CounterResult *foundDiagonals, size_t hitCounter) {
    size_t overflowHitCount = 0;
#define MERGE_CASE(x) \
    case x: overflowHitCount = (diagonalScoring == true) ? \
                                cachedOperation##x->mergeElementsByDiagonal(foundDiagonals,hitCounter) : \
                                cachedOperation##x->mergeElementsByScore(foundDiagonals,hitCounter); \
    break;

    switch (activeCounter){
        FOR_EACH(MERGE_CASE,2,4,8,16,32,64,128,256,512,1024,2048)
    }
#undef MERGE_CASE
    return overflowHitCount;
}

size_t QueryMatcher::keepMaxScoreElementOnly(CounterResult *foundDiagonals, size_t resultSize) {
    size_t retSize = 0;
#define MAX_CASE(x) case x: retSize = cachedOperation##x->keepMaxScoreElementOnly(foundDiagonals, resultSize); break;
    switch (activeCounter){
        FOR_EACH(MAX_CASE,2,4,8,16,32,64,128,256,512,1024,2048)
    }
#undef MAX_CASE
    return retSize;
}
#undef FOR_EACH
#undef GET_MACRO
#undef FE_11
#undef FE_10
#undef FE_9
#undef FE_8
#undef FE_7
#undef FE_6
#undef FE_5
#undef FE_4
#undef FE_3
#undef FE_2
#undef FE_1