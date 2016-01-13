//
// Created by mad on 5/26/15.
//
#include <SubstitutionMatrix.h>
#include "QueryTemplateLocalFast.h"
#include "QueryScoreLocal.h"
#include "QueryScore.h"

QueryTemplateLocalFast::QueryTemplateLocalFast(BaseMatrix *m, IndexTable *indexTable,
                                               unsigned int *seqLens, short kmerThr,
                                               double kmerMatchProb, int kmerSize, size_t dbSize,
                                               unsigned int maxSeqLen, unsigned int effectiveKmerSize,
                                               size_t maxHitsPerQuery, bool aaBiasCorrection,
                                               bool diagonalScoring)
        : QueryTemplateMatcher(m, indexTable, seqLens,
                               kmerThr, kmerMatchProb, kmerSize,
                               dbSize, aaBiasCorrection, maxSeqLen) {
    // assure that the whole database can be matched (extreme case)
    // this array will need 500 MB for 50 Mio. sequences ( dbSize * 2 * 5byte)
    this->maxDbMatches = dbSize * 2;
    this->dbSize = dbSize;
    this->resList = (hit_t *) mem_align(ALIGN_INT, QueryScore::MAX_RES_LIST_LEN * sizeof(hit_t) );
    this->databaseHits = new IndexEntryLocal[maxDbMatches];
    this->foundDiagonals = new CounterResult[dbSize];
    this->lastSequenceHit = this->databaseHits + maxDbMatches;
    this->indexPointer = new IndexEntryLocal*[maxSeqLen + 1];
    this->diagonalScoring = diagonalScoring;
    // data for histogram of score distribution
    this->scoreSizes = new unsigned int[QueryScoreLocal::SCORE_RANGE];
    memset(scoreSizes, 0, QueryScoreLocal::SCORE_RANGE * sizeof(unsigned int));
    this->maxHitsPerQuery = maxHitsPerQuery;
    // this array will need 128 * (maxDbMatches / 128) * 5byte ~ 500MB for 50 Mio. Sequences
    this->counter = new CountInt32Array(dbSize, maxDbMatches / 128 );
    // needed for p-value calc.
    this->mu = kmerMatchProb;
    this->logMatchProb = log(kmerMatchProb);
    this->logScoreFactorial = new double[QueryScoreLocal::SCORE_RANGE];
    QueryScoreLocal::computeFactorial(logScoreFactorial, QueryScoreLocal::SCORE_RANGE);

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
    diagonalMatcher = NULL;
    if(this->diagonalScoring == true) {
        diagonalMatcher = new DiagonalMatcher(maxSeqLen, m, indexTable->getSequenceLookup());
    }
}

QueryTemplateLocalFast::~QueryTemplateLocalFast(){
    free(resList);
    delete [] scoreSizes;
    delete [] databaseHits;
    delete [] indexPointer;
    delete [] foundDiagonals;
    delete counter;
    delete [] logScoreFactorial;
    delete [] seqLens;
    delete [] compositionBias;
    if(diagonalMatcher != NULL){
        delete diagonalMatcher;
    }
}

size_t QueryTemplateLocalFast::evaluateBins(IndexEntryLocal **hitsByIndex, CounterResult *output, size_t outputSize,
                                            unsigned short indexFrom, unsigned short indexTo, bool diagonalScoring) {
    size_t localResultSize = 0;
    if(diagonalScoring == true){
        localResultSize += counter->extractDiagonals(hitsByIndex, output, outputSize, indexFrom, indexTo);
    } else{
        localResultSize += counter->scoreDiagonals(hitsByIndex, output, outputSize, indexFrom, indexTo);
    }
    return localResultSize;
}

std::pair<hit_t *, size_t> QueryTemplateLocalFast::matchQuery (Sequence * seq, unsigned int identityId){
    seq->resetCurrPos();
//    std::cout << "Id: " << seq->getId() << std::endl;
    memset(scoreSizes, 0, QueryScoreLocal::SCORE_RANGE * sizeof(unsigned int));

    // bias correction
    if(aaBiasCorrection == true){
        SubstitutionMatrix::calcLocalAaBiasCorrection(m, seq->int_sequence, seq->L, compositionBias);
    } else {
        memset(compositionBias, 0, sizeof(float) * seq->L);
    }

    size_t resultSize = match(seq, compositionBias);
    std::pair<hit_t *, size_t > queryResult;
    if(diagonalScoring == true) {
        // write diagonal scores in count value
        diagonalMatcher->processQuery(seq, compositionBias, foundDiagonals, resultSize, 0);
        memset(scoreSizes, 0, QueryScoreLocal::SCORE_RANGE * sizeof(unsigned int));
        resultSize = counter->keepMaxScoreElementOnly(foundDiagonals, resultSize);
        updateScoreBins(foundDiagonals, resultSize);
        unsigned int diagonalThr = QueryScoreLocal::computeScoreThreshold(scoreSizes, this->maxHitsPerQuery);
        diagonalThr = std::max((unsigned int)30, diagonalThr);
        queryResult = getResult(foundDiagonals, resultSize, seq->L, identityId, diagonalThr, true);
    }else{
        unsigned int thr = QueryScoreLocal::computeScoreThreshold(scoreSizes, this->maxHitsPerQuery);
        queryResult = getResult(foundDiagonals, resultSize, seq->L, identityId, thr, false);
    }
    if(queryResult.second > 1){
        if (identityId != UINT_MAX){
            if(diagonalScoring == true) {
                std::sort(resList + 1, resList + queryResult.second, QueryScore::compareHitsByDiagonalScore);
            }else {
                std::sort(resList + 1, resList + queryResult.second, QueryScore::compareHitsByPValue);
            }
        } else{
            if(diagonalScoring == true) {
                std::sort(resList, resList + queryResult.second, QueryScore::compareHitsByDiagonalScore);
            } else {
                std::sort(resList, resList + queryResult.second, QueryScore::compareHitsByPValue);

            }
        }
    }
    return queryResult;
}

size_t QueryTemplateLocalFast::match(Sequence *seq, float *compositionBias) {
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
                const size_t hitCount = evaluateBins(indexPointer, foundDiagonals + overflowHitCount,
                                                     dbSize - overflowHitCount, indexStart, current_i, diagonalScoring);
                if(overflowHitCount != 0){ //merge lists
                    // hitCount is max. dbSize so there can be no overflow in mergeElemens
                    if(diagonalScoring == true) {
                        overflowHitCount = counter->mergeElementsByDiagonal(foundDiagonals,
                                                                            overflowHitCount + hitCount);
                    }else{
                        overflowHitCount = counter->mergeElementsByScore(foundDiagonals,
                                                                         overflowHitCount + hitCount);
                    }
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
    size_t hitCount = evaluateBins(indexPointer, foundDiagonals + overflowHitCount, dbSize - overflowHitCount,
                                   indexStart, indexTo, diagonalScoring);
    //fill the output
    if(overflowHitCount != 0){ // overflow occurred
        if(diagonalScoring == true){
            hitCount = counter->mergeElementsByDiagonal(foundDiagonals,
                                                        overflowHitCount + hitCount);

        }else {
            hitCount = counter->mergeElementsByScore(foundDiagonals,
                                                     overflowHitCount + hitCount);

        }
    }
    updateScoreBins(foundDiagonals, hitCount);
    stats->doubleMatches = getDoubleDiagonalMatches();
    stats->kmersPerPos   = ((double)kmerListLen/(double)seq->L);
    stats->querySeqLen   = seq->L;
    stats->dbMatches     = overflowNumMatches + numMatches;
    return hitCount;
}

size_t QueryTemplateLocalFast::getDoubleDiagonalMatches(){
    size_t retValue = 0;
    for(size_t i = 1; i < QueryScoreLocal::SCORE_RANGE; i++){
        retValue += scoreSizes[i] * i;
        //std::cout << scoreSizes[i] * i << std::endl;
    }
    return retValue;
}

void QueryTemplateLocalFast::updateScoreBins(CounterResult *result, size_t elementCount) {
    for(size_t i = 0; i < elementCount; i++){
        scoreSizes[result[i].count]++;
    }
}

std::pair<hit_t *, size_t>  QueryTemplateLocalFast::getResult(CounterResult * results,
                                                              size_t resultSize, const int l,
                                                              const unsigned int id,
                                                              const unsigned short thr,
                                                              const bool diagonalScoring) {
    size_t elementCounter = 0;
    if (id != UINT_MAX){
        hit_t * result = (resList + 0);
        const unsigned short rawScore  = std::min(QueryScoreLocal::SCORE_RANGE-1, (size_t) l);
        result->seqId = id;
        result->prefScore = rawScore;
        result->diagonal = 0;
        //result->pScore = (((float)rawScore) - mu)/ sqrtMu;
        result->pScore = (diagonalScoring == true) ? 0.0 : -QueryScoreLocal::computeLogProbability(rawScore, seqLens[id],
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
            result->pScore =  (diagonalScoring == true) ? 0.0 :  -QueryScoreLocal::computeLogProbability(scoreCurr, seqLens[seqIdCurr],
                                                                                                 mu, logMatchProb, logScoreFactorial[scoreCurr]);
            elementCounter++;
            if (elementCounter >= QueryScore::MAX_RES_LIST_LEN)
                break;
        }
    }
    std::pair<hit_t *, size_t>  pair = std::make_pair(this->resList, elementCounter);
    return pair;
}
