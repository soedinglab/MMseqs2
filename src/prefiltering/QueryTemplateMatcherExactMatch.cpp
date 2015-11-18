//
// Created by mad on 5/26/15.
//
#include <SubstitutionMatrix.h>
#include "QueryTemplateMatcherExactMatch.h"
#include "QueryScoreLocal.h"
#include "QueryScore.h"

QueryTemplateMatcherExactMatch::QueryTemplateMatcherExactMatch(BaseMatrix *m, IndexTable *indexTable,
                                                               unsigned int *seqLens, short kmerThr,
                                                               double kmerMatchProb, int kmerSize, size_t dbSize,
                                                               unsigned int maxSeqLen, int effectiveKmerSize,
                                                               size_t maxHitsPerQuery, bool aaBiasCorrection)
        : QueryTemplateMatcher(m, indexTable, seqLens, kmerThr, kmerMatchProb, kmerSize, dbSize, aaBiasCorrection, maxSeqLen) {
    // assure that the whole database can be matched (extreme case)
    // this array will need 500 MB for 50 Mio. sequences ( dbSize * 2 * 5byte)
    this->maxDbMatches = dbSize * 4;
    this->binSize =  this->maxDbMatches/BINCOUNT;
    this->resList = (hit_t *) mem_align(ALIGN_INT, QueryScore::MAX_RES_LIST_LEN * sizeof(hit_t) );
    this->binData = new CounterResult[maxDbMatches];
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
            this->seqLens[i] = (float) (seqLens[i] - effectiveKmerSize + 1);
        else
            this->seqLens[i] = 1.0f;
    }
    compositionBias = new float[maxSeqLen];
}

QueryTemplateMatcherExactMatch::~QueryTemplateMatcherExactMatch(){
    free(resList);
    delete [] scoreSizes;
    delete [] binData;
    delete counter;
    delete [] logScoreFactorial;
    delete [] seqLens;
    delete [] compositionBias;
}

size_t QueryTemplateMatcherExactMatch::evaluateBins(CounterResult *inputOutput, size_t N) {
    size_t localResultSize = 0;
    localResultSize += counter->findDuplicates(binData, binSize, diagonalBins, BINCOUNT, binData);
    return localResultSize;
}


std::pair<hit_t *, size_t> QueryTemplateMatcherExactMatch::matchQuery (Sequence * seq, unsigned int identityId){
    seq->resetCurrPos();
    setupBinPointer();
    memset(scoreSizes, 0, QueryScoreLocal::SCORE_RANGE * sizeof(unsigned int));
    size_t resultSize = match(seq);
    unsigned int thr = QueryScoreLocal::computeScoreThreshold(scoreSizes, this->maxHitsPerQuery);
    return getResult(resultSize, seq->L, identityId, thr);
}


void QueryTemplateMatcherExactMatch::setupBinPointer() {
    // Example binCount = 3
    // bin start             |-----------------------|-----------------------| bin end
    //    segments[bin_step][0]
    //                            segments[bin_step][1]
    //                                                    segments[bin_step][2]
    size_t curr_pos = 0;
    for(size_t bin = 0; bin < BINCOUNT; bin++){
        diagonalBins[bin] = binData + curr_pos;
        curr_pos += binSize;
    }
}



size_t QueryTemplateMatcherExactMatch::match(Sequence *seq){
    // go through the query sequence
    size_t kmerListLen = 0;
    size_t numMatches = 0;
    size_t overflowNumMatches = 0;
    size_t overflowHitCount = 0;
    //size_t pos = 0;
    stats->diagonalOverflow = false;
    CounterResult* sequenceHits = binData;
    CounterResult* lastSequenceHit = binData + this->maxDbMatches;
    // bias correction
    if(aaBiasCorrection == true){
        SubstitutionMatrix::calcLocalAaBiasCorrection(m, seq->int_sequence, seq->L, compositionBias);
    } else {
        memset(compositionBias, 0, sizeof(float) * seq->L);
    }
    size_t seqListSize;
    const CounterResult * lastPosition = (binData + BINCOUNT * binSize);
    while(seq->hasNextKmer()){
        const int* kmer = seq->nextKmer();
        const int * pos = seq->getKmerPositons();
        float biasCorrection = 0;
        for (int i = 0; i < kmerSize; i++){
            biasCorrection += compositionBias[pos[i]];
        }
        // round bias to next higher or lower value
        short bias = (short) (biasCorrection < 0.0) ? biasCorrection - 0.5: biasCorrection + 0.5;
        short kmerMatchScore = std::max(kmerThr - bias, 0);
        // adjust kmer threshold based on composition bias
        kmerGenerator->setThreshold(kmerMatchScore);
        const ScoreMatrix kmerList = kmerGenerator->generateKmerList(kmer);
        kmerListLen += kmerList.elementSize;
        const unsigned short i =  seq->getCurrentPosition();
        // match the index table
        for (unsigned int kmerPos = 0; kmerPos < kmerList.elementSize; kmerPos++) {
            // generate k-mer list
            const IndexEntryLocal *entries = indexTable->getDBSeqList<IndexEntryLocal>(kmerList.index[kmerPos],
                                                                                       &seqListSize);
            // detected overflow while matching
            if (sequenceHits + seqListSize >= lastSequenceHit) {
                stats->diagonalOverflow = true;
                sequenceHits = binData + overflowHitCount;
                const size_t hitCount = evaluateBins(sequenceHits, numMatches);
                if(overflowHitCount != 0){ //merge lists
                    // hitCount is max. dbSize so there can be no overflow in mergeElemens
                    overflowHitCount = counter->mergeElements(binData, overflowHitCount + hitCount);
                }else{
                    overflowHitCount = hitCount;
                }
                sequenceHits = binData + overflowHitCount;
                overflowNumMatches += numMatches;
                numMatches = 0;
            };

//            seqListSize = (sequenceHits + seqListSize >= lastSequenceHit) ? 0 : seqListSize;

            for (unsigned int seqIdx = 0; LIKELY(seqIdx < seqListSize); seqIdx++) {
                IndexEntryLocal entry = entries[seqIdx];
                const unsigned char j = entry.position_j;
                const unsigned int binId = entry.seqId & CountInt32Array::MASK_0_5;
                const unsigned int seqId = entry.seqId ;
                const unsigned char diagonal = (i - j);
                diagonalBins[binId]->id = seqId;
                diagonalBins[binId]->count = diagonal;
                diagonalBins[binId] += (diagonalBins[binId] < lastPosition) ? 1 : 0;
            }
            numMatches += seqListSize;
        }
    }
    CounterResult* inputOutputHits = binData + overflowHitCount;
    size_t hitCount = evaluateBins(inputOutputHits, numMatches);
    //fill the output
    if(overflowHitCount != 0){ // overflow occurred
        hitCount = counter->mergeElements(binData, overflowHitCount + hitCount);
    }
    updateScoreBins(binData, hitCount);
    stats->doubleMatches = getDoubleDiagonalMatches();
    stats->kmersPerPos   = ((double)kmerListLen/(double)seq->L);
    stats->querySeqLen   = seq->L;
    stats->dbMatches     = overflowNumMatches + numMatches;
    return hitCount;
}

size_t QueryTemplateMatcherExactMatch::getDoubleDiagonalMatches(){
    size_t retValue = 0;
    for(size_t i = 1; i < QueryScoreLocal::SCORE_RANGE; i++){
        retValue += scoreSizes[i] * i;
        //std::cout << scoreSizes[i] * i << std::endl;
    }
    return retValue;
}

void QueryTemplateMatcherExactMatch::updateScoreBins(CounterResult *result, size_t elementCount) {
    for(size_t i = 0; i < elementCount; i++){
        scoreSizes[result[i].count]++;
    }
}

std::pair<hit_t *, size_t>  QueryTemplateMatcherExactMatch::getResult(size_t resultSize, const int l,
                                                                      const unsigned int id,
                                                                      const unsigned short thr) {
    size_t elementCounter = 0;
    if (id != UINT_MAX){
        hit_t * result = (resList + 0);
        const unsigned short rawScore  = std::min(QueryScoreLocal::SCORE_RANGE-1, (size_t) l);
        result->seqId = id;
        result->prefScore = rawScore;
        //result->pScore = (((float)rawScore) - mu)/ sqrtMu;
        result->pScore = -QueryScoreLocal::computeLogProbability(rawScore, seqLens[id],
                                                                 mu, logMatchProb, logScoreFactorial[rawScore]);
        elementCounter++;
    }

    for (size_t i = 0; i < resultSize; i++) {
        const unsigned int seqIdCurr = binData[i].id;
        const unsigned int scoreCurr = binData[i].count;
        // write result to list
        if(scoreCurr >= thr && id != seqIdCurr){
            hit_t *result = (resList + elementCounter);
            result->seqId = seqIdCurr;
            result->prefScore = scoreCurr;
            //printf("%d\t%d\t%f\t%f\t%f\t%f\t%f\n", result->seqId, scoreCurr, seqLens[seqIdCurr], mu, logMatchProb, logScoreFactorial[scoreCurr]);
            result->pScore = -QueryScoreLocal::computeLogProbability(scoreCurr, seqLens[seqIdCurr],
                                                                     mu, logMatchProb, logScoreFactorial[scoreCurr]);

            elementCounter++;
            if (elementCounter >= QueryScore::MAX_RES_LIST_LEN)
                break;
        }
    }

    if(elementCounter > 1){
        if (id != UINT_MAX){
            std::sort(resList + 1, resList + elementCounter, QueryScore::compareHits);
        }
        else{
            std::sort(resList, resList + elementCounter, QueryScore::compareHits);
        }
    }
    std::pair<hit_t *, size_t>  pair = std::make_pair(this->resList, elementCounter);
    return pair;
}
