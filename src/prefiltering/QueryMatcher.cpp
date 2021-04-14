#include "SubstitutionMatrix.h"
#include "QueryMatcher.h"
#include "FastSort.h"
#include "Util.h"

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

QueryMatcher::QueryMatcher(IndexTable *indexTable, SequenceLookup *sequenceLookup,
                           BaseMatrix *kmerSubMat, BaseMatrix *ungappedAlignmentSubMat,
                           short kmerThr, int kmerSize, size_t dbSize,
                           unsigned int maxSeqLen, size_t maxHitsPerQuery, bool aaBiasCorrection,
                           bool diagonalScoring, unsigned int minDiagScoreThr, bool takeOnlyBestKmer)
                            : idx(indexTable->getAlphabetSize(), kmerSize)
{
    this->kmerSubMat = kmerSubMat;
    this->ungappedAlignmentSubMat = ungappedAlignmentSubMat;
    this->indexTable = indexTable;
    this->kmerSize = kmerSize;
    this->kmerThr = kmerThr;
    this->kmerGenerator = new KmerGenerator(kmerSize, indexTable->getAlphabetSize(), kmerThr);
    this->aaBiasCorrection = aaBiasCorrection;
    this->takeOnlyBestKmer = takeOnlyBestKmer;
    this->stats = new statistics_t();
    // assure that the whole database can be matched (extreme case)
    // this array will need 500 MB for 50 Mio. sequences ( dbSize * 2 * 5byte)
    this->dbSize = dbSize;
    this->foundDiagonalsSize = std::max((size_t)1000000, dbSize);
    this->maxDbMatches = std::max((size_t)1000000, dbSize) * 2;
    // we can never find more hits than dbSize
    this->maxHitsPerQuery = std::min(maxHitsPerQuery, dbSize);
    this->resList = (hit_t *) mem_align(ALIGN_INT, maxHitsPerQuery * sizeof(hit_t) );
    this->databaseHits = new(std::nothrow) IndexEntryLocal[maxDbMatches];
    Util::checkAllocation(databaseHits, "Can not allocate databaseHits memory in QueryMatcher");
    this->foundDiagonals = (CounterResult*)calloc(foundDiagonalsSize, sizeof(CounterResult));
    Util::checkAllocation(foundDiagonals, "Can not allocate foundDiagonals memory in QueryMatcher");
    this->lastSequenceHit = this->databaseHits + maxDbMatches;
    this->indexPointer = new(std::nothrow) IndexEntryLocal*[maxSeqLen + 1];
    Util::checkAllocation(indexPointer, "Can not allocate indexPointer memory in QueryMatcher");
    this->diagonalScoring = diagonalScoring;
    this->minDiagScoreThr = minDiagScoreThr;
    // data for histogram of score distribution
    this->scoreSizes = new unsigned int[SCORE_RANGE];
    memset(scoreSizes, 0, SCORE_RANGE * sizeof(unsigned int));
    // this array will need 128 * (maxDbMatches / 128) * 5byte ~ 500MB for 50 Mio. Sequences
    initDiagonalMatcher(dbSize, maxDbMatches);
//    this->diagonalMatcher = new CacheFriendlyOperations(dbSize, maxDbMatches / 128 );
    // needed for p-value calc.
    ungappedAlignment = NULL;
    if (diagonalScoring) {
        ungappedAlignment = new UngappedAlignment(maxSeqLen, ungappedAlignmentSubMat, sequenceLookup);
    }
    compositionBias = new float[maxSeqLen];
}

QueryMatcher::~QueryMatcher(){
    deleteDiagonalMatcher(activeCounter);
    free(resList);
    delete[] scoreSizes;
    delete[] databaseHits;
    delete[] indexPointer;
    free(foundDiagonals);
    delete[] compositionBias;
    if(ungappedAlignment != NULL){
        delete ungappedAlignment;
    }
    delete stats;
    delete kmerGenerator;
}

std::pair<hit_t*, size_t> QueryMatcher::matchQuery(Sequence *querySeq, unsigned int identityId) {
    querySeq->resetCurrPos();
//    std::cout << "Id: " << querySeq->getId() << std::endl;
    memset(scoreSizes, 0, SCORE_RANGE * sizeof(unsigned int));

    // bias correction
    if(aaBiasCorrection == true){
        if(Parameters::isEqualDbtype(querySeq->getSeqType(), Parameters::DBTYPE_AMINO_ACIDS)) {
            SubstitutionMatrix::calcLocalAaBiasCorrection(kmerSubMat, querySeq->numSequence, querySeq->L, compositionBias);
        }else{
            memset(compositionBias, 0, sizeof(float) * querySeq->L);
        }
    } else {
        memset(compositionBias, 0, sizeof(float) * querySeq->L);
    }

    size_t resultSize = match(querySeq, compositionBias);
    std::pair<hit_t *, size_t> queryResult;
    if (diagonalScoring) {
        // write diagonal scores in count value
        ungappedAlignment->processQuery(querySeq, compositionBias, foundDiagonals, resultSize);
        memset(scoreSizes, 0, SCORE_RANGE * sizeof(unsigned int));

        resultSize = keepMaxScoreElementOnly(foundDiagonals, resultSize);

        updateScoreBins(foundDiagonals, resultSize);
        unsigned int diagonalThr = computeScoreThreshold(scoreSizes, this->maxHitsPerQuery);
        diagonalThr = std::max(minDiagScoreThr, diagonalThr);

        // sort to not lose highest scoring hits if > 150.000 hits are searched
        if(resultSize < foundDiagonalsSize / 2){
            unsigned int maxDiagonalScoreThr = (UCHAR_MAX - ungappedAlignment->getQueryBias());
            bool scoreIsTruncated = (diagonalThr >= maxDiagonalScoreThr) ? true : false;
            size_t elementsCntAboveDiagonalThr = radixSortByScoreSize(scoreSizes, foundDiagonals + resultSize, diagonalThr, foundDiagonals, resultSize);
            if (scoreIsTruncated == true) {
                memset(scoreSizes, 0, SCORE_RANGE * sizeof(unsigned int));
                std::pair<size_t, unsigned int> rescoreResult = rescoreHits(querySeq, scoreSizes, foundDiagonals + resultSize, elementsCntAboveDiagonalThr, ungappedAlignment, maxDiagonalScoreThr);
                size_t newResultSize = rescoreResult.first;
                unsigned int maxSelfScoreMinusDiag = rescoreResult.second;
                elementsCntAboveDiagonalThr = radixSortByScoreSize(scoreSizes, foundDiagonals, 0, foundDiagonals + resultSize, newResultSize);
                queryResult = getResult<UNGAPPED_DIAGONAL_SCORE>(foundDiagonals, elementsCntAboveDiagonalThr, identityId, 0, ungappedAlignment, maxSelfScoreMinusDiag);
            }else{
                queryResult = getResult<UNGAPPED_DIAGONAL_SCORE>(foundDiagonals + resultSize, elementsCntAboveDiagonalThr, identityId, diagonalThr, ungappedAlignment, false);
            }
            stats->truncated = 0;
        }else{
            //Debug(Debug::WARNING) << "Sequence " << querySeq->getDbKey() << " produces too many hits. Results might be truncated\n";
            queryResult = getResult<UNGAPPED_DIAGONAL_SCORE>(foundDiagonals, resultSize, identityId, diagonalThr, ungappedAlignment, false);
            stats->truncated = 1;
        }
    }else{
        unsigned int thr = computeScoreThreshold(scoreSizes, this->maxHitsPerQuery);
        thr = std::max(minDiagScoreThr, thr);
        if(resultSize < foundDiagonalsSize / 2) {
            int elementsCntAboveDiagonalThr = radixSortByScoreSize(scoreSizes, foundDiagonals + resultSize, thr, foundDiagonals, resultSize);
            queryResult = getResult<KMER_SCORE>(foundDiagonals + resultSize, elementsCntAboveDiagonalThr, identityId, thr, ungappedAlignment, false);
            stats->truncated = 0;
        }else{
//            Debug(Debug::WARNING) << "Sequence " << querySeq->getDbKey() << " produces too many hits. Results might be truncated\n";
            queryResult = getResult<KMER_SCORE>(foundDiagonals, resultSize, identityId, thr, ungappedAlignment, false);
            stats->truncated = 1;
        }
    }
    if(queryResult.second > 1){
        if (identityId != UINT_MAX){
            SORT_SERIAL(resList + 1, resList + queryResult.second, hit_t::compareHitsByScoreAndId);
        } else{
            SORT_SERIAL(resList, resList + queryResult.second, hit_t::compareHitsByScoreAndId);
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
    while (seq->hasNextKmer()) {
        const unsigned char *kmer = seq->nextKmer();
        const unsigned char *pos = seq->getAAPosInSpacedPattern();
        const unsigned short current_i = seq->getCurrentPosition();

        float biasCorrection = 0;
        for (int i = 0; i < kmerSize; i++){
            biasCorrection += compositionBias[current_i + static_cast<short>(pos[i])];
        }
        if (seq->kmerContainsX()) {
            indexTo = current_i;
            indexPointer[current_i] = sequenceHits;
            continue;
        }
        // round bias to next higher or lower value
        short bias = static_cast<short>((biasCorrection < 0.0) ? biasCorrection - 0.5: biasCorrection + 0.5);
        short kmerMatchScore = std::max(kmerThr - bias, 0);

        // adjust kmer threshold based on composition bias
        kmerGenerator->setThreshold(kmerMatchScore);

        const size_t *index;
        size_t exactKmer;
        size_t kmerElementSize;
        if (takeOnlyBestKmer) {
            kmerElementSize = 1;
            exactKmer = idx.int2index(kmer);
            index = &exactKmer;
        } else {
            std::pair<size_t*, size_t> kmerList = kmerGenerator->generateKmerList(kmer);
            kmerElementSize = kmerList.second;
            index = kmerList.first;
        }
        //std::cout << kmer << std::endl;
        indexPointer[current_i] = sequenceHits;
        // match the index table

        //idx.printKmer(kmerList.index[0], kmerSize, m->num2aa);
        //std::cout << "\t" << kmerMatchScore << std::endl;
        kmerListLen += kmerElementSize;

        for (unsigned int kmerPos = 0; kmerPos < kmerElementSize; kmerPos++) {
            const IndexEntryLocal *entries = indexTable->getDBSeqList(index[kmerPos], &seqListSize);
            // DEBUG
            //std::cout << seq->getDbKey() << std::endl;
            //idx.printKmer(index[kmerPos], kmerSize, kmerSubMat->num2aa);
            //std::cout << "\t" << current_i << "\t"<< index[kmerPos] << std::endl;
            //for (size_t i = 0; i < seqListSize; i++) {
            //    char diag = entries[i].position_j - current_i;
            //    std::cout << "(" << entries[i].seqId << " " << (int) diag << ")\t";
            //}
            //std::cout << std::endl;

            // detected overflow while matching
            if ((sequenceHits + seqListSize) >= lastSequenceHit) {
                stats->diagonalOverflow = true;
                // last pointer
                indexPointer[current_i + 1] = sequenceHits;
                //std::cout << "Overflow in i=" << indexStart << std::endl;
                const size_t hitCount = findDuplicates(indexPointer,
                                                       foundDiagonals + overflowHitCount,
                                                       foundDiagonalsSize - overflowHitCount,
                                                       indexStart, current_i, (diagonalScoring == false));

                if (overflowHitCount != 0) {
                    // merge lists, hitCount is max. dbSize so there can be no overflow in mergeElements
                    overflowHitCount = mergeElements(foundDiagonals, hitCount + overflowHitCount);
                } else {
                    overflowHitCount = hitCount;
                }
                // reset pointer position
                sequenceHits = databaseHits;
                indexPointer[current_i] = sequenceHits;
                indexStart = current_i;
                overflowNumMatches += numMatches;
                numMatches = 0;
                // TODO might delete this?
                if ((sequenceHits + seqListSize) >= lastSequenceHit){
                    goto outer;
                }
            }
            memcpy(sequenceHits, entries, sizeof(IndexEntryLocal) * seqListSize);
            sequenceHits += seqListSize;
            numMatches += seqListSize;
        }
        indexTo = current_i;
    }
    outer:
    indexPointer[indexTo + 1] = databaseHits + numMatches;
    // fill the output
    size_t hitCount = findDuplicates(indexPointer, foundDiagonals + overflowHitCount,
                                     foundDiagonalsSize - overflowHitCount, indexStart, indexTo, (diagonalScoring == false));
    if (overflowHitCount != 0) {
        // overflow occurred
        hitCount = mergeElements(foundDiagonals, overflowHitCount + hitCount);
    }
    stats->doubleMatches = 0;
    if (diagonalScoring == false) {
        // remove double entries
        updateScoreBins(foundDiagonals, hitCount);
        stats->doubleMatches = getDoubleDiagonalMatches();
    }
    stats->kmersPerPos = ((double)kmerListLen/(double)seq->L);
    stats->querySeqLen = seq->L;
    stats->dbMatches   = overflowNumMatches + numMatches;

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

template <int TYPE>
std::pair<hit_t*, size_t> QueryMatcher::getResult(CounterResult * results,
                                                  size_t resultSize,
                                                  const unsigned int id,
                                                  const unsigned short thr,
                                                  UngappedAlignment *align,
                                                  const int rescaleScore) {
    size_t currentHits = 0;
    if (id != UINT_MAX) {
        hit_t *result = (resList + 0);
        unsigned short rawScore;
        if (TYPE == KMER_SCORE) {
            rawScore = UCHAR_MAX;
        } else {
            rawScore = USHRT_MAX;
        }
        result->seqId = id;
        result->prefScore = rawScore;
        result->diagonal = 0;
        //result->pScore = (((float)rawScore) - mu)/ sqrtMu;
        currentHits++;
    }

    for (size_t i = 0; i < resultSize && currentHits < maxHitsPerQuery; i++) {
        const unsigned int seqIdCurr = results[i].id;
        const unsigned int scoreCurr = results[i].count;
        const unsigned int diagCurr = results[i].diagonal;

        bool aboveThreshold = scoreCurr >= thr;
        bool isNotQueryId = (id != seqIdCurr);
        // write result to list
        //std::cout << i << "\t" << results[i].id << "\t" << (int)results[i].count << "\t" << results[i].diagonal << std::endl;
        if (aboveThreshold && isNotQueryId) {
            hit_t *result = (resList + currentHits);
            result->seqId = seqIdCurr;
            result->prefScore = scoreCurr;
            result->diagonal = diagCurr;
            //printf("%d\t%d\t%f\t%f\t%f\t%f\t%f\n", result->seqId, scoreCurr, seqLens[seqIdCurr], mu, logMatchProb, logScoreFactorial[scoreCurr]);
            if (TYPE == KMER_SCORE) {
                // make sure score is not great > 255
                result->prefScore = scoreCurr;
            } else {
                //need to get the real score
                if (rescaleScore != 0) {
                    unsigned int newScore = (UCHAR_MAX - align->getQueryBias());
                    newScore += (scoreCurr * rescaleScore / 255);
                    result->prefScore = newScore;
                } else if (static_cast<int>(scoreCurr) >= (UCHAR_MAX - align->getQueryBias())) {
                    unsigned int newScore = align->scoreSingelSequenceByCounterResult(results[i]);
                    result->prefScore = newScore;
                }
            }
            currentHits++;
        }
    }

    return std::make_pair(resList, currentHits);
}

void QueryMatcher::initDiagonalMatcher(size_t dbsize, unsigned int maxDbMatches) {
    uint64_t l2CacheSize = Util::getL2CacheSize();
#define INIT(x) cachedOperation##x = new CacheFriendlyOperations<x>(dbsize, maxDbMatches/x); \
                activeCounter = x;
    if(dbsize/2 < l2CacheSize){
        INIT(2)
    }else if(dbsize/4 < l2CacheSize){
        INIT(4)
    }else if(dbsize/8 < l2CacheSize){
        INIT(8)
    }else if(dbsize/16 < l2CacheSize){
        INIT(16)
    }else if(dbsize/32 < l2CacheSize){
        INIT(32)
    }else if(dbsize/64 < l2CacheSize){
        INIT(64)
    }else if(dbsize/128 < l2CacheSize){
        INIT(128)
    }else if(dbsize/256 < l2CacheSize){
        INIT(256)
    }else if(dbsize/512 < l2CacheSize){
        INIT(512)
    }else if(dbsize/1024 < l2CacheSize){
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

size_t QueryMatcher::findDuplicates(IndexEntryLocal **hitsByIndex,
                                   CounterResult *output, size_t outputSize,
                                   unsigned short indexFrom, unsigned short indexTo,
                                   bool computeTotalScore) {
    size_t localResultSize = 0;
#define COUNT_CASE(x) case x: localResultSize += cachedOperation##x->findDuplicates(hitsByIndex, output, outputSize, indexFrom, indexTo, computeTotalScore); break;
    switch (activeCounter){
        FOR_EACH(COUNT_CASE,2,4,8,16,32,64,128,256,512,1024,2048)
    }
#undef COUNT_CASE
    return localResultSize;
}

size_t QueryMatcher::mergeElements(CounterResult *foundDiagonals, size_t hitCounter) {
    size_t overflowHitCount = 0;
#define MERGE_CASE(x) \
    case x: overflowHitCount = diagonalScoring ? \
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

size_t QueryMatcher::radixSortByScoreSize(const unsigned int * scoreSizes,
                                        CounterResult *writePos,
                                        const unsigned int scoreThreshold,
                                        const CounterResult *results,
                                        const size_t resultSize) {
    CounterResult * ptr[SCORE_RANGE];
    ptr[0] = writePos+resultSize;
    CounterResult * ptr_prev=ptr[0];
    for(unsigned int i = 0; i < SCORE_RANGE; i++){
        ptr[i] = ptr_prev - scoreSizes[i];
        ptr_prev = ptr[i];
    }
    size_t aboveThresholdCnt = 0;
    for (size_t i = 0; i < resultSize; i++) {
        const unsigned int scoreCurr = results[i].count;
        if(scoreCurr >= scoreThreshold) {
            aboveThresholdCnt++;
            CounterResult*res = ptr[scoreCurr];
            res->id = results[i].id;
            res->count = results[i].count;
            res->diagonal = results[i].diagonal;
            ptr[scoreCurr]++;
        }
    }
    return aboveThresholdCnt;
}

std::pair<size_t, unsigned int> QueryMatcher::rescoreHits(Sequence * querySeq, unsigned int * scoreSizes, CounterResult *results,
        size_t resultSize, UngappedAlignment *align, int lowerBoundScore) {
    size_t elements = 0;
    const unsigned char * query = querySeq->numSequence;
    int maxSelfScore = align->scoreSingleSequence(std::make_pair(query, querySeq->L), 0,0);

    maxSelfScore = std::min(maxSelfScore, USHRT_MAX);
    maxSelfScore = (maxSelfScore-lowerBoundScore);
    maxSelfScore = std::max(1, maxSelfScore);
    float fltMaxSelfScore = static_cast<float>(maxSelfScore);
    for (size_t i = 0; i < resultSize && results[i].count >= lowerBoundScore; i++) {
        unsigned int newScore = align->scoreSingelSequenceByCounterResult(results[i]);
        newScore -= lowerBoundScore;
        float score = static_cast<float>(std::min(newScore, static_cast<unsigned int>(USHRT_MAX)));
        results[i].count = static_cast<unsigned char>((score/fltMaxSelfScore) * static_cast<float>(UCHAR_MAX) + 0.5);
        scoreSizes[results[i].count] += 1;
        elements++;
    }
    return std::make_pair(elements, maxSelfScore);
}

template std::pair<hit_t *, size_t>  QueryMatcher::getResult<0>(CounterResult * results, size_t resultSize,
                                                    const unsigned int id, const unsigned short thr,
                                                    UngappedAlignment * align, const int rescaleScore);
template std::pair<hit_t *, size_t>  QueryMatcher::getResult<1>(CounterResult * results, size_t resultSize,
                                                                const unsigned int id, const unsigned short thr,
                                                                UngappedAlignment * align, const int rescaleScore);

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
