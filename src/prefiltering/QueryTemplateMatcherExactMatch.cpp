//
// Created by mad on 5/26/15.
//

#include "QueryTemplateMatcherExactMatch.h"

QueryTemplateMatcherExactMatch::QueryTemplateMatcherExactMatch(BaseMatrix *m, IndexTable *indexTable,
                                                                   unsigned int *seqLens, short kmerThr,
                                                                   double kmerMatchProb, int kmerSize, size_t dbSize,
                                                                   unsigned int maxSeqLen)
        : QueryTemplateMatcher(m, indexTable, seqLens, kmerThr, kmerMatchProb, kmerSize, dbSize, false, maxSeqLen) {
    //this->idxer = new Indexer(m->alphabetSize, kmerSize);
    this->resList = (hit_t *) mem_align(ALIGN_INT, QueryScore::MAX_RES_LIST_LEN * sizeof(hit_t) );
    this->binData = new unsigned int[BIN_COUNT * BIN_SIZE];
    //this->BIN_SIZE = QueryScore::MAX_RES_LIST_LEN / BIN_COUNT;
//    this->foundSequences = new unsigned int[MAX_DB_MATCHES];
//    memset(foundSequences, 0, sizeof(unsigned int) * MAX_DB_MATCHES);
    // needed as buffer
    this->counterOutput = new unsigned int[MAX_DB_MATCHES];
    memset(counterOutput, 0, sizeof(unsigned int) * MAX_DB_MATCHES);
    this->counter = new CountInt32Array(dbSize, BIN_SIZE ); //TODO what is the right choice?!?
}

QueryTemplateMatcherExactMatch::~QueryTemplateMatcherExactMatch(){
    delete counter;
    delete [] counterOutput;
    delete [] binData;
    free(resList);
}

size_t QueryTemplateMatcherExactMatch::evaluateBins(unsigned int *output) {
    size_t localResultSize = 0;
//    const size_t N =  (diagonalBins[0] - binData);
//    memcpy(output, binData, sizeof(unsigned int) * N);
//    localResultSize += N;
    const unsigned int * lastPointer = output + MAX_DB_MATCHES;
    for(size_t bin = 0; bin < BIN_COUNT; bin++){
        unsigned int *binStartPos = (binData + bin * BIN_SIZE);
        const size_t N =  (diagonalBins[bin] - binStartPos);
//        std::sort(binStartPos, binStartPos + N);
//        unsigned int prefId = binStartPos[0];
//        for(size_t a = 1; a < N; a++){
//            unsigned int currId = binStartPos[a];
//            if(prefId == currId){
//                output[localResultSize] = binStartPos[a];
//                localResultSize++;
//            }
//            prefId = currId;
//        }
        if(output + localResultSize >= lastPointer)
            break;

        localResultSize += counter->countElements(binStartPos, N, output + localResultSize, lastPointer);

    }
    return localResultSize;
}

void QueryTemplateMatcherExactMatch::reallocBinMemory(const unsigned int binCount, const size_t binSize) {
    delete [] binData;
    binData = new unsigned int[binCount * binSize];
}

bool QueryTemplateMatcherExactMatch::checkForOverflowAndResizeArray() {
//    const unsigned int * bin_ref_pointer = binData;
//    unsigned int * lastPosition = (binData + BIN_COUNT * BIN_SIZE) - 1;
//    for (size_t bin = 0; bin < BIN_COUNT; bin++) {
//        const unsigned int *binStartPos = (bin_ref_pointer + bin * BIN_SIZE);
//        const size_t n = (diagonalBins[bin] - binStartPos);
//        // if one bin has more elements than BIN_SIZE
//        // or the current bin pointer is at the end of the binDataFrame
//        // reallocate new memory
//        if( n > BIN_SIZE || (diagonalBins[bin] - lastPosition) == 0) {
//            // overflow detected
//            // find nearest upper power of 2^(x)
//            std::cout << "Diagonal Found overlow" << std::endl;
//            this->BIN_SIZE = pow(2, ceil(log(BIN_SIZE + 1)/log(2)));
//            reallocBinMemory(BIN_COUNT, this->BIN_SIZE);
//            return true;
//        }
//    }
    return false;
}

void QueryTemplateMatcherExactMatch::setupBinPointer() {
    // Example binCount = 3
    // bin start             |-----------------------|-----------------------| bin end
    //    segments[bin_step][0]
    //                            segments[bin_step][1]
    //                                                    segments[bin_step][2]
    size_t curr_pos = 0;
    for(size_t bin = 0; bin < BIN_COUNT; bin++){
        diagonalBins[bin] = binData + curr_pos;
        curr_pos += BIN_SIZE;
    }
}

std::pair<hit_t *, size_t> QueryTemplateMatcherExactMatch::matchQuery (Sequence * seq, unsigned int identityId){
    seq->resetCurrPos();
    setupBinPointer();

    match(seq);

    return getResult(seq->L, identityId, 0);
}


void QueryTemplateMatcherExactMatch::match(Sequence* seq){
    // go through the query sequence
    size_t kmerListLen = 0;
    size_t numMatches = 0;
    //size_t pos = 0;
    size_t seqListSize = 0;

    while(seq->hasNextKmer()){
        const int* kmer = seq->nextKmer();

        const ScoreMatrix kmerList = kmerGenerator->generateKmerList(kmer);
        kmerListLen += kmerList.elementSize;
        const unsigned short i =  seq->getCurrentPosition();
        // match the index table
//        int pos_matched = 0;
        for (unsigned int kmerPos = 0; kmerPos < kmerList.elementSize; kmerPos++) {
            // generate k-mer list
            const IndexEntryLocal * entries = indexTable->getDBSeqList<IndexEntryLocal>(kmerList.index[kmerPos], &seqListSize);

            const unsigned int * lastPosition = (binData + BIN_COUNT * BIN_SIZE) - 1;
            for (unsigned int seqIdx = 0; LIKELY(seqIdx < seqListSize); seqIdx++){
                IndexEntryLocal entry = entries[seqIdx];
                const unsigned char j = entry.position_j;
                const unsigned int seqId = entry.seqId;
                const unsigned char diagonal = (i - j);
                const unsigned char diagonal2 = diagonal/static_cast<unsigned char>(64);
                *diagonalBins[diagonal2] = seqId;
                diagonalBins[diagonal2] += ((diagonalBins[diagonal2] - lastPosition) != 0) ? 1 : 0; //TODO what to do with it?
                //if((diagonalBins[diagonal] - lastPosition) == 0)
               // std::cout << (int)i << " " << (int) j << " " << (int) diagonal << std::endl;
            }
            numMatches += seqListSize;
        }
    }
    //filles the output
    size_t doubleMatches = evaluateBins(counterOutput);

    // needed to call here to get the LocalResultSize
    //Debug(Debug::WARNING) << "QUERY: " << seq->getDbKey();
    //Debug(Debug::WARNING) << " score = " << overall_score;
    //Debug(Debug::WARNING) << " matched at " << match_pos << " positions. ";
    //Debug(Debug::WARNING) << match_num << " times.\n";
    // write statistics
   // std::cout << doubleMatches << std::endl;
    stats->doubleMatches = doubleMatches;
    stats->kmersPerPos   = ((double)kmerListLen/(double)seq->L);
    stats->querySeqLen   = seq->L;
    stats->dbMatches     = numMatches;
    //std::cout << seq->getDbKey() << " " <<  seq->getId() << " " << stats->doubleMatches << " "
    //          << stats->kmersPerPos << " " << stats->dbMatches << std::endl;

    //    delete indexer;
}

std::pair<hit_t *, size_t>  QueryTemplateMatcherExactMatch::getResult(const int l,
                                                                      const unsigned int id,
                                                                      const unsigned short thr) {
    size_t elementCounter = 0;
    size_t resultSize = stats->doubleMatches;
    //std::sort(foundSequences, foundSequences + resultSize);
    std::sort(counterOutput, counterOutput + resultSize);
    if (id != UINT_MAX){
        hit_t * result = (resList + 0);
        const unsigned short rawScore  = l;
        result->seqId = id;
        result->prefScore = rawScore;
        result->zScore = rawScore;
        elementCounter++;
    }

    unsigned int seqIdPrev = counterOutput[0];
    unsigned short scoreMax = 1;

    if(resultSize == 1){
        hit_t *result = (resList + elementCounter);
        result->seqId = seqIdPrev;
        result->zScore = scoreMax;
        result->prefScore = scoreMax;
        elementCounter++;
    }
    for (unsigned int i = 1; i < resultSize; i++) {
        const unsigned int seqIdCurr = counterOutput[i];
        // if new sequence occurs or end of data write the result back
        scoreMax += 1;
        if (seqIdCurr != seqIdPrev || i == (resultSize - 1)) {
            // write result to list
            if(scoreMax > thr && id != seqIdPrev){
                hit_t *result = (resList + elementCounter);
                result->seqId = seqIdPrev;
                result->zScore = scoreMax;
                result->prefScore = scoreMax;
                elementCounter++;
                if (elementCounter >= QueryScore::MAX_RES_LIST_LEN)
                    break;
            }
            seqIdPrev = seqIdCurr;
            // reset values
            scoreMax = 0;  // current best score
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
