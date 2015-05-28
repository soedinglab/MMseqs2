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
    this->foundSequences = new unsigned int[MAX_DB_MATCHES];
    memset(foundSequences, 0, sizeof(unsigned int) * MAX_DB_MATCHES);
    // needed as buffer
    this->counterOutpot = new unsigned int[MAX_DB_MATCHES];
    memset(counterOutpot, 0, sizeof(unsigned int) * MAX_DB_MATCHES);
    this->counter = new CountInt32Array(dbSize, MAX_DB_MATCHES / 512);
}


QueryTemplateMatcherExactMatch::~QueryTemplateMatcherExactMatch(){
    delete counter;
    delete [] foundSequences;
    free(resList);
}

std::pair<hit_t *, size_t> QueryTemplateMatcherExactMatch::matchQuery (Sequence * seq, unsigned int identityId){
    seq->resetCurrPos();

    match(seq);

    return getResult(seq->L, identityId, 1);
}


void QueryTemplateMatcherExactMatch::match(Sequence* seq){
    unsigned int * entries;
    size_t indexTabListSize = 0;
    // go through the query sequence
    size_t kmerListLen = 0;

    size_t pos = 0;
    while(seq->hasNextKmer()){
        const int* kmer = seq->nextKmer();

        ScoreMatrix kmerList = kmerGenerator->generateKmerList(kmer);
        kmerListLen += kmerList.elementSize;
        // match the index table
//        int pos_matched = 0;
        for (unsigned int i = 0; i < kmerList.elementSize; i++) {
            // generate k-mer list
            entries = indexTable->getDBSeqList<unsigned int>(kmerList.index[i], &indexTabListSize);

            if (pos + indexTabListSize >= MAX_DB_MATCHES)
                break;
            for (unsigned int seqIdx = 0; LIKELY(seqIdx < indexTabListSize); seqIdx++) {
                foundSequences[pos] = entries[seqIdx];
                pos++;
            }
        }
    }


    // needed to call here to get the LocalResultSize
    //Debug(Debug::WARNING) << "QUERY: " << seq->getDbKey();
    //Debug(Debug::WARNING) << " score = " << overall_score;
    //Debug(Debug::WARNING) << " matched at " << match_pos << " positions. ";
    //Debug(Debug::WARNING) << match_num << " times.\n";
    // write statistics
    stats->doubleMatches = 0;
    stats->kmersPerPos   = ((double)kmerListLen/(double)seq->L);
    stats->querySeqLen   = seq->L;
    stats->dbMatches     = pos;
    //std::cout << seq->getDbKey() << " " <<  seq->getId() << " " << stats->doubleMatches << " "
    //          << stats->kmersPerPos << " " << stats->dbMatches << std::endl;

    //    delete indexer;
}

std::pair<hit_t *, size_t>  QueryTemplateMatcherExactMatch::getResult(const int l, const unsigned int id,
                                                                      const unsigned short thr) {
    size_t elementCounter = 0;
    //std::sort(foundSequences, foundSequences + resultSize);
    const size_t resultSize = counter->countElements(foundSequences, stats->dbMatches, counterOutpot);
    std::sort(counterOutpot, counterOutpot + resultSize);
    if (id != UINT_MAX){
        hit_t * result = (resList + 0);
        const unsigned short rawScore  = l;
        result->seqId = id;
        result->prefScore = rawScore;
        result->zScore = rawScore;
        elementCounter++;
    }

    unsigned int seqIdPrev = counterOutpot[0];
    unsigned short scoreMax = 1;
    for (unsigned int i = 1; i < resultSize; i++) {
        const unsigned int seqIdCurr = counterOutpot[i];
        // if new sequence occurs or end of data write the result back
        scoreMax += 1;
        if (seqIdCurr != seqIdPrev || i == (resultSize - 1)) {
            // write result to list
            if(scoreMax > thr){
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
    std::pair<hit_t *, size_t>  pair = std::make_pair(this->resList, elementCounter);
    return pair;
}
