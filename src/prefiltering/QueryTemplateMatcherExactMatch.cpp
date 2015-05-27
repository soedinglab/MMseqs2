//
// Created by mad on 5/26/15.
//

#include <CountInt32Array.h>
#include "QueryTemplateMatcherExactMatch.h"
#include "QueryTemplateMatcher.h"

QueryTemplateMatcherExactMatch::QueryTemplateMatcherExactMatch(BaseMatrix *m, IndexTable *indexTable,
                                                               unsigned int *seqLens, short kmerThr,
                                                               double kmerMatchProb, int kmerSize, size_t dbSize,
                                                               unsigned int maxSeqLen, size_t maxHitsPerQuery) : QueryTemplateMatcher(m, indexTable, seqLens, kmerThr, kmerMatchProb,
                                                                                                                                      kmerSize, dbSize, false, maxSeqLen) {
    this->maxHitsPerQuery = maxHitsPerQuery;
    this->idxer = new Indexer(m->alphabetSize, kmerSize);
    this->resList = (hit_t *) mem_align(ALIGN_INT, QueryScore::MAX_RES_LIST_LEN * sizeof(hit_t) );
    this->foundSequences = new unsigned int[100000];
    //this->counter = new CountInt32Array(dbSize, 2048);
}


QueryTemplateMatcherExactMatch::~QueryTemplateMatcherExactMatch(){
    delete idxer;
    delete [] foundSequences;
    free(resList);
    //delete counter;
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
        // generate k-mer list
        unsigned int kmerCode = idxer->getNextKmerIndex(kmer, kmerSize);
        entries = indexTable->getDBSeqList<unsigned int>(kmerCode, &indexTabListSize);

        // add the scores for the k-mer to the overall score for this query sequence
        // for the overall score, bit/2 is a sufficient sensitivity and we can use the capacity of unsigned short max score in QueryScore better
        // std::cout << "i: " << seq->getCurrentPosition() << std::endl;
        for (unsigned int seqIdx = 0; LIKELY(seqIdx < indexTabListSize); seqIdx++){
            foundSequences[pos] = entries[seqIdx];
            pos++;
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
    const size_t resultSize = stats->dbMatches;
    std::sort(foundSequences, foundSequences + resultSize);

    if (id != UINT_MAX){
        hit_t * result = (resList + 0);
        const unsigned short rawScore  = l;
        result->seqId = id;
        result->prefScore = rawScore;
        result->zScore = rawScore;
        elementCounter++;
    }

    unsigned int seqIdPrev = foundSequences[0];
    unsigned short scoreMax = 1;
    for (unsigned int i = 1; i < resultSize; i++) {
        const unsigned int seqIdCurr = foundSequences[i];
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
