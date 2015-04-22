#include "QueryTemplateMatcherLocal.h"
#include "QueryScoreLocal.h"
#include "QueryTemplateMatcher.h"

QueryTemplateMatcherLocal::QueryTemplateMatcherLocal(BaseMatrix *m,
                                                     IndexTable *indexTable,
                                                     unsigned int *seqLens,
                                                     short kmerThr,
                                                     double kmerMatchProb,
                                                     int kmerSize,
                                                     size_t effectiveKmerSize,
                                                     size_t dbSize,
                                                     bool fastMode,
                                                     unsigned int maxSeqLen,
                                                     size_t maxHitsPerQuery) : QueryTemplateMatcher(m, indexTable, seqLens, kmerThr, kmerMatchProb,
                                                       kmerSize, dbSize, fastMode, maxSeqLen) {
    this->queryScore = new QueryScoreLocal(dbSize, seqLens, effectiveKmerSize, kmerThr, kmerMatchProb);
    this->maxHitsPerQuery = maxHitsPerQuery;
    this->fastMode = fastMode;
}


QueryTemplateMatcherLocal::~QueryTemplateMatcherLocal(){
    delete queryScore;
}

std::pair<hit_t *, size_t> QueryTemplateMatcherLocal::matchQuery (Sequence * seq, unsigned int identityId){
    queryScore->reset();
    seq->resetCurrPos();
    
    match(seq);

    unsigned int scoreThreshold;
    if(fastMode == true){
        scoreThreshold = 1;
    }else{
        scoreThreshold = queryScore->computeScoreThreshold(this->maxHitsPerQuery);
    }

    return queryScore->getResult(seq->L, scoreThreshold, identityId);
}


void QueryTemplateMatcherLocal::match(Sequence* seq){
    IndexEntryLocal* entries;
    size_t indexTabListSize = 0;
    // go through the query sequence
    size_t kmerListLen = 0;

    unsigned int pos = 0;
    while(seq->hasNextKmer()){
        const int* kmer = seq->nextKmer();
        // generate k-mer list
        ScoreMatrix kmerList = kmerGenerator->generateKmerList(kmer);
        kmerListLen += kmerList.elementSize;
        // match the index table
        //        int pos_matched = 0;
        for (unsigned int i = 0; i < kmerList.elementSize; i++){

            entries = indexTable->getDBSeqList<IndexEntryLocal>(kmerList.index[i], &indexTabListSize);

            /*            if (seq->getId() == 1 && pos == 2 ){
             std::cout << "\t\t";
             indexer->printKmer(retList[i].second, kmerSize, m->int2aa);
             std::cout << " " << retList[i].first << "\n";
             }*/
            //            for (int j = 0; j < indexTabListSize; j++){
            //                    std::cout << "Similar k-mer list pos: " << i << ", score: " << kmerList.index[i] << ", kmer idx: " << kmerList.score[i] << "\n";
            //                    pos_matched = 1;
            //                    std::cout << pos << " ";
            //                    indexer->printKmer(kmer, kmerSize, m->int2aa);
            //                    std::cout << "->";
            //                    indexer->printKmer(kmerList.index[i], kmerSize, m->int2aa);
            //                    std::cout << "\n";
            //                    std::cout << "\t" << kmerMatchScore << "\n";
            //                    overall_score+=kmerMatchScore;
            //                    match_num++;
            //            }

            // add the scores for the k-mer to the overall score for this query sequence
            // for the overall score, bit/2 is a sufficient sensitivity and we can use the capacity of unsigned short max score in QueryScore better
            // std::cout << "i: " << seq->getCurrentPosition() << std::endl;
            queryScore->addScoresLocal(entries, seq->getCurrentPosition(), indexTabListSize);
        }
        pos++;
    }
    if(fastMode == false) {
        queryScore->updateScoreBins();
    }
    // needed to call here to get the LocalResultSize
    //Debug(Debug::WARNING) << "QUERY: " << seq->getDbKey();
    //Debug(Debug::WARNING) << " score = " << overall_score;
    //Debug(Debug::WARNING) << " matched at " << match_pos << " positions. ";
    //Debug(Debug::WARNING) << match_num << " times.\n";
    // write statistics
    stats->doubleMatches = queryScore->getLocalResultSize();
    stats->kmersPerPos   = ((double)kmerListLen/(double)seq->L);
    stats->querySeqLen   = seq->L;
    stats->dbMatches     = queryScore->getNumMatches();
    //std::cout << seq->getDbKey() << " " <<  seq->getId() << " " << stats->doubleMatches << " "
    //          << stats->kmersPerPos << " " << stats->dbMatches << std::endl;

    //    delete indexer;
}
