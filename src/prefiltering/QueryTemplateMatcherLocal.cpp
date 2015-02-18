#include "QueryTemplateMatcherLocal.h"
#include "QueryScoreLocal.h"
#include "QueryTemplateMatcher.h"

QueryTemplateMatcherLocal::QueryTemplateMatcherLocal(BaseMatrix *m,
                                                     IndexTable *indexTable,
                                                     unsigned int *seqLens,
                                                     short kmerThr,
                                                     double kmerMatchProb,
                                                     int kmerSize,
                                                     int dbSize,
                                                     bool aaBiasCorrection,
                                                     int maxSeqLen,
                                                     float zscoreThr) : QueryTemplateMatcher(m, indexTable, seqLens, kmerThr, kmerMatchProb,
                                                                                                kmerSize, dbSize, aaBiasCorrection, maxSeqLen, zscoreThr) {
    this->queryScore = new QueryScoreLocal(dbSize, seqLens, kmerSize, kmerThr, kmerMatchProb, zscoreThr, 262144);
}


QueryTemplateMatcherLocal::~QueryTemplateMatcherLocal(){
    delete queryScore;
}

std::pair<hit_t *, size_t> QueryTemplateMatcherLocal::matchQuery (Sequence * seq, unsigned int identityId){
    queryScore->reset();
    seq->resetCurrPos();
    
    match(seq);


//    queryScore->setPrefilteringThresholds();
    
    return queryScore->getResult(seq->L,0);
}


void QueryTemplateMatcherLocal::match(Sequence* seq){
    IndexEntryLocal* entries;
    size_t indexTabListSize = 0;
    // go through the query sequence
    int kmerListLen = 0;
    
    float biasCorrection = 0;
    for (int i = 0; i < kmerSize && i < seq->L; i++)
        biasCorrection += deltaS[i];
    
    
    unsigned int pos = 0;
    short zero = 0;
    while(seq->hasNextKmer()){
        const int* kmer = seq->nextKmer();
        // generate k-mer list
        ScoreMatrix kmerList = kmerGenerator->generateKmerList(kmer);
        kmerListLen += kmerList.elementSize;
        // match the index table
        //        int pos_matched = 0;
        short biasCorrection_short = (short) biasCorrection;
        for (unsigned int i = 0; i < kmerList.elementSize; i++){
            // avoid unsigned short overflow
            short kmerMatchScore = kmerList.score[i] + biasCorrection_short;
            // avoid unsigned short overflow
            kmerMatchScore = std::max(kmerMatchScore, zero);
            
            
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
            queryScore->addScoresLocal(entries, seq->getCurrentPosition(), indexTabListSize, kmerMatchScore/4);
        }
        pos++;
    }
    // needed to call here to get the LocalResultSize
    //Debug(Debug::WARNING) << "QUERY: " << seq->getDbKey();
    //Debug(Debug::WARNING) << " score = " << overall_score;
    //Debug(Debug::WARNING) << " matched at " << match_pos << " positions. ";
    //Debug(Debug::WARNING) << match_num << " times.\n";
    // write statistics
    stats->doubleMatches = queryScore->getLocalResultSize();
    stats->kmersPerPos   = ((float)kmerListLen/(float)seq->L);
    stats->dbMatches     = queryScore->getNumMatches();
    //    delete indexer;
}
