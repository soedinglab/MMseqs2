#include "QueryTemplateMatcher.h"
QueryTemplateMatcher::QueryTemplateMatcher ( ExtendedSubstitutionMatrix* _2merSubMatrix,
                                             ExtendedSubstitutionMatrix* _3merSubMatrix,
                                             IndexTable * indexTable,
                                             short kmerThr,
                                             short prefThr,
                                             int kmerSize, 
                                             int dbSize,
                                             int alphabetSize){
    this->indexTable = indexTable;
    this->kmerSize = kmerSize;
    this->alphabetSize = alphabetSize;
    this->kmerGenerator = new KmerGenerator(kmerSize, alphabetSize, kmerThr, _3merSubMatrix, _2merSubMatrix);
    this->queryScore = new QueryScore(dbSize, prefThr);
}

QueryTemplateMatcher::~QueryTemplateMatcher (){
    delete indexTable;
    delete kmerGenerator;
    delete queryScore;
}

std::list<hit_t>* QueryTemplateMatcher::matchQuery (Sequence * seq){
    queryScore->reset();
    seq->resetCurrPos();

    // DEBUGGING
//    Indexer* idxer = new Indexer(alphabetSize, kmerSize);
 
    int* seqList;
    unsigned int kmerIdx;
    int listSize = 0;
    // go through the query sequence
    int* testKmer = new int[kmerSize];
    int kmerListLen = 0;
    int numMatches = 0;
    while(seq->hasNextKmer(kmerSize)){
        const int* kmer = seq->nextKmer(kmerSize);
        
        // generate k-mer list
        kmer_list kmerList = kmerGenerator->generateKmerList(kmer);
        kmerListLen += kmerList.count;
        std::vector<std::pair<short,size_t> > * retList = kmerList.score_kmer_list;
        
        // match the index table
        for (int i = 0; i < kmerList.count; i++){
            std::pair<short,size_t> kmerMatch = retList->at(i);
            seqList = indexTable->getDBSeqList(kmerMatch.second, &listSize);
            numMatches += listSize;

            // add the scores for the k-mer to the overall score for this query sequence
            queryScore->addScores(seqList, listSize, kmerMatch.first);
        }
    } 
    seq->stats->kmersPerPos = ((float)kmerListLen/(float)seq->L);
    seq->stats->dbMatches = numMatches;

//    delete idxer;
    return queryScore->getResult();
}
