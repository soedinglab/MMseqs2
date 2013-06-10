#include "QueryTemplateMatcher.h"
#include "QueryScoreGlobal.h"
QueryTemplateMatcher::QueryTemplateMatcher ( ExtendedSubstitutionMatrix* _2merSubMatrix,
                                             ExtendedSubstitutionMatrix* _3merSubMatrix,
                                             IndexTable * indexTable,
                                             unsigned short * seqLens,
                                             short kmerThr,
                                             float prefThr,
                                             int kmerSize, 
                                             int dbSize,
                                             int alphabetSize){
    this->indexTable = indexTable;
    this->kmerSize = kmerSize;
    this->alphabetSize = alphabetSize;
    this->kmerGenerator = new KmerGenerator(kmerSize, alphabetSize, kmerThr, _3merSubMatrix, _2merSubMatrix);
    this->queryScore    = new QueryScoreGlobal(dbSize, seqLens, prefThr, kmerSize);
}

QueryTemplateMatcher::~QueryTemplateMatcher (){
    delete indexTable;
    delete kmerGenerator;
    delete queryScore;
}

std::list<hit_t>* QueryTemplateMatcher::matchQuery (Sequence * seq){
    queryScore->reset();
    seq->resetCurrPos();

/* DEBUGGING
 *    int* testKmer = new int[kmerSize];
 *    unsigned int kmerIdx;
 *    Indexer* idxer = new Indexer(alphabetSize, kmerSize);
 */
 
    int* seqList;
    int listSize = 0;
    // go through the query sequence
    int kmerListLen = 0;
    int numMatches = 0;
    while(seq->hasNextKmer(kmerSize)){
        const int* kmer = seq->nextKmer(kmerSize);
        /* DEBUGGING
         * idxer->printKmer(kmer, kmerSize, seq->int2aa);
         * std::cout << "\n";
         */
        
        // generate k-mer list
        KmerGeneratorResult kmerList = kmerGenerator->generateKmerList(kmer);
        kmerListLen += kmerList.count;
        std::pair<short,unsigned int> * retList = kmerList.scoreKmerList;
        
        // match the index table
        for (unsigned int i = 0; i < kmerList.count; i++){
            std::pair<short,unsigned int> kmerMatch = retList[i];
            /* DEBUGGING
             * std::cout << "\t";
             * idxer->printKmer(testKmer, kmerMatch.second, kmerSize, seq->int2aa);
             */
             seqList = indexTable->getDBSeqList(kmerMatch.second, &listSize);
             /* DEBUGGING
              * std::cout << " (" << listSize << ")\n";
              */
            numMatches += listSize;

            // add the scores for the k-mer to the overall score for this query sequence
            queryScore->addScores(seqList, listSize, (unsigned short)kmerMatch.first);
        }
    } 
    seq->stats->kmersPerPos = ((float)kmerListLen/(float)seq->L);
    seq->stats->dbMatches = numMatches;

    /* DEBUGGING
     * delete idxer;
     */
    return queryScore->getResult(seq->L);
}

void QueryTemplateMatcher::printStats(){
    this->queryScore->printStats();
}
