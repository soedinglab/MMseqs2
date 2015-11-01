#include "QueryTemplateMatcherGlobal.h"
#include "QueryScoreGlobal.h"

QueryTemplateMatcherGlobal::QueryTemplateMatcherGlobal(BaseMatrix *m,
        IndexTable *indexTable,
        unsigned int *seqLens,
        short kmerThr,
        double kmerMatchProb,
        int kmerSize,
        size_t effectiveKmerSize,
        size_t dbSize,
        bool aaBiasCorrection,
        int maxSeqLen,
        float zscoreThr) : QueryTemplateMatcher(m, indexTable, seqLens, kmerThr, kmerMatchProb, kmerSize, dbSize, aaBiasCorrection, maxSeqLen) {
    this->queryScore = new QueryScoreGlobal(dbSize, seqLens, effectiveKmerSize, kmerThr, kmerMatchProb, zscoreThr);
    this->deltaS = new float[maxSeqLen];
    memset(this->deltaS, 0, maxSeqLen * sizeof(float));

}


QueryTemplateMatcherGlobal::~QueryTemplateMatcherGlobal(){
    delete queryScore;
    delete [] deltaS;
}

// calculate local amino acid bias correction score for each position in the sequence
void QueryTemplateMatcherGlobal::calcLocalAaBiasCorrection(Sequence* seq){
    const int windowSize = 40;
    if (seq->L < windowSize + 1)
        memset(this->deltaS, 0, seq->L * sizeof(float));
    else{
        // calculate local amino acid bias
        for (int i = 0; i < seq->L; i++){
            const int minPos = std::max(0, (i - windowSize/2));
            const int maxPos = std::min(seq->L, (i + windowSize/2));
            const int _2d = maxPos - minPos;
            // negative score for the amino acids in the neighborhood of i
            int sumSubScores = 0;
            short * subMat = m->subMatrix[seq->int_sequence[i]];
            for (int j = minPos; j < maxPos; j++){
                sumSubScores += (j != i) ? subMat[seq->int_sequence[j]] : 0;
            }
            float deltaS_i = (float) sumSubScores;
            deltaS_i /= -1.0 * _2d;
            // positive score for the background score distribution for i
            for (int a = 0; a < m->alphabetSize; a++){
                deltaS_i += m->pBack[a] * subMat[a];
            }

            deltaS[i] = deltaS_i;
        }
    }
}



std::pair<hit_t *, size_t> QueryTemplateMatcherGlobal::matchQuery (Sequence * seq, unsigned int identityId){
    queryScore->reset();
    seq->resetCurrPos();
    
    if (this->aaBiasCorrection)
        this->calcLocalAaBiasCorrection(seq);

    match(seq);

    queryScore->setPrefilteringThresholds();

    return queryScore->getResult(seq->L,identityId);
}

void QueryTemplateMatcherGlobal::match(Sequence* seq){

    unsigned int* seqList;
    size_t indexTabListSize = 0;
    // go through the query sequence
    int kmerListLen = 0;

    float biasCorrection = 0;
    for (int i = 0; i < kmerSize && i < seq->L; i++)
        biasCorrection += deltaS[i];

    int pos = 0;
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
            
            
            seqList = indexTable->getDBSeqList<unsigned int>(kmerList.index[i], &indexTabListSize);
            
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
            queryScore->addScores(seqList, indexTabListSize, (kmerMatchScore/4));
        }
        biasCorrection -= deltaS[pos];
        biasCorrection += deltaS[pos + kmerSize];
        pos++;
    }
    //Debug(Debug::WARNING) << "QUERY: " << seq->getDbKey();
    //Debug(Debug::WARNING) << " score = " << overall_score;
    //Debug(Debug::WARNING) << " matched at " << match_pos << " positions. ";
    //Debug(Debug::WARNING) << match_num << " times.\n";
    // write statistics
    stats->doubleMatches = 0;
    stats->kmersPerPos   = ((double)kmerListLen/(double)seq->L);
    stats->querySeqLen   = seq->L;
    stats->dbMatches     = queryScore->getNumMatches();

//    delete indexer;

}
