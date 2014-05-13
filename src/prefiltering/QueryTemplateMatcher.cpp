#include "QueryTemplateMatcher.h"
#include "QueryScoreGlobal.h"
QueryTemplateMatcher::QueryTemplateMatcher ( BaseMatrix* m,
        ExtendedSubstitutionMatrix* _2merSubMatrix,
        ExtendedSubstitutionMatrix* _3merSubMatrix,
        IndexTable * indexTable,
        unsigned short * seqLens,
        short kmerThr,
        double kmerMatchProb,
        int kmerSize, 
        int dbSize,
        bool aaBiasCorrection,
        int maxSeqLen,
        float zscoreThr){
    this->m = m;
    this->indexTable = indexTable;
    this->kmerSize = kmerSize;
    this->kmerGenerator = new KmerGenerator(kmerSize, m->alphabetSize, kmerThr, _3merSubMatrix, _2merSubMatrix);
    this->queryScore    = new QueryScoreGlobal(dbSize, seqLens, kmerSize, kmerThr, kmerMatchProb, zscoreThr);
    this->aaBiasCorrection = aaBiasCorrection;

    this->deltaS = new float[maxSeqLen];
    memset(this->deltaS, 0, maxSeqLen * sizeof(float));
}

QueryTemplateMatcher::~QueryTemplateMatcher (){
    delete[] deltaS;
    delete kmerGenerator;
    delete queryScore;
}

void QueryTemplateMatcher::calcLocalAaBiasCorrection(Sequence* seq){
    int windowSize = 40;
    if (seq->L < windowSize + 1)
        memset(this->deltaS, 0, seq->L * sizeof(float));
    else{
        float deltaS_i;
        int minPos;
        int maxPos;
        int _2d;
        // calculate local amino acid bias 
        for (int i = 0; i < seq->L; i++){
            deltaS_i = 0.0;
            minPos = std::max(0, (i - windowSize/2));
            maxPos = std::min(seq->L, (i + windowSize/2));
            _2d = maxPos - minPos;

            // negative score for the amino acids in the neighborhood of i
            for (int j = minPos; j < maxPos; j++){
                if (j != i){
                    deltaS_i += m->subMatrix[seq->int_sequence[i]][seq->int_sequence[j]];
                }
            }
            deltaS_i /= -1.0 * _2d;
            // positive score for the background score distribution for i
            for (int a = 0; a < m->alphabetSize; a++)
                deltaS_i += m->pBack[a] * m->subMatrix[seq->int_sequence[i]][a];

            deltaS[i] = deltaS_i;
        }
    }
}

std::pair<hit_t *, size_t> QueryTemplateMatcher::matchQuery (Sequence * seq, unsigned int identityId){
    queryScore->reset();
    seq->resetCurrPos();

    match(seq);

    queryScore->setPrefilteringThresholds();

    return queryScore->getResult(seq->L, identityId);
}

void QueryTemplateMatcher::match(Sequence* seq){

    seq->resetCurrPos();
    
    if (this->aaBiasCorrection)
        calcLocalAaBiasCorrection(seq);

    int* seqList;
    int indexTabListSize = 0;
    // go through the query sequence
    int kmerListLen = 0;
    int numMatches = 0;

    float biasCorrection = 0;
    for (int i = 0; i < kmerSize && i < seq->L; i++)
        biasCorrection += deltaS[i];

    int pos = 0;
    short zero = 0;
    while(seq->hasNextKmer(kmerSize)){
        const int* kmer = seq->nextKmer(kmerSize);
        // generate k-mer list
        KmerGeneratorResult kmerList = kmerGenerator->generateKmerList(kmer);
        kmerListLen += kmerList.count;
        std::pair<short,unsigned int> * retList = kmerList.scoreKmerList;

        // match the index table
        for (unsigned int i = 0; i < kmerList.count; i++){
            std::pair<short,unsigned int> kmerMatch = retList[i];
            short kmerMatchScore = kmerMatch.first + (short) biasCorrection;
            // avoid unsigned short overflow
            kmerMatchScore = std::max(kmerMatchScore, zero);

            seqList = indexTable->getDBSeqList(kmerMatch.second, &indexTabListSize);
            numMatches += indexTabListSize;

            // add the scores for the k-mer to the overall score for this query sequence
            // for the overall score, bit/2 is a sufficient sensitivity and we can use the capacity of unsigned short max score in QueryScore better
            queryScore->addScores(seqList, indexTabListSize, (kmerMatchScore/4));
        }
        biasCorrection -= deltaS[pos];
        biasCorrection += deltaS[pos + kmerSize];
        pos++;
    }
    // write statistics
    seq->stats->kmersPerPos = ((float)kmerListLen/(float)seq->L);
    seq->stats->dbMatches = numMatches;

}
