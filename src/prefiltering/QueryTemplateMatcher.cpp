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
    this->kmerGenerator = new KmerGenerator(kmerSize, m->alphabetSize, kmerThr);
    this->kmerGenerator->setDivideStrategy(_3merSubMatrix->scoreMatrix, _2merSubMatrix->scoreMatrix );
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

std::list<hit_t>* QueryTemplateMatcher::matchQuery (Sequence * seq){
    queryScore->reset();
    seq->resetCurrPos();

    match(seq);

    queryScore->setPrefilteringThresholds();

    return queryScore->getResult(seq->L);
}

void QueryTemplateMatcher::match(Sequence* seq){

//    Indexer* indexer = new Indexer(m->alphabetSize, kmerSize);
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

/*    int overall_score = 0;
    int match_num = 0;
    int match_pos = 0;*/

    int pos = 0;
//    std::cout << "\nQUERY: " << seq->getDbKey() << "\n";
    while(seq->hasNextKmer(kmerSize)){
        const int* kmer = seq->nextKmer(kmerSize);
        // generate k-mer list
        ScoreMatrix kmerList = kmerGenerator->generateKmerList(kmer);
        kmerListLen += kmerList.elementSize;

        // match the index table
//        int pos_matched = 0;
        for (unsigned int i = 0; i < kmerList.elementSize; i++){
            short kmerMatchScore = kmerList.score[i] + (short) biasCorrection;
            seqList = indexTable->getDBSeqList(kmerList.index[i], &indexTabListSize);
            numMatches += indexTabListSize;

/*            if (seq->getId() == 1 && pos == 2 ){
                std::cout << "\t\t";
                indexer->printKmer(retList[i].second, kmerSize, m->int2aa);
                std::cout << " " << retList[i].first << "\n";
            }*/
/*            for (int j = 0; j < indexTabListSize; j++){
                if (seqList[j] == 1){
                    std::cout << "Similar k-mer list pos: " << i << ", score: " << retList[i].first << ", kmer idx: " << retList[i].second << "\n";
                    pos_matched = 1;
                    std::cout << pos << " ";
                    indexer->printKmer(kmer, kmerSize, m->int2aa);
                    std::cout << "->";
                    indexer->printKmer(kmerMatch.second, kmerSize, m->int2aa);
                    std::cout << "\n";
                    std::cout << "\t" << kmerMatchScore << "\n";
                    overall_score+=kmerMatchScore;
                    match_num++;
                }
            }*/

            // add the scores for the k-mer to the overall score for this query sequence
            // for the overall score, bit/2 is a sufficient sensitivity and we can use the capacity of unsigned short max score in QueryScore better
            queryScore->addScores(seqList, indexTabListSize, (kmerMatchScore/4));
        }
        biasCorrection -= deltaS[pos];
        biasCorrection += deltaS[pos + kmerSize];
//        if(pos_matched == 1)
//            match_pos++;
        pos++;
    }
    //Debug(Debug::WARNING) << "QUERY: " << seq->getDbKey();
    //Debug(Debug::WARNING) << " score = " << overall_score;
    //Debug(Debug::WARNING) << " matched at " << match_pos << " positions. ";
    //Debug(Debug::WARNING) << match_num << " times.\n";
    // write statistics
    seq->stats->kmersPerPos = ((float)kmerListLen/(float)seq->L);
    seq->stats->dbMatches = numMatches;

//    delete indexer;
}
