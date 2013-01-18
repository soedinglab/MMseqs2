#include "QueryTemplateMatcher.h"
QueryTemplateMatcher::QueryTemplateMatcher ( short kmerThreshold,
                                             short queryScoreThreshold, 
                                             size_t kmer_size,
                                             SubstitutionMatrix * substitutionMatrix,
                                             IndexTable * indexTable){
    this->kmer_size = kmer_size;
    this->indexTable = indexTable;
    ReducedMatrix  ReducedMatrix(substitutionMatrix->probMatrix,
                               substitutionMatrix->aa2int,
                               substitutionMatrix->int2aa,
                               substitutionMatrix->ALPHABET_SIZE,
                               substitutionMatrix->ALPHABET_SIZE-16);
    ExtendedSubstitutionMatrix threeExtendedSubstitutionMatrix(ReducedMatrix.reduced_Matrix, 3,
                                                               ReducedMatrix.reduced_alphabet_size);
    ExtendedSubstitutionMatrix twoExtendedSubstitutionMatrix(ReducedMatrix.reduced_Matrix, 2,
                                                             ReducedMatrix.reduced_alphabet_size);
    this->kmerGenerator = new KmerGenerator(kmer_size, 
                                            ReducedMatrix.reduced_alphabet_size, kmerThreshold, 
                                            &threeExtendedSubstitutionMatrix, &twoExtendedSubstitutionMatrix);
}

QueryTemplateMatcher::~QueryTemplateMatcher (){
    
}



std::list<hit_t>* QueryTemplateMatcher::matchQuery (Sequence * seq){
    queryScore->reset();
    
    while(seq->hasNextKmer(this->kmer_size)){
        const int * curr_pos= seq->nextKmer(kmer_size);
        kmer_list kmerList=kmerGenerator->generateKmerList(curr_pos);
        std::vector<std::pair<short,size_t> > * kmer_vec=kmerList.score_kmer_list;
        for(int i = 0; i < kmerList.count;i++){
            const std::pair<short,size_t> curr_kmer_pair=kmer_vec->at(i);
            const size_t kmerIdx = curr_kmer_pair.second;
            const short kmerScore =  curr_kmer_pair.first;
            size_t listSize = 0;
            size_t* seqList = indexTable->getDBSeqList(kmerIdx, &listSize);
            queryScore->addScores(seqList,listSize,kmerScore);
        }
    }
    return queryScore->getResult();    
}
