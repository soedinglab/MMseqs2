#include "KmerGenerator.h"
#include <emmintrin.h>
#include <mmintrin.h>
#include <smmintrin.h>

KmerGenerator::KmerGenerator(size_t kmerSize,size_t alphabetSize, short threshold,
                             ExtendedSubstitutionMatrix * three,ExtendedSubstitutionMatrix * two ){
    this->threshold = threshold;
    this->kmerSize = kmerSize;
    this->three = three;
    this->two = two;
    this->indexer = new Indexer((int) alphabetSize, (int)kmerSize);
    calcDivideStrategy();
}

KmerGenerator::~KmerGenerator(){
    delete [] this->stepMultiplicator;
    delete [] this->highestScorePerArray;
    delete [] this->possibleRest;
    delete [] this->kmerIndex;
    delete [] this->divideStep;
    delete [] this->matrixLookup;
    for(size_t i = 0 ; i < this->divideStepCount - 1; i++){
        delete[] outputScoreArray[i];
        delete[] outputIndexArray[i];

    }
    delete [] outputScoreArray;
    delete [] outputIndexArray;
    delete indexer;
}

void KmerGenerator::calcDivideStrategy(){
    const size_t threeDivideCount = this->kmerSize / 3;

    switch(kmerSize%3){
        case 0:
            this->divideStepCount=threeDivideCount;
            this->matrixLookup= new ExtendedSubstitutionMatrix*[divideStepCount];
            this->divideStep = new unsigned int[divideStepCount];
            for(size_t i = 0; i < threeDivideCount; i++){
                this->divideStep[i] = 3;
                this->matrixLookup[i] = three;
            }
            break;
        case 1: 
            this->divideStepCount=threeDivideCount+1;
            this->matrixLookup= new ExtendedSubstitutionMatrix*[divideStepCount];
            
            this->divideStep = new unsigned int[divideStepCount];
            for(size_t i = 0; i < threeDivideCount-1; i++){
                this->divideStep[i] = 3;
                this->matrixLookup[i] = three;
            }
            this->divideStep[threeDivideCount-1]=2;
            this->matrixLookup[threeDivideCount-1] = two;
            
            this->divideStep[threeDivideCount]=2;
            this->matrixLookup[threeDivideCount] = two;
            
            break;
        case 2:
            this->divideStepCount=threeDivideCount+1;
            this->matrixLookup= new ExtendedSubstitutionMatrix*[divideStepCount];
            this->divideStep = new unsigned int[divideStepCount];
            for(size_t i = 0; i < threeDivideCount; i++){
                this->divideStep[i] = 3;
                this->matrixLookup[i] = three;
            }
            this->divideStep[threeDivideCount]=2;
            this->matrixLookup[threeDivideCount] = two;
            
            break;
    }
    this->stepMultiplicator = new unsigned int[divideStepCount];
    this->highestScorePerArray= new short[divideStepCount];
    // init possibleRest 
    this->possibleRest= new short[divideStepCount];
    this->possibleRest[divideStepCount-1]=0;
    this->kmerIndex = new unsigned int[divideStepCount];
    initResultList(divideStepCount);
//    std::reverse(this->matrixLookup, &this->matrixLookup[divideStepCount]);
//    std::reverse(this->divideStep, &this->divideStep[divideStepCount]);
    Debug(Debug::WARNING) << "Divide step for kmer = ";
    for(int i = 0; i < divideStepCount;i++){
        Debug(Debug::WARNING) << divideStep[i] << " ";
    }
    Debug(Debug::WARNING) << "\n";
}


void KmerGenerator::initResultList(size_t divide_steps){
    outputScoreArray = new short *[divide_steps];
    outputIndexArray = new unsigned int *[divide_steps];

    for(size_t i = 0 ; i < divide_steps - 1; i++){
        outputScoreArray[i] = (short *)        Util::mem_align(16,MAX_KMER_RESULT_SIZE*sizeof(short));
        outputIndexArray[i] = (unsigned int *) Util::mem_align(16,MAX_KMER_RESULT_SIZE*sizeof(unsigned int));
    }
}


KmerGeneratorResult KmerGenerator::generateKmerList(const int * int_seq){
    int dividerBefore=0;
    KmerGeneratorResult retList;

    // pre compute phase
    // find first threshold
    for(size_t i = 0; i < this->divideStepCount; i++){
        int divider=divideStep[i];
        
        const unsigned int index= this->indexer->int2index(int_seq,dividerBefore,dividerBefore+divider);
        this->kmerIndex[i]=index;

        stepMultiplicator[i]=this->indexer->powers[dividerBefore];

        ExtendedSubstitutionMatrix * extMatrix= this->matrixLookup[i];
        const ExtendedSubstitutionMatrix::ScoreMatrix * scoreMatrix = extMatrix->scoreMatrix;
        // get highest element in array for index
        const short score = scoreMatrix->score[index*scoreMatrix->rowSize];
        this->highestScorePerArray[i] = score; //highest score
        dividerBefore+=divider;
        
    }
    for(size_t i = this->divideStepCount -1; i >= 1 ; i--){
        this->possibleRest[i-1] =this->highestScorePerArray[i]+ possibleRest[i];
    }
    
    // create kmer list
    short cutoff1=this->threshold - this->possibleRest[0];
    size_t index=this->kmerIndex[0];
    ExtendedSubstitutionMatrix * extMatrix= this->matrixLookup[0];
    const ExtendedSubstitutionMatrix::ScoreMatrix * inputScoreMatrix = extMatrix->scoreMatrix;
    size_t sizeInputMatrix= extMatrix->size;
    short        * inputScoreArray = &inputScoreMatrix->score[index*inputScoreMatrix->rowSize];
    unsigned int * inputIndexArray = &inputScoreMatrix->index[index*inputScoreMatrix->rowSize];
    size_t i;
    for(i = 0; i < this->divideStepCount-1; i++){
        const size_t index=this->kmerIndex[i+1];
        extMatrix= this->matrixLookup[i+1];
        const ExtendedSubstitutionMatrix::ScoreMatrix * nextScoreMatrix = extMatrix->scoreMatrix;
        short        * nextScoreArray = &nextScoreMatrix->score[index*nextScoreMatrix->rowSize];
        unsigned int * nextIndexArray = &nextScoreMatrix->index[index*nextScoreMatrix->rowSize];

        int lastElm=calculateArrayProduct(inputScoreArray,
                                          inputIndexArray,
                                          sizeInputMatrix,
                                          nextScoreArray,
                                          nextIndexArray,
                                          extMatrix->size,
                                          outputScoreArray[i],
                                          outputIndexArray[i],
                                          cutoff1,
                                          possibleRest[i+1],
                                          stepMultiplicator[i+1]);
        if(lastElm==-1){
            retList.count=0;
            retList.score=NULL;
            retList.index=NULL;
            return retList;
        }
        inputScoreArray=this->outputScoreArray[i];
        inputIndexArray=this->outputIndexArray[i];
        cutoff1 = this->threshold - this->outputScoreArray[i][lastElm]; //because old data can be under it
        retList.count = lastElm;
        sizeInputMatrix = retList.count;
    }
    retList.score=outputScoreArray[i-1];
    retList.index=outputIndexArray[i-1];

    return retList;
}




int KmerGenerator::calculateArrayProduct(const short        * scoreArray1,
                                         const unsigned int * indexArray1,
                                         const size_t array1Size,
                                         const short        * scoreArray2,
                                         const unsigned int * indexArray2,
                                         const size_t array2Size,
                                         short              * outputScoreArray,
                                         unsigned int       * outputIndexArray,
                                         const short cutoff1,
                                         const short possibleRest,
                                         const unsigned int pow){
    int counter=0;
    const __m128i * scoreArray2_simd = (const __m128i *) scoreArray2;
    const __m128i * indexArray2_simd = (const __m128i *) indexArray2;
    __m128i pow_simd     = _mm_set1_epi32(pow);

    for(size_t i = 0 ; i < array1Size; i++){
        const short score_i = scoreArray1[i];
        const unsigned int kmer_i = indexArray1[i];
        
        if(score_i < cutoff1 )
            break;
        const short cutoff2=this->threshold-score_i-possibleRest;
        

        __m128i cutoff2_simd = _mm_set1_epi16(cutoff2);
        __m128i score_i_simd = _mm_set1_epi16(score_i);
        __m128i kmer_i_simd  = _mm_set1_epi32(kmer_i);
        const size_t SIMD_SIZE = 8;
        const size_t array2SizeSIMD = (array2Size/SIMD_SIZE)+1;
        for(size_t j = 0; j < array2SizeSIMD; j++){
            if(counter >= (int) MAX_KMER_RESULT_SIZE) //TODO
                return counter;
            
            __m128i score_j_simd   = _mm_load_si128(scoreArray2_simd + j);
            __m128i kmer_j_1_simd  = _mm_load_si128(indexArray2_simd + (j*2));
            __m128i kmer_j_2_simd  = _mm_load_si128(indexArray2_simd + (j*2+1));
            // score_j < cutoff2 -> fffff, score_j > cutoff2 -> 0000
            __m128i cmp = _mm_cmplt_epi16 (score_j_simd, cutoff2_simd);
            const unsigned int score_j_lt_cutoff = _mm_movemask_epi8(cmp);
            
            
            __m128i * scoreOutput_simd = (__m128i *) &outputScoreArray[counter];
            __m128i * indexOutput_simd = (__m128i *) &outputIndexArray[counter];
            _mm_storeu_si128(scoreOutput_simd,   _mm_add_epi16(score_i_simd,score_j_simd));
            _mm_storeu_si128(indexOutput_simd,   _mm_add_epi32(kmer_i_simd,_mm_mullo_epi32(kmer_j_1_simd, pow_simd)));
            _mm_storeu_si128(indexOutput_simd+1, _mm_add_epi32(kmer_i_simd,_mm_mullo_epi32(kmer_j_2_simd, pow_simd)));
            counter += std::min(SIMD_SIZE,array2Size - (j*SIMD_SIZE)); //protect from running to far
             // if(score_j < cutoff2)
            if (score_j_lt_cutoff > 0){
                for(int vec_index = 0; vec_index < SIMD_SIZE; vec_index++){
                    if(CHECK_BIT(score_j_lt_cutoff,vec_index*2)){ // all with 1 is not a result
                        counter--;
                    }
                }
                break;
            }
        }
    }
    return counter;
}

