#include "KmerGenerator.h"

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
    for(size_t i = 0 ; i < this->divideStepCount; i++){
        delete[] outputScoreArray[i];
        delete[] outputIndexArray[i];

    }
    delete [] outputScoreArray;
    delete [] outputIndexArray;
    delete indexer;
}

void KmerGenerator::calcDivideStrategy(){
    const size_t threeDivideCount = this->kmerSize /3;
    
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
}


void KmerGenerator::initResultList(size_t divide_steps){
    outputScoreArray = new short *[divide_steps];
    outputIndexArray = new unsigned int *[divide_steps];

    for(size_t i = 0 ; i < divide_steps; i++){
        outputScoreArray[i] = (short *)        Util::mem_align(16,VEC_LIMIT*sizeof(short));
        outputIndexArray[i] = (unsigned int *) Util::mem_align(16,VEC_LIMIT*sizeof(unsigned int));
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
    const ExtendedSubstitutionMatrix::ScoreMatrix * scoreMatrix = extMatrix->scoreMatrix;

    const size_t sizeExtendedMatrix= extMatrix->size;
    
    short        * inputScoreArray = &scoreMatrix->score[index*scoreMatrix->rowSize];
    unsigned int * inputIndexArray = &scoreMatrix->index[index*scoreMatrix->rowSize];
    size_t i;
    for(i = 0; i < this->divideStepCount-1; i++){
        const size_t index=this->kmerIndex[i+1];
        extMatrix= this->matrixLookup[i+1];
        const ExtendedSubstitutionMatrix::ScoreMatrix * scoreMatrix = extMatrix->scoreMatrix;

        short        * nextScoreArray = &scoreMatrix->score[index*scoreMatrix->rowSize];
        unsigned int * nextIndexArray = &scoreMatrix->index[index*scoreMatrix->rowSize];

        int lastElm=calculateArrayProduct(inputScoreArray,
                                          inputIndexArray,
                                          sizeExtendedMatrix,
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
        retList.count=lastElm+1;
        
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
    int counter=-1;
    for(size_t i = 0 ; i< array1Size;i++){
        const short score_i = scoreArray1[i];
        const unsigned int kmer_i = indexArray1[i];
        if(score_i < cutoff1 )
            break;
        const short cutoff2=this->threshold-score_i-possibleRest;
        for(size_t j = 0; j < array2Size;j++){
            if(counter+1 >= (int) VEC_LIMIT)
                return counter;
            
            const short score_j = scoreArray2[j];
            const unsigned int kmer_j = indexArray2[j];

            if(score_j < cutoff2)
                break;
            counter++;
            outputScoreArray[counter]=score_i+score_j;
            outputIndexArray[counter]=kmer_i+(kmer_j*pow);

        }
    }
    return counter;
}

