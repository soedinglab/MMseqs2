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
    for(size_t i = 0 ; i < this->divideStepCount - 1; i++){
        delete[] outputArray[i];
    }
    delete [] outputArray;
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
}


void KmerGenerator::initResultList(size_t divide_steps){
    outputArray = new std::pair<short, unsigned int> *[divide_steps];
    for(size_t i = 0 ; i < divide_steps - 1; i++){
        outputArray[i] = new std::pair<short, unsigned int> [VEC_LIMIT];
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
        // get highest element in array for index
        const std::pair<short,unsigned int>  score=(const std::pair<short,unsigned int> ) extMatrix->scoreMatrix[index][0];
        this->highestScorePerArray[i]=score.first; //highest score
        dividerBefore+=divider;
        
    }
    for(size_t i = this->divideStepCount -1; i >= 1 ; i--){
        this->possibleRest[i-1] =this->highestScorePerArray[i]+ possibleRest[i];
    }
    
    // create kmer list
    short cutoff1=this->threshold - this->possibleRest[0];
    size_t index=this->kmerIndex[0];
    ExtendedSubstitutionMatrix * extMatrix= this->matrixLookup[0];
    size_t sizeInputMatrix= extMatrix->size;
    const std::pair<short,unsigned int> * inputArray=extMatrix->scoreMatrix[index];
    
    size_t i;
    for(i = 0; i < this->divideStepCount-1; i++){
        const size_t index=this->kmerIndex[i+1];
        extMatrix= this->matrixLookup[i+1];
        const std::pair<short, unsigned int >  * nextIndexScoreArray=extMatrix->scoreMatrix[index];
        
        int lastElm=calculateArrayProduct(inputArray,
                                   sizeInputMatrix,
                                   nextIndexScoreArray,
                                   extMatrix->size,
                                   outputArray[i],
                                   cutoff1,
                                   possibleRest[i+1],
                                   this->stepMultiplicator[i+1]);
        if(lastElm==-1){
            retList.count=0;
            retList.scoreKmerList=NULL;
            return retList;
        }
            
        inputArray=(const std::pair<short,unsigned int> *) this->outputArray[i];
        cutoff1 = -1000; // we need all that came through 
        retList.count = lastElm + 1;
        sizeInputMatrix = retList.count; // because old data can be under it
    }
    retList.scoreKmerList=outputArray[i-1];
    return retList;
}




int KmerGenerator::calculateArrayProduct( const std::pair<short,unsigned int> * array1,
                                          const size_t array1Size,
                                          const std::pair<short,unsigned int> * array2,
                                          const size_t array2Size,
                                          std::pair<short,unsigned int> * outputvec,
                                          const short cutoff1,const short possibleRest,
                                          const unsigned int pow){
    int counter=-1;
    for(size_t i = 0 ; i< array1Size;i++){
        const std::pair<short, unsigned int> array1Pair=array1[i];
        const short score_i = array1Pair.first;
        const unsigned int kmer_i = array1Pair.second;
        if(score_i < cutoff1 )
            break;
        const short cutoff2=this->threshold-score_i-possibleRest;
        for(size_t j = 0; j < array2Size;j++){
            if(counter+1 >= (int) VEC_LIMIT)
                return counter;
            
            const std::pair<short, unsigned int> array2Pair=array2[j];
            const short score_j = array2Pair.first;
            const unsigned int kmer_j = array2Pair.second;

            if(score_j < cutoff2)
                break;
            counter++;
            std::pair<short, unsigned int> * result=&outputvec[counter];
            result->first=score_i+score_j;
            result->second=kmer_i+(kmer_j*pow);

        }
    }
    return counter;
}

