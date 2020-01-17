#include "KmerGenerator.h"
#include <algorithm>    // std::reverse
#include <MathUtil.h>
#include "simd.h"


KmerGenerator::KmerGenerator(size_t kmerSize, size_t alphabetSize, short threshold ){
    this->threshold = threshold;
    this->kmerSize = kmerSize;
    this->indexer = new Indexer((int) alphabetSize, (int)kmerSize);
//    calcDivideStrategy();
}

void KmerGenerator::setThreshold(short threshold){
    this->threshold = threshold;
}
KmerGenerator::~KmerGenerator(){
    delete [] this->stepMultiplicator;
    delete [] this->highestScorePerArray;
    delete [] this->possibleRest;
    delete [] this->kmerIndex;
    delete [] this->divideStep;
    delete [] this->matrixLookup;
    for(size_t i = 0 ; i < 2; i++){
        free( outputScoreArray[i] );
        free( outputIndexArray[i] );
    }
    delete [] outputScoreArray;
    delete [] outputIndexArray;
    delete indexer;
}

void KmerGenerator::setDivideStrategy(ScoreMatrix ** one){
    this->divideStepCount = kmerSize;
    this->matrixLookup = new ScoreMatrix*[divideStepCount];
    this->divideStep   = new unsigned int[divideStepCount];
    for(size_t i = 0; i < kmerSize; i++){
        this->divideStep[i] = 1;
        this->matrixLookup[i] = one[i];
    }
    initDataStructure();
}

void KmerGenerator::setDivideStrategy(ScoreMatrix * three, ScoreMatrix * two){
    const size_t threeDivideCount = this->kmerSize / 3;

    switch(kmerSize%3){
        case 0:
            this->divideStepCount=threeDivideCount;
            this->matrixLookup= new ScoreMatrix*[divideStepCount];
            this->divideStep = new unsigned int[divideStepCount];
            for(size_t i = 0; i < threeDivideCount; i++){
                this->divideStep[i] = 3;
                this->matrixLookup[i] = three;
            }
            break;
        case 1:
            this->divideStepCount=threeDivideCount+1;
            this->matrixLookup= new ScoreMatrix*[divideStepCount];
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
            this->matrixLookup= new ScoreMatrix*[divideStepCount];
            this->divideStep = new unsigned int[divideStepCount];
            for(size_t i = 0; i < threeDivideCount; i++){
                this->divideStep[i] = 3;
                this->matrixLookup[i] = three;
            }
            this->divideStep[threeDivideCount]=2;
            this->matrixLookup[threeDivideCount] = two;

            break;
    }

    initDataStructure();
    std::reverse(this->matrixLookup, &this->matrixLookup[divideStepCount]);
    std::reverse(this->divideStep, &this->divideStep[divideStepCount]);
}


void KmerGenerator::initDataStructure(){
    this->stepMultiplicator = new size_t[divideStepCount];
    this->highestScorePerArray = new short[divideStepCount];
    // init possibleRest
    this->possibleRest = new short[divideStepCount];
    this->possibleRest[divideStepCount-1] = 0;
    this->kmerIndex = new size_t[divideStepCount];

    outputScoreArray = new short *[2];
    outputIndexArray = new size_t *[2];

    for(size_t i = 0 ; i < 2; i++){
        outputScoreArray[i] = (short *)  mem_align(ALIGN_INT, MAX_KMER_RESULT_SIZE * sizeof(short));
        outputIndexArray[i] = (size_t *) mem_align(ALIGN_INT, MAX_KMER_RESULT_SIZE * sizeof(size_t));
    }
}


std::pair<size_t *, size_t> KmerGenerator::generateKmerList(const unsigned char * int_seq, bool addIdentity){
    int dividerBefore=0;
    // pre compute phase
    // find first threshold
    for(size_t i = 0; i < this->divideStepCount; i++){
        const int divider=divideStep[i];

        const unsigned int index = this->indexer->int2index(int_seq,dividerBefore,dividerBefore+divider);
        this->kmerIndex[i] = index;

        stepMultiplicator[i] = this->indexer->powers[dividerBefore];

        const ScoreMatrix * scoreMatrix = this->matrixLookup[i];
        // get highest element in array for index
        const short score = scoreMatrix->score[index*scoreMatrix->rowSize];
        this->highestScorePerArray[i] = score; //highest score
        dividerBefore+=divider;

    }
    for(size_t i = this->divideStepCount -1; i >= 1 ; i--){
        this->possibleRest[i-1] = this->highestScorePerArray[i] + possibleRest[i];
    }

    // create kmer list
    short cutoff1 = this->threshold - this->possibleRest[0];
    const size_t index  = this->kmerIndex[0];
    const ScoreMatrix * inputScoreMatrix = this->matrixLookup[0];
    size_t sizeInputMatrix = inputScoreMatrix->elementSize;
    const short  * inputScoreArray = &inputScoreMatrix->score[index*inputScoreMatrix->rowSize];
    for(size_t pos = 0; pos < inputScoreMatrix->rowSize && inputScoreArray[pos] >= cutoff1; pos++){
        outputIndexArray[1][pos] = inputScoreMatrix->index[index*inputScoreMatrix->rowSize + pos];
    }
    const size_t * inputIndexArray = outputIndexArray[1];
    size_t i;

    for(i = 0; i < this->divideStepCount-1; i++){
        const size_t index = this->kmerIndex[i+1];
        const ScoreMatrix * nextScoreMatrix = this->matrixLookup[i+1];
        const short        * nextScoreArray = &nextScoreMatrix->score[index*nextScoreMatrix->rowSize];
        const unsigned int * nextIndexArray = &nextScoreMatrix->index[index*nextScoreMatrix->rowSize];

        const size_t lastElm=calculateArrayProduct(inputScoreArray,
                                                   inputIndexArray,
                                                   sizeInputMatrix,
                                                   nextScoreArray,
                                                   nextIndexArray,
                                                   nextScoreMatrix->elementSize,
                                                   outputScoreArray[i%2],
                                                   outputIndexArray[i%2],
                                                   cutoff1,
                                                   possibleRest[i+1],
                                                   stepMultiplicator[i+1]);

        inputScoreArray = this->outputScoreArray[i%2];
        inputIndexArray = this->outputIndexArray[i%2];
        cutoff1 = -1000; //all must be inspected
        sizeInputMatrix = lastElm;
    }

    // add identity of input kmer
    if(addIdentity  && sizeInputMatrix == 0){
        outputScoreArray[0][0] = 0;
        outputIndexArray[0][0] = 0;
        // create first kmer
        for(unsigned int z = 0; z < this->divideStepCount; z++){
            const size_t index = this->kmerIndex[z];
            const ScoreMatrix * nextScoreMatrix = this->matrixLookup[z];
            const short        * nextScoreArray = &nextScoreMatrix->score[index*nextScoreMatrix->rowSize];
            const unsigned int * nextIndexArray = &nextScoreMatrix->index[index*nextScoreMatrix->rowSize];
            outputScoreArray[0][0] += nextScoreArray[0];
            outputIndexArray[0][0] += static_cast<size_t>(nextIndexArray[0]) * stepMultiplicator[z];
        }

        return std::make_pair(outputIndexArray[0], 1);
    }
    return std::make_pair(outputIndexArray[(i-1)%2], sizeInputMatrix);
}


size_t KmerGenerator::calculateArrayProduct(const short        * __restrict scoreArray1,
                                            const size_t       * __restrict indexArray1,
                                            const size_t array1Size,
                                            const short        * __restrict scoreArray2,
                                            const unsigned int * __restrict indexArray2,
                                            const size_t array2Size,
                                            short              * __restrict outputScoreArray,
                                            size_t             * __restrict outputIndexArray,
                                            const short cutoff1,
                                            const short possibleRest,
                                            const size_t pow){
    size_t counter=0;
    for(size_t i = 0 ; i< array1Size;i++){
        const short score_i = scoreArray1[i];
        const size_t kmer_i = indexArray1[i];
        if(score_i < cutoff1 )
            break;
        const short cutoff2=this->threshold-score_i-possibleRest;
        for(size_t j = 0; j < array2Size && (counter+1 < (int) MAX_KMER_RESULT_SIZE) && (scoreArray2[j] >= cutoff2); j++){
            const short score_j = scoreArray2[j];
            const size_t kmer_j = indexArray2[j];
            outputScoreArray[counter]=score_i+score_j;
            outputIndexArray[counter]=kmer_i+(kmer_j*pow);
            counter++;
        }
        if(counter+1 >= (int) MAX_KMER_RESULT_SIZE){
            return counter;
        }
    }
    return counter;
}

