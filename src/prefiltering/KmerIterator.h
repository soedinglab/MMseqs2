class KmerIterator {

    public:
        KmerIterator(unsigned int stepSize,
                     const int8_t * pattern,
                     unsigned int patternSize) :
                stepSize(stepSize), patternSize(patternSize), pattern(pattern)
        {
            kmerSize = 0;
            for(unsigned int i = 0; i < patternSize; i++){
                if(pattern[i]){
                    kmerSize++;
                }
            }
            buffer = new unsigned int[kmerSize];
            reset();
        }

        ~KmerIterator(){
            delete buffer;
        }
    
        inline bool hasNext(){
            return (currentStep < sequenceSizeMinusPattern);
        }

        inline void setSequenceSize(unsigned int seqSize){
            this->sequenceSizeMinusPattern = (seqSize - patternSize) + 1;
        }
        
        inline unsigned int getPosition(){
            return idx;

        }
        // Generates lookup indexes for rotating kmer with spaced pattern
        // Example with stepSize 17
        // i0 = (i0 + 0.7*17) % 17
        // i0 = 0, i = i0 = 0, 17, 34, 51, … 170, 187.
        // i = i0 = 11, 28, 49, … 181, 198.
        // i = i0 = 5, 22, 39, … 175, 192.
        inline const unsigned int* nextKmer(){
            if(hasNext()){
                currentStep++;
                idx = i0+run*stepSize;
                if( idx >= sequenceSizeMinusPattern){
                    i0 = (i0 + int(0.7*stepSize)) % stepSize;
                    run = 0;
                    idx = i0;
                }
                run += 1;
                unsigned int pos = 0;
                for(unsigned int i = 0; i < this->patternSize; i++) {
                    if(pattern[i]) {
                        buffer[pos++] = (idx + i);
                    }
                }
                return buffer;
            }
            return 0;
        }

        inline void reset(){
            idx = 0;
            run = 0;
            i0  = 0;
            currentStep = 0;
        }


    private:
        const unsigned int stepSize;
        const unsigned int patternSize;
        const int8_t * pattern;
        unsigned int idx;
        unsigned int run;
        unsigned int i0;
        unsigned int currentStep;
        int kmerSize;
        unsigned int * buffer;
        unsigned int sequenceSizeMinusPattern;
};

