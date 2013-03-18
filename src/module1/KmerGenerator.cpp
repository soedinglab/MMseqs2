#include "KmerGenerator.h"

    KmerGenerator::KmerGenerator(size_t kmer_size,size_t alphabet_size, short threshold, 
                                 ExtendedSubstitutionMatrix * three,ExtendedSubstitutionMatrix * two ){
        this->threshold = threshold;
        this->kmer_size = kmer_size;
        this->three = three;
        this->two = two;
        this->indexer = new Indexer(alphabet_size, kmer_size);
        calc_divide_strategy();
    }

    KmerGenerator::~KmerGenerator(){
        delete [] this->pow_per_step;
        delete [] this->max_score_per_vec;
        delete [] this->possible_rest;
        delete [] this->kmer_index;
        delete [] this->divide_steps;
        delete [] this->matrixLookup;
        for(size_t i = 0 ; i < this->divide_steps_count; i++){
            delete outputvec[i];
        }     
        delete [] outputvec;
    }

    void KmerGenerator::calc_divide_strategy(){
        size_t three_dividecount = this->kmer_size /3;
        
        switch(kmer_size%3){
            case 0:
                this->divide_steps_count=three_dividecount;
                this->matrixLookup= new ExtendedSubstitutionMatrix*[divide_steps_count];
                this->divide_steps = new size_t[divide_steps_count];
                for(int i = 0; i < three_dividecount; i++){
                    this->divide_steps[i] = 3;
                    this->matrixLookup[i] = three;
                }
                break;
            case 1: 
                this->divide_steps_count=three_dividecount+1;
                this->matrixLookup= new ExtendedSubstitutionMatrix*[divide_steps_count];
                
                this->divide_steps = new size_t[divide_steps_count];
                for(int i = 0; i < three_dividecount-1; i++){
                    this->divide_steps[i] = 3;
                    this->matrixLookup[i] = three;
                }
                this->divide_steps[three_dividecount-1]=2;
                this->matrixLookup[three_dividecount-1] = two;
                
                this->divide_steps[three_dividecount]=2;
                this->matrixLookup[three_dividecount] = two;
                
                break;
            case 2:
                this->divide_steps_count=three_dividecount+1;
                this->matrixLookup= new ExtendedSubstitutionMatrix*[divide_steps_count];
                this->divide_steps = new size_t[divide_steps_count];
                for(int i = 0; i < three_dividecount; i++){
                    this->divide_steps[i] = 3;
                    this->matrixLookup[i] = three;
                }
                this->divide_steps[three_dividecount]=2;
                this->matrixLookup[three_dividecount] = two;
                
                break;
        }
        
        this->pow_per_step = new size_t[divide_steps_count];
        this->max_score_per_vec= new short[divide_steps_count];
        this->possible_rest= new short[divide_steps_count];
        this->possible_rest[divide_steps_count-1]=0;
        this->kmer_index = new unsigned int[divide_steps_count];
        init_result_lists(divide_steps_count);
    }


    void KmerGenerator::init_result_lists(size_t divide_steps){
        outputvec = new std::vector<std::pair<short, size_t> >*[divide_steps];
        for(size_t i = 0 ; i < divide_steps; i++){
            outputvec[i] = new std::vector<std::pair<short, size_t> >(VEC_LIMIT);
        }
    }
    kmer_list KmerGenerator::generateKmerList(const int * kmer){
        size_t divider_before=0;
        kmer_list retList;
        //find first threshold
        for(int i = 0; i < this->divide_steps_count; i++){
            size_t divider=divide_steps[i];
            const unsigned int index= this->indexer->int2index(kmer,divider_before,divider_before+divider);
            pow_per_step[i]=this->indexer->powers[divider_before];
            this->kmer_index[i]=index;
            ExtendedSubstitutionMatrix * extMatrix= this->matrixLookup[i];
            std::pair<short,size_t> score=extMatrix->scoreMatrix[index]->at(0);
            this->max_score_per_vec[i]=score.first; //highest score
            divider_before+=divider;
        }
        for(size_t i = this->divide_steps_count -1; i >= 1 ; i--){
            this->possible_rest[i-1] =this->max_score_per_vec[i]+ possible_rest[i];
        }
        // create kmer list
        short cutoff1=this->threshold - this->possible_rest[0];
        
        
        size_t index=this->kmer_index[0];
        ExtendedSubstitutionMatrix * extMatrix= this->matrixLookup[0];
        std::vector<std::pair<short, size_t > > * vec_input=extMatrix->scoreMatrix[index];
        size_t i;
        for(i = 0; i < this->divide_steps_count-1; i++){
            size_t index=this->kmer_index[i+1];
            extMatrix= this->matrixLookup[i+1];
            std::vector<std::pair<short, size_t > > * next_vec=extMatrix->scoreMatrix[index];
            
            int lastElm=calcProduct(vec_input,next_vec,outputvec[i], 
                                       cutoff1,possible_rest[i+1],
                                      this->pow_per_step[i+1]);
            if(lastElm==-1){
                retList.count=0;
                retList.score_kmer_list=NULL;
                return retList;
            }
                
            vec_input=this->outputvec[i];
            cutoff1 = this->threshold - this->outputvec[i]->at(lastElm).first; //because old data can be under it   
            retList.count=lastElm+1;
            
        }
        retList.score_kmer_list=outputvec[i-1];
        return retList;
    }




    int KmerGenerator::calcProduct(const std::vector<std::pair<short, size_t> > * vec1,
                                      const std::vector<std::pair<short, size_t> > * vec2,
                                      std::vector<std::pair<short, size_t> > * outputvec, 
                                      const short cutoff1,const short possible_rest,const size_t pow){
        int counter=-1;
        for(size_t i = 0 ; i< vec1->size();i++){
            const std::pair<short, size_t> vec1_pair=vec1->at(i);
            const short score_i = vec1_pair.first;
            const size_t kmer_i = vec1_pair.second;
        if(score_i < cutoff1 )
            break;
        const short cutoff2=this->threshold-score_i-possible_rest;
        const size_t vec2_size=vec2->size();
        for(int j = 0; j < vec2_size;j++){
            if(counter+1 >= VEC_LIMIT)
                return counter;
            
            const std::pair<short, size_t> vec2_pair=vec2->at(j);
            const short score_j = vec2_pair.first;
            const size_t kmer_j = vec2_pair.second;

            if(score_j < cutoff2)
                break;
            counter++;
            std::pair<short, size_t> * result=&outputvec->at(counter);
            result->first=score_i+score_j;
            result->second=kmer_i+(kmer_j*pow);

        }
    }
    return counter;
}
