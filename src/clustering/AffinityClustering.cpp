#include <glob.h>
#include <stddef.h>
#include <float.h>
#include "AffinityClustering.h"

AffinityClustering::AffinityClustering(size_t set_count, size_t unique_element_count, size_t all_element_count,unsigned int *element_size_lookup, double **similarities, unsigned int **setids,  size_t iterationnumber, double input_lambda) {
    this->set_count=set_count;
    this->unique_element_count=unique_element_count;
    this->all_element_count=all_element_count;
    this->element_size_lookup=element_size_lookup;
    this->similarities=similarities;
    this->setids=setids;
    this->input_lambda=input_lambda;
    this->iterationnumber= iterationnumber;

}

AffinityClustering::~AffinityClustering(){
    //delete set_data;
}



std::list<set *> AffinityClustering::execute(){
    //dampening factor (default value is 0.5)
    double lambda=0;
    std::list<set *> result;
    double * currentsimilarities;
    //data structures
    double ** availabilities=new double *[set_count];
    double * availabilitiesData = new double[all_element_count];
    double ** responsibilities=new double *[set_count];
    double * responsibilitiesData=new double [all_element_count];
    double * selfavailabilities=new double [set_count];


    std::fill_n(availabilitiesData, all_element_count , 0);
    size_t curr_pos  = 0;
    for(size_t i = 0; i < set_count; i++){
        availabilities[i] = &availabilitiesData[curr_pos];
        responsibilities[i] = &responsibilitiesData[curr_pos];
        similarities[i][0]=0;//input preference minimal to get a lot of clusters TODO set by parameter
        curr_pos  += element_size_lookup[i];
    }

    for (size_t j = 0; j < iterationnumber; j++) {
    //responsibilities
    for(size_t i = 0; i < set_count; i++) {
        const unsigned int *currentset = setids[i];
        currentsimilarities = similarities[i];
        const size_t currentsetsize=element_size_lookup[i];

        double maxresp1 = -DBL_MAX;
        double maxresp2 = -DBL_MAX;
        size_t maxposition1 = -1;
        // detemine 2 max values for responsibility formula
        for (size_t k = 0; k < currentsetsize ; k++) {
          //  std::cout << currentset[k] <<"\t" << currentsimilarities[k]<< "\n";
            if(availabilities[i][k]+currentsimilarities[k]>maxresp1){
                maxposition1=currentset[k];
                maxresp2=maxresp1;
                maxresp1=availabilities[i][k]+currentsimilarities[k];

            }else if (availabilities[i][k]+currentsimilarities[k]>maxresp2){
                maxresp2=availabilities[i][k]+currentsimilarities[k];
            }
        }

        //compute responsibilities
        for (size_t k = 0; k < currentsetsize ; k++) {
            if(currentset[k]==maxposition1){
                responsibilities[i][k]=responsibilities[i][k]*lambda+(1-lambda)*currentsimilarities[k]-(maxresp2);
            }else{
                responsibilities[i][k]=responsibilities[i][k]*lambda+(1-lambda)*currentsimilarities[k]-(maxresp1);
            }
        }
    }
    //availability update
        //determine self availabilities
        std::fill_n(selfavailabilities, set_count, 0);
        for(size_t i = 0; i < set_count; i++) {
            const unsigned int *currentset = setids[i];
            for (size_t k = 1; k < element_size_lookup[i] ; k++) {
                //random memory access
                selfavailabilities[currentset[k]]+=std::max(0.0, responsibilities[i][k]);
            }
        }
        //determine other availabilities
        for(size_t i = 0; i < set_count; i++) {
                const unsigned int *currentset = setids[i];
                availabilities[i][0]=availabilities[i][0]*lambda+(1-lambda)*selfavailabilities[i];
            for (size_t k = 1; k < element_size_lookup[i] ; k++) {
                //random memory access
                 availabilities[i][k]=availabilities[i][k]*lambda+(1-lambda)*std::min(0.0,responsibilities[currentset[k]][0]+(selfavailabilities[currentset[k]]-std::max(0.0,responsibilities[i][k])));
            }
        }
        //after initialisation, set lambda input value
        lambda=input_lambda;
        /*
        for(size_t i = 0; i < set_count; i++) {
            currentset = setids[i];
            for (int k = 0; k < element_size_lookup[i]; k++) {
                std::cout << j <<"\t"<< i << "\t" << currentset[k] << "\t" << similarities[i][k] << "\t" << availabilities[i][k] << "\t" << responsibilities[i][k] << "\n";
            }
        }
        */
    }
    //compute result
    set* sets = new set[set_count];
    memset(sets, 0, sizeof(set *)*(set_count+1));

    for(size_t i = 0; i < set_count; i++) {
        const unsigned int *currentset = setids[i];
        int maxk=0;
        double maxvalue=-DBL_MAX;
        for (size_t k = 0; k < element_size_lookup[i] ; k++) {
            //std::cout << i<<"\t"<< currentset[k]<<"\t"<< similarities[i][k]<<"\t"<< availabilities[i][k] <<"\t"<< responsibilities[i][k]<<"\n";
           // std::cout << i << "\t" << currentset[k] << "\t" << similarities[i][k] << "\n";
            if(maxvalue<availabilities[i][k]+responsibilities[i][k]){
                maxvalue=availabilities[i][k]+responsibilities[i][k];
                maxk=currentset[k];
            }
            //debugging
           // add_to_set(currentset[k],&sets[i],i);
        }
        //add i to set k
        add_to_set(i,&sets[maxk],maxk);
    }



    for(size_t i = 0; i < this->set_count; i++) {
        set * max_set = &sets[i];
        if (max_set->elements == NULL)
            continue;
        result.push_back(max_set); // O(1)
    }

    return result;
}


void AffinityClustering::add_to_set(const unsigned int element_id, set * curr_set, const unsigned int set_id){
       // set up doubled linked list + fill element_lookup

        // init element with id
        set::element * curr_element_ptr=new set::element();
        curr_element_ptr->element_id=element_id;
        curr_element_ptr->last=NULL;
        if(curr_set->elements == NULL) {// first ptr is not yet set
            curr_element_ptr->next == NULL;
        }
        else {
            set::element *element_first_ptr = curr_set->elements;
            element_first_ptr->last=curr_element_ptr;
            curr_element_ptr->next = element_first_ptr;
        }
        // navigation (double linked list, parent)
        curr_element_ptr->parent_set = curr_set;
    curr_set->elements=curr_element_ptr;
    curr_set->set_id=set_id;




}