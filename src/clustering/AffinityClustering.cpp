#include <glob.h>
#include <stddef.h>
#include <float.h>
#include "AffinityClustering.h"

AffinityClustering::AffinityClustering(size_t set_count, size_t unique_element_count, int all_element_count,unsigned int *element_size_lookup, double **similarities, unsigned int **setids,  size_t iterationnumber, double input_lambda) {
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
    unsigned int * currentset;
    double * currentsimilarities;
    bool init=true;
    //data structures
    double ** availabilities=new double *[all_element_count];
    double ** responsibilities=new double *[all_element_count];
    double * selfavailabilities=new double [set_count];

    for (int j = 0; j < iterationnumber; j++) {
    //responsibilities
    for(size_t i = 0; i < set_count; i++) {
        currentset = setids[i];
        currentsimilarities = similarities[i];
        //init availabilities and responsibility array, set preference score
        if (init) {
            availabilities[i] = new double[element_size_lookup[i]];
            std::fill_n(availabilities[i], element_size_lookup[i], 0);
            responsibilities[i] = new double[element_size_lookup[i]];
            currentsimilarities[0]=0; //input preference minimal to get a lot of clusters TODO set by parameter
            }

        double maxresp1 = -DBL_MAX;
        double maxresp2 = -DBL_MAX;
        int maxposition1 = -1;
        // detemine 2 max values for responsibility formula
        for (int k = 0; k < element_size_lookup[i] ; k++) {
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
        for (int k = 0; k < element_size_lookup[i] ; k++) {
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
            currentset = setids[i];
            for (int k = 1; k < element_size_lookup[i] ; k++) {
                selfavailabilities[currentset[k]]+=std::max(0.0, responsibilities[i][k]);
            }
        }
        //determine other availabilities
        for(size_t i = 0; i < set_count; i++) {
                currentset = setids[i];
                availabilities[i][0]=availabilities[i][0]*lambda+(1-lambda)*selfavailabilities[i];
            for (int k = 1; k < element_size_lookup[i] ; k++) {
                 availabilities[i][k]=availabilities[i][k]*lambda+(1-lambda)*std::min(0.0,responsibilities[currentset[k]][0]+(selfavailabilities[currentset[k]]-std::max(0.0,responsibilities[i][k])));
            }
        }
        init=false;
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
        int maxk=0;
        double maxvalue=-DBL_MAX;
        currentset = setids[i];
        for (int k = 0; k < element_size_lookup[i] ; k++) {
            //std::cout << i<<"\t"<< currentset[k]<<"\t"<< similarities[i][k]<<"\t"<< availabilities[i][k] <<"\t"<< responsibilities[i][k]<<"\n";
           // std::cout << i << "\t" << currentset[k] << "\t" << similarities[i][k] << "\n";
            if(maxvalue<availabilities[i][k]+responsibilities[i][k]){
                maxvalue=availabilities[i][k]+responsibilities[i][k];
                maxk=currentset[k];
            }
        }
        //add i to set k
        add_to_set(i,&sets[maxk],i);
    }



    for(int i = 0; i < this->set_count; i++) {
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