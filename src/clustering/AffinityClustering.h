#ifndef AFFINITY_CLUSTERING_H
#define AFFINITY_CLUSTERING_H

#include <iostream>
#include <list>
#include "SetElement.h"
#include "LinearMultiArray.h"

class AffinityClustering {
    public:
        AffinityClustering(size_t set_count, size_t unique_element_count, size_t all_element_count,
                unsigned int *element_size_lookup, float **similarities, unsigned int **setids,  size_t iterationnumber, unsigned int convergenceIterations,float input_lambda, float preference,std::list<int>*validids);
        ~AffinityClustering();


        std::list<set *> execute();
    static void add_to_set(const unsigned int element_id,  set * curr_set,const unsigned int set_id);

    private:
        unsigned int iterationnumber;

        size_t set_count;
        size_t unique_element_count;
        size_t all_element_count;
        unsigned int convergenceIterations;
        unsigned int *element_size_lookup;
        float **similarities;
        unsigned int **setids;
        float input_lambda;
        float preference;
        std::list<int> * validids;
};

#endif
