#ifndef AFFINITY_CLUSTERING_H
#define AFFINITY_CLUSTERING_H

#include <iostream>
#include <list>
#include <string.h>
#include "SetElement.h"
#include "LinearMultiArray.h"

class AffinityClustering {
    public:
        AffinityClustering(size_t set_count, size_t unique_element_count, int all_element_count,
                unsigned int *element_size_lookup, double **similarities, unsigned int **setids,  size_t iterationnumber, double input_lambda);
        ~AffinityClustering();


        std::list<set *> execute();
    void add_to_set(const unsigned int element_id,  set * curr_set,const unsigned int set_id);

    private:
        unsigned int iterationnumber;

        size_t set_count;
        size_t unique_element_count;
        int all_element_count;
        unsigned int *element_size_lookup;
        double **similarities;
        unsigned int **setids;
        double input_lambda;
};

#endif
