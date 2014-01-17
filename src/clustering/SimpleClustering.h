#ifndef SIMPLE_CLUSTERING_H
#define SIMPLE_CLUSTERING_H

#include <iostream>
#include <list>
#include <string.h>
#include "SetElement.h"
#include "LinearMultiArray.h"

class SimpleClustering {
    public:
        SimpleClustering(unsigned int set_count,
                unsigned int unique_element_size,
                unsigned int all_element_count,
                unsigned int * element_size_lookup);
        ~SimpleClustering();

        void add_set(const unsigned int * element_ids,
                const int element_size);
        void removeSet(set * s);
        set::element * unplug_element(set::element * element_to_unplug,set::element * first_element);
        std::list<set *> execute();

    private:
        unsigned int add_position;
        int unique_element_size;
        int all_element_count;
        int set_count;

        // pointer to a corresponding set for each representative sequence
        set* sets;

        int curr_pos;

        linear_multi_array<set::element *> * element_lookup;

        set::element * set_elements;
};

#endif
