//
//  SetCover.h
//  graphcluster
//
//  Created by Martin Steinegger on 21.05.13.
//  Copyright (c) 2013 Martin Steinegger. All rights reserved.
//

#ifndef __graphcluster__SetCover__
#define __graphcluster__SetCover__

#include <iostream>
#include <list>
#include "SetElement.h"
#include "LinearMultiArray.h"
class SetCover {
public:
    SetCover(unsigned int set_size,
              unsigned int element_size,
              unsigned int weight_range,
              unsigned int all_element_count,
              unsigned int * element_size_lookup);
    ~SetCover();

    void add_set(const int set_id, const int set_weight,
                 const unsigned int * element_ids,
                 const unsigned short * weights,
                 const int element_size);
    std::list<set *> execute_set_cover();
/*
    get_highest_weighted_set 
    input 
            int start_pos is the search start position 
    output
            returns the highest sets from a certaint position
*/
    set * get_highest_weighted_set(int start_pos);
private:

    unsigned int add_position;
    int element_size;
    int set_size;
    int weight_range;
    set ** ordered_by_score_set;
    set::element * set_elements;
    set * sets;

    linear_multi_array<set::element *> * element_lookup;
    // methodes
    void removeSet(set * s);
    set::element * unplug_element(set::element * element_to_unplug,set::element * first_element);
    void unplug_set(set * set_to_remove);
    set * create_set_at_weight_position(int weight,set * set_to_add);

};
#endif /* defined(__graphcluster__SetCover__) */
