//
//  element.h
//  graphcluster
//
//  Created by Martin Steinegger on 21.05.13.
//  Copyright (c) 2013 Martin Steinegger. All rights reserved.
//

#ifndef GRAPHCLUSTER_ELEMENT
#define GRAPHCLUSTER_ELEMENT

struct set {
    struct element {
        set * parent_set;
        element * last;
        element * next;
        unsigned int element_id;
        unsigned short weight;

    };
    element * elements;
    set * last;
    set * next;
    unsigned int set_id;
    unsigned short weight;
};



#endif
