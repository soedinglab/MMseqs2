#include "SimpleClustering.h"

SimpleClustering::SimpleClustering(unsigned int set_count,
                     unsigned int unique_element_size,
                     unsigned int all_element_count,
                     unsigned int * element_size_lookup){
    this->set_count = set_count;
    this->unique_element_size = unique_element_size;

    this->sets = new set[set_count];
    memset(this->sets, 0, sizeof(set *)*(set_count+1));

    this->element_lookup     = new linear_multi_array<set::element*>(element_size_lookup, unique_element_size, all_element_count);
    this->set_elements = (set::element *) malloc(sizeof(set::element)*all_element_count);

    this->curr_pos = 0;
    this->add_position = 0;

}

SimpleClustering::~SimpleClustering(){
    delete element_lookup;
    delete[] sets;
    free(set_elements);
}

void SimpleClustering::add_set(const unsigned int * element_ids, const int element_size){

    set::element * element_last_ptr = NULL;
    set::element * element_first_ptr = NULL;
    set * curr_set = &this->sets[curr_pos];
    curr_set->next = NULL;
    curr_set->last = NULL;
    curr_set->elements = NULL;
    curr_set->set_id = curr_pos;

    // set up doubled linked list + fill element_lookup
    for(int i = 0 ; i < element_size; i++) {
        // init element with id
        set::element * curr_element_ptr=&set_elements[add_position+i];
        curr_element_ptr->element_id=element_ids[i];
        if(element_first_ptr == NULL) // first ptr is not yet set
            element_first_ptr = curr_element_ptr;
        // navigation (double linked list, parent)
        curr_element_ptr->parent_set = curr_set;
        curr_element_ptr->last = element_last_ptr;
        if(element_last_ptr != NULL) // not in first iteration
            element_last_ptr->next = curr_element_ptr;
        element_last_ptr = curr_element_ptr;
        // element_lookup fill up
        int element_id=curr_element_ptr->element_id;
        element_lookup->add_value_at(element_id,curr_element_ptr);
    }
    add_position += element_size;
    if(element_last_ptr!=NULL)
        element_last_ptr->next = NULL; // last element points to NULL
    curr_set->elements = element_first_ptr; // set element pointer to current_set

    this->curr_pos++;
}

void SimpleClustering::removeSet(set * s){
    set::element * element=s->elements;
    int s_set_id = s->set_id;
    do{ // for(Element element in elements
        set::element * element_to_remove;
        int element_id=element->element_id;

        std::pair<set::element**,int> element_lookup_structure=element_lookup->get_array(element_id);
        set::element ** element_lookup_array = element_lookup_structure.first;
        int array_size = element_lookup_structure.second;
        for(int i =0; i < array_size;i++){
            element_to_remove = element_lookup_array[i];
            set * parent_set = element_to_remove->parent_set;
            if(parent_set!=NULL && parent_set->set_id!=s_set_id){
                parent_set->elements = unplug_element(element_to_remove,parent_set->elements);
            }
            //delete element_to_remove;
        }

    }while((element=element->next) != NULL);
}

set::element * SimpleClustering::unplug_element(set::element * element_to_unplug,set::element * first_element) {
    set::element * last_element=element_to_unplug->last;
    set::element * next_element=element_to_unplug->next;
    element_to_unplug->last       = NULL;
    element_to_unplug->next       = NULL;
    element_to_unplug->parent_set = NULL;

    if(last_element == NULL && next_element==NULL){
        return NULL;
    }if(last_element == NULL){ // first element
        next_element->last = NULL;
        return next_element;
    } else if (next_element == NULL) { // end of list
        last_element->next = NULL;
    } else { // middle of the list
        last_element->next = next_element;
        next_element->last = last_element;
    }
    return first_element;
}

std::list<set *> SimpleClustering::execute(){
    std::list<set *> result;
    for(int i = 0; i < this->set_count; i++) {
        set * max_set = &this->sets[i];
        if (max_set->elements == NULL)
            continue;
        removeSet(max_set);
        result.push_back(max_set); // O(1)
    }
    return result;
}
