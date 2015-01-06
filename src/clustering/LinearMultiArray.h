#ifndef __graphcluster__LinearMultiArray__
#define __graphcluster__LinearMultiArray__
#include <algorithm>
#include <cstdlib>
#include "Debug.h"
#include "Util.h"


template<class T> class linear_multi_array
{
public:
    linear_multi_array(unsigned int *element_size_lookup,
            size_t element_size,
            size_t size)
    {
        if(size == 0){
            Debug(Debug::ERROR) << "ERROR: Size of linear_multi_array is 0. " << "\n";
            EXIT(EXIT_FAILURE);
        }
        this->size = size; //at least one
        multi_array = (T*) malloc( sizeof(T) * size);
        next_pos_to_write = (T**) malloc( sizeof(T *) * (element_size + 1));
        size_t last_element = 0;
        size_t last_element_2 = 0;

        for(size_t i = 0; i <= element_size; i++){
            std::swap(element_size_lookup[0], element_size_lookup[i+1]);
            element_size_lookup[i+1] = last_element + element_size_lookup[i+1];
            last_element = element_size_lookup[i+1];
            if(element_size_lookup[i+1] > 2147483648)
                std::cout << "FUck my life " << last_element  << " " <<  element_size_lookup[i+1] << std::endl;
        }

        this->element_offset = element_size_lookup;
        for(size_t i = 0; i <= element_size; i++){
            next_pos_to_write[i] = &multi_array[element_offset[i]];
        }
    }

    ~linear_multi_array() {free(multi_array); free(next_pos_to_write);}

    const std::pair<T*, unsigned int> get_array(size_t i); // get value
    void add_value_at(size_t i, T element); // push it
    T   *multi_array;                   // the array
    size_t size;                  // size of the whole array
    
private:
    unsigned int *element_offset;    // element lookuptable
    T   **next_pos_to_write;
};

template<class T> const std::pair<T*, unsigned int> linear_multi_array<T>::get_array(size_t i)
{
    std::pair<T*,unsigned int> ret_value;
    size_t from = element_offset[i];
    size_t to   = element_offset[i+1];
    size_t size = to - from;
    ret_value.second = size;
    ret_value.first  = &multi_array[from];
    return ret_value;
}

template<class T> void linear_multi_array<T>::add_value_at(size_t i, T element)
{
//    if(element_offset[i] == 0){
//        Debug(Debug::ERROR) << "ERROR: Array position " << i << " is has more elements than expected. " << "\n";
//        EXIT(EXIT_FAILURE);
//    }
//    element_offset[i]--;
    *(next_pos_to_write[i]) = element;
    next_pos_to_write[i]++;
}

#endif
