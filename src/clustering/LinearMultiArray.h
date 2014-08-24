#ifndef __graphcluster__LinearMultiArray__
#define __graphcluster__LinearMultiArray__
#include <algorithm>
#include <cstdlib>




template<class T> class linear_multi_array
{
public:
    linear_multi_array(unsigned int * element_size_lookup,
                       unsigned int element_size,
                       unsigned int size)
    {
        this->size = std::max((unsigned int)1,size); //at least one
        multi_array = (T*) malloc( sizeof(T)*size);
        next_pos_to_write = (T**) malloc( sizeof(T *)*(element_size+1));
        int last_element=0;
        for(unsigned int i = 0; i <= element_size; i++){
            this->xorSwap(&element_size_lookup[0], &element_size_lookup[i+1]);
            element_size_lookup[i+1] = last_element+element_size_lookup[i+1];
            last_element = element_size_lookup[i+1];
        }
        this->element_offset = element_size_lookup;
        
        
        for(unsigned int i = 0; i <= element_size; i++){
            next_pos_to_write[i] = &multi_array[element_offset[i]];
        }

    }
        
        
    
    ~linear_multi_array() {free(multi_array); free(next_pos_to_write);}

    const std::pair<T*, int> get_array(int i); // get value
    void add_value_at(int i,T element); // push it
    T   *multi_array;                   // the array
    unsigned int size;                  // size of the whole array
    
private:
    unsigned int *element_offset;    // element lookuptable
    T   **next_pos_to_write;
    void xorSwap (unsigned int * x, unsigned int * y);
};

template<class T> void linear_multi_array<T>::xorSwap (unsigned int * x, unsigned int * y) {
    *x ^= *y;
    *y ^= *x;
    *x ^= *y;
}

template<class T> const std::pair<T*, int> linear_multi_array<T>::get_array(int i)
{
    std::pair<T*,int> ret_value;
    int from = element_offset[i];
    int to   = element_offset[i+1];
    int size = to - from;
    ret_value.second = size;
    ret_value.first  = &multi_array[from];
    return ret_value;
}

template<class T> void linear_multi_array<T>::add_value_at(int i,T element)
{
    *(next_pos_to_write[i]) = element;
    next_pos_to_write[i]++;
}

#endif
