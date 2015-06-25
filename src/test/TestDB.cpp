#include <iostream>
#include <list>
#include <algorithm>
#include <math.h>
#include <SetCover2.h>


#include "Clustering.h"
#include "SetElement.h"

#include "DBReader.h"
#include "DBWriter.h"

int randBetween(size_t start, size_t end){
    return rand() % (end - start) + start;
}


Clustering::set_data create_random_data(size_t elemSize, size_t setSize) {

    Clustering::set_data ret_struct;


    ret_struct.uniqu_element_count = elemSize;
    ret_struct.set_count = setSize;

    // set element size lookup
    unsigned int * element_buffer = new unsigned int[elemSize ];
    unsigned int * element_size   = new unsigned int[elemSize + 2];
    memset(element_size, 0, sizeof(unsigned int) * (elemSize + 2));
    ret_struct.element_size_lookup = element_size;
    // return set size
    unsigned int * set_size=new unsigned int[setSize];
    ret_struct.set_sizes = set_size;
    // set sets
    unsigned int ** sets = new unsigned int*[setSize];
    ret_struct.sets = sets;
    // set weights
    unsigned short ** weights = new unsigned short*[setSize];
    ret_struct.weights = weights;
    ret_struct.max_weight = 0;
    ret_struct.all_element_count = 0;

    int empty = 0;
    unsigned int ELEMENTS_IN_RECORD = 2;
    char * words[ELEMENTS_IN_RECORD];
    char * dbKey = new char[32+1];
    unsigned int * distribution = new unsigned int[setSize];
    size_t elementCount = 0;
    for(unsigned int i = 0; i < setSize; i++){
        distribution[i] = randBetween(1, setSize-(setSize/3));
        elementCount += distribution[i];
    }
    unsigned int * elements = new unsigned int[elementCount];
    unsigned short * weight = new unsigned short[elementCount];
    std::fill_n(weight, elementCount, 1);
    size_t curr_start_pos = 0;

    // the reference id of the elements is always their id in the sequence database
    for(unsigned int i = 0; i < setSize; i++) {
        unsigned int element_counter = 0;
        size_t cnt = 0;
        element_buffer[element_counter++] = i;
        element_size[i]++;
        ret_struct.all_element_count++;
        cnt++;
        for(unsigned int j = 1; j < distribution[i]; j++){
            unsigned int curr_element = randBetween(0, elemSize);
            cnt++;
            element_buffer[element_counter++] = curr_element;
            element_size[curr_element]++;
            ret_struct.all_element_count++;
        }

        if (cnt == 0){
            empty++;
        }
        // max_weight can not be gibber than 2^16
        if(element_counter > SHRT_MAX){
            Debug(Debug::ERROR)  << "ERROR: Set has too much elements. Set name is "
                    << dbKey << " and has has the weight " << element_counter <<".\n";
        }
        ret_struct.max_weight = std::max((size_t)element_counter, ret_struct.max_weight);
        memcpy(elements + curr_start_pos, element_buffer, sizeof(unsigned int) * element_counter);
        weights[i] = (weight + curr_start_pos);
        sets[i] = (elements + curr_start_pos);
        set_size[i] = element_counter;
        curr_start_pos += element_counter;
    }
    ret_struct.startElementsArray = elements;
    ret_struct.startWeightsArray = weight;
    if (empty > 0)
        Debug(Debug::WARNING) << empty << " input sets were empty!\n";
    delete [] element_buffer;
    delete [] dbKey;
    delete [] distribution;
    return ret_struct;
}


int main(int argc, char **argv)
{
    // DBReader test
    // argv[1] = ffindex_data_file, argv[2] = ffindex_index_file
    Clustering::set_data set_data;
    std::list<SetCover2::Set *> ret;
    // memory test

    Debug(Debug::INFO) << "Clustering mode: SET COVER\n";
    Debug(Debug::INFO) << "Reading the data...\n";
    for(unsigned int i = 1; i < 10000; i++){
        set_data = create_random_data(1000, 1000);

//        Debug(Debug::INFO) << "\nInit set cover...\n";
        SetCover2 setcover2(set_data.set_count, (const unsigned int *) set_data.set_sizes,
                            set_data.sets, (const unsigned short **) set_data.weights);
//        for(size_t i = 0; i < set_data.set_count; i++) {
//            std::cout << i << " ( ";
//
//            for (size_t j = 0; j < set_data.set_sizes[i]; j++) {
//                std::cout << set_data.sets[i][j] << " ";
//
//            }
//            std::cout << ")" << std::endl;
//
//        }
//        Debug(Debug::INFO) << "Adding sets...\n";
//        for(size_t i = 0; i < set_data.set_count; i++){
//            setcover.add_set(i+1, set_data.set_sizes[i],
//                    (const unsigned int*)   set_data.sets[i],
//                    (const unsigned short*) set_data.weights[i],
//                    set_data.set_sizes[i]);
//        }

//        Debug(Debug::WARNING) << "Calculating the clustering.\n";
        ret = setcover2.execute_set_cover();
//        Debug(Debug::INFO) << "done.\n";
//
//        Debug(Debug::INFO) << "Writing results...\n";
        //writeData(ret);
        std::cout << ret.size() << std::endl;
        std::list<SetCover2::Set *>::const_iterator iterator;
//        for (iterator = ret.begin(); iterator != ret.end(); ++iterator) {
//            std::cout << "Set: " << (*iterator)->setId << " ( ";
//
//            // first entry is the representative sequence
//
//            for(size_t i = 0; i < set_data.set_sizes[ (*iterator)->setId ]; i++){
//                if((*iterator)->elements[i] != UINT_MAX)
//                std::cout << (*iterator)->elements[i] <<  ", " ;
//            }
//            std::cout << " )" << std::endl;
//
//        }
//        Debug(Debug::INFO) << "...done.\n";
        // DBWriter test
        delete[] set_data.startWeightsArray;
        delete[] set_data.startElementsArray;
        delete[] set_data.weights;
        delete[] set_data.sets;
        delete[] set_data.set_sizes;
        delete[] set_data.element_size_lookup;
    }
}
