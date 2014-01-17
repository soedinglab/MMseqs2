#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <iostream>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <istream>
#include <cstring>
#include <stdlib.h>

#include "LinearMultiArray.h"
#include "SetCover.h"
#include "SimpleClustering.h"
#include "../commons/DBReader.h"
#include "../commons/DBWriter.h"

class Clustering {

    public:

        Clustering (std::string seqDBBase, std::string alnResultsDBBase, std::string outDBBAse);

        struct set_data {
            // one set contains pointers to the cluster member ids
            unsigned int ** sets;
            unsigned short ** weights;
            unsigned int * set_sizes;
            unsigned int * element_size_lookup;
            unsigned int set_count;
            unsigned int uniqu_element_count;
            unsigned int all_element_count;
            unsigned int max_weight;
        };

        void run(int mode);

        static const int SET_COVER = 0;
        static const int GREEDY = 1;

    private:
        // check if every element is member in only one cluster
        bool validate_result(std::list<set *> * ret,unsigned int uniqu_element_count);

        // read data for set cover
        set_data read_in_set_data_set_cover();

        // read data for the clustering by the sequence length
        set_data read_in_set_data_simple_clustering();

        void writeData(std::list<set *> ret);

        int mode;

        DBReader* seqDbr;

        DBReader* alnDbr;

        DBWriter* dbw;
};
#endif
