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
#include <time.h>
#include <sys/time.h>

#include "LinearMultiArray.h"
#include "SetCover.h"
#include "SimpleClustering.h"
#include "../commons/DBReader.h"
#include "../commons/DBWriter.h"
#include "../commons/Log.h"
#include "../commons/Debug.h"

class Clustering {

    public:

        Clustering (std::string seqDB, std::string seqDBIndex,
                std::string alnResultsDB, std::string alnResultsDBIndex,
                std::string outDB, std::string outDBIndex, 
                float seqIdThr, int validateClustering, int maxListLen);

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
        set_data read_in_set_data();

        void writeData(std::list<set *> ret);

        int mode;

        DBReader* seqDbr;

        DBReader* alnDbr;

        DBWriter* dbw;

        float seqIdThr;

        int validate;

        // maximum length of the alignment lists 
        // after the maximum length, the reading of the alignment list for a query is aborted
        int maxListLen;
};
#endif
