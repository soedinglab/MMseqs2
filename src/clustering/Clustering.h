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
#include <cstdlib> //the standard C library header

#include "LinearMultiArray.h"
#include "SetCover.h"
#include "SimpleClustering.h"
#include "AffinityClustering.h"
#include "../commons/DBReader.h"
#include "../commons/DBWriter.h"
#include "../commons/Log.h"
#include "../commons/Debug.h"
#include "Parameters.h"


class Clustering {

    public:

        Clustering (std::string seqDB, std::string seqDBIndex,
                std::string alnResultsDB, std::string alnResultsDBIndex,
                std::string outDB, std::string outDBIndex, 
                int validateClustering, float seqId, int maxListLen, unsigned int maxIteration,
                unsigned int convergenceIterations,float dampingFactor,
                int similarityScoreType, float preference, int threads);

        struct set_data {
            // one set contains pointers to the cluster member ids
            unsigned int ** sets;
            float ** similarities;
            unsigned short ** weights;
            unsigned int * set_sizes;
            unsigned int * element_size_lookup;
            size_t set_count;
            size_t unique_element_count;
            size_t all_element_count;
            size_t max_weight;
            unsigned int * startElementsArray;
            unsigned short * startWeightsArray;
            std::list<int>* validids;
        };

        void run(int mode);



    private:
        // check if every element is member in only one cluster
        bool validate_result(std::list<set *> * ret,unsigned int uniqu_element_count);

        // read data for set cover
        set_data read_in_set_data(int mode);


        void writeData(DBWriter *dbw, std::list<set *> ret);

        DBReader* seqDbr;
        DBReader* alnDbr;

        int validate;

        // maximum length of the alignment lists 
        // after the maximum length, the reading of the alignment list for a query is aborted
        unsigned int maxListLen;

    //values for affinity clustering
        unsigned int maxIteration;
        unsigned int convergenceIterations;
        float dampingFactor;
        int similarityScoreType;
        float preference;
        // sequence id. thr.
        float seqIdThr;
        int threads;
        std::string outDBIndex;
        std::string outDB;
};
#endif
