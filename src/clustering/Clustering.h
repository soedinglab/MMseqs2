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
#include "ClusteringAlgorithms.h"
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
                int validateClustering,  unsigned int maxIteration,
                int similarityScoreType, int threads);

               void run(int mode);



    private:
        // check if every element is member in only one cluster
        bool validate_result(std::list<set *> * ret,unsigned int uniqu_element_count);



        void writeData(DBWriter *dbw, std::list<set *> ret);

        DBReader<unsigned int>* seqDbr;
        DBReader<unsigned int>* alnDbr;

        int validate;

        //values for affinity clustering
        unsigned int maxIteration;
        int similarityScoreType;
        int threads;
        std::string outDBIndex;
        std::string outDB;
};
#endif
