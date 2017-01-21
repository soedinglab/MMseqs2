#include <iostream>
#include <string>
#include <sys/time.h>

#include "Prefiltering.h"
#include "Util.h"
#include "Parameters.h"

#include "MMseqsMPI.h"
#ifdef OPENMP
#include <omp.h>
#endif

int prefilter(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 3, true, false, MMseqsParameter::COMMAND_PREFILTER);

    MMseqsMPI::init(argc, argv);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    struct timeval start, end;
    gettimeofday(&start, NULL);

    Prefiltering pref(par.db2, par.db2Index, par);
    Debug(Debug::INFO) << "Initialising data structures...\n";
    gettimeofday(&end, NULL);
    time_t sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "Time for init: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n\n";
    gettimeofday(&start, NULL);
#ifdef HAVE_MPI
    std::pair<std::string, std::string> filenamePair = Util::createTmpFileNames(par.db3, par.db3Index, MMseqsMPI::rank);

    //TODO check again :(
    if(pref.getSplit() > MMseqsMPI::numProc){
        // if split size is great than nodes than we have to
        // distribute all splits equally over all nodes
        unsigned int * splitCntPerProc = new unsigned int[MMseqsMPI::numProc];
        memset(splitCntPerProc, 0, sizeof(unsigned int) * MMseqsMPI::numProc);
        for(int i = 0; i < pref.getSplit(); i++){
            splitCntPerProc[i % MMseqsMPI::numProc] += 1;
        }
        int fromSplit = 0;
        for(int i = 0; i < MMseqsMPI::rank; i++){
            fromSplit += splitCntPerProc[i];
        }
        pref.runSplits(par.db1,par.db1Index, filenamePair.first.c_str(), filenamePair.second.c_str(), fromSplit, splitCntPerProc[MMseqsMPI::rank]);
        delete[] splitCntPerProc;
    } else {
        // if more nodes exists than splits are needed set split to the amount of nodes
        // each node needs to compute 1 target split
        // OR
        // if database fits into the memory of a single node split by query.
        // each node should just compute 1 query split
        pref.split = MMseqsMPI::numProc;
        pref.runSplits(par.db1,par.db1Index, filenamePair.first.c_str(), filenamePair.second.c_str(), MMseqsMPI::rank, 1);
    }
#else
    pref.runAllSplits(par.db1, par.db1Index, par.db3, par.db3Index);
#endif

#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    if(MMseqsMPI::rank == 0){
        std::vector<std::pair<std::string, std::string> > splitFiles;
        for(int procs = 0; procs < MMseqsMPI::numProc; procs++){
            splitFiles.push_back(Util::createTmpFileNames(par.db3, par.db3Index, procs));
        }
        // merge output ffindex databases
        pref.mergeFiles(splitFiles, par.db3, par.db3Index);
    }
#endif
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nOverall time for prefiltering run: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
