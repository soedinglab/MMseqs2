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
    MMseqsMPI::init(argc, argv);

    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 3, true, 0, MMseqsParameter::COMMAND_PREFILTER);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    struct timeval start, end;
    gettimeofday(&start, NULL);

    Debug(Debug::INFO) << "Initialising data structures...\n";

    int queryDbType = DBReader<unsigned int>::parseDbType(par.db1.c_str());
    int targetDbType = DBReader<unsigned int>::parseDbType(par.db2.c_str());
    if(queryDbType == -1 || targetDbType == -1){
        Debug(Debug::ERROR) << "Please recreate your database or add a .dbtype file to your sequence/profile database.\n";
        EXIT(EXIT_FAILURE);
    }
//    if(queryDbType == DBReader<unsigned int>::DBTYPE_NUC || targetDbType == DBReader<unsigned int>::DBTYPE_NUC){
//        Debug(Debug::ERROR) << "The prefilter does not support nucleotide sequences.\n";
//        EXIT(EXIT_FAILURE);
//    }
    if(queryDbType == DBReader<unsigned int>::DBTYPE_PROFILE && targetDbType == DBReader<unsigned int>::DBTYPE_PROFILE ){
        Debug(Debug::ERROR) << "Only the query OR the target database can be a profile database.\n";
        EXIT(EXIT_FAILURE);
    }
    if(queryDbType != DBReader<unsigned int>::DBTYPE_PROFILE && targetDbType == DBReader<unsigned int>::DBTYPE_PROFILE_STATE ){
        Debug(Debug::ERROR) << "The query has to be a profile when using a target profile state database.\n";
        EXIT(EXIT_FAILURE);
    }else if (queryDbType == DBReader<unsigned int>::DBTYPE_PROFILE && targetDbType == DBReader<unsigned int>::DBTYPE_PROFILE_STATE){
        queryDbType = Sequence::PROFILE_STATE_PROFILE;
    }

    Prefiltering pref(par.db2, par.db2Index, queryDbType, targetDbType, par);
    gettimeofday(&end, NULL);
    time_t sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "Time for init: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n\n";

#ifdef HAVE_MPI
    pref.runMpiSplits(par.db1, par.db1Index, par.db3, par.db3Index);
#else
    pref.runAllSplits(par.db1, par.db1Index, par.db3, par.db3Index);
#endif

    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nOverall time for prefiltering run: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";

    EXIT(EXIT_SUCCESS);
}
