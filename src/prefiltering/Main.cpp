
#include "Prefiltering.h"
#include "Util.h"
#include "Parameters.h"
#include "MMseqsMPI.h"
#include "Timer.h"

#include <iostream>
#include <string>

#ifdef OPENMP
#include <omp.h>
#endif

int prefilter(int argc, const char **argv, const Command& command) {
    MMseqsMPI::init(argc, argv);

    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 3, true, 0, MMseqsParameter::COMMAND_PREFILTER);

    Timer timer;

    int queryDbType = DBReader<unsigned int>::parseDbType(par.db1.c_str());

    std::string indexStr = PrefilteringIndexReader::searchForIndex(par.db2);
    int targetDbType = DBReader<unsigned int>::parseDbType(par.db2.c_str());
    if(Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_INDEX_DB) == true){
        DBReader<unsigned int> dbr(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        dbr.open(DBReader<unsigned int>::NOSORT);
        PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(&dbr);
        targetDbType = data.seqType;
        dbr.close();
    }
    if (queryDbType == -1 || targetDbType == -1) {
        Debug(Debug::ERROR) << "Please recreate your database or add a .dbtype file to your sequence/profile database.\n";
        return EXIT_FAILURE;
    }
    if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) && Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_HMM_PROFILE)) {
        Debug(Debug::ERROR) << "Only the query OR the target database can be a profile database.\n";
        return EXIT_FAILURE;
    }

    if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_AMINO_ACIDS) && Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_NUCLEOTIDES)) {
        Debug(Debug::ERROR) << "The prefilter can not search amino acids against nucleotides. Something might got wrong while createdb or createindex.\n";
        return EXIT_FAILURE;
    }
    if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES) && Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_AMINO_ACIDS)) {
        Debug(Debug::ERROR) << "The prefilter can not search nucleotides against amino acids. Something might got wrong while createdb or createindex.\n";
        return EXIT_FAILURE;
    }
    if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) == false && Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_PROFILE_STATE_SEQ)) {
        Debug(Debug::ERROR) << "The query has to be a profile when using a target profile state database.\n";
        return EXIT_FAILURE;
    } else if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) && Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_PROFILE_STATE_SEQ)) {
        queryDbType = Parameters::DBTYPE_PROFILE_STATE_PROFILE;
    }
    Prefiltering pref(par.db2, par.db2Index, queryDbType, targetDbType, par);
    //Debug(Debug::INFO) << "Time for init: " << timer.lap() << "\n";

#ifdef HAVE_MPI
    pref.runMpiSplits(par.db1, par.db1Index, par.db3, par.db3Index, par.localTmp);
#else
    pref.runAllSplits(par.db1, par.db1Index, par.db3, par.db3Index);
#endif

    return EXIT_SUCCESS;
}
