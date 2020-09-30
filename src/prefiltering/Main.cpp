#include "Prefiltering.h"
#include "Util.h"
#include "Parameters.h"
#include "MMseqsMPI.h"
#include "DBReader.h"
#include "Timer.h"
#include "FileUtil.h"

#ifdef OPENMP
#include <omp.h>
#endif

int prefilter(int argc, const char **argv, const Command& command) {
    MMseqsMPI::init(argc, argv);

    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_PREFILTER);

    Timer timer;
    int queryDbType = FileUtil::parseDbType(par.db1.c_str());
    int targetDbType = FileUtil::parseDbType(par.db2.c_str());
    if(Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_INDEX_DB) == true) {
        DBReader<unsigned int> dbr(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
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

    Prefiltering pref(par.db1, par.db1Index, par.db2, par.db2Index, queryDbType, targetDbType, par);

#ifdef HAVE_MPI
    int runRandomId = 0;
    if (par.localTmp != "") {
        std::srand(std::time(nullptr)); // use current time as seed for random generator
        runRandomId = std::rand();
        runRandomId = runRandomId / 2; // to avoid the unlikely case of overflowing later
    }
    pref.runMpiSplits(par.db3, par.db3Index, par.localTmp, runRandomId);
#else
    pref.runAllSplits(par.db3, par.db3Index);
#endif

    return EXIT_SUCCESS;
}
