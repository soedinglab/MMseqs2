#include "Parameters.h"
#include "Util.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "HeaderSummarizer.h"

#ifdef OPENMP
#include <omp.h>
#endif

int setextendeddbtype(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    int dbtype = FileUtil::parseDbType(par.db1.c_str());
    // check if dbtype uses isCompressed flag
    bool isCompressed = (dbtype & (1 << 31));
    dbtype = DBReader<unsigned int>::setExtendedDbtype(dbtype,  par.extendedDbtype);
    DBWriter::writeDbtypeFile(par.db1.c_str(), dbtype, isCompressed);
    return EXIT_SUCCESS;
}
