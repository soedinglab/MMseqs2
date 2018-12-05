#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include <sstream>

int touchdb(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 1);
    DBReader<unsigned int> dbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    dbr.open(DBReader<unsigned int>::NOSORT);
    dbr.readMmapedDataInMemory();
    dbr.close();
    return EXIT_SUCCESS;
}
