#include "Util.h"
#include "Parameters.h"
#include "DBReader.h"


int mvdb(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);
    DBReader<unsigned int>::moveDb(par.db1.c_str(), par.db2.c_str());
    return EXIT_SUCCESS;
}
