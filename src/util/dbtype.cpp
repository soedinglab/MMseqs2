#include "Debug.h"
#include "DBReader.h"
#include "Parameters.h"
#include "Util.h"

int dbtype(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 1, false);
    Debug(Debug::INFO) << DBReader<unsigned int>::getDbTypeName(DBReader<unsigned int>::parseDbType(par.db1.c_str()));
    EXIT(EXIT_SUCCESS);
}
