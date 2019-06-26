#include "Debug.h"
#include "DBReader.h"
#include "Parameters.h"
#include "Util.h"
#include "FileUtil.h"

int dbtype(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, false, 0, 0);
    Debug(Debug::INFO) << Parameters::getDbTypeName(FileUtil::parseDbType(par.db1.c_str()));
    EXIT(EXIT_SUCCESS);
}
