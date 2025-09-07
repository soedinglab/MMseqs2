#include "FileUtil.h"
#include "Parameters.h"
#include "DBReader.h"

int rmdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    DBReader<IdType>::removeDb(par.db1);
    return EXIT_SUCCESS;
}

int mvdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    DBReader<IdType>::moveDb(par.db1.c_str(), par.db2.c_str());
    return EXIT_SUCCESS;
}

int cpdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    DBReader<IdType>::copyDb(par.db1.c_str(), par.db2.c_str());
    return EXIT_SUCCESS;
}

int lndb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    DBReader<IdType>::softlinkDb(par.db1.c_str(), par.db2.c_str());
    return EXIT_SUCCESS;
}

int aliasdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    std::string alias = FileUtil::baseName(par.db2.c_str());
    DBReader<IdType>::aliasDb(par.db1.c_str(), alias);
    return EXIT_SUCCESS;
}
