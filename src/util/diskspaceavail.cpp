#include "Command.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "Parameters.h"

extern const char* version;

int diskspaceavail(int, const char**, const Command&) {
    Parameters &par = Parameters::getInstance();
    size_t diskLimit = FileUtil::getFreeSpace(FileUtil::dirName(par.db1).c_str());
    Debug(Debug::INFO) << diskLimit << "\n"; // in bytes
    EXIT(EXIT_SUCCESS);
}
