#include "Parameters.h"
#include "Util.h"
#include "PrefilteringIndexReader.h"
#include "MemoryMapped.h"

int touchdb(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string db = par.db1;

    std::string indexDB = PrefilteringIndexReader::searchForIndex(db);
    if (indexDB.empty() == false) {
        db = indexDB;
    }

    MemoryMapped map(db, MemoryMapped::WholeFile, MemoryMapped::CacheHint::SequentialScan);
    Util::touchMemory(reinterpret_cast<const char*>(map.getData()), map.mappedSize());

    return EXIT_SUCCESS;
}
