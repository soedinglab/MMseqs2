#include "Debug.h"
#include "Parameters.h"
#include "FileUtil.h"
#include "MappingReader.h"

int createbintaxmapping(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    MappingReader reader(par.db1, false);
    std::pair<char*, size_t> serialized = MappingReader::serialize(reader);
    FILE* handle = fopen(par.db2.c_str(), "w");
    if (handle == NULL) {
        Debug(Debug::ERROR) << "Could not open " << par.db2 << " for writing\n";
        return EXIT_FAILURE;
    }
    size_t written = fwrite(serialized.first, serialized.second * sizeof(char), 1, handle);
    free(serialized.first);
    if (written != 1) {
        Debug(Debug::ERROR) << "Could not write to " << par.db2 << "\n";
        return EXIT_FAILURE;
    }
    if (fclose(handle) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << par.db2 << "\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
