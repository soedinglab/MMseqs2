#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "FileUtil.h"
#include "Debug.h"

#ifdef OPENMP
#include <omp.h>
#endif

int unpackdb(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    int mode = DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA;
    if (par.unpackNameMode == Parameters::UNPACK_NAME_ACCESSION) {
        mode |= DBReader<unsigned int>::USE_LOOKUP;
    }
    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, mode);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    if (FileUtil::directoryExists(par.db2.c_str()) == false && FileUtil::makeDir(par.db2.c_str()) == false) {
        Debug(Debug::ERROR) << "Cannot create output folder " << par.db2 << "\n";
        EXIT(EXIT_FAILURE);
    }

    size_t entries = reader.getSize();
    Debug::Progress progress(entries);

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

#pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < entries; ++i) {
            progress.updateProgress();
            unsigned int key = reader.getDbKey(i);
            std::string name = par.db2;
            if (name.back() != '/') {
                name.append(1, '/');
            }
            if (par.unpackNameMode == Parameters::UNPACK_NAME_ACCESSION) {
                size_t lookupId = reader.getLookupIdByKey(key);
                name.append(FileUtil::sanitizeFilename(reader.getLookupEntryName(lookupId)));
            } else {
                name.append(SSTR(key));
            }
            name.append(par.unpackSuffix);
            FILE* handle = FileUtil::openAndDelete(name.c_str(), "w");
            fwrite(reader.getData(i, thread_idx), sizeof(char), reader.getEntryLen(i) - 1, handle);
            fclose(handle);
        }
    }
    reader.close();
    return EXIT_SUCCESS;
}
