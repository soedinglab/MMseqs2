#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "FileUtil.h"
#include "Debug.h"
#ifdef HAVE_ZLIB
#include <zlib.h>
#endif

#ifdef OPENMP
#include <omp.h>
#endif

int unpackdb(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string lookupFile = par.db1 + ".lookup";
    if (par.unpackNameMode == Parameters::UNPACK_NAME_ACCESSION && FileUtil::fileExists(lookupFile.c_str()) == false) {
        Debug(Debug::INFO) << "No lookup file for " << FileUtil::baseName(par.db1) << " found, using key-based file naming\n";
        par.unpackNameMode = Parameters::UNPACK_NAME_KEY;
    }

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

    size_t localThreads = 1;
#ifdef OPENMP
    localThreads = std::max(std::min((size_t)par.threads, reader.getSize()), (size_t)1);
#endif

    size_t entries = reader.getSize();
    Debug::Progress progress(entries);
#pragma omp parallel num_threads(localThreads)
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

            const char* cname = name.c_str();

            if (FileUtil::fileExists(cname) == true) {
                if (FileUtil::directoryExists(cname) == true) {
                    Debug(Debug::ERROR) << "Cannot open directory " << name << " for writing\n";
                    continue;
                }
                FileUtil::remove(cname);
            }

            if (Util::endsWith(".gz", name.c_str()) == true) {
#ifdef HAVE_ZLIB
                gzFile handle = gzopen(cname, "w");
                if (handle == NULL) {
                    Debug(Debug::ERROR) << "Cannot not open " << name << " for writing\n";
                    continue;
                }
                size_t len = reader.getEntryLen(i) - 1;
                int n = gzwrite(handle ,reader.getData(i, thread_idx), len * sizeof(char));
                if ((size_t)n != len) {
                    Debug(Debug::ERROR) << "Cannot not write " << name << "\n";
                    continue;
                }
                if (gzclose(handle) != 0) {
                    Debug(Debug::ERROR) << "Cannot not close " << name << "\n";
                    continue;
                }
#else
                Debug(Debug::ERROR) << "MMseqs2 was not compiled with zlib support. Cannot write compressed output\n";
                EXIT(EXIT_FAILURE);
#endif
            } else {
                FILE* handle = fopen(cname, "w");
                if (handle == NULL) {
                    Debug(Debug::ERROR) << "Cannot not open " << name << " for writing\n";
                    continue;
                }
                size_t len = reader.getEntryLen(i) - 1;
                int n = fwrite(reader.getData(i, thread_idx), sizeof(char), len, handle);
                if ((size_t)n != len) {
                    Debug(Debug::ERROR) << "Cannot not write " << name << "\n";
                    continue;
                }
                if (fclose(handle) != 0) {
                    Debug(Debug::ERROR) << "Cannot not close " << name << "\n";
                    continue;
                }
            }
        }
    }
    reader.close();
    return EXIT_SUCCESS;
}
