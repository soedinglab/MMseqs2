#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "Debug.h"

#include <fstream>

#ifdef OPENMP
#include <omp.h>
#endif

int mapresult(int argn, const char **argv) {
    std::string usage;
    usage.append("Maps each line of the input database with the lookup file.\n");
    usage.append("USAGE: <ffindexInDB> <mappingFile> <ffindexOutDB>\n");
    usage.append("\nDesigned and implemented by Milot Mirdita <milot@mirdita.de>.\n");

    Parameters par;
    par.parseParameters(argn, argv, usage, par.onlyverbosity, 3);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<std::string> reader(par.db1.c_str(), par.db1Index.c_str());
    reader.open(DBReader<std::string>::LINEAR_ACCCESS);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), par.threads);
    writer.open();

    size_t entries = reader.getSize();

    std::map<unsigned int, std::string> mapping = Util::readLookup(par.db2);

#pragma omp for schedule(dynamic, 100)
    for (size_t i = 0; i < entries; ++i) {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        std::string key = reader.getDbKey(i);
        std::istringstream data(reader.getData(i));
        std::ostringstream ss;

        std::string line;
        while (std::getline(data, line)) {
            char* rest;
            unsigned int k = static_cast<unsigned int>(strtoul(line.c_str(), &rest, 10));
            if ((rest != line.c_str() && *rest != '\0') || errno == ERANGE) {
                Debug(Debug::WARNING) << "Invalid entry in line " << line << "!\n";
                continue;
            }

            ss << mapping[k] << "\n";
        }

        std::string result = ss.str();
        writer.write(result.c_str(), result.length(), key.c_str(), thread_idx);
    }

    writer.close();
    reader.close();

    return EXIT_SUCCESS;
}

