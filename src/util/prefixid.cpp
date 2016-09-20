#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "Debug.h"

#include <fstream>

#ifdef OPENMP
#include <omp.h>
#endif

int prefixid(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<std::string> reader(par.db1.c_str(), par.db1Index.c_str());
    reader.open(DBReader<std::string>::NOSORT);

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads);
    writer.open();

    size_t entries = reader.getSize();

    std::map<unsigned int, std::string> mapping = Util::readLookup(par.mappingFile);

#pragma omp parallel for schedule(dynamic, 100)
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
            if (par.mappingFile.length() > 0) {
                unsigned int k = static_cast<unsigned int>(strtoul(key.c_str(), NULL, 10));
                ss << mapping[k] << "\t" << line << "\n";
            } else {
                ss << key << "\t" << line << "\n";
            }
        }

        std::string result = ss.str();
        writer.writeData(result.c_str(), result.length(), key.c_str(), thread_idx);
    }

    writer.close();
    reader.close();

    return EXIT_SUCCESS;
}

