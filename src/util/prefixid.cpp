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

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str());
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads);
    writer.open();

    size_t entries = reader.getSize();

    const bool useUserPrefix = par.prefix != "";

    Debug(Debug::INFO) << "Start prefixing database.\n";

    std::map<unsigned int, std::string> mapping = Util::readLookup(par.mappingFile);

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

#pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < entries; ++i) {
            Debug::printProgress(i);

            unsigned int key = reader.getDbKey(i);
            std::istringstream data(reader.getData(i));
            std::ostringstream ss;

            std::string line;
            while (std::getline(data, line)) {
                if (useUserPrefix) {
                    ss << par.prefix << "\t" << line << "\n";
                } else if (par.mappingFile.length() > 0) {
                    ss << mapping[key] << "\t" << line << "\n";
                } else {
                    ss << key << "\t" << line << "\n";
                }
            }

            std::string result = ss.str();
            writer.writeData(result.c_str(), result.length(), key, thread_idx);
        }
    }
    writer.close();

    Debug(Debug::INFO) << "\nDone.\n";

    reader.close();

    return EXIT_SUCCESS;
}

