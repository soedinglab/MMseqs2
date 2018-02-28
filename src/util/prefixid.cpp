#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "Debug.h"

#include <fstream>

#ifdef OPENMP
#include <omp.h>
#endif


int addid(const std::string &db1, const std::string &db1Index, const std::string &db2, const std::string &db2Index, 
const bool tsvOut, const std::string &mappingFile, const std::string &userStrToAdd, const bool isPrefix, const int threads) {
#ifdef OPENMP
    omp_set_num_threads(threads);
#endif

    DBReader<unsigned int> reader(db1.c_str(), db1Index.c_str());
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(db2.c_str(), db2Index.c_str(), threads);
    writer.open();
    bool shouldWriteNullByte = !tsvOut;

    size_t entries = reader.getSize();

    Debug(Debug::INFO) << "Start adding to database.\n";

    std::map<unsigned int, std::string> mapping = Util::readLookup(mappingFile);

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
                std::string strToAdd = "";
                if (userStrToAdd != "") {
                    strToAdd = userStrToAdd;
                } else if (mappingFile.length() > 0) {
                    strToAdd = mapping[key];
                } else {
                    std::ostringstream sstmp;
                    sstmp << key;
                    strToAdd = sstmp.str();
                }

                if (isPrefix) {
                    ss << strToAdd << "\t" << line << "\n";
                } else {
                    ss << line << "\t" << strToAdd << "\n";
                }
            }

            std::string result = ss.str();
            writer.writeData(result.c_str(), result.length(), key, thread_idx, shouldWriteNullByte);
        }
    }
    writer.close();

    Debug(Debug::INFO) << "\nDone.\n";

    reader.close();

    return EXIT_SUCCESS;
}

int prefixid(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);
    return(addid(par.db1, par.db1Index, par.db2, par.db2Index, par.tsvOut, par.mappingFile, par.prefix, true, par.threads));
}

int suffixid(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);
    return(addid(par.db1, par.db1Index, par.db2, par.db2Index, par.tsvOut, par.mappingFile, par.prefix, false, par.threads));
}

