#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "Debug.h"

#ifdef OPENMP
#include <omp.h>
#endif


int addid(const std::string &db1, const std::string &db1Index, const std::string &db2, const std::string &db2Index, 
const bool tsvOut, const std::string &mappingFile, const std::string &userStrToAdd, const bool isPrefix, const int threads, const int compressed) {
    DBReader<unsigned int> reader(db1.c_str(), db1Index.c_str(), threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    const bool shouldCompress = tsvOut == false && compressed == true;
    // TODO: does generic db make more sense than copying db type here?
    const int dbType = tsvOut == true ? Parameters::DBTYPE_OMIT_FILE : reader.getDbtype();
    DBWriter writer(db2.c_str(), db2Index.c_str(), threads, shouldCompress, dbType);
    writer.open();
    const bool shouldWriteNullByte = !tsvOut;

    size_t entries = reader.getSize();
    Debug::Progress progress(entries);

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
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            std::istringstream data(reader.getData(i, thread_idx));
            std::ostringstream ss;

            std::string line;
            while (std::getline(data, line)) {
                std::string strToAdd = "";
                if (userStrToAdd != "") {
                    strToAdd = userStrToAdd;
                } else if (mappingFile.length() > 0) {
                    strToAdd = mapping[key];
                } else {
                    strToAdd = SSTR(key);
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
    writer.close(tsvOut);
    reader.close();

    return EXIT_SUCCESS;
}

int prefixid(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);
    return(addid(par.db1, par.db1Index, par.db2, par.db2Index, par.tsvOut, par.mappingFile, par.prefix, true, par.threads, par.compressed));
}

int suffixid(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);
    return(addid(par.db1, par.db1Index, par.db2, par.db2Index, par.tsvOut, par.mappingFile, par.prefix, false, par.threads, par.compressed));
}

