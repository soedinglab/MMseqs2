#include "Parameters.h"
#include "Util.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "HeaderSummarizer.h"

#ifdef OPENMP
#include <omp.h>
#endif

int summarizeheaders(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<std::string> queryReader(par.db1.c_str(), par.db1Index.c_str());
    queryReader.open(DBReader<std::string>::NOSORT);

    DBReader<std::string> targetReader(par.db2.c_str(), par.db2Index.c_str());
    targetReader.open(DBReader<std::string>::NOSORT);

    DBReader<std::string> reader(par.db3.c_str(), par.db3Index.c_str());
    reader.open(DBReader<std::string>::NOSORT);

    DBWriter writer(par.db4.c_str(), par.db4Index.c_str(), par.threads);
    writer.open();

    UniprotHeaderSummarizer summarizer;

    Debug(Debug::INFO) << "Start writing to file " << par.db4 << "\n";

    #pragma omp for schedule(dynamic, 100)
    for (size_t i = 0; i < reader.getSize(); ++i) {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        std::string id = reader.getDbKey(i);
        char* data = reader.getData(i);

        std::vector<std::string> headers;

        std::istringstream inStream(data);
        std::string line;
        size_t entry = 0;
        std::string representative;
        while (std::getline(inStream, line))
        {
            char* header;
            if(entry == 0) {
                header = queryReader.getDataByDBKey(line);
                representative = line;
            } else {
                header = targetReader.getDataByDBKey(line);
            }

            headers.emplace_back(header);
            entry++;
        }

        std::ostringstream oss;
        oss << par.summaryPrefix << "-" << representative << "|" << summarizer.summarize(headers);

        std::string summary = oss.str();
        writer.writeData(summary.c_str(), summary.length(), id.c_str(), thread_idx);
    }
    writer.close();
    reader.close();
    targetReader.close();
    queryReader.close();

    return EXIT_SUCCESS;
}
