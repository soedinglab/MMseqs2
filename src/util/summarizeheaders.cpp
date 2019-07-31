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
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> queryReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    queryReader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> targetReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    targetReader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> reader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::NOSORT);

    DBWriter writer(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    writer.open();

    HeaderSummarizer * summarizer;
    if(par.headerType == Parameters::HEADER_TYPE_METACLUST){
        summarizer = new MetaclustHeaderSummarizer;
    }else if(par.headerType == Parameters::HEADER_TYPE_UNICLUST) {
        summarizer = new UniprotHeaderSummarizer;
    }else {
        Debug(Debug::ERROR) << "Header type is not supported\n";
        return EXIT_FAILURE;
    }
    Debug::Progress progress(reader.getSize());

#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
#pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();

            unsigned int id = reader.getDbKey(i);
            char *data = reader.getData(i, thread_idx);

            std::vector<std::string> headers;

            std::istringstream inStream(data);
            std::string line;
            size_t entry = 0;
            std::string representative;
            while (std::getline(inStream, line)) {
                char *header;
                if (entry == 0) {
                    header = queryReader.getDataByDBKey((unsigned int) strtoul(line.c_str(), NULL, 10), thread_idx);

                    representative = line;
                } else {
                    header = targetReader.getDataByDBKey((unsigned int) strtoul(line.c_str(), NULL, 10), thread_idx);
                }
                headers.emplace_back(header);
                entry++;
            }

            std::ostringstream oss;
            oss << par.summaryPrefix << "-" << representative << "|" << summarizer->summarize(headers);

            std::string summary = oss.str();
            writer.writeData(summary.c_str(), summary.length(), id, thread_idx);
        }
    }
    writer.close();
    reader.close();
    targetReader.close();
    queryReader.close();
    delete summarizer;
    return EXIT_SUCCESS;
}
