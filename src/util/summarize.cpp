#include "Parameters.h"
#include "Util.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "HeaderSummarizer.h"

#ifdef OPENMP
#include <omp.h>
#endif

int summarize(int argc, const char** argv) {
    std::string usage("Summarizes all the headers in a fasta entry and prepend the summary to the fasta entry.\n");
    usage.append("Written by Milot Mirdita (milot@mirdita.de)\n");
    usage.append("USAGE: prefilter <inDB> <outDB>\n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.onlyverbosity, 2);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<std::string> reader(par.db2.c_str(), par.db2Index.c_str());
    reader.open(DBReader<std::string>::NOSORT);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), par.threads);
    writer.open();

    UniprotHeaderSummarizer summarizer;

    Debug(Debug::INFO) << "Start writing to file " << par.db3 << "\n";

    #pragma omp for schedule(dynamic, 100)
    for (int i = 0; i < reader.getSize(); ++i) {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        std::string id = reader.getDbKey(i);
        char* data = reader.getData(i);

        std::vector<std::string> headers;

        std::istringstream inStream(data);
        std::string line;
        while (std::getline(inStream, line))
        {
            if(line[0] == '>') {
                headers.emplace_back(line);
            }
        }

        std::string summary = summarizer.summarize(headers, par.summaryPrefix);

        std::ostringstream outStream;
        outStream << summary << data;
        std::string out = outStream.str();
        writer.write(out.c_str(), out.length(), id.c_str(), thread_idx);
    }
    writer.close();
    reader.close();

    return EXIT_SUCCESS;
}
