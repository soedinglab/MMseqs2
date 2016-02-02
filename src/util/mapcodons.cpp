#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "Log.h"
#include "Debug.h"

extern "C" {
#include "ext/fmemopen.h"
}

#include <unistd.h>
#include "kseq.h"
KSEQ_INIT(int, read)

#ifdef OPENMP
#include <omp.h>
#endif


int mapcodons(int argn, const char **argv)
{
    std::string usage;
    usage.append("Maps back the original nucleotides to a MSA produced by MMseqs\n");
    usage.append("USAGE: <msaDB> <ncDB> <outDB>\n");
    usage.append("\nDesigned and implemented by Milot Mirdita <milot@mirdita.de>.\n");

    Parameters par;
    par.parseParameters(argn, argv, usage, par.onlyverbosity, 4);
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<std::string> reader(DBReader<std::string>(par.db1.c_str(), par.db1Index.c_str()));
    reader.open(DBReader<std::string>::LINEAR_ACCCESS);

    DBReader<std::string> ncReader(DBReader<std::string>(par.db3.c_str(), par.db3Index.c_str()));
    ncReader.open(DBReader<std::string>::NOSORT);

    //DBWriter writer(par.db4.c_str(), par.db4Index.c_str(), par.threads);
    //writer.open();

    const size_t LINE_BUFFER_SIZE = 1000000;

#pragma omp parallel
    {
        std::string buffer = "";
        buffer.reserve(LINE_BUFFER_SIZE);


#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < reader.getSize(); id++) {
            Log::printProgress(id);
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif

            char *data = reader.getData(id);
            size_t dataLength = reader.getSeqLens(id);

            FILE* file = fmemopen(data, dataLength, "r");
            kseq_t *seq = kseq_init(fileno(file));
            while (kseq_read(seq) >= 0) {
                Debug(Debug::INFO) << seq->comment.s << "\n";
                //char * ncSeq = ncReader.getData(id);
                //size_t ncSeqLength = ncReader.getSeqLens(id);
            }
            kseq_destroy(seq);
            fclose(file);

            //writer.write(buffer.c_str(), buffer.length(), (char*) reader.getDbKey(id).c_str(), thread_idx);
        }
    }

    ncReader.close();
    reader.close();
    //writer.close();
    return EXIT_SUCCESS;
}
