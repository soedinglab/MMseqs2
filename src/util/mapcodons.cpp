
#include <Parameters.h>
#include <DBReader.h>
#include <DBWriter.h>
#include <Util.h>
#include <Log.h>

#ifdef OPENMP
#include <omp.h>
#endif

int mapcodons(std::string inputDb, std::string outputDb, int threads, int column, std::string regexStr) {
    DBReader<std::string> reader(DBReader<std::string>(inputDb.c_str(), (std::string(inputDb).append(".index")).c_str()));
    reader.open(DBReader<std::string>::LINEAR_ACCCESS);

    DBWriter writer(outputDb.c_str(), (std::string(outputDb).append(".index")).c_str(), threads);
    writer.open();

    const size_t LINE_BUFFER_SIZE = 1000000;

#pragma omp parallel
    {
        char *lineBuffer = new char[LINE_BUFFER_SIZE];
        char *columnValue = new char[LINE_BUFFER_SIZE];
        char **columnPointer = new char*[column + 1];
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

            while (*data != '\0') {
                //Util::getLine(data, lineBuffer);
                size_t foundElements = Util::getWordsOfLine(lineBuffer, columnPointer, column + 1);

                data = Util::skipLine(data);
            }

            writer.write(buffer.c_str(), buffer.length(), (char*) reader.getDbKey(id).c_str(), thread_idx);
            buffer.clear();
        }
        delete [] lineBuffer;
        delete [] columnValue;
        delete [] columnPointer;
    }

    reader.close();
    writer.close();
    return 0;
}

int filterdb(int argn, const char **argv)
{
    std::string usage;
    usage.append("Filter a database by column regex\n");
    usage.append("USAGE:  <ffindexDB> <outDB>\n");
    usage.append("\nDesigned and implemented by Martin Steinegger <martin.steinegger@mpibpc.mpg.de>.\n");

    Parameters par;
    par.parseParameters(argn, argv, usage, par.filterDb, 2);
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    return dofilter(par.db1,
                    par.db2,
                    par.threads,
                    par.filterColumn,
                    par.filterColumnRegex);
}
