//
// Created by mad on 1/27/16.
//


#include <Parameters.h>
#include <DBReader.h>
#include <DBWriter.h>
#include <Util.h>
#include <Log.h>
#ifdef OPENMP
#include <omp.h>
#endif

int dofilter(std::string inputDb, std::string outputDb, int threads, int column, std::string regexStr) {
    DBReader<std::string>* dataDb=new DBReader<std::string>(inputDb.c_str(),(std::string(inputDb).append(".index")).c_str());
    dataDb->open(DBReader<std::string>::LINEAR_ACCCESS);
    DBWriter* dbw = new DBWriter(outputDb.c_str(), (std::string(outputDb).append(".index")).c_str(), threads);
    dbw->open();

    regex_t regex;
    int status = regcomp(&regex, regexStr.c_str(), REG_EXTENDED | REG_NEWLINE);
    if (status != 0 ){
        Debug(Debug::INFO) << "Error in regex " << regexStr << "\n";
        EXIT(EXIT_FAILURE);
    }
#pragma omp parallel
    {
        size_t LINE_BUFFER_SIZE = 1000000;
        char *lineBuffer = new char[LINE_BUFFER_SIZE];
        char *columnValue = new char[LINE_BUFFER_SIZE];
        char **columnPointer = new char*[column + 1];
        std::string buffer = "";
        buffer.reserve(LINE_BUFFER_SIZE);
#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < dataDb->getSize(); id++) {
            Log::printProgress(id);
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            char *data = dataDb->getData(id);

            while (*data != '\0') {
                Util::getLine(data, lineBuffer);
                size_t foundElements = Util::getWordsOfLine(lineBuffer, columnPointer, column + 1);
                if(foundElements < column  ){
                    Debug(Debug::ERROR) << "Column=" << column << " does not exist in line " << lineBuffer << "\n";
                    EXIT(EXIT_FAILURE);
                }
                // if column is last column
                if(column == foundElements){
                    const ptrdiff_t entrySize = Util::skipLine(data) - columnPointer[column];
                    memcpy(columnValue, columnPointer[column], entrySize);
                    columnValue[entrySize] = '\0';
                }else{
                    const ptrdiff_t entrySize = columnPointer[column + 1] - columnPointer[column];
                    memcpy(columnValue, columnPointer[column], entrySize);
                    columnValue[entrySize] = '\0';
                }
                int nomatch = regexec(&regex, columnValue, 0, NULL, 0);

                if(!(nomatch)){
                    buffer.append(lineBuffer);
                }
                data = Util::skipLine(data);
            }

            dbw->write(buffer.c_str(), buffer.length(), (char*) dataDb->getDbKey(id).c_str(), thread_idx);
            buffer.clear();
        }
        delete [] lineBuffer;
        delete [] columnValue;
        delete [] columnPointer;
    }
    regfree(&regex);
    dataDb->close();
    dbw->close();
    delete dataDb;
    delete dbw;
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
