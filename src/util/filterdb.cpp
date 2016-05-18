#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "Log.h"
#include "Debug.h"
#include "filterdb.h"


#ifdef OPENMP
#include <omp.h>
#endif

int ffindexFilter::initFiles() {
    dataDb=new DBReader<unsigned int>(inDB.c_str(),(std::string(inDB).append(".index")).c_str());
    dataDb->open(DBReader<std::string>::LINEAR_ACCCESS);
	
    dbw = new DBWriter(outDB.c_str(), (std::string(outDB).append(".index")).c_str(), threads);
    dbw->open();
	return 0;
}

ffindexFilter::ffindexFilter(std::string inDB,	std::string outDB, int threads, size_t column, std::string regexStr):
inDB(inDB),outDB(outDB),threads(threads),column(column),regexStr(regexStr) {

	initFiles();
	
	mode = REGEX_FILTERING;
	
    int status = regcomp(&regex, regexStr.c_str(), REG_EXTENDED | REG_NEWLINE);
    if (status != 0 ){
        Debug(Debug::INFO) << "Error in regex " << regexStr << "\n";
        EXIT(EXIT_FAILURE);
    }
}

ffindexFilter::ffindexFilter(std::string inDB, std::string outDB, std::string filterFile, int threads, size_t column,bool positiveFiltering=true):
inDB(inDB),outDB(outDB),threads(threads),column(column), filterFile(filterFile),positiveFiltering(positiveFiltering){

	initFiles();
	
	
	mode = FILE_FILTERING;
	
	// Fill the filter with the data contained in the file
	std::ifstream filterFileStream;
	filterFileStream.open(filterFile);
	std::string line;
	while (std::getline(filterFileStream,line))
	{
		filter.push_back(line);
	}
	
	std::stable_sort(filter.begin(), filter.end(), compareString());
	
}


ffindexFilter::~ffindexFilter() {
	
	if (mode == REGEX_FILTERING)
		regfree(&regex);
		
	dataDb->close();
	dbw->close();
	delete dataDb;
	delete dbw;
}



int ffindexFilter::runFilter(){
	
    const size_t LINE_BUFFER_SIZE = 1000000;
#pragma omp parallel
    {
        char *lineBuffer = new char[LINE_BUFFER_SIZE];
        char *columnValue = new char[LINE_BUFFER_SIZE];
        char **columnPointer = new char*[column + 1];
        std::string buffer = "";
        buffer.reserve(LINE_BUFFER_SIZE);
#pragma omp for schedule(static)
        for (size_t id = 0; id < dataDb->getSize(); id++) {

            Log::printProgress(id);
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            char *data = dataDb->getData(id);
            size_t dataLength = dataDb->getSeqLens(id);

            while (*data != '\0') {
                if(!Util::getLine(data, dataLength, lineBuffer, LINE_BUFFER_SIZE)) {
                    Debug(Debug::WARNING) << "Warning: Identifier was too long and was cut off!\n";
                    continue;
                }
                size_t foundElements = Util::getWordsOfLine(lineBuffer, columnPointer, column + 1);
                if(foundElements < column  ){
                    Debug(Debug::ERROR) << "Column=" << column << " does not exist in line " << lineBuffer << "\n";
                    EXIT(EXIT_FAILURE);
                }
				  
				  size_t colStrLen;
                // if column is last column
                if(column == foundElements){
                    const ptrdiff_t entrySize = Util::skipLine(data) - columnPointer[(column - 1)];
                    memcpy(columnValue, columnPointer[column - 1], entrySize);
                    columnValue[entrySize] = '\0';
					  colStrLen = entrySize;
                }else{
                    const ptrdiff_t entrySize = columnPointer[column] - columnPointer[(column - 1)];
                    memcpy(columnValue, columnPointer[column - 1], entrySize);
                    columnValue[entrySize] = '\0';
					  colStrLen = entrySize;
                }
					int nomatch;
					if (mode == REGEX_FILTERING)
						nomatch = regexec(&regex, columnValue, 0, NULL, 0);
					else if(mode == FILE_FILTERING)
					{
						columnValue[Util::getLastNonWhitespace(columnValue,colStrLen)] = '\0'; // remove the whitespaces at the end
						std::string toSearch(columnValue);
						
						std::vector<std::string>::iterator foundInFilter = std::upper_bound(filter.begin(), filter.end(), toSearch, compareString());
						if (foundInFilter != filter.end() && toSearch.compare(*foundInFilter) == 0)
						{ // Found in filter
							if (positiveFiltering)
								nomatch = 0; // add to the output
							else
								nomatch = 1;
						} else {
							// not found in the filter
							if (positiveFiltering)
								nomatch = 1; // do NOT add to the output
							else
								nomatch = 0;
						}
			
					} else // Unknown filtering mode, keep all entries
						nomatch = 0;
                 
                if(!(nomatch)){
                    buffer.append(lineBuffer);
                    buffer.append("\n");
                }
                data = Util::skipLine(data);
            }

            dbw->write(buffer.c_str(), buffer.length(), (char*) SSTR(dataDb->getDbKey(id)).c_str(), thread_idx);
            buffer.clear();
        }
        delete [] lineBuffer;
        delete [] columnValue;
        delete [] columnPointer;
    }

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

	if (par.filteringFile == "")
	{
		ffindexFilter filter(par.db1,
						par.db2,
						par.threads,
						static_cast<size_t>(par.filterColumn),
						par.filterColumnRegex);
						
		return filter.runFilter();
	} else {
		std::cout<<"Filtering by file "<<par.filteringFile << std::endl;
		ffindexFilter filter(par.db1,
						par.db2,
						par.filteringFile,
						par.threads,
						static_cast<size_t>(par.filterColumn),
						par.positiveFilter);
		return filter.runFilter();
	}
	
	return -1;
}
