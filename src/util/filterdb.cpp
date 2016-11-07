#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "Debug.h"
#include "filterdb.h"

#include <fstream>
#include <iostream>
#include <algorithm>


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

ffindexFilter::ffindexFilter(std::string inDB, std::string outDB, int threads, size_t column, float compValue,
                 std::string compOperator) :
		inDB(inDB), outDB(outDB), threads(threads), column(column), trimToOneColumn(false),compValue(compValue),compOperator(compOperator) {

	initFiles();
	mode = NUMERIC_COMPARISON;
    
}



ffindexFilter::ffindexFilter(std::string inDB, std::string outDB, int threads, size_t column, int numberOfLines) :
		inDB(inDB), outDB(outDB), threads(threads), column(column), trimToOneColumn(false),numberOfLines(numberOfLines) {

	initFiles();
	mode = GET_FIRST_LINES;
}

ffindexFilter::ffindexFilter(std::string inDB, std::string outDB, int threads, size_t column, std::string regexStr,
							 bool trimToOneColumn) :
		inDB(inDB), outDB(outDB), threads(threads), column(column), regexStr(regexStr), trimToOneColumn(trimToOneColumn) {

	initFiles();

	mode = REGEX_FILTERING;

	int status = regcomp(&regex, regexStr.c_str(), REG_EXTENDED | REG_NEWLINE);
	if (status != 0 ){
		Debug(Debug::INFO) << "Error in regex " << regexStr << "\n";
		EXIT(EXIT_FAILURE);
	}
}

ffindexFilter::ffindexFilter(std::string inDB, std::string outDB, std::string filterFile, int threads, size_t column,
							 bool positiveFiltering) :
		inDB(inDB), outDB(outDB), filterFile(filterFile), threads(threads), column(column), trimToOneColumn(false),
		positiveFiltering(positiveFiltering) {

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

ffindexFilter::ffindexFilter(std::string inDB, std::string outDB, std::string filterFile, int threads, size_t column) :
		inDB(inDB), outDB(outDB), filterFile(filterFile), threads(threads), column(column), trimToOneColumn(false) {

	initFiles();

	mode = FILE_MAPPING;

	// Fill the filter with the data contained in the file
	std::ifstream filterFileStream;
	filterFileStream.open(filterFile);
	std::string line;
	while (std::getline(filterFileStream,line))
	{
		std::string keyOld,keyNew;
		std::istringstream lineToSplit(line);
		std::getline(lineToSplit,keyOld,'\t');
		std::getline(lineToSplit,keyNew,'\t');


		mapping.push_back(std::make_pair(keyOld, keyNew));
	}

	std::stable_sort(mapping.begin(), mapping.end(), compareFirstString());

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

			Debug::printProgress(id);
			int thread_idx = 0;
#ifdef OPENMP
			thread_idx = omp_get_thread_num();
#endif
			char *data = dataDb->getData(id);
			size_t dataLength = dataDb->getSeqLens(id);
			int counter = 0;
			while (*data != '\0') {
				if(!Util::getLine(data, dataLength, lineBuffer, LINE_BUFFER_SIZE)) {
					Debug(Debug::WARNING) << "Warning: Identifier was too long and was cut off!\n";
					data = Util::skipLine(data);
					continue;
				}
				size_t foundElements = Util::getWordsOfLine(lineBuffer, columnPointer, column + 1);
				if(foundElements < column  ){
					Debug(Debug::ERROR) << "Column=" << column << " does not exist in line " << lineBuffer << "\n";
					EXIT(EXIT_FAILURE);
				}
                
				counter++;
				size_t colStrLen;
				// if column is last column
				if(column == foundElements){
					const size_t entrySize = Util::skipNoneWhitespace(columnPointer[(column - 1)]); //Util::skipLine(data)
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
				if(mode == GET_FIRST_LINES){
					nomatch = 0; // output the line
					if(counter > numberOfLines){
						nomatch = 1; // hide the line in the output
					}
				} else if (mode == NUMERIC_COMPARISON) {
                    float toCompare = atof(columnValue);
                    if (compOperator == GREATER_OR_EQUAL) {
                        nomatch = !(toCompare >= compValue); // keep if the comparison is true
                    } else if(compOperator == LOWER_OR_EQUAL) {
                        nomatch = !(toCompare <= compValue); // keep if the comparison is true
                    } else if(compOperator == EQUAL) {
                        nomatch = !(toCompare == compValue); // keep if the comparison is true
                    } else {
                        nomatch = 0;
                    } 

                } else if (mode == REGEX_FILTERING){
					nomatch = regexec(&regex, columnValue, 0, NULL, 0);
				}
				else // i.e. (mode == FILE_FILTERING || mode == FILE_MAPPING)
				{
					columnValue[Util::getLastNonWhitespace(columnValue,colStrLen)] = '\0'; // remove the whitespaces at the end
					std::string toSearch(columnValue);

					if (mode == FILE_FILTERING)
					{
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
					} else if(mode == FILE_MAPPING) {
						std::vector<std::pair<std::string,std::string>>::iterator foundInFilter = std::lower_bound(mapping.begin(), mapping.end(), toSearch, compareToFirstString());
                      
                      nomatch = 1; // by default, do NOT add to the output
                      
                      char *newLineBuffer = new char[LINE_BUFFER_SIZE];
                      size_t newLineBufferIndex = 0;
                      char *endLine = lineBuffer + dataLength;
                      *newLineBuffer = '\0';
                      
                      for (size_t i = 0;i<dataLength;i++)
                          if (lineBuffer[i] == '\n' || lineBuffer[i] == '\0')
                          {
                              endLine = lineBuffer+i;
                              break;
                          }
                      size_t fieldLength = Util::skipNoneWhitespace(columnPointer[column-1]);
                      
		      // Output all the possible mapping value
                      while (foundInFilter != mapping.end() && toSearch.compare(foundInFilter->first) == 0)
                        {
                            nomatch = 0; // add to the output
                            
                            // copy the previous columns
                            memcpy(newLineBuffer + newLineBufferIndex,lineBuffer,columnPointer[column-1] - columnPointer[0]);
                            newLineBufferIndex += columnPointer[column-1] - columnPointer[0];
                            
                            // map the current column value
                            memcpy(newLineBuffer + newLineBufferIndex,(foundInFilter->second).c_str(),(foundInFilter->second).length());
                            newLineBufferIndex += (foundInFilter->second).length();
                            
                            
                            // copy the next columns                            
                            if (foundElements > column)
                            {
                                memcpy(newLineBuffer + newLineBufferIndex,columnPointer[column-1]+fieldLength,endLine - (columnPointer[column-1]+fieldLength));
                                newLineBufferIndex += endLine - (columnPointer[column-1]+fieldLength);
                            } else {
                                newLineBuffer[newLineBufferIndex++] = '\n';
                            }
                            newLineBuffer[newLineBufferIndex] = '\0';
                            
                            foundInFilter++;
                        }
                      if(!nomatch)
                        memcpy(lineBuffer,newLineBuffer,newLineBufferIndex+1);

                      delete [] newLineBuffer;

					} else // Unknown filtering mode, keep all entries
						nomatch = 0;

				}

				if(!(nomatch)){
					if (trimToOneColumn)
						buffer.append(columnValue);
					else
						buffer.append(lineBuffer);
                  if (buffer.back() != '\n')
                        buffer.append("\n");
				}
				data = Util::skipLine(data);
			}

            dbw->writeData(buffer.c_str(), buffer.length(), (char *) SSTR(dataDb->getDbKey(id)).c_str(), thread_idx);
			buffer.clear();
		}
		delete [] lineBuffer;
		delete [] columnValue;
		delete [] columnPointer;
	}

	return 0;
}

int filterdb(int argc, const char **argv, const Command& command) {
	Parameters& par = Parameters::getInstance();
	par.parseParameters(argc, argv, command, 2);

#ifdef OPENMP
	omp_set_num_threads(par.threads);
#endif

	if (par.filteringFile != "")
	{
		std::cout<<"Filtering by file "<<par.filteringFile << std::endl;
		ffindexFilter filter(par.db1,
							 par.db2,
							 par.filteringFile,
							 par.threads,
							 static_cast<size_t>(par.filterColumn),
							 par.positiveFilter);
		return filter.runFilter();
	} else if(par.mappingFile != ""){
		std::cout<<"Mapping keys by file "<<par.mappingFile << std::endl;
		ffindexFilter filter(par.db1,
							 par.db2,
							 par.mappingFile,
							 par.threads,
							 static_cast<size_t>(par.filterColumn));
		return filter.runFilter();
	} else if(par.extractLines > 0){ // GET_FIRST_LINES mode
		ffindexFilter filter(par.db1,
							 par.db2,
							 par.threads,
							 static_cast<size_t>(par.filterColumn),
                           par.extractLines);
		return filter.runFilter();
	} else if(par.compOperator != ""){ // NUMERIC_COMPARISON mode
		ffindexFilter filter(par.db1,
							 par.db2,
							 par.threads,
							 static_cast<size_t>(par.filterColumn),
                           par.compValue,
                           par.compOperator);
		return filter.runFilter();
	} else { // RegEx filter
		ffindexFilter filter(par.db1,
							 par.db2,
							 par.threads,
							 static_cast<size_t>(par.filterColumn),
							 par.filterColumnRegex,par.trimToOneColumn);
		return filter.runFilter();
	}

	return -1;
}
