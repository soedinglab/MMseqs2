#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "Debug.h"
#include "filterdb.h"

#include <fstream>
#include <iostream>
#include <algorithm>

#include <omptl/omptl_algorithm>

#ifdef OPENMP
#include <omp.h>
#endif


struct compareFirstEntry {
    bool operator()(const std::pair <double,std::string> &lhs,
                    const std::pair <double,std::string> &rhs) const {
        return (lhs.first < rhs.first);
    }
};


struct compareFirstEntryDecreasing {
    bool operator()(const std::pair <double,std::string> &lhs,
                    const std::pair <double,std::string> &rhs) const {
        return (lhs.first > rhs.first);
    }
};



int ffindexFilter::initFiles() {
	dataDb=new DBReader<unsigned int>(inDB.c_str(),(std::string(inDB).append(".index")).c_str());
	dataDb->open(DBReader<unsigned int>::LINEAR_ACCCESS);

	dbw = new DBWriter(outDB.c_str(), (std::string(outDB).append(".index")).c_str(), threads);
	dbw->open();
	return 0;
}




ffindexFilter::ffindexFilter(Parameters &par) {
    inDB = std::string(par.db1);
    outDB = std::string(par.db2);
    threads = par.threads;
    column  = static_cast<size_t>(par.filterColumn);
    trimToOneColumn = par.trimToOneColumn;
    positiveFiltering = par.positiveFilter;
    shouldAddSelfMatch = par.includeIdentity;
    
	initFiles();


    if (par.sortEntries)
    {
        mode = SORT_ENTRIES;
        std::cout<<"Filtering by sorting entries."<<std::endl;
        sortingMode = par.sortEntries;
    } else if (par.filteringFile != "")
	{
        mode = FILE_FILTERING;
        std::cout<<"Filtering with a filter files."<<std::endl;
        filterFile = par.filteringFile;
        // Fill the filter with the data contained in the file
        std::ifstream filterFileStream;
        filterFileStream.open(filterFile);
        std::string line;
        while (std::getline(filterFileStream,line))
        {
            filter.push_back(line);
        }

        std::stable_sort(filter.begin(), filter.end(), compareString());
    } else if(par.mappingFile != "")
    {

        mode = FILE_MAPPING;
        std::cout<<"Filtering by mapping values."<<std::endl;
        filterFile = par.mappingFile;

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
    } else if(par.extractLines > 0){ // GET_FIRST_LINES mode
        mode = GET_FIRST_LINES;
        numberOfLines = par.extractLines;
        std::cout << "Filtering by extracting the first " << numberOfLines << " lines.\n";
    } else if (par.beatsFirst){
        mode = BEATS_FIRST;
        std::cout << "Filter by numerical comparison to first row.\n";
        compOperator = par.compOperator;
    } else if(par.compOperator != "") {
        mode = NUMERIC_COMPARISON;
        std::cout << "Filtering by numerical comparison.\n";
        compValue = par.compValue;
        compOperator = par.compOperator;
    } else {
        mode = REGEX_FILTERING;
        std::cout << "Filtering by RegEx.\n";
        regexStr = par.filterColumnRegex;
        int status = regcomp(&regex, regexStr.c_str(), REG_EXTENDED | REG_NEWLINE);
        if (status != 0 ){
            Debug(Debug::INFO) << "Error in regex " << regexStr << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
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
            unsigned int queryKey = dataDb->getDbKey(id);
			size_t dataLength = dataDb->getSeqLens(id);
			int counter = 0;
            
            std::vector<std::pair<double, std::string>> toSort;
            bool addSelfMatch = false;

			while (*data != '\0') {
                if (shouldAddSelfMatch)
                {
                    char dbKeyBuffer[255 + 1];
                    Util::parseKey(data, dbKeyBuffer);
                    const unsigned int curKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                    addSelfMatch = (queryKey == curKey);
                    
                }
                    
				if(!Util::getLine(data, dataLength, lineBuffer, LINE_BUFFER_SIZE)) {
					Debug(Debug::WARNING) << "Warning: Identifier was too long and was cut off!\n";
					data = Util::skipLine(data);
					continue;
				}

                counter++;
                size_t foundElements = 1;
                if (mode != GET_FIRST_LINES) {
                    foundElements = Util::getWordsOfLine(lineBuffer, columnPointer, column + 1);
                    if(foundElements < column  ){
                        Debug(Debug::ERROR) << "Column=" << column << " does not exist in line " << lineBuffer << "\n";
                        EXIT(EXIT_FAILURE);
                    }

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

                    columnValue[Util::getLastNonWhitespace(columnValue,colStrLen)] = '\0'; // remove the whitespaces at the end
                }

				int nomatch = 0;
				if(mode == GET_FIRST_LINES){
					nomatch = 0; // output the line
					if(counter > numberOfLines){
						nomatch = 1; // hide the line in the output
					}
				} else if (mode == NUMERIC_COMPARISON) {
                    double toCompare = strtod(columnValue, NULL);
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
				} else if (mode == BEATS_FIRST) {
                    if (counter == 1) {
                        compValue = strtod(columnValue, NULL);
                    } else {
                        double toCompare = strtod(columnValue, NULL);
                        if (compOperator == GREATER_OR_EQUAL) {
                            nomatch = !(toCompare >= compValue); // keep if the comparison is true
                        } else if(compOperator == LOWER_OR_EQUAL) {
                            nomatch = !(toCompare <= compValue); // keep if the comparison is true
                        } else if(compOperator == EQUAL) {
                            nomatch = !(toCompare == compValue); // keep if the comparison is true
                        } else {
                            nomatch = 0;
                        }
                    }
                } else {
                    // i.e. (mode == FILE_FILTERING || mode == FILE_MAPPING)
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

					}  else if(mode == SORT_ENTRIES) {
                        toSort.push_back(std::make_pair<double, std::string>(std::strtod(columnValue, NULL), lineBuffer));
                        nomatch = 1; // do not put anything in the output buffer
                  }
                    else // Unknown filtering mode, keep all entries
						nomatch = 0;

				}
                
                if (addSelfMatch)
                    nomatch = 0;
                    
                    
				if(!(nomatch)){
                    if (trimToOneColumn)
                    {
                        buffer.append(columnValue);
                    }
                    else
                    {
						buffer.append(lineBuffer);
                    }
                    
                    if (buffer.back() != '\n')
                        buffer.append("\n");
				}
				data = Util::skipLine(data);
			}

            if(mode == SORT_ENTRIES)
            {
                if (sortingMode ==INCREASING)
                    omptl::sort(toSort.begin(),toSort.end(),compareFirstEntry());
                else if (sortingMode == DECREASING)
                    omptl::sort(toSort.begin(),toSort.end(),compareFirstEntryDecreasing());
                else if (sortingMode == SHUFFLE)
                {
                    srand ( unsigned ( time(0) ) );
                    std::random_shuffle(toSort.begin(),toSort.end());
                }

                
                for (size_t i = 0; i< toSort.size(); i++)
                {
                    buffer.append(toSort[i].second);    
                    if (buffer.back() != '\n')
                      buffer.append("\n");
                }
                
            }
            dbw->writeData(buffer.c_str(), buffer.length(), dataDb->getDbKey(id), thread_idx);
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

    ffindexFilter filter(par);
    return filter.runFilter();
}
