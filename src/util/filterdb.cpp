#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "Debug.h"
#include "filterdb.h"
#include "FileUtil.h"

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
	dataDb=new DBReader<unsigned int>(inDB.c_str(),(std::string(inDB).append(".index")).c_str(), threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
	dataDb->open(DBReader<unsigned int>::LINEAR_ACCCESS);

	dbw = new DBWriter(outDB.c_str(), (std::string(outDB).append(".index")).c_str(), threads, compressed, dataDb->getDbtype());
	dbw->open();
	return 0;
}




ffindexFilter::ffindexFilter(Parameters &par) {
    inDB = std::string(par.db1);
    outDB = std::string(par.db2);
    threads = par.threads;
    compressed = par.compressed;
    column  = static_cast<size_t>(par.filterColumn);
    columnToTake = par.columnToTake;
    trimToOneColumn = par.trimToOneColumn;
    positiveFiltering = par.positiveFilter;
    shouldAddSelfMatch = par.includeIdentity;
    parser = NULL;
    
	initFiles();


    if (par.sortEntries)
    {
        mode = SORT_ENTRIES;
        std::cout<<"Filtering by sorting entries."<<std::endl;
        sortingMode = par.sortEntries;
    } else if (par.filteringFile != "") {
        mode = FILE_FILTERING;
//        filter.reserve(1000000);
        std::cout << "Filtering with a filter files." << std::endl;
        filterFile = par.filteringFile;
        // Fill the filter with the data contained in the file
        std::vector<std::string> filenames;
        if (FileUtil::fileExists(filterFile.c_str())) {
            filenames.push_back(filterFile);
        } else if (FileUtil::fileExists((filterFile + ".dbtype").c_str())) {
            filenames = FileUtil::findDatafiles(filterFile.c_str());
        } else {
            Debug(Debug::ERROR) << "File " << filterFile << " does not exist.\n";
            EXIT(EXIT_FAILURE);
        }
        char *line = new char[65536];
        size_t len = 0;
        char * key=new char[65536];
        for(size_t i = 0; i < filenames.size(); i++) {
            FILE *orderFile = fopen(filenames[i].c_str(), "r");
            while (getline(&line, &len, orderFile) != -1) {
                size_t offset = 0;
                // ignore \0 in data files
                // to support datafiles as input
                while (offset < len && line[offset] == '\0') {
                    offset++;
                }
                if (offset >= len) {
                    break;
                }
                Util::parseKey(line + offset, key);
                filter.emplace_back(key);
            }
            fclose(orderFile);
        }
        delete [] key;
        delete [] line;
        omptl::sort(filter.begin(), filter.end());
        std::vector<std::string>::iterator last = std::unique(filter.begin(), filter.end());
        filter.erase(last, filter.end());
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
    } else if(!par.joinDB.empty()){
        mode = JOIN_DB;
        std::string joinIndex(par.joinDB);
        joinIndex.append(".index");
        joinDB = new DBReader<unsigned int>(par.joinDB.c_str(), joinIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        joinDB->open(DBReader<unsigned int>::NOSORT);
        std::cout << "Joining targets to query database.\n";
    } else if (!par.compPos.empty()) {
        mode = COMPUTE_POSITIONS;
        std::string swapIndex = par.compPos + ".index";
        swapDB = new DBReader<unsigned int>(par.compPos.c_str(), swapIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        swapDB->open(DBReader<unsigned int>::NOSORT);
        std::string A = par.db3;
        std::cout << "Swapping fields\n";
    } else if (!par.clusterFile.empty()){
        mode = TRANSITIVE_REPLACE;
        std::string clusterIndex = par.clusterFile + ".index";
        clusterDB = new DBReader<unsigned int>(par.clusterFile.c_str(), clusterIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        clusterDB->open(DBReader<unsigned int>::NOSORT);
        std::cout << "Replacing target Field by clusters Genes\n";
    } else if (par.beatsFirst){
        mode = BEATS_FIRST;
        std::cout << "Filter by numerical comparison to first row.\n";
        compOperator = par.compOperator;
    } else if(par.compOperator != "") {
        mode = NUMERIC_COMPARISON;
        std::cout << "Filtering by numerical comparison.\n";
        compValue = par.compValue;
        compOperator = par.compOperator;
    } else if (par.filterExpression != "") {
        mode = EXPRESSION_FILTERING;
        parser = new ExpressionParser(par.filterExpression.c_str());
        if (parser->isOk() == false){
            Debug(Debug::INFO) << "Error in expression " << par.filterExpression << "\n";
            EXIT(EXIT_FAILURE);
        }
        bindableParserColumns = parser->findBindableIndices();
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
    Debug::Progress progress(dataDb->getSize());

#pragma omp parallel
	{
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

		char *lineBuffer = new char[LINE_BUFFER_SIZE];
		char *columnValue = new char[LINE_BUFFER_SIZE];
		const char **columnPointer = new const char*[column + 1];

		std::string buffer = "";
		buffer.reserve(LINE_BUFFER_SIZE);

#pragma omp for schedule(dynamic, 10)
		for (size_t id = 0; id < dataDb->getSize(); id++) {
			progress.updateProgress();

			char *data = dataDb->getData(id,  thread_idx);
            unsigned int queryKey = dataDb->getDbKey(id);
			size_t dataLength = dataDb->getSeqLens(id);
			int counter = 0;
            
            std::vector<std::pair<double, std::string>> toSort;
            bool addSelfMatch = false;

			while (*data != '\0') {
                if (shouldAddSelfMatch) {
                    char dbKeyBuffer[255 + 1];
                    Util::parseKey(data, dbKeyBuffer);
                    const unsigned int curKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                    addSelfMatch = (queryKey == curKey);
                }
                    
				if(!Util::getLine(data, dataLength, lineBuffer, LINE_BUFFER_SIZE)) {
					Debug(Debug::WARNING) << "Identifier was too long and was cut off!\n";
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
                    } else if (compOperator == LOWER_OR_EQUAL) {
                        nomatch = !(toCompare <= compValue); // keep if the comparison is true
                    } else if (compOperator == EQUAL) {
                        nomatch = !(toCompare == compValue); // keep if the comparison is true
                    } else {
                        nomatch = 0;
                    }
                } else if (mode == EXPRESSION_FILTERING) {
				    const char* columnPointers[128];
                    Util::getWordsOfLine(lineBuffer, columnPointers, 128);
				    for (size_t i = 0; i < bindableParserColumns.size(); ++i) {
				        size_t columnToBind = bindableParserColumns[i];
				        char* rest;
				        const double value = strtod(columnPointers[columnToBind], &rest);
                        if ((rest == columnPointers[columnToBind]) || errno == ERANGE) {
                            Debug(Debug::WARNING) << "Can not parse column " << columnToBind << "!\n";
                            continue;
                        }
				        parser->bind(columnToBind, value);
				    }
				    const double result = parser->evaluate();
                    nomatch = (result == 0);
                } else if (mode == REGEX_FILTERING){
					nomatch = regexec(&regex, columnValue, 0, NULL, 0);
                } else if (mode == JOIN_DB){
                    size_t newId = joinDB->getId(static_cast<unsigned int>(strtoul(columnValue, NULL, 10)));
                    size_t originalLength = strlen(lineBuffer);
                    // Replace the last \n
                    lineBuffer[originalLength - 1] = '\t';
                    char* fullLine = joinDB->getData(newId, thread_idx);
                    // either append the full line (default mode):
                    if (columnToTake == -1) {
                        size_t fullLineLength = joinDB->getSeqLens(newId);
                        // Appending join database entry to query database entry
                        memcpy(lineBuffer + originalLength, fullLine, fullLineLength);
                    }
                    // or a specified column:
                    else {
                        if(*fullLine != '\0'){
                            std::vector<std::string> splittedLine = Util::split(fullLine, "\t") ;
                            char* newValue = const_cast<char *>(splittedLine[columnToTake].c_str());
                            size_t valueLength = joinDB->getSeqLens(newId);
                            // Appending join database entry to query database entry
                            memcpy(lineBuffer + originalLength, newValue, valueLength);
                        }
                    }
                }
                else if (mode == TRANSITIVE_REPLACE) {
                    std::string singleGene;

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

                    char *clusterGenes = clusterDB->getDataByDBKey(Util::fast_atoi<unsigned int>(columnValue), thread_idx);
                    std::stringstream stringstreamClusterGenes(clusterGenes);


                    while (std::getline(stringstreamClusterGenes, singleGene)) {

                        if (!singleGene.empty()) {

                            // copy the previous columns
                            memcpy(newLineBuffer + newLineBufferIndex,lineBuffer,columnPointer[column-1] - columnPointer[0]);
                            newLineBufferIndex += columnPointer[column-1] - columnPointer[0];

                            // map the current column value
                            memcpy(newLineBuffer + newLineBufferIndex,singleGene.c_str(),singleGene.length());
                            newLineBufferIndex += singleGene.length();


                            // copy the next columns
                            if (foundElements > column)
                            {
                                memcpy(newLineBuffer + newLineBufferIndex,columnPointer[column-1]+fieldLength,endLine - (columnPointer[column-1]+fieldLength));
                                newLineBufferIndex += endLine - (columnPointer[column-1]+fieldLength);
                            } else {
                                newLineBuffer[newLineBufferIndex++] = '\n';
                            }

                            if( newLineBuffer[newLineBufferIndex-1] != '\n') {
                                newLineBuffer[newLineBufferIndex++] = '\n';
                            }

                            newLineBuffer[newLineBufferIndex] = '\0';

                        }

                    }

                    if(!nomatch)
                        memcpy(lineBuffer,newLineBuffer,newLineBufferIndex+1);
                    delete [] newLineBuffer;

                }
                else if (mode == COMPUTE_POSITIONS) {
				    // Optimise it
                    std::vector<std::string> splittedOriginalLine = Util::split(lineBuffer, "\t");
                    char *lineWithNewFields = swapDB->getDataByDBKey(
                            static_cast<unsigned int>(strtoul(columnValue, NULL, 10)), thread_idx) ;
                    std::vector<std::string> splittedLineWithNewFields = Util::split(lineWithNewFields, "\t");

                    unsigned long posStart = std::stoul(splittedOriginalLine[7].c_str()) * 3;
                    unsigned long posStop = std::stoul(splittedOriginalLine[8].c_str()) * 3;

                    if (posStart < posStop) {
                        unsigned long startPosHitOnGenome =
                                std::stoul(splittedLineWithNewFields[7].c_str(), NULL) + posStart;
                        unsigned long endPosHitOnGenome =
                                std::stoul(splittedLineWithNewFields[7].c_str(), NULL) + posStop;
                        splittedOriginalLine.insert(splittedOriginalLine.begin() + 8,
                                                    std::to_string(startPosHitOnGenome) + "\t");
                        splittedOriginalLine.insert(splittedOriginalLine.begin() + 10,
                                                    std::to_string(endPosHitOnGenome) + "\t");
                    } else {
                        unsigned long startPosHitOnGenome =
                                std::stoul(splittedLineWithNewFields[7].c_str(), NULL) - posStart;
                        unsigned long endPosHitOnGenome = std::stoul(splittedLineWithNewFields[7].c_str(), NULL) - posStop;
                        splittedOriginalLine.insert(splittedOriginalLine.begin() + 8, std::to_string(startPosHitOnGenome) + "\t");
                        splittedOriginalLine.insert(splittedOriginalLine.begin() + 10, std::to_string(endPosHitOnGenome) + "\t");
                    }
                    std::string tempBuffer = "" ;
                    for (std::vector<std::string>::const_iterator i = splittedOriginalLine.begin(); i != splittedOriginalLine.end(); ++i) {
                        tempBuffer.append(*i);
                        tempBuffer.append("\t");
                    }
                    tempBuffer.append("\n") ;
                    strcpy(lineBuffer, tempBuffer.c_str()) ;
                }
                else if (mode == BEATS_FIRST){
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
                    if (mode == FILE_FILTERING) {
                        std::vector<std::string>::iterator foundInFilter = std::upper_bound(filter.begin(),
                                                                                            filter.end(), toSearch,
                                                                                            compareString());
                        if (foundInFilter != filter.end() && toSearch.compare(*foundInFilter) == 0) {
                            // Found in filter
                            if (positiveFiltering)
                                nomatch = 0;
                            else
                                nomatch = 1;
                        } else {
                            // not found in the filter
                            if (positiveFiltering)
                                nomatch = 1;
                            else
                                nomatch = 0;
                        }
                    } else if (mode == FILE_MAPPING) {
                        std::vector<std::pair<std::string, std::string>>::iterator foundInFilter = std::lower_bound(
                                mapping.begin(), mapping.end(), toSearch, compareToFirstString());

                        // by default, do NOT add to the output
                        nomatch = 1;

                        char *newLineBuffer = new char[LINE_BUFFER_SIZE];
                        size_t newLineBufferIndex = 0;
                        char *endLine = lineBuffer + dataLength;
                        *newLineBuffer = '\0';

                        for (size_t i = 0; i < dataLength; i++) {
                            if (lineBuffer[i] == '\n' || lineBuffer[i] == '\0') {
                                endLine = lineBuffer + i;
                                break;
                            }
                        }
                        size_t fieldLength = Util::skipNoneWhitespace(columnPointer[column - 1]);

                        // Output all the possible mapping value
                        while (foundInFilter != mapping.end() && toSearch.compare(foundInFilter->first) == 0) {
                            nomatch = 0;

                            // copy the previous columns
                            memcpy(newLineBuffer + newLineBufferIndex, lineBuffer,
                                   columnPointer[column - 1] - columnPointer[0]);
                            newLineBufferIndex += columnPointer[column - 1] - columnPointer[0];

                            // map the current column value
                            memcpy(newLineBuffer + newLineBufferIndex, (foundInFilter->second).c_str(),
                                   (foundInFilter->second).length());
                            newLineBufferIndex += (foundInFilter->second).length();

                            // copy the next columns
                            if (foundElements > column) {
                                memcpy(newLineBuffer + newLineBufferIndex, columnPointer[column - 1] + fieldLength,
                                       endLine - (columnPointer[column - 1] + fieldLength));
                                newLineBufferIndex += endLine - (columnPointer[column - 1] + fieldLength);
                            } else {
                                newLineBuffer[newLineBufferIndex++] = '\n';
                            }
                            newLineBuffer[newLineBufferIndex] = '\0';

                            foundInFilter++;
                        }
                        if (nomatch == 0) {
                            memcpy(lineBuffer, newLineBuffer, newLineBufferIndex + 1);
                        }
                        delete[] newLineBuffer;
                    } else if (mode == SORT_ENTRIES) {
                        toSort.push_back(
                                std::make_pair<double, std::string>(std::strtod(columnValue, NULL), lineBuffer));
                        nomatch = 1; // do not put anything in the output buffer
                    }
                    else // Unknown filtering mode, keep all entries
						nomatch = 0;

				}
                
                if (addSelfMatch) {
                    nomatch = 0;
                }  
                    
				if(!(nomatch)) {
                    if (trimToOneColumn) {
                        buffer.append(columnValue);
                    }
                    else {
						buffer.append(lineBuffer);
                    }
                    
                    if (buffer.back() != '\n') {
                        buffer.append("\n");
                    }
				}
				data = Util::skipLine(data);
			}

            if(mode == SORT_ENTRIES)
            {
                if (sortingMode ==INCREASING)
                    std::stable_sort(toSort.begin(),toSort.end(),compareFirstEntry());
                else if (sortingMode == DECREASING)
                    std::stable_sort(toSort.begin(),toSort.end(),compareFirstEntryDecreasing());
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

    ffindexFilter filter(par);
    return filter.runFilter();
}
