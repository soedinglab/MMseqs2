#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "Debug.h"
#include "FileUtil.h"
#include "ExpressionParser.h"
#include "FastSort.h"
#include <fstream>
#include <random>
#include <iostream>

#include <regex.h>

#ifdef OPENMP
#include <omp.h>
#endif

#define REGEX_FILTERING       0
#define FILE_FILTERING        1
#define FILE_MAPPING          2
#define GET_FIRST_LINES       3
#define NUMERIC_COMPARISON    4
#define SORT_ENTRIES          5
#define BEATS_FIRST           6
#define JOIN_DB               7
#define EXPRESSION_FILTERING 10


enum ComparisonOperator {
    OP_GEQ,
    OP_LEQ,
    OP_EQ,
    OP_IN_P,
    OP_OUT_P,
    OP_EQ_P,

    OP_INVALID
};

ComparisonOperator mapOperator(const std::string& op) {
    if (op == "ge") return OP_GEQ;
    if (op == "le") return OP_LEQ;
    if (op == "e")  return OP_EQ;
    if (op == "ip") return OP_IN_P;
    if (op == "op") return OP_OUT_P;
    if (op == "ep") return OP_EQ_P;
    return OP_INVALID;
}

#define INCREASING 1
#define DECREASING 2
#define SHUFFLE    3

struct compareString {
    bool operator() (const std::string& lhs, const std::string& rhs) const{
        return (lhs.compare(rhs) <= 0);
    }
};

struct compareFirstString {
    bool operator() (const std::pair<std::string, std::string>& lhs, const std::pair<std::string,std::string>& rhs) const{
        return (lhs.first.compare(rhs.first) <= 0);
    }
};

struct compareToFirstString {
    bool operator() (const std::pair<std::string,std::string>& lhs,const std::string& rhs) const{
        return (lhs.first.compare(rhs) < 0);
    }
};

struct compareFirstEntry {
    bool operator()(const std::pair<double, std::string> &lhs, const std::pair<double, std::string> &rhs) const {
        return (lhs.first < rhs.first);
    }
};

struct compareFirstEntryDecreasing {
    bool operator()(const std::pair<double, std::string> &lhs, const std::pair<double, std::string> &rhs) const {
        return (lhs.first > rhs.first);
    }
};

int filterdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    const size_t column = static_cast<size_t>(par.filterColumn);
    const int columnToTake = par.columnToTake;
    const bool trimToOneColumn = par.trimToOneColumn;
    // positiveFilter = true => outDB = inDB \intersect filter ; othw : outDB = inDB - filter
    const bool positiveFiltering = par.positiveFilter;
    const bool shouldAddSelfMatch = par.includeIdentity;
    const ComparisonOperator compOperator = mapOperator(par.compOperator);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, reader.getDbtype());
    writer.open();

    // FILE_FILTERING
    std::vector<std::string> filter;

    // FILE_MAPPING
    std::vector<std::pair<std::string, std::string>> mapping;

    // JOIN_DB
    DBReader<unsigned int>* helper = NULL;

    // REGEX_FILTERING
    regex_t regex;
    std::random_device rng;
    std::mt19937 urng(rng());
    int mode;
    if (par.sortEntries != 0) {
        mode = SORT_ENTRIES;
        Debug(Debug::INFO) << "Filtering by sorting entries\n";
    } else if (par.filteringFile.empty() == false) {
        mode = FILE_FILTERING;
        Debug(Debug::INFO) << "Filtering using file(s)\n";
        // Fill the filter with the data contained in the file
        std::vector<std::string> filenames;
        if (FileUtil::fileExists(par.filteringFile.c_str())) {
            filenames.push_back(par.filteringFile);
        } else if (FileUtil::fileExists((par.filteringFile + ".dbtype").c_str())) {
            filenames = FileUtil::findDatafiles(par.filteringFile.c_str());
        } else {
            Debug(Debug::ERROR) << "File " << par.filteringFile << " does not exist\n";
            EXIT(EXIT_FAILURE);
        }
        char key[65536];
        for (size_t i = 0; i < filenames.size(); i++) {
            FILE * orderFile = fopen(filenames[i].c_str(), "r");
            int c;
            size_t offset = 0;
            bool inKey = true;
            // parse first column in each line without tripping over additional null bytes
            // as we allow database data files as input
            while ((c = fgetc(orderFile)) != EOF) {
                if (c == '\n') {
                    if (offset > 0) {
                        key[offset] = '\0';
                        offset = 0;
                        filter.emplace_back(key);
                    }
                    inKey = true;
                    continue;
                }
                if (c == ' ' || c == '\t') {
                    inKey = false;
                    continue;
                }
                if (c == '\0' || inKey == false) {
                    continue;
                }

                key[offset] = c;
                offset++;

                if (offset == 65536) {
                    Debug(Debug::ERROR) << "Input in file " << filenames[i] << " too long\n";
                    EXIT(EXIT_FAILURE);
                }
            }
            if (inKey == true && offset > 0) {
                key[offset] = '\0';
                filter.emplace_back(key);
            }
            if (fclose(orderFile) != 0) {
                Debug(Debug::ERROR) << "Cannot close file " << filenames[i] << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        SORT_PARALLEL(filter.begin(), filter.end());
        std::vector<std::string>::iterator last = std::unique(filter.begin(), filter.end());
        filter.erase(last, filter.end());
    } else if (par.mappingFile.empty() == false) {
        mode = FILE_MAPPING;
        Debug(Debug::INFO) << "Filtering by mapping values\n";
        std::ifstream ss(par.mappingFile);
        std::string line;
        std::string keyOld, keyNew;
        while (std::getline(ss, line)) {
            std::istringstream lineToSplit(line);
            std::getline(lineToSplit, keyOld, '\t');
            std::getline(lineToSplit, keyNew, '\t');
            mapping.emplace_back(keyOld, keyNew);
        }
        std::stable_sort(mapping.begin(), mapping.end(), compareFirstString());
    } else if (par.extractLines > 0) {
        mode = GET_FIRST_LINES;
        Debug(Debug::INFO) << "Filtering by extracting the first " << par.extractLines << " lines\n";
    } else if (par.joinDB.empty() == false) {
        mode = JOIN_DB;
        std::string joinIndex = par.joinDB + ".index";
        helper = new DBReader<unsigned int>(par.joinDB.c_str(), joinIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        helper->open(DBReader<unsigned int>::NOSORT);
        Debug(Debug::INFO) << "Joining databases by column value\n";
    } else if (par.beatsFirst == true) {
        mode = BEATS_FIRST;
        Debug(Debug::INFO) << "Filtering by numerical comparison to first row\n";
    } else if (par.compOperator.empty() == false) {
        mode = NUMERIC_COMPARISON;
        Debug(Debug::INFO) << "Filtering by numerical comparison\n";
    } else if (par.filterExpression.empty() == false) {
        mode = EXPRESSION_FILTERING;
    } else {
        mode = REGEX_FILTERING;
        Debug(Debug::INFO) << "Filtering using regular expression\n";
        int status = regcomp(&regex, par.filterColumnRegex.c_str(), REG_EXTENDED | REG_NEWLINE);
        if (status != 0) {
            Debug(Debug::ERROR) << "Error in regex " << par.filterColumnRegex << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

    const size_t LINE_BUFFER_SIZE = 1000000;
    Debug::Progress progress(reader.getSize());
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        char *lineBuffer = new char[LINE_BUFFER_SIZE];
        char *columnValue = new char[LINE_BUFFER_SIZE];
        const char **columnPointer = new const char *[column + 1];

        char *newLineBuffer = new char[LINE_BUFFER_SIZE];

        double referenceValue = 0;

        std::string buffer = "";
        buffer.reserve(LINE_BUFFER_SIZE);

        std::vector<std::pair<double, std::string>> toSort;

        char dbKeyBuffer[255 + 1];

        // EXPRESSION_FILTERING
        ExpressionParser* parser = NULL;
        std::vector<int> bindableParserColumns;

        if (mode == EXPRESSION_FILTERING) {
            parser = new ExpressionParser(par.filterExpression.c_str());
            if (parser->isOk() == false) {
                Debug(Debug::INFO) << "Error in expression " << par.filterExpression << "\n";
                EXIT(EXIT_FAILURE);
            }
            bindableParserColumns = parser->findBindableIndices();
        }

#pragma omp for schedule(dynamic, 10)
        for (size_t id = 0; id < reader.getSize(); ++id) {
            progress.updateProgress();

            char *data = reader.getData(id, thread_idx);
            unsigned int queryKey = reader.getDbKey(id);
            size_t dataLength = reader.getEntryLen(id);
            int counter = 0;

            bool addSelfMatch = false;

            while (*data != '\0') {
                if (shouldAddSelfMatch) {
                    Util::parseKey(data, dbKeyBuffer);
                    const unsigned int curKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                    addSelfMatch = (queryKey == curKey);
                }

                if (!Util::getLine(data, dataLength, lineBuffer, LINE_BUFFER_SIZE)) {
                    Debug(Debug::WARNING) << "Identifier was too long and was cut off!\n";
                    data = Util::skipLine(data);
                    continue;
                }

                counter++;
                size_t foundElements = 1;
                if (mode != GET_FIRST_LINES || trimToOneColumn) {
                    foundElements = Util::getWordsOfLine(lineBuffer, columnPointer, column + 1);
                    if (foundElements < column) {
                        Debug(Debug::ERROR) << "Column=" << column << " does not exist in line " << lineBuffer << "\n";
                        EXIT(EXIT_FAILURE);
                    }

                    size_t colStrLen;
                    // if column is last column
                    if (column == foundElements) {
                        const size_t entrySize = Util::skipNoneWhitespace(columnPointer[(column - 1)]);
                        memcpy(columnValue, columnPointer[column - 1], entrySize);
                        columnValue[entrySize] = '\0';
                        colStrLen = entrySize;
                    } else {
                        const ptrdiff_t entrySize = columnPointer[column] - columnPointer[(column - 1)];
                        memcpy(columnValue, columnPointer[column - 1], entrySize);
                        columnValue[entrySize] = '\0';
                        colStrLen = entrySize;
                    }

                    // remove the whitespaces at the end
                    columnValue[Util::getLastNonWhitespace(columnValue, colStrLen)] = '\0';
                }

                int nomatch = 0;
                if (mode == GET_FIRST_LINES) {
                    // output the line
                    nomatch = 0;
                    if (counter > par.extractLines) {
                        // hide the line in the output
                        nomatch = 1;
                    }
                } else if (mode == NUMERIC_COMPARISON) {
                    double toCompare = strtod(columnValue, NULL);
                    if (compOperator == OP_GEQ) {
                        nomatch = !(toCompare >= par.compValue);
                    } else if (compOperator == OP_LEQ) {
                        nomatch = !(toCompare <= par.compValue);
                    } else if (compOperator == OP_EQ) {
                        nomatch = !(toCompare == par.compValue);
                    } else {
                        nomatch = 0;
                    }
                } else if (mode == EXPRESSION_FILTERING) {
                    const char *columnPointers[128];
                    Util::getWordsOfLine(lineBuffer, columnPointers, 128);
                    for (size_t i = 0; i < bindableParserColumns.size(); ++i) {
                        size_t columnToBind = bindableParserColumns[i];
                        char *rest;
                        errno = 0;
                        const double value = strtod(columnPointers[columnToBind], &rest);
                        if ((rest == columnPointers[columnToBind]) || errno == ERANGE) {
                            Debug(Debug::WARNING) << "Can not parse column " << columnToBind << "!\n";
                            continue;
                        }
                        parser->bind(columnToBind, value);
                    }
                    const double result = parser->evaluate();
                    nomatch = (result == 0);
                } else if (mode == REGEX_FILTERING) {
                    nomatch = regexec(&regex, columnValue, 0, NULL, 0);
                } else if (mode == JOIN_DB) {
                    size_t newId = helper->getId(static_cast<unsigned int>(strtoul(columnValue, NULL, 10)));
                    if (newId != UINT_MAX) {
                        size_t originalLength = strlen(lineBuffer);
                        // Continue the string by replacing the null byte
                        lineBuffer[originalLength] = '\t';
                        originalLength++;
                        char *fullLine = helper->getData(newId, thread_idx);
                        if (columnToTake == -1) {
                            // either append the full line (default mode)
                            size_t fullLineLength = helper->getEntryLen(newId);
                            // Appending join database entry to query database entry
                            memcpy(lineBuffer + originalLength, fullLine, fullLineLength);
                        } else if (*fullLine != '\0') {
                            // or a specified column
                            std::vector<std::string> splittedLine = Util::split(fullLine, "\t");
                            char *newValue = const_cast<char *>(splittedLine[columnToTake].c_str());
                            size_t valueLength = helper->getEntryLen(newId);
                            // Appending join database entry to query database entry
                            memcpy(lineBuffer + originalLength, newValue, valueLength);
                        }
                        nomatch = 0;
                    } else {
                        nomatch = 1;
                    }
                } else if (mode == BEATS_FIRST) {
                    if (counter == 1) {
                        referenceValue = strtod(columnValue, NULL);
                    } else {
                        double toCompare = strtod(columnValue, NULL);
                        if (compOperator == OP_GEQ) {
                            nomatch = !(toCompare >= referenceValue);
                        } else if (compOperator == OP_LEQ) {
                            nomatch = !(toCompare <= referenceValue);
                        } else if (compOperator == OP_EQ) {
                            nomatch = !(toCompare == referenceValue);
                        } else if (compOperator == OP_IN_P) {
                            nomatch = !(toCompare >= (referenceValue * par.compValue));
                        } else if (compOperator == OP_OUT_P) {
                            nomatch = !(toCompare <= (referenceValue * par.compValue));
                        } else if (compOperator == OP_EQ_P) {
                            nomatch = !(toCompare == referenceValue * par.compValue);
                        } else {
                            nomatch = 0;
                        }
                    }
                } else if (mode == FILE_FILTERING) {
                    std::string toSearch(columnValue);
                    std::vector<std::string>::iterator it = std::upper_bound(filter.begin(), filter.end(), toSearch, compareString());
                    if (it != filter.end() && toSearch.compare(*it) == 0) {
                        // Found in filter
                        if (positiveFiltering) {
                            nomatch = 0;
                        } else {
                            nomatch = 1;
                        }
                    } else {
                        // not found in the filter
                        if (positiveFiltering) {
                            nomatch = 1;
                        } else {
                            nomatch = 0;
                        }
                    }
                } else if (mode == FILE_MAPPING) {
                    std::string toSearch(columnValue);
                    std::vector<std::pair<std::string, std::string>>::iterator it
                        = std::lower_bound(mapping.begin(), mapping.end(), toSearch, compareToFirstString());

                    // by default, do NOT add to the output
                    nomatch = 1;

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
                    while (it != mapping.end() && toSearch.compare(it->first) == 0) {
                        nomatch = 0;

                        // copy the previous columns
                        memcpy(newLineBuffer + newLineBufferIndex, lineBuffer, columnPointer[column - 1] - columnPointer[0]);
                        newLineBufferIndex += columnPointer[column - 1] - columnPointer[0];

                        // map the current column value
                        memcpy(newLineBuffer + newLineBufferIndex, (it->second).c_str(), (it->second).length());
                        newLineBufferIndex += (it->second).length();

                        // copy the next columns
                        if (foundElements > column) {
                            memcpy(newLineBuffer + newLineBufferIndex, columnPointer[column - 1] + fieldLength,
                                   endLine - (columnPointer[column - 1] + fieldLength));
                            newLineBufferIndex += endLine - (columnPointer[column - 1] + fieldLength);
                        } else {
                            newLineBuffer[newLineBufferIndex++] = '\n';
                        }
                        newLineBuffer[newLineBufferIndex] = '\0';

                        ++it;
                    }
                    if (nomatch == 0) {
                        memcpy(lineBuffer, newLineBuffer, newLineBufferIndex + 1);
                    }
                } else if (mode == SORT_ENTRIES) {
                    toSort.emplace_back(std::strtod(columnValue, NULL), lineBuffer);
                    // do not put anything in the output buffer
                    nomatch = 1;
                } else {
                    // Unknown filtering mode, keep all entries
                    nomatch = 0;
                }

                if (addSelfMatch) {
                    nomatch = 0;
                }

                if (nomatch == false) {
                    if (trimToOneColumn) {
                        buffer.append(columnValue);
                    } else {
                        buffer.append(lineBuffer);
                    }

                    if (buffer.back() != '\n') {
                        buffer.append(1, '\n');
                    }
                }
                data = Util::skipLine(data);
            }

            if (mode == SORT_ENTRIES) {
                if (par.sortEntries == INCREASING) {
                    std::stable_sort(toSort.begin(), toSort.end(), compareFirstEntry());
                } else if (par.sortEntries == DECREASING) {
                    std::stable_sort(toSort.begin(), toSort.end(), compareFirstEntryDecreasing());
                } else if (par.sortEntries == SHUFFLE) {
                    std::shuffle(toSort.begin(), toSort.end(), urng);
                }

                for (size_t i = 0; i < toSort.size(); i++) {
                    buffer.append(toSort[i].second);
                    if (buffer.back() != '\n') {
                        buffer.append(1, '\n');
                    }
                }
                toSort.clear();
            }
            writer.writeData(buffer.c_str(), buffer.length(), queryKey, thread_idx);
            buffer.clear();
        }

        if (parser != NULL) {
            delete parser;
        }

        delete[] lineBuffer;
        delete[] columnValue;
        delete[] columnPointer;
        delete[] newLineBuffer;
    }
    writer.close();
    reader.close();

    if (helper != NULL) {
        helper->close();
        delete helper;
    }

    if (mode == REGEX_FILTERING) {
        regfree(&regex);
    }

    return EXIT_SUCCESS;
}
