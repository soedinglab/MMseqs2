#include "Parameters.h"
#include "DBWriter.h"
#include "Util.h"
#include "FileUtil.h"
#include "Debug.h"
#include "UniprotKB.h"

#ifdef HAVE_ZLIB
#include "gzstream.h"
#endif

#include <fstream>
#include <set>

std::string getPrimaryAccession(const std::string &accession) {
    size_t end = accession.find_first_of(';');
    if(UNLIKELY(end == std::string::npos)) {
        Debug(Debug::ERROR) << "Could not extract primary accession!\n";
        EXIT(EXIT_FAILURE);
    }
    return accession.substr(0, end);
}

std::vector<unsigned int> getEnabledColumns(const std::string& columns,
                                            const std::string* columnNames, size_t columnCount) {
    std::vector<std::string> kbColumns = Util::split(columns, ",");
    std::set<unsigned int> enabledColumns;
    for (std::vector<std::string>::const_iterator it = kbColumns.begin(); it != kbColumns.end(); ++it) {
        if (Util::isNumber(*it)) {
            char* rest;
            unsigned int col = static_cast<unsigned int>(strtoul((*it).c_str(), &rest, 10));
            if ((rest != (*it).c_str() && *rest != '\0') || errno == ERANGE) {
                Debug(Debug::ERROR) << "Invalid selected column: " << (*it) << "!\n";
                EXIT(EXIT_FAILURE);
            }

            if (col >= columnCount) {
                Debug(Debug::ERROR) << "Invalid selected column: " << col << "!\n";
                EXIT(EXIT_FAILURE);
            }
            enabledColumns.insert(col);
        } else {
            for (size_t i = 0; i < columnCount; ++i) {
                if (columnNames[i] == (*it)) {
                    enabledColumns.emplace(i);
                    break;
                }
            }
        }
    }

    return std::vector<unsigned int>(enabledColumns.begin(), enabledColumns.end());
}

void setConvertKbDefaults(Parameters* par, unsigned int maxColumns) {
    std::ostringstream ss;
    for (unsigned int i = 0; i < maxColumns - 1; ++i) {
        ss << i << ",";
    }
    ss << maxColumns - 1;

    par->kbColumns = ss.str();
}

int convertkb(int argc, const char **argv, const Command& command) {
    UniprotKB kb;
    size_t columns = static_cast<unsigned int>(kb.getColumnCount());

    Parameters& par = Parameters::getInstance();
    setConvertKbDefaults(&par, columns);
    par.parseParameters(argc, argv, command, 2);

    std::vector<unsigned int> enabledColumns = getEnabledColumns(par.kbColumns, kb.columnNames, kb.getColumnCount());

    std::istream *kbIn;
    if (Util::endsWith(".gz", par.db1)) {
#ifdef HAVE_ZLIB
        kbIn = new igzstream(par.db1.c_str());
#else
        Debug(Debug::ERROR) << "MMseqs was not compiled with zlib support. Can not read compressed input!\n";
        EXIT(EXIT_FAILURE);
#endif
    } else {
        kbIn = new std::ifstream(par.db1);
    }


    if (kbIn->fail()) {
        Debug(Debug::ERROR) << "File " << par.db1 << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    DBWriter **writers = new DBWriter*[columns];
    for (std::vector<unsigned int>::const_iterator it = enabledColumns.begin(); it != enabledColumns.end(); ++it) {
        std::string dataFile = par.db2 + "_" + kb.columnNames[*it];
        std::string indexFile = par.db2 + "_" + kb.columnNames[*it] + ".index";
        writers[*it] = new DBWriter(dataFile.c_str(), indexFile.c_str(), 1);
        writers[*it]->open();
    }

    std::map<std::string, unsigned int> mapping = Util::readLookupReverse(par.mappingFile, true);
    const bool doMapping = mapping.size() > 0;

    std::string lookupFile = par.db2;
    lookupFile.append(".lookup");
    std::ofstream* lookupStream = NULL;
    if (!doMapping) {
        lookupStream = new std::ofstream(lookupFile);

        if(lookupStream->fail()) {
            Debug(Debug::ERROR) << "Could not open " << lookupFile << " for writing.";
            EXIT(EXIT_FAILURE);
        }
    }


    Debug(Debug::INFO) << "Extracting data...\n";
    std::string line;
    unsigned int i = 0;
    while (std::getline(*kbIn, line)) {
        if (line.length() < 2) {
            Debug(Debug::WARNING) << "Invalid line" << "\n";
            continue;
        }

        if (kb.readLine(line.c_str())) {
            Debug::printProgress(i);
            std::string accession = getPrimaryAccession(kb.getColumn(UniprotKB::COL_KB_AC));

            for (std::vector<unsigned int>::const_iterator it = enabledColumns.begin(); it != enabledColumns.end(); ++it) {
                std::string column = kb.getColumn(*it);

                unsigned int key = i;
                if (doMapping) {
                    std::map<std::string, unsigned int>::iterator it = mapping.find(accession);
                    if (it == mapping.end()) {
                        Debug(Debug::WARNING) << "Could not find accession " << accession << " in mapping!\n";
                        continue;
                    }
                    key = (*it).second;
                }

                writers[*it]->writeData(column.c_str(), column.length(), key);
            }

            if (!doMapping) {
                *lookupStream << i << "\t" << accession << "\n";
            }

            i++;
        }
    }
    delete kbIn;
    if (!doMapping) {
        lookupStream->close();
        delete lookupStream;
    }
    Debug(Debug::INFO) << "\n";

    for (std::vector<unsigned int>::const_iterator it = enabledColumns.begin(); it != enabledColumns.end(); ++it) {
        writers[*it]->close();
        delete writers[*it];
    }
    delete[] writers;

    return EXIT_SUCCESS;
}
