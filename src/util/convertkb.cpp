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
    if (UNLIKELY(end == std::string::npos)) {
        Debug(Debug::ERROR) << "Could not extract primary accession!\n";
        EXIT(EXIT_FAILURE);
    }
    return accession.substr(0, end);
}

std::vector<unsigned int> getEnabledColumns(const std::string &columns,
                                            const std::string *columnNames, size_t columnCount) {
    std::vector<std::string> kbColumns = Util::split(columns, ",");
    std::set<unsigned int> enabledColumns;
    for (std::vector<std::string>::const_iterator it = kbColumns.begin(); it != kbColumns.end(); ++it) {
        if (Util::isNumber(*it)) {
            char *rest;
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

void setConvertKbDefaults(Parameters *par, unsigned int maxColumns) {
    std::ostringstream ss;
    for (unsigned int i = 0; i < maxColumns - 1; ++i) {
        ss << i << ",";
    }
    ss << maxColumns - 1;

    par->kbColumns = ss.str();
}

int convertkb(int argc, const char **argv, const Command &command) {
    UniprotKB kb;
    size_t columns = static_cast<unsigned int>(kb.getColumnCount());

    Parameters &par = Parameters::getInstance();
    setConvertKbDefaults(&par, columns);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    std::string outputBase = par.filenames.back();
    par.filenames.pop_back();

    std::vector<unsigned int> enabledColumns = getEnabledColumns(par.kbColumns, kb.columnNames, kb.getColumnCount());

    DBWriter **writers = new DBWriter*[columns];
    for (std::vector<unsigned int>::const_iterator it = enabledColumns.begin(); it != enabledColumns.end(); ++it) {
        std::string dataFile = outputBase + "_" + kb.columnNames[*it];
        std::string indexFile = outputBase + "_" + kb.columnNames[*it] + ".index";
        writers[*it] = new DBWriter(dataFile.c_str(), indexFile.c_str(), 1, par.compressed, Parameters::DBTYPE_GENERIC_DB);
        writers[*it]->open();
    }

    DBReader<unsigned int>* reader = NULL;
    std::ofstream *lookupStream = NULL;

    const bool doMapping = FileUtil::fileExists(par.mappingFile.c_str());
    if (!doMapping) {
        std::string lookupFile = outputBase + ".lookup";
        lookupStream = new std::ofstream(lookupFile);
        if (lookupStream->fail()) {
            Debug(Debug::ERROR) << "Could not open " << lookupFile << " for writing.";
            EXIT(EXIT_FAILURE);
        }
    } else {
        reader = new DBReader<unsigned int>(par.mappingFile.c_str(), par.mappingFile.c_str(), 1, DBReader<unsigned int>::USE_LOOKUP_REV);
        reader->open(DBReader<unsigned int>::NOSORT);
    }

    Debug::Progress progress;
    for (std::vector<std::string>::const_iterator it = par.filenames.begin(); it != par.filenames.end(); ++it) {
        std::istream *kbIn;
        if (Util::endsWith(".gz", *it)) {
#ifdef HAVE_ZLIB
            kbIn = new igzstream((*it).c_str());
#else
            Debug(Debug::ERROR) << "MMseqs2 was not compiled with zlib support. Can not read compressed input\n";
            EXIT(EXIT_FAILURE);
#endif
        } else {
            kbIn = new std::ifstream(*it);
        }

        if (kbIn->fail()) {
            Debug(Debug::ERROR) << "File " << (*it) << " not found\n";
            EXIT(EXIT_FAILURE);
        }

        Debug(Debug::INFO) << "Extracting data from " << (*it) << "\n";
        std::string line;
        unsigned int i = 0;
        while (std::getline(*kbIn, line)) {
            if (line.length() < 2) {
                Debug(Debug::WARNING) << "Invalid entry\n";
                continue;
            }

            if (kb.readLine(line.c_str())) {
                progress.updateProgress();
                std::string accession = getPrimaryAccession(kb.getColumn(UniprotKB::COL_KB_AC));

                for (std::vector<unsigned int>::const_iterator it = enabledColumns.begin();
                     it != enabledColumns.end(); ++it) {
                    std::string column = kb.getColumn(*it);

                    unsigned int key = i;
                    if (doMapping) {
                        size_t lookupId = reader->getLookupIdByAccession(accession);
                        if (lookupId == SIZE_MAX) {
                            Debug(Debug::WARNING) << "Could not find accession " << accession << " in lookup\n";
                            continue;
                        }
                        key = reader->getLookupKey(lookupId);
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
    }

    for (std::vector<unsigned int>::const_iterator it = enabledColumns.begin(); it != enabledColumns.end(); ++it) {
        writers[*it]->close();
        delete writers[*it];
    }
    delete[] writers;

    if (doMapping) {
        reader->close();
        delete reader;
    } else {
        lookupStream->close();
        delete lookupStream;
    }

    return EXIT_SUCCESS;
}
