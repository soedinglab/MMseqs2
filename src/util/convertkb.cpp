#include "Parameters.h"
#include "DBWriter.h"
#include "Util.h"
#include "FileUtil.h"
#include "Debug.h"
#include "UniprotKB.h"

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

std::vector<unsigned int> getEnabledColumns(const std::string& columns, size_t maxColumn) {
    std::vector<std::string> kbColumns = Util::split(columns, ",");
    std::set<unsigned int> enabledColumns;
    for (std::vector<std::string>::const_iterator it = kbColumns.begin(); it != kbColumns.end(); ++it) {
        char* rest;
        unsigned int col = static_cast<unsigned int>(strtoul((*it).c_str(), &rest, 10));
        if ((rest != (*it).c_str() && *rest != '\0') || errno == ERANGE) {
            Debug(Debug::ERROR) << "Invalid selected column: " << (*it) << "!\n";
            EXIT(EXIT_FAILURE);
        }

        if (col >= columns) {
            Debug(Debug::ERROR) << "Invalid selected column: " << col << "!\n";
            EXIT(EXIT_FAILURE);
        }
        enabledColumns.insert(col);
    }

    return std::vector<unsigned int>(enabledColumns.begin(), enabledColumns.end());
}

void setConvertKbDefaults(Parameters* par, size_t maxColumns) {
    std::ostringstream ss;
    for (int i = 0; i < maxColumns - 1; ++i) {
        ss << i << ",";
    }
    ss << maxColumns - 1;

    par->kbColumns = ss.str();
}

int convertkb(int argn, const char **argv) {
    std::string usage;
    usage.append("Turns an UniprotKB file into multiple ffindex database for every KB column.\n");
    usage.append("USAGE: <uniprotKB> <outDbPrefix>\n");
    usage.append("\nDesigned and implemented by Milot Mirdita <milot@mirdita.de>.\n");

    UniprotKB kb;
    size_t columns = kb.getColumnCount();

    Parameters par;
    setConvertKbDefaults(&par, columns);
    par.parseParameters(argn, argv, usage, par.convertkb, 2);

    std::vector<unsigned int> enabledColumns = getEnabledColumns(par.kbColumns, columns);


    std::ifstream kbIn(par.db1);
    if (kbIn.fail()) {
        Debug(Debug::ERROR) << "File " << par.db1 << " not found!\n";
        EXIT(EXIT_FAILURE);
    }


    DBWriter **writers = new DBWriter*[columns];
    for (size_t i = 0; i < columns; ++i) {
        std::string dataFile = par.db2 + "_" + kb.columnNames[i];
        std::string indexFile = par.db2 + "_" + kb.columnNames[i] + ".index";
        writers[i] = new DBWriter(dataFile.c_str(), indexFile.c_str(), 1);
        writers[i]->open();
    }

    std::string line;
    while (std::getline(kbIn, line)) {
        if (line.length() < 2) {
            Debug(Debug::WARNING) << "Invalid line" << "\n";
            continue;
        }

        if (kb.readLine(line.c_str())) {
            std::string accession = getPrimaryAccession(kb.getColumn(UniprotKB::COL_KB_AC));
            for (std::vector<unsigned int>::const_iterator it = enabledColumns.begin(); it != enabledColumns.end(); ++it) {
                std::string column = kb.getColumn(*it);
                writers[*it]->write(column.c_str(), column.length(), accession.c_str());
            }
        }
    }

    for (size_t i = 0; i < columns; ++i) {
        writers[i]->close();
        delete writers[i];
    }
    delete[] writers;

    return EXIT_SUCCESS;
}
