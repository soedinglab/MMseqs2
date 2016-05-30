#include "Parameters.h"
#include "DBWriter.h"
#include "Util.h"
#include "FileUtil.h"
#include "Debug.h"
#include "UniprotKB.h"

#include <fstream>

int convertkb(int argn, const char **argv) {
    std::string usage;
    usage.append("Turns an UniprotKB file into multiple ffindex database for every KB column.\n");
    usage.append("USAGE: <uniprotKB> <outDbPrefix>\n");
    usage.append("\nDesigned and implemented by Milot Mirdita <milot@mirdita.de>.\n");

    Parameters par;
    par.parseParameters(argn, argv, usage, par.onlyverbosity, 2);


    std::ifstream kbIn(par.db1);
    if (kbIn.fail()) {
        Debug(Debug::ERROR) << "File " << par.db1 << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    UniprotKB kb;

    size_t columns = kb.getColumnCount();
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
            std::string accession = kb.getColumn(UniprotKB::COL_KB_ID);
            for (size_t i = 1; i < columns; ++i) {
                std::string column = kb.getColumn(i);
                writers[i]->write(column.c_str(), column.length(), accession.c_str());
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
