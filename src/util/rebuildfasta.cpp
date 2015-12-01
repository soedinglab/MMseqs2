/*
 * rebuildfasta
 * written by Milot Mirdita <milot@mirdita.de>
 */

#define _GNU_SOURCE 1
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include <cstring>
#include <cstdio>
#include <Parameters.h>

#include "DBReader.h"
#include "Debug.h"
#include "Util.h"

const char header_start[] = {'>'};
const char newline[] = {'\n'};

int rebuildfasta(int argc, const char * argv[])
{
    std::string usage("Converts an mmseqs ffindex database back to a fasta file. \n");
    usage.append("Written by Milot Mirdita <milot@mirdita.de>.\n\n");
    usage.append("USAGE: <dbIn> <fastaOut>\n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.rebuildfasta, 2);

    Debug::setDebugLevel(par.verbosity);

    std::string data_filename = par.db1;
    std::string index_filename = par.db1Index;

    std::string data_filename_hdr(data_filename);
    data_filename_hdr.append("_h");

    std::string index_filename_hdr(data_filename);
    index_filename_hdr.append("_h.index");

    DBReader db(data_filename.c_str(), index_filename.c_str());
    db.open(DBReader::NOSORT);

    DBReader db_header(data_filename_hdr.c_str(), index_filename_hdr.c_str());
    db_header.open(DBReader::NOSORT);

    FILE *fastaFP =  Util::openFileOrDie(par.db2.c_str(), "w", false);

    DBReader& from = db;
    if(par.useHeaderFile) {
        from = db_header;
    }

    Debug(Debug::INFO) << "Start writing file to " << par.db2 << "\n";
    for(size_t i = 0; i < from.getSize(); i++){
        std::string key = from.getDbKey(i);

        const char* header_data = db_header.getDataByDBKey(key.c_str());

        fwrite(header_start, sizeof(char), 1, fastaFP);
        fwrite(header_data, sizeof(char), strlen(header_data) - 1, fastaFP);
        fwrite(newline, sizeof(char), 1, fastaFP);

        const char* body_data = db.getDataByDBKey(key.c_str());
        fwrite(body_data, sizeof(char), strlen(body_data) - 1, fastaFP);
        fwrite(newline, sizeof(char), 1, fastaFP);

    }

    fclose(fastaFP);
    db_header.close();
    db.close();

    return EXIT_SUCCESS;
}
