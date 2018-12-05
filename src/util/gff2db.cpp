/*
 * gff2db
 * written by Milot Mirdita <milot@mirdita.de>
 */

#include <string>
#include <fstream>
#include <sstream>
#include <climits>

#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

int gff2db(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 3);

    DBReader<std::string> reader(par.db2.c_str(), par.db2Index.c_str(), 1, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    DBReader<std::string> headerReader(par.hdr2.c_str(), par.hdr2Index.c_str(), 1, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);

    reader.open(DBReader<std::string>::NOSORT);
    headerReader.open(DBReader<std::string>::NOSORT);

    DBWriter out_writer(par.db3.c_str(), par.db3Index.c_str(), 1, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    DBWriter out_hdr_writer(par.hdr3.c_str(), par.hdr3Index.c_str(), 1, par.compressed, Parameters::DBTYPE_GENERIC_DB);

    out_writer.open();
    out_hdr_writer.open();

    bool shouldCompareType = par.gffType.length() > 0;

    unsigned int entries_num = 0;

    std::ifstream  file_in(par.db1);
    std::string    gff_line;
    while(std::getline(file_in, gff_line)) {
        entries_num++;

        // line is a comment
        if(gff_line.find_first_of("#", 0, 1) != std::string::npos) {
            continue;
        }

        std::vector<std::string> fields = Util::split(gff_line, "\t");

        // gff has always 9 fields
        if(fields.size() != 9) {
            Debug(Debug::WARNING) << "Invalid GFF format in line " << entries_num << "!";
            continue;
        }

        std::string name(fields[0]);

        std::string type(fields[2]);

        if(shouldCompareType && type.compare(par.gffType) != 0) {
            continue;
        }

        char* rest;
        size_t start = strtoull(fields[3].c_str(), &rest, 10);
        if ((rest != fields[3].c_str() && *rest != '\0') || errno == ERANGE) {
            Debug(Debug::WARNING) << "Invalid start position format in line " << entries_num << "!\n";
            continue;
        }
        size_t end = strtoull(fields[4].c_str(), &rest, 10);
        if ((rest != fields[4].c_str() && *rest != '\0') || errno == ERANGE) {
            Debug(Debug::WARNING) << "Invalid end position format in line " << entries_num << "!\n";
            continue;
        }

        if (start == end) {
            Debug(Debug::WARNING) << "Invalid sequence length in line " << entries_num << "!\n";
            continue;
        }

        size_t length = end - start;

        size_t headerId = headerReader.getId(name);

        if(headerId == UINT_MAX) {
            Debug(Debug::ERROR) << "GFF entry not found in database: " << name << "!\n";
            return EXIT_FAILURE;
        }

        char* header = headerReader.getData(headerId, 0);
        size_t headerLength = headerReader.getSeqLens(headerId);

        char* body = reader.getDataByDBKey(name, 0);

        if(!header || !body) {
            Debug(Debug::ERROR) << "GFF entry not found in database: " << name << "!\n";
            return EXIT_FAILURE;
        }


        // header
        char* buffer = new char[headerLength + 128];
        if(shouldCompareType) {
            snprintf(buffer, headerLength + 128, "%s %s:%zu-%zu\n", header, type.c_str(), start, end);
        } else {
            snprintf(buffer, headerLength + 128, "%s %zu-%zu\n", header, start, end);
        }

        unsigned int id = par.identifierOffset + entries_num;

        // hack: header contains a new line, lets just overwrite the new line with a space
        buffer[headerLength - 2] = ' ';
        out_hdr_writer.writeData(buffer, strlen(buffer), id);
        delete[] buffer;

        // sequence
        char* bodyBuffer = new char[length + 1];
        strncpy(bodyBuffer, body + start, length);
        bodyBuffer[length] = '\n';
        out_writer.writeData(bodyBuffer, length + 1, id);
        delete[] bodyBuffer;
    }
    file_in.close();

    out_hdr_writer.close();
    out_writer.close();

    return EXIT_SUCCESS;
}



