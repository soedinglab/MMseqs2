#include <string>
#include <fstream>
#include <climits>

#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Parameters.h"
#include "Util.h"

int maskbygff(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 3);

    DBReader<std::string> ffindexReader(par.db2.c_str(), par.db2Index.c_str(),
                                        DBReader<std::string>::USE_DATA | DBReader<std::string>::USE_WRITABLE);
    ffindexReader.open(DBReader<std::string>::NOSORT);

    bool shouldCompareType = par.gffType.length() > 0;

    size_t entries_num = 0;

    std::ifstream  gffFile(par.db1);
    std::string gffLine;
    DBReader<std::string>::Index* index = ffindexReader.getIndex();
    char * data = (char *)ffindexReader.getData();
    while(std::getline(gffFile, gffLine)) {
        entries_num++;

        // line is a comment
        if(gffLine.find_first_of("#", 0, 1) != std::string::npos) {
            continue;
        }

        std::vector<std::string> fields = Util::split(gffLine, "\t");

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

        // gff start and end are supposed to be 1-indexed
        if (end <= start || end == 0 || start == 0) {
            Debug(Debug::WARNING) << "Invalid sequence length in line " << entries_num << "!\n";
            continue;
        }

        // so 0-index them now
        start -= 1;
        end -= 1;

        size_t id = ffindexReader.getId(name);
        if(id == UINT_MAX) {
            Debug(Debug::ERROR) << "GFF entry not found in fasta ffindex: " << name << "!\n";
            return EXIT_FAILURE;
        }

        char* body = data + index[id].offset;
        for (char* i = body + start; i <= body + end; ++i)
        {
            *i = 'X';
        }

    }
    gffFile.close();

    unsigned int* seqLengths = ffindexReader.getSeqLens();

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str());
    writer.open();

    std::string headerFilename(par.db2);
    headerFilename.append("_h");

    std::string headerIndexFilename(par.db2);
    headerIndexFilename.append("_h.index");

    DBReader<std::string> headerReader(headerFilename.c_str(), headerIndexFilename.c_str());
    headerReader.open(DBReader<std::string>::NOSORT);

    DBReader<std::string>::Index* headerIndex = headerReader.getIndex();
    unsigned int* headerLengths = headerReader.getSeqLens();
    char * headerData = (char *)headerReader.getData();

    std::string headerOutFilename(par.db3);
    headerOutFilename.append("_h");

    std::string headerIndexOutFilename(par.db3);
    headerIndexOutFilename.append("_h.index");

    DBWriter headerWriter(headerOutFilename.c_str(), headerIndexOutFilename.c_str());
    headerWriter.open();

    for(size_t i = 0; i < ffindexReader.getSize(); ++i ) {
        std::string id;
        if (par.useHeader) {
            id = index[i].id;
        } else {
            id = SSTR(par.identifierOffset + i);
        }

        // ignore nulls
        writer.writeData(data + index[i].offset, seqLengths[i] - 1, id.c_str());
        headerWriter.writeData(headerData + headerIndex[i].offset, headerLengths[i] - 1, id.c_str());
    }
    headerWriter.close();
    writer.close();

    return EXIT_SUCCESS;
}
