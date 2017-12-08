#include <cstdio>
#include "Parameters.h"

#include "DBReader.h"
#include "Debug.h"
#include "Util.h"

int createtsv(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 3, true, Parameters::PARSE_VARIADIC);

    const bool hasTargetDB = par.filenames.size() > 3;
    size_t i = 0;

    const std::string& query = par.filenames[i++];
    const std::string& target = hasTargetDB ? par.filenames[i++] : "";
    const std::string& result = par.filenames[i++];
    const std::string& tsv = par.filenames[i++];


    Debug(Debug::INFO) << "Query file is " << query << "\n";
    DBReader<unsigned int> qHeader(std::string(query + "_h").c_str(), std::string(query + "_h.index").c_str());
    qHeader.open(DBReader<unsigned int>::NOSORT);
    qHeader.readMmapedDataInMemory();

    bool sameDatabase = false;
    DBReader<unsigned int> *tHeader = NULL;
    if (hasTargetDB) {
        if (query == target) {
            sameDatabase = true;
            tHeader = &qHeader;
        } else {
            Debug(Debug::INFO) << "Target file is " << target << "\n";
            tHeader = new DBReader<unsigned int>(std::string(target + "_h").c_str(),
                                                 std::string(target + "_h.index").c_str());
            tHeader->open(DBReader<unsigned int>::NOSORT);
            tHeader->readMmapedDataInMemory();
        }
    }

    Debug(Debug::INFO) << "Data file is " << result << "\n";
    DBReader<unsigned int> reader(result.c_str(), std::string(result + ".index").c_str());
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    FILE *file = fopen(tsv.c_str(), "w");

    Debug(Debug::INFO) << "Start writing file to " << tsv << "\n";
    char *dbKey = new char[par.maxSeqLen + 1];
    for (size_t i = 0; i < reader.getSize(); i++) {
        unsigned int queryKey = reader.getDbKey(i);
        std::string queryHeader = Util::parseFastaHeader(qHeader.getDataByDBKey(queryKey));

        char *data = reader.getData(i);
        size_t entryIndex = 0;
        while (*data != '\0') {
            Util::parseKey(data, dbKey);
            size_t keyLen = strlen(dbKey);

            std::string targetAccession;
            if (tHeader != NULL) {
                unsigned int targetKey = (unsigned int) strtoul(dbKey, NULL, 10);
                targetAccession = Util::parseFastaHeader(tHeader->getDataByDBKey(targetKey));
            } else {
                targetAccession = dbKey;
            }

            if (par.firstSeqRepr && !entryIndex) {
                queryHeader = targetAccession;
            }

            char *nextLine = Util::skipLine(data);

            // write to file
            fwrite(queryHeader.c_str(), sizeof(char), queryHeader.length(), file);
            fwrite("\t", sizeof(char), 1, file);
            fwrite(targetAccession.c_str(), sizeof(char), targetAccession.length(), file);
            fwrite("\t", sizeof(char), 1, file);
            fwrite(data + keyLen, sizeof(char), (nextLine - (data + keyLen)) - 1, file);
            fwrite("\n", sizeof(char), 1, file);

            data = nextLine;
            entryIndex++;
        }
    }
    delete[] dbKey;

    Debug(Debug::INFO) << "Done.\n";

    fclose(file);
    if (sameDatabase == false && tHeader != NULL) {
        tHeader->close();
        delete tHeader;
    }

    qHeader.close();
    reader.close();

    return EXIT_SUCCESS;
}
