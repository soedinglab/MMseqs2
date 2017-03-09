#include <sstream>
#include <fstream>

#include "Parameters.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

int tsv2db(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str());
    writer.open();

    std::ifstream  tsv(par.db1);
    if (tsv.fail()) {
        Debug(Debug::ERROR) << "File " << par.db1 << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    std::ostringstream ss;
    char keyData[255];
    bool skippedFirst = false;
    std::string lastKey;
    std::string line;
    while(std::getline(tsv, line)) {
        char* current = (char*) line.c_str();
        Util::parseKey(current, keyData);
        std::string key(keyData);

        if (key != lastKey && skippedFirst == true) {
            std::string result = ss.str();
            unsigned int keyId = strtoull(lastKey.c_str(), NULL, 10);
            writer.writeData(result.c_str(), result.length(), keyId);
            ss.str("");
            ss.clear();

        }

        char *restStart = current + key.length();
        restStart = restStart + Util::skipWhitespace(restStart);
        char *restEnd = restStart;
        restEnd = Util::seekToNextEntry(restEnd) - 1;

        std::string rest(restStart, restEnd - restStart);

        skippedFirst = true;
        ss << rest << "\n";
        lastKey = key;
    }

    std::string result = ss.str();
    if (result != "") {
        unsigned int keyId = strtoull(lastKey.c_str(), NULL, 10);
        writer.writeData(result.c_str(), result.length(), keyId);
    }

    writer.close();

    return EXIT_SUCCESS;
}
