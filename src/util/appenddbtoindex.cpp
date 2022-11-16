#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

int appenddbtoindex(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    std::string outDb = par.filenames.back();
    par.filenames.pop_back();

    // read in database keys for the new database entries and validate that we have enough
    std::vector<unsigned int> keys;
    {
        std::vector<std::string> ids = Util::split(par.idList, ",");
        keys.reserve(ids.size());
        for (size_t i = 0; i < ids.size(); ++i) {
            char *rest;
            errno = 0;
            unsigned int key = strtoul(ids[i].c_str(), &rest, 10);
            if ((rest != ids[i].c_str() && *rest != '\0') || errno == ERANGE) {
                Debug(Debug::ERROR) << "Could not read key " << ids[i] << "\n";
                return EXIT_FAILURE;
            }
            keys.emplace_back(key);
        }
        if (keys.size() != par.filenames.size()) {
            Debug(Debug::ERROR) << "Same number of databases and keys are needed\n";
            return EXIT_FAILURE;
        }
        // fail early if duplicates are found
        std::vector<unsigned int> check(keys.begin(), keys.end());
        std::sort(check.begin(), check.end());
        for (size_t i = 1; i < check.size(); ++i) {
            if (check[i - 1] == check[i] || (check[i - 1] + 1) == check[i]) {
                Debug(Debug::ERROR) << "Duplicate ID given. Each database takes two consecutive IDs.\n";
                return EXIT_FAILURE;
            }
        }
    }

    // if we have a split database, make one new split where we append the new files
    FILE* outDataHandle = NULL;
    {
        std::string checkName = outDb + ".0";
        size_t cnt = 0;
        while (FileUtil::fileExists(checkName.c_str()) == true) {
            cnt++;
            checkName = outDb + "." + SSTR(cnt);
        }
        if (cnt == 0) {
            outDataHandle = FileUtil::openFileOrDie(outDb.c_str(), "ab", true);
        } else {
            outDataHandle = FileUtil::openFileOrDie(checkName.c_str(), "wb", false);
        }
    }

    std::string outIndexName = outDb + ".index";
    size_t offset = 0;
    {
        DBReader<unsigned int> outReader(outDb.c_str(), outIndexName.c_str(), 1, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        outReader.open(DBReader<unsigned int>::NOSORT);
        // validate that given keys dont exist already
        for (size_t i = 0; i < keys.size(); ++i) {
            if (outReader.getId(keys[i]) != UINT_MAX) {
                Debug(Debug::ERROR) << "Key " << keys[i] << " already exists in database\n";
                return EXIT_FAILURE;
            }
            if (outReader.getId(keys[i]+1) != UINT_MAX) {
                Debug(Debug::ERROR) << "Key " << (keys[i]+1) << " already exists in database\n";
                return EXIT_FAILURE;
            }
        }
        offset = outReader.getTotalDataSize();
        outReader.close();
    }

    const char nullbyte = '\0';
    char buffer[8192];
    FILE* outIndexHandle = FileUtil::openFileOrDie(outIndexName.c_str(), "a", true);
    for (size_t i = 0; i < par.filenames.size(); ++i) {
        const unsigned int key = keys[i];
        const std::string& inDb = par.filenames[i];
        const std::string inIndexName = inDb + ".index";

        DBReader<unsigned int> reader(inDb.c_str(), inIndexName.c_str(), 1, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        reader.open(DBReader<unsigned int>::HARDNOSORT);

        char* data = DBReader<unsigned int>::serialize(reader);
        size_t inSize = DBReader<unsigned int>::indexMemorySize(reader);
        size_t written = fwrite(data, 1, inSize, outDataHandle);
        free(data);
        if (written != inSize) {
            Debug(Debug::ERROR) << "Cannot write to data file " << outDb << "\n";
            EXIT(EXIT_FAILURE);
        }
        written = fwrite(&nullbyte, sizeof(char), 1, outDataHandle);
        if (written != 1) {
            Debug(Debug::ERROR) << "Cannot write to data file " << outDb << "\n";
            EXIT(EXIT_FAILURE);
        }
        inSize += 1;
        size_t len = DBWriter::indexToBuffer(buffer, key, offset, inSize);
        written = fwrite(buffer, sizeof(char), len, outIndexHandle);
        if (written != len) {
            Debug(Debug::ERROR) << "Cannot write to index file " << outIndexName << "\n";
            EXIT(EXIT_FAILURE);
        }
        offset += inSize;

        inSize = reader.getTotalDataSize();
        for (size_t idx = 0; idx < reader.getDataFileCnt(); idx++) {
            char* data = reader.getDataForFile(idx);
            size_t size = reader.getDataSizeForFile(idx);
            written = fwrite(data, sizeof(char), size, outDataHandle);
            if (written != size) {
                Debug(Debug::ERROR) << "Cannot write to data file " << outDb << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        reader.close();

        written = fwrite(&nullbyte, sizeof(char), 1, outDataHandle);
        if (written != 1) {
            Debug(Debug::ERROR) << "Cannot write to data file " << outDb << "\n";
            EXIT(EXIT_FAILURE);
        }
        inSize += 1;
        len = DBWriter::indexToBuffer(buffer, key + 1, offset, inSize);
        written = fwrite(buffer, sizeof(char), len, outIndexHandle);
        if (written != len) {
            Debug(Debug::ERROR) << "Cannot write to index file " << outIndexName << "\n";
            EXIT(EXIT_FAILURE);
        }
        offset += inSize;
    }

    if (fclose(outDataHandle) != 0) {
        Debug(Debug::ERROR) << "Cannot close data file " << outDb << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fclose(outIndexHandle) != 0) {
        Debug(Debug::ERROR) << "Cannot close index file " << outIndexName << "\n";
        EXIT(EXIT_FAILURE);
    }

    DBWriter::sortIndex(outIndexName.c_str(), outIndexName.c_str(), false);

    return EXIT_SUCCESS;
}
