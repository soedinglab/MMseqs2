#include "Parameters.h"
#include "FileUtil.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

#include <climits>

int createsubdb(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    bool isIndex = false;
    FILE *orderFile = NULL;
    if (FileUtil::fileExists(par.db1Index.c_str())) {
        orderFile = fopen(par.db1Index.c_str(), "r");
        isIndex = true;
    } else {
        if(FileUtil::fileExists(par.db1.c_str())){
            orderFile = fopen(par.db1.c_str(), "r");
        }else{
            Debug(Debug::ERROR) << "File " << par.db1 << " does not exist.\n";
            EXIT(EXIT_FAILURE);
        }
    }

    const bool lookupMode = par.dbIdMode == Parameters::ID_MODE_LOOKUP;
    int dbMode = DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA;
    if (lookupMode) {
        dbMode |= DBReader<unsigned int>::USE_LOOKUP_REV;
    }
    DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str(), 1, dbMode);
    reader.open(DBReader<unsigned int>::NOSORT);
    const bool isCompressed = reader.isCompressed();

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), 1, 0, Parameters::DBTYPE_OMIT_FILE);
    writer.open();
    // getline reallocs automatic
    char *line = NULL;
    size_t len = 0;
    char dbKey[256];
    unsigned int prevKey = 0;
    bool isOrdered = true;
    char* result;
    char newLine = '\n';
    std::vector<std::string> arr;
    while (getline(&line, &len, orderFile) != -1) {
        Util::parseKey(line, dbKey);
        arr = Util::split(line, "\t");
        unsigned int key;
        if (lookupMode) {
            size_t lookupId = reader.getLookupIdByAccession(dbKey);
            if (lookupId == SIZE_MAX) {
                Debug(Debug::WARNING) << "Could not find name " << dbKey << " in lookup\n";
                continue;
            }
            key = reader.getLookupKey(lookupId);
        } else {
            key = Util::fast_atoi<unsigned int>(dbKey);
        }

        isOrdered &= (prevKey <= key);
        prevKey = key;
        const size_t id = reader.getId(key);
        if (id >= UINT_MAX) {
            Debug(Debug::WARNING) << "Key " << dbKey << " not found in database\n";
            continue;
        }
        if (par.subDbMode == Parameters::SUBDB_MODE_SOFT) {
            writer.writeIndexEntry(key, reader.getOffset(id), reader.getEntryLen(id), 0);
        } else if (isIndex == true || arr.size() == 1) { 
            char* data = reader.getDataUncompressed(id);
            size_t originalLength = reader.getEntryLen(id);
            size_t entryLength = std::max(originalLength, static_cast<size_t>(1)) - 1;

            if (isCompressed) {
                // copy also the null byte since it contains the information if compressed or not
                entryLength = *(reinterpret_cast<unsigned int *>(data)) + sizeof(unsigned int) + 1;
                writer.writeData(data, entryLength, key, 0, false, false);
            } else {
                writer.writeData(data, entryLength, key, 0, true, false);
            }
            // do not write null byte since
            writer.writeIndexEntry(key, writer.getStart(0), originalLength, 0);
        } else {
            char* data = reader.getDataUncompressed(id);
            size_t originalLength = reader.getEntryLen(id);
            size_t entryLength = std::max(originalLength, static_cast<size_t>(1)) - 1;
            int totalLength = 0;
            if (isCompressed) {
                entryLength = *(reinterpret_cast<unsigned int *>(data)) + sizeof(unsigned int) + 1;
            }
            if (arr.size()%2 == 0) {
                Debug(Debug::ERROR) << "Input list not in format\n";
            } else {
                result = new char[entryLength];
                for (int ord = 0 ; ord < int((arr.size()-1)/2); ord ++) {
                    int currLength = std::stoi(arr[ord * 2 + 2])  - std::stoi(arr[ord * 2 + 1]) + 1;
                    strncpy(result + totalLength, data + std::stoi(arr[ord * 2 + 1]), currLength);
                    totalLength += currLength;
                }
                if (isCompressed) {
                    //TODO
                } else {
                    writer.writeData(result, totalLength, key, 0, false, false);
                }
                writer.writeAdd(&newLine, 1, 0);
                writer.writeIndexEntry(key, writer.getStart(0), totalLength + 2, 0);
                delete [] result;
                result = nullptr;
            }
        }
    }
    // merge any kind of sequence database
    const bool shouldMerge = Parameters::isEqualDbtype(reader.getDbtype(), Parameters::DBTYPE_HMM_PROFILE)
                             || Parameters::isEqualDbtype(reader.getDbtype(), Parameters::DBTYPE_AMINO_ACIDS)
                             || Parameters::isEqualDbtype(reader.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES);
    writer.close(shouldMerge, !isOrdered);
    if (par.subDbMode == Parameters::SUBDB_MODE_SOFT) {
        DBReader<unsigned int>::softlinkDb(par.db2, par.db3, DBFiles::DATA);
    }
    DBWriter::writeDbtypeFile(par.db3.c_str(), reader.getDbtype(), isCompressed);
    DBReader<unsigned int>::softlinkDb(par.db2, par.db3, DBFiles::SEQUENCE_ANCILLARY);
    free(line);
    reader.close();
    if (fclose(orderFile) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << par.db1 << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
