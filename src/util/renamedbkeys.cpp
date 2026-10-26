#include "Parameters.h"
#include "FileUtil.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "FastSort.h"

#include <climits>

static bool compareToFirst(const std::pair<DBKeyType, DBKeyType>& lhs, const std::pair<DBKeyType, DBKeyType>& rhs){
    return (lhs.first <= rhs.first);
}

void copyEntry(DBKeyType oldKey, DBKeyType newKey, DBReader<DBKeyType>& reader, DBWriter& writer, bool isCompressed, int subDbMode) {
    const size_t id = reader.getId(oldKey);
    if (id == DB_ENTRY_NOT_FOUND) {
        Debug(Debug::ERROR) << "Key " << oldKey << " not found in database\n";
        EXIT(EXIT_FAILURE);
    }
    if (subDbMode == Parameters::SUBDB_MODE_SOFT) {
        writer.writeIndexEntry(newKey, reader.getOffset(id), reader.getEntryLen(id), 0);
    } else {
        char *data = reader.getDataUncompressed(id);
        size_t originalLength = reader.getEntryLen(id);
        size_t entryLength = std::max(originalLength, static_cast<size_t>(1)) - 1;

        if (isCompressed) {
            // copy also the null byte since it contains the information if compressed or not
            entryLength = *(reinterpret_cast<unsigned int *>(data)) + sizeof(unsigned int) + 1;
            writer.writeData(data, entryLength, newKey, 0, false, false);
        } else {
            writer.writeData(data, entryLength, newKey, 0, true, false);
        }
        // do not write null byte since
        writer.writeIndexEntry(newKey, writer.getStart(0), originalLength, 0);
    }
}

int renamedbkeys(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    FILE *orderFile = NULL;
    if (FileUtil::fileExists(par.db1.c_str())) {
        orderFile = fopen(par.db1.c_str(), "r");
    } else {
        Debug(Debug::ERROR) << "File " << par.db1 << " does not exist.\n";
        EXIT(EXIT_FAILURE);
    }

    FILE* newLookupFile = NULL;
    unsigned int mode = DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA;
    if (FileUtil::fileExists((par.db2 + ".lookup").c_str())) {
        mode |= DBReader<DBKeyType>::USE_LOOKUP;
        newLookupFile = FileUtil::openAndDelete((par.db3 + ".lookup").c_str(), "w");
    }
    DBReader<DBKeyType> reader(par.db2.c_str(), par.db2Index.c_str(), 1, mode);
    reader.open(DBReader<DBKeyType>::NOSORT);
    const bool isCompressed = reader.isCompressed();

    FILE* newMappingFile = NULL;
    std::vector<std::pair<DBKeyType, DBKeyType>> mapping;
    std::vector<std::pair<DBKeyType, DBKeyType>> newMapping;
    if (FileUtil::fileExists((par.db2 + "_mapping").c_str())) {
        mapping.reserve(reader.getSize());
        newMapping.reserve(reader.getSize());
        bool isSorted = Util::readMappingDBKey(par.db2 + "_mapping", mapping);
        if (isSorted == false) {
            std::stable_sort(mapping.begin(), mapping.end(), compareToFirst);
        }
        newMappingFile = FileUtil::openAndDelete((par.db3 + "_mapping").c_str(), "w");
    }

    bool isHeaderCompressed = false;
    DBReader<DBKeyType>* headerReader = NULL;
    if (FileUtil::fileExists(par.hdr2dbtype.c_str())) {
        headerReader = new DBReader<DBKeyType>(par.hdr2.c_str(), par.hdr2Index.c_str(), 1, DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA);
        headerReader->open(DBReader<DBKeyType>::NOSORT);
        isHeaderCompressed = headerReader->isCompressed();
    }

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), 1, 0, Parameters::DBTYPE_OMIT_FILE);
    writer.open();

    DBWriter* headerWriter = NULL;
    if (headerReader != NULL) {
        headerWriter = new DBWriter(par.hdr3.c_str(), par.hdr3Index.c_str(), 1, 0, Parameters::DBTYPE_OMIT_FILE);
        headerWriter->open();
    }

    DBReader<DBKeyType>::LookupEntry* lookup = NULL;
    std::vector<DBReader<DBKeyType>::LookupEntry> newLookup;
    if (newLookupFile != NULL) {
        lookup = reader.getLookup();
        newLookup.reserve(reader.getLookupSize());
    }

    char *line = NULL;
    size_t len = 0;
    const char *fields[2];
    // getline malloc/reallocs automatically
    while (getline(&line, &len, orderFile) != -1) {
        const size_t columns = Util::getWordsOfLine(line, fields, 2);
        if (columns < 2) {
            Debug(Debug::WARNING) << "Not enough columns in mapping file\n";
            continue;
        }
        const DBKeyType oldKey = Util::fast_atoi<DBKeyType>(fields[0]);
        const DBKeyType newKey = Util::fast_atoi<DBKeyType>(fields[1]);

        copyEntry(oldKey, newKey, reader, writer, isCompressed, par.subDbMode);
        if (lookup != NULL) {
            size_t lookupId = reader.getLookupIdByKey(oldKey);
            DBReader<DBKeyType>::LookupEntry entry = lookup[lookupId];
            entry.id = newKey;
            newLookup.emplace_back(entry);
        }

        if (mapping.size() > 0) {
            std::pair<DBKeyType, DBKeyType> val;
            val.first = oldKey;
            std::vector<std::pair<DBKeyType, DBKeyType>>::iterator mappingIt;
            mappingIt = std::upper_bound(mapping.begin(), mapping.end(), val, compareToFirst);
            if (mappingIt != mapping.end() && mappingIt->first == val.first) {
                val.first = newKey;
                val.second = mappingIt->second;
                newMapping.emplace_back(val);
            }
        }

        if (headerReader != NULL && headerWriter != NULL) {
            copyEntry(oldKey, newKey, *headerReader, *headerWriter, isHeaderCompressed, par.subDbMode);
        }
    }
    // merge any kind of sequence database
    writer.close(headerWriter != NULL);
    DBWriter::writeDbtypeFile(par.db3.c_str(), reader.getDbtype(), isCompressed);
    if (par.subDbMode == Parameters::SUBDB_MODE_SOFT) {
        DBReader<DBKeyType>::softlinkDb(par.db2, par.db3, DBFiles::DATA);
    }
    if (newMappingFile != NULL) {
        SORT_PARALLEL(newMapping.begin(), newMapping.end(), compareToFirst);
        std::string buffer;
        for (size_t i = 0; i < newMapping.size(); ++i) {
            buffer.append(SSTR(newMapping[i].first));
            buffer.append(1, '\t');
            buffer.append(SSTR(newMapping[i].second));
            buffer.append(1, '\n');
            fwrite(buffer.c_str(), sizeof(char), buffer.size(), newMappingFile);
            buffer.clear();
        }
        fclose(newMappingFile);
    }

    if (newLookupFile != NULL) {
        SORT_PARALLEL(newLookup.begin(), newLookup.end(), DBReader<DBKeyType>::LookupEntry::compareById);
        std::string lookupBuffer;
        lookupBuffer.reserve(2048);
        for (size_t i = 0; i < newLookup.size(); ++i) {
            reader.lookupEntryToBuffer(lookupBuffer, newLookup[i]);
            fwrite(lookupBuffer.c_str(), sizeof(char), lookupBuffer.size(), newLookupFile);
            lookupBuffer.clear();
        }
        fclose(newLookupFile);
    }

    if (headerWriter != NULL) {
        headerWriter->close(true);
        delete headerWriter;
        DBWriter::writeDbtypeFile(par.hdr3.c_str(), headerReader->getDbtype(), isHeaderCompressed);
        if (par.subDbMode == Parameters::SUBDB_MODE_SOFT) {
            DBReader<DBKeyType>::softlinkDb(par.db2, par.db3, DBFiles::HEADER);
        }
    }
    if (par.subDbMode == Parameters::SUBDB_MODE_SOFT) {
        DBReader<DBKeyType>::softlinkDb(par.db2, par.db3, (DBFiles::Files) (DBFiles::SOURCE | DBFiles::TAX_MERGED | DBFiles::TAX_NAMES | DBFiles::TAX_NODES | DBFiles::TAX_BINARY));
    } else {
        DBReader<DBKeyType>::copyDb(par.db2, par.db3, (DBFiles::Files) (DBFiles::SOURCE | DBFiles::TAX_MERGED | DBFiles::TAX_NAMES | DBFiles::TAX_NODES | DBFiles::TAX_BINARY));
    }

    free(line);
    if (headerReader != NULL) {
        headerReader->close();
        delete headerReader;
    }
    reader.close();
    fclose(orderFile);

    return EXIT_SUCCESS;
}
