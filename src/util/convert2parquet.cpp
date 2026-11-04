/*
 * convert2parquet
 * Writes an MMseqs2 database as parquet: one file, or a directory of part files for large databases.
 * The layout (columns) follows the database type from the .dbtype file, see ParquetDBWriter::dbTypeToFields.
 */

#include "Parameters.h"
#include "DBReader.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "ParquetDBWriter.h"

int convert2parquet(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    // key sorted, so keyorder and the key range planning work on contiguous ids
    DBReader<DBKeyType> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<DBKeyType>::USE_INDEX|DBReader<DBKeyType>::USE_DATA);
    reader.open(DBReader<DBKeyType>::SORT_BY_ID);

    const int dbtype = reader.getDbtype();
    // compressed profiles are binary, their decompressed length is not known from the index
    const bool compressedProfile = Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_HMM_PROFILE)
                                   && reader.isCompressed() == DBReader<DBKeyType>::COMPRESSED;
    if (ParquetDBWriter::isSupportedDbType(dbtype) == false || compressedProfile) {
        Debug(Debug::ERROR) << "DB currently not supported as parquet (" << Parameters::getDbTypeName(dbtype) << ")\n";
        reader.close();
        return EXIT_FAILURE;
    }

    // --parquet-schema: required for generic DBs (unknown line layout), a column selection for the others
    std::vector<ParquetDBWriter::ColumnSpec> schemaColumns;
    const bool isGeneric = Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_GENERIC_DB);
    if (par.parquetSchema.empty() == false) {
        std::string error;
        if (ParquetDBWriter::parseSchema(par.parquetSchema, schemaColumns, error) == false) {
            Debug(Debug::ERROR) << "Invalid --parquet-schema: " << error << "\n";
            reader.close();
            return EXIT_FAILURE;
        }
    } else if (isGeneric) {
        Debug(Debug::ERROR) << "Generic DBs need --parquet-schema '{name:type, ...}' describing the tab separated columns\n";
        reader.close();
        return EXIT_FAILURE;
    }

    const ParquetDBWriter::WriterOptions options = ParquetDBWriter::writerOptionsFromParameters(par);

    // key -> accession (written by createdb and others): <DB>.lookup, or --lookup-db (a .lookup file or a
    // DB name) when the lookup belongs to a DB of another name
    const std::string suffix = ".lookup";
    std::string lookupFile = par.lookupDb.empty() ? par.db1 : par.lookupDb;
    if (lookupFile.size() < suffix.size() || lookupFile.compare(lookupFile.size() - suffix.size(), suffix.size(), suffix) != 0) {
        lookupFile += suffix;
    }
    DBReader<DBKeyType>* lookupReader = NULL;
    if (FileUtil::fileExists(lookupFile.c_str())) {
        // DBReader reads <name>.lookup and always an index; the input DB's index has the keys that are looked up
        const std::string lookupName = lookupFile.substr(0, lookupFile.size() - suffix.size());
        lookupReader = new DBReader<DBKeyType>(lookupName.c_str(), par.db1Index.c_str(), 1, DBReader<DBKeyType>::USE_LOOKUP);
        lookupReader->open(DBReader<DBKeyType>::NOSORT);
        // keys and lookup are both sorted by key: count the keys that have an accession
        const DBReader<DBKeyType>::LookupEntry* lookup = lookupReader->getLookup();
        const size_t lookupSize = lookupReader->getLookupSize();
        size_t matched = 0;
        for (size_t id = 0, l = 0; id < reader.getSize() && l < lookupSize;) {
            const DBKeyType key = reader.getDbKey(id);
            if (lookup[l].id < key) {
                l++;
            } else {
                matched += (lookup[l].id == key);
                id++;
            }
        }
        Debug(Debug::INFO) << matched << " of " << reader.getSize() << " keys have an accession in " << lookupFile << "\n";
        if (matched == 0) {
            Debug(Debug::WARNING) << "No key of " << par.db1 << " is in " << lookupFile << ", is it the lookup of this DB?\n";
        }
    } else if (par.lookupDb.empty() == false) {
        Debug(Debug::ERROR) << "Lookup file " << lookupFile << " not found\n";
        reader.close();
        return EXIT_FAILURE;
    } else {
        Debug(Debug::INFO) << "No lookup file " << lookupFile << ", keys are written without accessions\n";
    }

    ParquetDBWriter::writeDbAsParquet(reader, par.db2, par.parquetSortMode, par.threads, options, lookupReader, schemaColumns);
    if (lookupReader != NULL) {
        lookupReader->close();
        delete lookupReader;
    }
    reader.close();
    return EXIT_SUCCESS;
}
