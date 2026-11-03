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

    // key sorted, so keyorder and the key range planning work on contiguous ids;
    // a .lookup (key -> accession, written by createdb and others) is converted as well
    int mode = DBReader<DBKeyType>::USE_INDEX|DBReader<DBKeyType>::USE_DATA;
    const bool withLookup = FileUtil::fileExists((par.db1 + ".lookup").c_str());
    if (withLookup) {
        mode |= DBReader<DBKeyType>::USE_LOOKUP;
    } else {
        Debug(Debug::INFO) << "No lookup file " << par.db1 << ".lookup, keys are written without accessions\n";
    }
    DBReader<DBKeyType> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, mode);
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

    ParquetDBWriter::writeDbAsParquet(reader, par.db2, par.parquetSortMode, par.threads, options, withLookup, schemaColumns);
    reader.close();
    return EXIT_SUCCESS;
}
