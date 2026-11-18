#ifndef PARQUETWRITER_H
#define PARQUETWRITER_H
#include "MemoryTracker.h"
#include "Parameters.h"
#include "TranslateNucl.h"
#include "Orf.h"
#include "DBReader.h"
#include <carquet/include/carquet/carquet.h>
#include <functional>

typedef std::pair<void*,carquet_physical_type> Record;

class ParquetDBWriter : public MemoryTracker{
public:
    // column codes for raw database export; start at 1000 so they never collide with Parameters::OUTFMT_*
    enum DbField {
        DBFIELD_KEY = 1000,
        DBFIELD_SEQUENCE,
        DBFIELD_POSITION,
        DBFIELD_PROFILE_AA,                                   // +0..19 -> A,C,D,...,Y
        DBFIELD_QUERY_RESIDUE = DBFIELD_PROFILE_AA + 20,
        DBFIELD_CONSENSUS_RESIDUE,
        DBFIELD_NEFF,
        DBFIELD_GAP_DEL,
        DBFIELD_GAP_INS,
        DBFIELD_QUERY,
        DBFIELD_TARGET,
        DBFIELD_SCORE,
        DBFIELD_SEQID,
        DBFIELD_EVALUE,
        DBFIELD_QSTART,
        DBFIELD_QEND,
        DBFIELD_QLEN,
        DBFIELD_TSTART,
        DBFIELD_TEND,
        DBFIELD_TLEN,
        DBFIELD_QORFSTART,
        DBFIELD_QORFEND,
        DBFIELD_TORFSTART,
        DBFIELD_TORFEND,
        DBFIELD_BACKTRACE,
        DBFIELD_REPRESENTATIVE,
        DBFIELD_MEMBER,
        DBFIELD_PREFSCORE,
        DBFIELD_DIAGONAL,
        DBFIELD_TAXID,
        DBFIELD_RANK,
        DBFIELD_NAME,
        DBFIELD_RANKS,
        DBFIELD_LINEAGE,
        DBFIELD_BLOB,
        DBFIELD_MSA,
        DBFIELD_DATA,
        DBFIELD_TAX_NUM_SEQS,       // aggregatetax: sequences of the set
        DBFIELD_TAX_ASSIGNED_SEQS,  // aggregatetax: sequences with a taxonomic assignment
        DBFIELD_TAX_AGREE_SEQS,     // aggregatetax: sequences agreeing with the selected taxon
        DBFIELD_TAX_SUPPORT,        // aggregatetax: fraction of support
        DBFIELD_ACCESSION,          // .lookup: accession of the key
        DBFIELD_FILE_NUMBER         // .lookup: input file / set number
    };
    // optional columns that cannot be derived from the dbtype alone
    enum DbFieldOptional {
        DBOPT_NONE        = 0,
        DBOPT_ORF         = 1,  // ALIGNMENT_RES: qOrfStart,qOrfEnd,tOrfStart,tOrfEnd
        DBOPT_BACKTRACE   = 2,  // ALIGNMENT_RES: backtrace
        DBOPT_TAX_RANKS   = 4,  // TAXONOMICAL_RESULT: ranks
        DBOPT_TAX_LINEAGE = 8,  // TAXONOMICAL_RESULT: lineage
        DBOPT_TAX_AGGREGATE = 16 // TAXONOMICAL_RESULT from aggregatetax: set counts and support
    };
    // row order of writeDbAsParquet, values of --parquet-sort-mode
    enum SortMode {
        SORT_MODE_AUTO = 0,      // keyorder if the DB fits into half of the available memory, else spillsort
        SORT_MODE_LINEAR = 1,    // DB file order, sequential reads, no key order
        SORT_MODE_KEYORDER = 2,  // key order, random reads
        SORT_MODE_SPILLSORT = 3  // sequential reads, key ranges sorted through temporary files
    };
    // column types of a layout or of --parquet-schema
    enum ColumnType { COL_INT8, COL_UINT8, COL_INT16, COL_INT32, COL_INT64, COL_FLOAT, COL_DOUBLE, COL_STRING, COL_BINARY };
    struct ColumnSpec {
        std::string name;
        ColumnType type;
    };
    // parse '{name:type, name:type, ...}'; false and a message in error on invalid input
    static bool parseSchema(const std::string& text, std::vector<ColumnSpec>& columns, std::string& error);
    static ColumnSpec fieldSpec(const int field);
    static const char* columnTypeName(const ColumnType type);
    // generic DB: one row per tab separated line, key + schema columns.
    // other DBs: the schema selects columns of the dbtype layout by name (types must match the layout)
    void initFromSchema(const char* resultFile, const int dbtype, const unsigned int optional,
                        const std::vector<ColumnSpec>& schemaColumns);
    struct WriterOptions {
        int compression = CARQUET_COMPRESSION_ZSTD;
        int compressionLevel = 3;
        size_t rowGroupSize = 128 * 1024 * 1024;
    };

    // --parquet-compression, --parquet-compression-level, --parquet-row-group-size
    static WriterOptions writerOptionsFromParameters(const Parameters& par);
    // customisation of writeDbAsParquet for other output formats (e.g. convertalis); empty = convert2parquet
    struct ConvertHooks {
        // file names of `files` parts, removing older results; default <output> or <output>/part-NNNNN.parquet
        std::function<std::vector<std::string>(size_t files)> partNames;
        // set up the columns of a part; default: DB layout or --parquet-schema
        std::function<void(ParquetDBWriter& writer, const std::string& file)> openPart;
        // rows of one entry, data is \0 terminated; default: writeEntry
        std::function<void(ParquetDBWriter& writer, DBKeyType key, char* data, size_t len, unsigned int thread_idx)> writeEntry;
        // at least this many files (one writer per file, so this is the parallelism for small DBs)
        size_t minFiles;
        ConvertHooks() : minFiles(0) {}
    };
    // convert a whole DB (reader opened with USE_INDEX|USE_DATA and SORT_BY_ID) into <output> when it
    // results in one file, or <output>/part-NNNNN.parquet otherwise. lookupReader (opened with USE_LOOKUP,
    // NULL = none): its lookup (key, accession, file_number) is written next to it, see lookupParquetName
    static void writeDbAsParquet(DBReader<DBKeyType>& reader, const std::string& output, const int sortMode,
                                 const int threads, const WriterOptions& options, DBReader<DBKeyType>* lookupReader,
                                 const std::vector<ColumnSpec>& schemaColumns = std::vector<ColumnSpec>(),
                                 const ConvertHooks& hooks = ConvertHooks());
    static bool isSupportedDbType(const int dbtype);
    // lookup table name: "_" prefix (skipped by Spark/Hive/Arrow datasets) and ".pq" extension
    // (not matched by *.parquet globs, e.g. in DuckDB)
    static std::string lookupParquetName(const std::string& output, const bool multiFile);
    // keys of entries that produced no rows (e.g. queries without hits): <dir>/_empty.pq or _<name>.empty.pq;
    // a reader needs them to know every key of the DB (ParquetDBReader)
    static std::string emptyParquetName(const std::string& output, const bool multiFile);
    // key-value metadata written by convert2parquet: the dbtype of the converted DB
    static constexpr const char* DBTYPE_METADATA_KEY = "mmseqs.dbtype";
    // remove <dir>/<prefix>_<number>.parquet(.tmp) of an earlier run (convertalis file naming)
    static void removeNumberedParquet(const std::string& dir, const std::string& prefix);
    // codec, level and row group size of the next init*() call
    void setWriterOptions(const WriterOptions& options);
    // parse one entry of this writer's dbtype into rows; key -> first column
    void writeEntry(const DBKeyType key, const char* data, const size_t len);

    ParquetDBWriter();
    void init(const char* resultFile,const std::vector<int>& outcodes);
    void initFromDbType(const char* resultFile, const int dbtype, const unsigned int optional = DBOPT_NONE);
    void initFromFields(const char* resultFile, const std::vector<int>& fieldList);
    // same as initFromDbType, but takes dbtype and optional columns from an opened reader (needs USE_DATA)
    void initFromDbReader(const char* resultFile, DBReader<DBKeyType>& reader, const int thread_idx = 0);
    static std::vector<int> dbTypeToFields(const int dbtype, const unsigned int optional);
    // returns DBOPT_* flags by counting the columns of the first non-empty entry
    static unsigned int detectOptionalColumns(DBReader<DBKeyType>& reader, const int thread_idx = 0);
    const std::vector<int>& getFields() const { return fields; }
    ~ParquetDBWriter();
    void checkStatus(const std::string message);
    void writeCell(const std::string& record_value, const int record_num); 
    void writeCell(const double& record_value, const int record_num); 
    void writeCell(const int& record_value, const int record_num);
    void writeCell(const int64_t& record_value, const int record_num); 
    void writeCell(const float& record_value, const int record_num);
    void addSeqBasedOnAln(const int record_num, const char *seq, unsigned int offset,
                        const std::string &bt, bool reverse, bool isReverseStrand,
                        bool translateSequence, const TranslateNucl &translateNucl);
    void writeBatchToFile(const std::vector<int>& outcodes);
    // marks the end of a row; flushes the batch to carquet once it holds flush_bytes of data
    void endRow();
    // rows ended so far (endRow calls), to see whether an entry produced rows
    size_t getRowsWritten() const { return rowsWritten; }
    // writes all pending rows to carquet, also called by close()
    void flush();
    void setFlushBytes(const size_t bytes);
    static const size_t DEFAULT_FLUSH_BYTES = 8 * 1024 * 1024;
    void close();
private:
    carquet_status_t status;
    carquet_status_t writer_close_status = CARQUET_ERROR_NOT_OPEN;

    carquet_error_t err;
    carquet_schema_t* schema;
    carquet_writer_options_t opts;
    carquet_writer_t* writer;
    // pending values per output column, indexed by column (sized by the init functions); a vector instead of a
    // std::map: writeCell looks the column up for every value
    std::vector<Record> record_batch;
    int record_counter;
    std::vector<int> fields;
    // data held in record_batch since the last flush
    size_t pending_bytes = 0;
    size_t pending_rows = 0;
    size_t flush_bytes = DEFAULT_FLUSH_BYTES;
    int dbtype = 0;
    // dbtype written as DBTYPE_METADATA_KEY on close; -1 = none (layouts not based on a DB, e.g. convertalis)
    int metadataDbtype = -1;
    size_t rowsWritten = 0;
    WriterOptions writerOptions;
    // per output column, for parsing and error messages
    std::vector<std::string> columnNames;
    std::vector<ColumnType> columnTypes;
    // columns of the DB entry (key + fields); projection maps them to output columns, -1 = not written
    // (always filled, identity when all columns are written in DB order: no extra branch per cell)
    int nativeColumns = 0;
    std::vector<int> projection;
    bool parseLines = false;    // generic DB with schema: tab separated lines
    DBKeyType currentKey = 0;   // entry being parsed, for error messages

    template <typename T> void put(const int nativeColumn, const T& value) {
        const int column = projection[nativeColumn];
        if (column >= 0) {
            writeCell(value, column);
        }
    }
    void initFromSpecs(const char* resultFile, const std::vector<ColumnSpec>& specs);
    void setIdentityProjection(const size_t columns);
    void addSpecColumn(const ColumnSpec& spec);
    void writeField(const int column, const char* s, const size_t n);
    static unsigned int detectTaxonomyColumns(const char* entry);

    void* typToPointer(carquet_physical_type type);
    void outcodeToField(const int outcode);
    void addColumn(const char* name, carquet_physical_type type, const carquet_logical_type_t* logical);
    void dbFieldToField(const int field);
    void openWriter(const char* resultFile, const size_t columnCount);
    void writeColumn(int column_number);
    void freeRecordBatch();
    void clearRecordBatch();
};

#endif