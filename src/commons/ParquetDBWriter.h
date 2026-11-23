#ifndef PARQUETWRITER_H
#define PARQUETWRITER_H
#include "MemoryTracker.h"
#include "Parameters.h"
#include "TranslateNucl.h"
#include "Orf.h"
#include "DBReader.h"
#include <carquet/include/carquet/carquet.h>
#include <functional>

class ParquetDBWriter : public MemoryTracker{
public:
    // column codes for raw database export; start at 1000 so they never collide with Parameters::OUTFMT_*
    enum DbField {
        DBFIELD_KEY = 1000,
        DBFIELD_SEQUENCE,
        DBFIELD_POSITION,
        DBFIELD_PROFILE_AA,                                   
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
        DBFIELD_TAX_NUM_SEQS,       
        DBFIELD_TAX_ASSIGNED_SEQS,  
        DBFIELD_TAX_AGREE_SEQS,     
        DBFIELD_TAX_SUPPORT,       
        DBFIELD_ACCESSION,          
        DBFIELD_FILE_NUMBER         
    };

    enum DbFieldOptional {
        DBOPT_NONE        = 0,
        DBOPT_ORF         = 1,  
        DBOPT_BACKTRACE   = 2,  
        DBOPT_TAX_RANKS   = 4,  
        DBOPT_TAX_LINEAGE = 8,  
        DBOPT_TAX_AGGREGATE = 16 
    };

    enum SortMode {
        SORT_MODE_AUTO = 0,      
        SORT_MODE_LINEAR = 1,    
        SORT_MODE_KEYORDER = 2,  
        SORT_MODE_SPILLSORT = 3  
    };

    enum ColumnType { COL_INT8, COL_UINT8, COL_INT16, COL_INT32, COL_INT64, COL_FLOAT, COL_DOUBLE, COL_STRING, COL_BINARY };
    struct ColumnSpec {
        std::string name;
        ColumnType type;
    };
    static bool parseSchema(const std::string& text, std::vector<ColumnSpec>& columns, std::string& error);
    static ColumnSpec fieldSpec(const int field);
    static const char* columnTypeName(const ColumnType type);

    void initFromSchema(const char* resultFile, const int dbtype, const unsigned int optional,
                        const std::vector<ColumnSpec>& schemaColumns);
    struct WriterOptions {
        int compression = CARQUET_COMPRESSION_ZSTD;
        int compressionLevel = 3;
        size_t rowGroupSize = 128 * 1024 * 1024;
    };

    static WriterOptions writerOptionsFromParameters(const Parameters& par);
    struct ConvertHooks {
        // <output> or <output>/part-NNNNN.parquet
        std::function<std::vector<std::string>(size_t files)> partNames;
        std::function<void(ParquetDBWriter& writer, const std::string& file)> openPart;
        std::function<void(ParquetDBWriter& writer, DBKeyType key, char* data, size_t len, unsigned int thread_idx)> writeEntry;
        size_t minFiles;
        ConvertHooks() : minFiles(0) {}
    };
    static void writeDbAsParquet(DBReader<DBKeyType>& reader, const std::string& output, const int sortMode,
                                 const int threads, const WriterOptions& options, DBReader<DBKeyType>* lookupReader,
                                 const std::vector<ColumnSpec>& schemaColumns = std::vector<ColumnSpec>(),
                                 const ConvertHooks& hooks = ConvertHooks());
    static bool isSupportedDbType(const int dbtype);
    static std::string lookupParquetName(const std::string& output, const bool multiFile);
    static std::string emptyParquetName(const std::string& output, const bool multiFile);
    static constexpr const char* DBTYPE_METADATA_KEY = "mmseqs.dbtype";
    static void removeNumberedParquet(const std::string& dir, const std::string& prefix);
    void setWriterOptions(const WriterOptions& options);
    void writeEntry(const DBKeyType key, const char* data, const size_t len);

    void init(const char* resultFile,const std::vector<int>& outcodes);
    void initFromDbType(const char* resultFile, const int dbtype, const unsigned int optional = DBOPT_NONE);
    void initFromFields(const char* resultFile, const std::vector<int>& fieldList);
    static std::vector<int> dbTypeToFields(const int dbtype, const unsigned int optional);
    static unsigned int detectOptionalColumns(DBReader<DBKeyType>& reader, const int thread_idx = 0);
    ~ParquetDBWriter();
    void checkStatus(const std::string& message);
    void writeCell(const char* s, const size_t n, const int record_num);
    void writeCell(const std::string& record_value, const int record_num);
    void writeCell(const double& record_value, const int record_num); 
    void writeCell(const int& record_value, const int record_num);
    void writeCell(const int64_t& record_value, const int record_num); 
    void writeCell(const float& record_value, const int record_num);
    void addSeqBasedOnAln(const int record_num, const char *seq, unsigned int offset,
                        const std::string &bt, bool reverse, bool isReverseStrand,
                        bool translateSequence, const TranslateNucl &translateNucl);
    void endRow();
    size_t getRowsWritten() const { return rowsWritten; }
    void flush();
    static const size_t DEFAULT_FLUSH_BYTES = 8 * 1024 * 1024;
    void close();
private:
    // all values of a batch are stored back to back in one arena, avoiding a heap allocation
    // and the std::string overhead per value; byte arrays are built at finalize, once the arena no longer moves
    struct StringColumn {
        std::vector<char> arena;
        std::vector<size_t> offsets;
        std::vector<carquet_byte_array_t> values;

        void add(const char* s, size_t len) {
            offsets.push_back(arena.size());
            if (arena.size() + len > arena.capacity() && len > arena.size()) {
                arena.reserve(arena.size() + len);
            }
            arena.insert(arena.end(), s, s + len);
        }

        // build a value char by char with append()
        void startValue() {
            offsets.push_back(arena.size());
        }

        void append(char c) {
            arena.push_back(c);
        }

        void finalize(){
            values.clear();
            values.reserve(offsets.size());
            uint8_t* base = reinterpret_cast<uint8_t*>(arena.data());
            for (size_t i = 0; i < offsets.size(); i++){
                const size_t end = (i + 1 < offsets.size()) ? offsets[i + 1] : arena.size();
                values.push_back(carquet_byte_array_t{base + offsets[i], static_cast<int32_t>(end - offsets[i])});
            }
        }

        void clear(size_t maxKeepBytes) {
            values.clear();
            offsets.clear();
            arena.clear();
            if (arena.capacity() > maxKeepBytes) {
                std::vector<char>().swap(arena);
                std::vector<carquet_byte_array_t>().swap(values);
                std::vector<size_t>().swap(offsets);
            }
        }
    };

    struct Column {
        carquet_physical_type type;
        std::vector<int32_t> i32;
        std::vector<int64_t> i64;
        std::vector<float> f32;
        std::vector<double> f64;
        StringColumn str;
        explicit Column(carquet_physical_type type) : type(type) {}
    };

    carquet_status_t status = CARQUET_OK;
    carquet_schema_t* schema = NULL;
    carquet_writer_t* writer = NULL;
    std::vector<Column> record_batch;
    size_t pending_bytes = 0;
    size_t pending_rows = 0;
    size_t flush_bytes = DEFAULT_FLUSH_BYTES;
    int dbtype = -1;
    size_t rowsWritten = 0;
    WriterOptions writerOptions;
    std::vector<std::string> columnNames;
    std::vector<ColumnType> columnTypes;
    int nativeColumns = 0;
    std::vector<int> projection;
    bool parseLines = false;
    DBKeyType currentKey = 0;

    template <typename T> void put(const int nativeColumn, const T& value) {
        const int column = projection[nativeColumn];
        if (column >= 0) {
            writeCell(value, column);
        }
    }
    void put(const int nativeColumn, const char* s, const size_t n) {
        const int column = projection[nativeColumn];
        if (column >= 0) {
            writeCell(s, n, column);
        }
    }
    void initFromSpecs(const char* resultFile, const std::vector<ColumnSpec>& specs);
    void setIdentityProjection(const size_t columns);
    void addSpecColumn(const ColumnSpec& spec);
    void writeField(const int column, const char* s, const size_t n);
    static unsigned int detectTaxonomyColumns(const char* entry);

    static ColumnSpec outcodeSpec(const int outcode);
    void addColumn(const char* name, carquet_physical_type type, const carquet_logical_type_t* logical);
    void openWriter(const char* resultFile);
    void writeColumn(int column_number);
};

#endif