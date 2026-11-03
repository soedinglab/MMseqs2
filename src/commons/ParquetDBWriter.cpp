#include "Debug.h"
#include "Util.h"
#include "Parameters.h"
#include "TranslateNucl.h"
#include "Orf.h"
#include "ParquetDBWriter.h"
#include "FileUtil.h"
#include "Sequence.h"

#include <fast_float/fast_float.h>
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <system_error>
#include <dirent.h>
#include <mutex>
#include <thread>
#include <unistd.h>

#ifdef OPENMP
#include <omp.h>
#endif


//#include <variant> //new
#include <carquet/include/carquet/carquet.h>
#define ZSTD_STATIC_LINKING_ONLY
#include <zstd.h>


// all values of a batch are stored back to back in one arena, avoiding a heap allocation
// and the std::string overhead per value; byte arrays are built at finalize, once the arena no longer moves
struct StringColumn {
    std::vector<char> arena;
    std::vector<size_t> offsets;
    std::vector<carquet_byte_array_t> values;

    void add(const char* s, size_t len) {
        offsets.push_back(arena.size());
        if (arena.size() + len > arena.capacity() && len > arena.size()) {
            // a value larger than the whole batch so far: grow to the exact size instead of doubling
            arena.reserve(arena.size() + len);
        }
        arena.insert(arena.end(), s, s + len);
    }

    void add(const std::string& s) {
        add(s.data(), s.size());
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

    // keeps up to maxKeepBytes of arena capacity for the next batch, so a single huge value
    // does not pin its memory for the rest of the run
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


void ParquetDBWriter::checkStatus(const std::string message){
    if (status != carquet_status_t::CARQUET_OK){
        Debug(Debug::ERROR) << message << "\n";
        EXIT(EXIT_FAILURE);
    }
}

ParquetDBWriter::ParquetDBWriter(){

};

ParquetDBWriter::~ParquetDBWriter(){
    if (writer_close_status == CARQUET_ERROR_NOT_OPEN){
        // the instance was not initialised
        return;
    }
    // close first, it flushes the pending rows from the record batch
    if (writer_close_status != CARQUET_ERROR_ALREADY_CLOSED){
        close();
    }
    freeRecordBatch();
}

void ParquetDBWriter::init(const char* resultFile,const std::vector<int>& outcodes){
    err = CARQUET_ERROR_INIT;
    schema = carquet_schema_create(&err);
    record_counter = 0;
    
    if (!schema) {
        status = CARQUET_ERROR_INTERNAL;
        checkStatus("Schema initialisation error");
    }

    for (size_t i = 0; i < outcodes.size(); i++){
        // creates schema column and returns void pointer
        outcodeToField(outcodes[i]);
        checkStatus("Adding column to schema failed");
    }
    openWriter(resultFile, outcodes.size());
};

void ParquetDBWriter::initFromDbType(const char* resultFile, const int dbtype, const unsigned int optional){
    this->dbtype = dbtype;
    initFromFields(resultFile, dbTypeToFields(dbtype, optional));
}

void ParquetDBWriter::initFromFields(const char* resultFile, const std::vector<int>& fieldList){
    std::vector<ColumnSpec> specs;
    for (size_t i = 0; i < fieldList.size(); i++){
        specs.push_back(fieldSpec(fieldList[i]));
    }
    fields = fieldList;
    nativeColumns = static_cast<int>(fieldList.size());
    projection.clear();
    initFromSpecs(resultFile, specs);
}

void ParquetDBWriter::initFromSpecs(const char* resultFile, const std::vector<ColumnSpec>& specs){
    err = CARQUET_ERROR_INIT;
    schema = carquet_schema_create(&err);
    record_counter = 0;

    if (!schema) {
        status = CARQUET_ERROR_INTERNAL;
        checkStatus("Schema initialisation error");
    }

    columnNames.clear();
    columnTypes.clear();
    for (size_t i = 0; i < specs.size(); i++){
        addSpecColumn(specs[i]);
        checkStatus("Adding column " + specs[i].name + " to schema failed");
    }
    openWriter(resultFile, specs.size());
}

void ParquetDBWriter::initFromDbReader(const char* resultFile, DBReader<DBKeyType>& reader, const int thread_idx){
    const int dbtype = reader.getDbtype();
    initFromDbType(resultFile, dbtype, detectOptionalColumns(reader, thread_idx));
}

unsigned int ParquetDBWriter::detectOptionalColumns(DBReader<DBKeyType>& reader, const int thread_idx){
    if (Parameters::isEqualDbtype(reader.getDbtype(), Parameters::DBTYPE_TAXONOMICAL_RESULT)) {
        // the first classified entry tells the columns apart; unclassified ones have empty optional fields
        unsigned int fallback = DBOPT_NONE;
        bool found = false;
        for (size_t id = 0; id < reader.getSize(); id++) {
            const char* data = reader.getData(id, thread_idx);
            if (data == NULL || *data == '\0') {
                continue;
            }
            if (found == false) {
                fallback = detectTaxonomyColumns(data);
                found = true;
            }
            if (Util::fast_atoi<int>(data) != 0) {
                return detectTaxonomyColumns(data);
            }
        }
        return fallback;
    }
    if (Parameters::isEqualDbtype(reader.getDbtype(), Parameters::DBTYPE_ALIGNMENT_RES) == false) {
        return DBOPT_NONE;
    }
    const char* words[Matcher::ALN_RES_WITH_ORF_AND_BT_COL_CNT + 1];
    for (size_t id = 0; id < reader.getSize(); id++) {
        const char* data = reader.getData(id, thread_idx);
        if (data == NULL || *data == '\0') {
            // queries without hits have empty entries, check the next one
            continue;
        }
        const size_t columns = Util::getFieldsOfLine(data, words, Matcher::ALN_RES_WITH_ORF_AND_BT_COL_CNT + 1);
        switch (columns) {
            case Matcher::ALN_RES_WITHOUT_BT_COL_CNT:
                return DBOPT_NONE;
            case Matcher::ALN_RES_WITH_BT_COL_CNT:
                return DBOPT_BACKTRACE;
            case Matcher::ALN_RES_WITH_ORF_POS_WITHOUT_BT_COL_CNT:
                return DBOPT_ORF;
            case Matcher::ALN_RES_WITH_ORF_AND_BT_COL_CNT:
                return DBOPT_ORF | DBOPT_BACKTRACE;
            default:
                Debug(Debug::ERROR) << "Unexpected number of alignment columns: " << columns << "\n";
                EXIT(EXIT_FAILURE);
        }
    }
    // all entries are empty, no optional columns to write
    return DBOPT_NONE;
}

// shared by init() and initFromDbType()
void ParquetDBWriter::openWriter(const char* resultFile, const size_t columnCount){
    carquet_writer_options_init(&opts);
    opts.compression = static_cast<carquet_compression_t>(writerOptions.compression);
    opts.compression_level = writerOptions.compressionLevel;
    opts.row_group_size = static_cast<int64_t>(writerOptions.rowGroupSize);
    opts.page_size = 1024 * 1024;
    opts.write_statistics = true;
    opts.write_crc = true;
    opts.write_page_index = true;
    opts.write_bloom_filters = true;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN; // might want to change, depends on performance
    opts.created_by = "MMseqs parquet writer";

    writer = carquet_writer_create(resultFile, schema, &opts, &err);
    carquet_schema_free(schema); // writer copies schema, can be freed early
    if (!writer) {
        freeRecordBatch();
        status = CARQUET_ERROR_INTERNAL;
        checkStatus("Writer initialisation error");
    }
    for (size_t i = 0; i < columnCount; i++){
        status = carquet_writer_set_column_encoding(writer, i, CARQUET_ENCODING_PLAIN);
        checkStatus("Column encoding overwrite");
    }
    writer_close_status = CARQUET_OK;
}

std::vector<int> ParquetDBWriter::dbTypeToFields(const int dbtype, const unsigned int optional){
    std::vector<int> f;
    const int type = dbtype & Parameters::DBTYPE_MASK; // strip extended/compressed bits
    switch (type) {
        case Parameters::DBTYPE_AMINO_ACIDS:
        case Parameters::DBTYPE_NUCLEOTIDES:
            f = {DBFIELD_KEY, DBFIELD_SEQUENCE};
            break;
        case Parameters::DBTYPE_HMM_PROFILE:
            f = {DBFIELD_KEY, DBFIELD_POSITION};
            for (int aa = 0; aa < 20; aa++){
                f.push_back(DBFIELD_PROFILE_AA + aa);
            }
            f.insert(f.end(), {DBFIELD_QUERY_RESIDUE, DBFIELD_CONSENSUS_RESIDUE, DBFIELD_NEFF,
                               DBFIELD_GAP_DEL, DBFIELD_GAP_INS});
            break;
        case Parameters::DBTYPE_ALIGNMENT_RES:
            f = {DBFIELD_QUERY, DBFIELD_TARGET, DBFIELD_SCORE, DBFIELD_SEQID, DBFIELD_EVALUE,
                 DBFIELD_QSTART, DBFIELD_QEND, DBFIELD_QLEN, DBFIELD_TSTART, DBFIELD_TEND, DBFIELD_TLEN};
            if (optional & DBOPT_ORF){
                f.insert(f.end(), {DBFIELD_QORFSTART, DBFIELD_QORFEND, DBFIELD_TORFSTART, DBFIELD_TORFEND});
            }
            if (optional & DBOPT_BACKTRACE){
                f.push_back(DBFIELD_BACKTRACE);
            }
            break;
        case Parameters::DBTYPE_CLUSTER_RES:
            f = {DBFIELD_REPRESENTATIVE, DBFIELD_MEMBER};
            break;
        case Parameters::DBTYPE_PREFILTER_RES:
        case Parameters::DBTYPE_PREFILTER_REV_RES:
            f = {DBFIELD_QUERY, DBFIELD_TARGET, DBFIELD_PREFSCORE, DBFIELD_DIAGONAL};
            break;
        case Parameters::DBTYPE_TAXONOMICAL_RESULT:
            f = {DBFIELD_QUERY, DBFIELD_TAXID, DBFIELD_RANK, DBFIELD_NAME};
            if (optional & DBOPT_TAX_AGGREGATE){
                f.insert(f.end(), {DBFIELD_TAX_NUM_SEQS, DBFIELD_TAX_ASSIGNED_SEQS, DBFIELD_TAX_AGREE_SEQS, DBFIELD_TAX_SUPPORT});
            }
            if (optional & DBOPT_TAX_RANKS){
                f.push_back(DBFIELD_RANKS);
            }
            if (optional & DBOPT_TAX_LINEAGE){
                f.push_back(DBFIELD_LINEAGE);
            }
            break;
        case Parameters::DBTYPE_INDEX_DB:
            f = {DBFIELD_KEY, DBFIELD_BLOB};
            break;
        case Parameters::DBTYPE_CA3M_DB:
        case Parameters::DBTYPE_MSA_DB:
            f = {DBFIELD_KEY, DBFIELD_MSA};
            break;
        case Parameters::DBTYPE_GENERIC_DB:
            f = {DBFIELD_KEY, DBFIELD_DATA};
            break;
        default:
            // OMIT_FILE, OFFSETDB, DIRECTORY, FLATFILE, SEQTAXDB, STDIN, URI have no record layout
            Debug(Debug::ERROR) << "Database type " << type << " has no parquet layout\n";
            EXIT(EXIT_FAILURE);
    }
    return f;
}

void ParquetDBWriter::addColumn(const char* name, carquet_physical_type type, const carquet_logical_type_t* logical){
    status = carquet_schema_add_column(schema, name, type, logical, CARQUET_REPETITION_REQUIRED, 0, 0);
    record_batch[record_counter] = Record(typToPointer(type), type);
    record_counter++;
}

ParquetDBWriter::ColumnSpec ParquetDBWriter::fieldSpec(const int field){
    if (field >= DBFIELD_PROFILE_AA && field < DBFIELD_PROFILE_AA + 20){
        static const char* aaNames[20] = {"A","C","D","E","F","G","H","I","K","L",
                                          "M","N","P","Q","R","S","T","V","W","Y"};
        return ColumnSpec{aaNames[field - DBFIELD_PROFILE_AA], COL_INT8};
    }
    switch (field) {
        // keys are always INT64 so uint32 keys > INT32_MAX stay positive and the schema is build independent
        case DBFIELD_KEY:               return ColumnSpec{"key",               COL_INT64};
        case DBFIELD_QUERY:             return ColumnSpec{"query",             COL_INT64};
        case DBFIELD_TARGET:            return ColumnSpec{"target",            COL_INT64};
        case DBFIELD_REPRESENTATIVE:    return ColumnSpec{"representative",    COL_INT64};
        case DBFIELD_MEMBER:            return ColumnSpec{"member",            COL_INT64};

        case DBFIELD_SEQUENCE:          return ColumnSpec{"sequence",          COL_STRING};
        case DBFIELD_BACKTRACE:         return ColumnSpec{"backtrace",         COL_STRING};
        case DBFIELD_RANK:              return ColumnSpec{"rank",              COL_STRING};
        case DBFIELD_NAME:              return ColumnSpec{"name",              COL_STRING};
        case DBFIELD_RANKS:             return ColumnSpec{"ranks",             COL_STRING};
        case DBFIELD_LINEAGE:           return ColumnSpec{"lineage",           COL_STRING};
        case DBFIELD_MSA:               return ColumnSpec{"msa",               COL_STRING};
        case DBFIELD_DATA:              return ColumnSpec{"data",              COL_STRING};
        case DBFIELD_ACCESSION:         return ColumnSpec{"accession",         COL_STRING};
        case DBFIELD_BLOB:              return ColumnSpec{"blob",              COL_BINARY};

        case DBFIELD_POSITION:          return ColumnSpec{"position",          COL_INT32};
        case DBFIELD_QUERY_RESIDUE:     return ColumnSpec{"query_residue",     COL_UINT8};
        case DBFIELD_CONSENSUS_RESIDUE: return ColumnSpec{"consensus_residue", COL_UINT8};
        case DBFIELD_NEFF:              return ColumnSpec{"neff",              COL_UINT8};
        case DBFIELD_GAP_DEL:           return ColumnSpec{"gap_del",           COL_UINT8};
        case DBFIELD_GAP_INS:           return ColumnSpec{"gap_ins",           COL_UINT8};

        case DBFIELD_SCORE:             return ColumnSpec{"score",             COL_INT32};
        case DBFIELD_SEQID:             return ColumnSpec{"seqId",             COL_FLOAT};
        case DBFIELD_EVALUE:            return ColumnSpec{"evalue",            COL_DOUBLE};
        case DBFIELD_QSTART:            return ColumnSpec{"qStart",            COL_INT32};
        case DBFIELD_QEND:              return ColumnSpec{"qEnd",              COL_INT32};
        case DBFIELD_QLEN:              return ColumnSpec{"qLen",              COL_INT32};
        case DBFIELD_TSTART:            return ColumnSpec{"tStart",            COL_INT32};
        case DBFIELD_TEND:              return ColumnSpec{"tEnd",              COL_INT32};
        case DBFIELD_TLEN:              return ColumnSpec{"tLen",              COL_INT32};
        case DBFIELD_QORFSTART:         return ColumnSpec{"qOrfStart",         COL_INT32};
        case DBFIELD_QORFEND:           return ColumnSpec{"qOrfEnd",           COL_INT32};
        case DBFIELD_TORFSTART:         return ColumnSpec{"tOrfStart",         COL_INT32};
        case DBFIELD_TORFEND:           return ColumnSpec{"tOrfEnd",           COL_INT32};
        case DBFIELD_PREFSCORE:         return ColumnSpec{"prefScore",         COL_INT32};
        case DBFIELD_DIAGONAL:          return ColumnSpec{"diagonal",          COL_INT16};
        case DBFIELD_TAXID:             return ColumnSpec{"taxid",             COL_INT32};
        case DBFIELD_TAX_NUM_SEQS:      return ColumnSpec{"num_seqs",          COL_INT32};
        case DBFIELD_TAX_ASSIGNED_SEQS: return ColumnSpec{"assigned_seqs",     COL_INT32};
        case DBFIELD_TAX_AGREE_SEQS:    return ColumnSpec{"agree_seqs",        COL_INT32};
        case DBFIELD_TAX_SUPPORT:       return ColumnSpec{"support",           COL_FLOAT};
        case DBFIELD_FILE_NUMBER:       return ColumnSpec{"file_number",       COL_INT64};
        default:
            Debug(Debug::ERROR) << "Tried to create schema for unkown database field, field number:" << field << "\n";
            EXIT(EXIT_FAILURE);
    }
}

void ParquetDBWriter::addSpecColumn(const ColumnSpec& spec){
    carquet_logical_type_t string_type = {};
    string_type.id = CARQUET_LOGICAL_STRING;
    carquet_logical_type_t int_type = {};
    int_type.id = CARQUET_LOGICAL_INTEGER;
    int_type.params.integer.is_signed = true;
    switch (spec.type) {
        case COL_INT8:   int_type.params.integer.bit_width = 8;  addColumn(spec.name.c_str(), CARQUET_PHYSICAL_INT32, &int_type); break;
        case COL_INT16:  int_type.params.integer.bit_width = 16; addColumn(spec.name.c_str(), CARQUET_PHYSICAL_INT32, &int_type); break;
        case COL_UINT8:
            int_type.params.integer.bit_width = 8;
            int_type.params.integer.is_signed = false;
            addColumn(spec.name.c_str(), CARQUET_PHYSICAL_INT32, &int_type);
            break;
        case COL_INT32:  addColumn(spec.name.c_str(), CARQUET_PHYSICAL_INT32, NULL); break;
        case COL_INT64:  addColumn(spec.name.c_str(), CARQUET_PHYSICAL_INT64, NULL); break;
        case COL_FLOAT:  addColumn(spec.name.c_str(), CARQUET_PHYSICAL_FLOAT, NULL); break;
        case COL_DOUBLE: addColumn(spec.name.c_str(), CARQUET_PHYSICAL_DOUBLE, NULL); break;
        case COL_STRING: addColumn(spec.name.c_str(), CARQUET_PHYSICAL_BYTE_ARRAY, &string_type); break;
        case COL_BINARY: addColumn(spec.name.c_str(), CARQUET_PHYSICAL_BYTE_ARRAY, NULL); break; // binary, not UTF-8
    }
    columnNames.push_back(spec.name);
    columnTypes.push_back(spec.type);
}

void ParquetDBWriter::dbFieldToField(const int field){
    addSpecColumn(fieldSpec(field));
}

void* ParquetDBWriter::typToPointer(carquet_physical_type type){
    switch (type){
        case carquet_physical_type::CARQUET_PHYSICAL_BYTE_ARRAY:
            return reinterpret_cast<void*>(new StringColumn());
        case carquet_physical_type::CARQUET_PHYSICAL_DOUBLE:
            return reinterpret_cast<void*>(new std::vector<double>());
        case carquet_physical_type::CARQUET_PHYSICAL_INT32:
            return reinterpret_cast<void*>(new std::vector<int32_t>());
        case carquet_physical_type::CARQUET_PHYSICAL_INT64:
            return reinterpret_cast<void*>(new std::vector<int64_t>());
        case carquet_physical_type::CARQUET_PHYSICAL_FLOAT:
            return reinterpret_cast<void*>(new std::vector<float>());
        default:
            status = CARQUET_ERROR_INTERNAL;
            checkStatus("Create pointer for unkown field type");
            return nullptr;
    }
}

void ParquetDBWriter::outcodeToField(const int outcode){
    carquet_logical_type_t string_type = {};
    string_type.id = CARQUET_LOGICAL_STRING;
    switch (outcode) {
        case Parameters::OUTFMT_QUERY:{
            status = carquet_schema_add_column(schema, "query", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* query = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(query,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TARGET:{
            status =carquet_schema_add_column(schema, "target", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* target = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(target,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_EVALUE:{
            status =carquet_schema_add_column(schema, "evalue", CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* evalue = typToPointer(CARQUET_PHYSICAL_DOUBLE);
            record_batch[record_counter] = Record(evalue,CARQUET_PHYSICAL_DOUBLE);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_GAPOPEN:{
            status =carquet_schema_add_column(schema, "gapopen", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* gapopen = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(gapopen,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_FIDENT:{
            status =carquet_schema_add_column(schema, "fident", CARQUET_PHYSICAL_FLOAT, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* findent = typToPointer(CARQUET_PHYSICAL_FLOAT);
            record_batch[record_counter] = Record(findent,CARQUET_PHYSICAL_FLOAT);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_PIDENT:{
            status =carquet_schema_add_column(schema, "pident", CARQUET_PHYSICAL_FLOAT, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* pindent = typToPointer(CARQUET_PHYSICAL_FLOAT);
            record_batch[record_counter] = Record(pindent,CARQUET_PHYSICAL_FLOAT);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_NIDENT:{
            status =carquet_schema_add_column(schema, "nident", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* nident = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(nident,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QSTART:{
            status =carquet_schema_add_column(schema, "Qstart", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qstart = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(qstart,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QEND:{
            status =carquet_schema_add_column(schema, "qend", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qend = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(qend,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QLEN:{
            status =carquet_schema_add_column(schema, "qlen", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qlen = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(qlen,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TSTART:{
            status = carquet_schema_add_column(schema, "tstart", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* tstart = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(tstart,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TEND:{
            status = carquet_schema_add_column(schema, "tend", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* tend = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(tend,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TLEN:{
            status = carquet_schema_add_column(schema, "tlen", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* tlen = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(tlen,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_ALNLEN:{
            status = carquet_schema_add_column(schema, "alnlen", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* alnlen = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(alnlen,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_RAW:{
            status = carquet_schema_add_column(schema, "raw", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* raw = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(raw,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_BITS:{
            status = carquet_schema_add_column(schema, "score", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* score = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(score,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_CIGAR:{
            status = carquet_schema_add_column(schema, "cigar", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* cigar = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(cigar,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QSEQ:{
            status = carquet_schema_add_column(schema, "qseq", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qseq = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(qseq,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TSEQ:{
            status = carquet_schema_add_column(schema, "tseq", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* tseq = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(tseq,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QHEADER:{
            status = carquet_schema_add_column(schema, "qheader", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qheader = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(qheader,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_THEADER:{
            status = carquet_schema_add_column(schema, "theader", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* theader = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(theader,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QALN:{
            status = carquet_schema_add_column(schema, "qaln", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qaln = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(qaln,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TALN:{
            status = carquet_schema_add_column(schema, "taln", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* taln = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(taln,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_MISMATCH:{
            status = carquet_schema_add_column(schema, "mismatch", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* mismatch = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(mismatch,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QCOV:{
            status = carquet_schema_add_column(schema, "qcov", CARQUET_PHYSICAL_FLOAT, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qcov = typToPointer(CARQUET_PHYSICAL_FLOAT);
            record_batch[record_counter] = Record(qcov,CARQUET_PHYSICAL_FLOAT);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TCOV:{
            status = carquet_schema_add_column(schema, "tcov", CARQUET_PHYSICAL_FLOAT, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* tcov = typToPointer(CARQUET_PHYSICAL_FLOAT);
            record_batch[record_counter] = Record(tcov,CARQUET_PHYSICAL_FLOAT);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QSET:{
            status = carquet_schema_add_column(schema, "qset", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qset = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(qset,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QSETID:{
#ifdef MMSEQS_INT64_IDS
            status = carquet_schema_add_column(schema, "qsetid", CARQUET_PHYSICAL_INT64, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qsetid = typToPointer(CARQUET_PHYSICAL_INT64);
            record_batch[record_counter] = Record(qsetid,CARQUET_PHYSICAL_INT64);
            record_counter++;
#else
            status = carquet_schema_add_column(schema, "qsetid", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qsetid = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(qsetid,CARQUET_PHYSICAL_INT32);
            record_counter++;
#endif
            break;
        }
        case Parameters::OUTFMT_TSET:{
            status = carquet_schema_add_column(schema, "tset", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qset = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(qset,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TSETID:{
#ifdef MMSEQS_INT64_IDS
            status = carquet_schema_add_column(schema, "tsetid", CARQUET_PHYSICAL_INT64, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* tsetid = typToPointer(CARQUET_PHYSICAL_INT64);
            record_batch[record_counter] = Record(tsetid,CARQUET_PHYSICAL_INT64);
            record_counter++;
#else
            status = carquet_schema_add_column(schema, "tsetid", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* tsetid = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(tsetid,CARQUET_PHYSICAL_INT32);
            record_counter++;
#endif
            break;
        }
        case Parameters::OUTFMT_TAXID:{
#ifdef MMSEQS_INT64_IDS
            status = carquet_schema_add_column(schema, "taxid", CARQUET_PHYSICAL_INT64, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* taxid = typToPointer(CARQUET_PHYSICAL_INT64);
            record_batch[record_counter] = Record(taxid,CARQUET_PHYSICAL_INT64);
            record_counter++;
#else
            status = carquet_schema_add_column(schema, "taxid", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* taxid = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(taxid,CARQUET_PHYSICAL_INT32);
            record_counter++;
#endif
            break;
        }
        case Parameters::OUTFMT_TAXNAME:{
            status = carquet_schema_add_column(schema, "taxname", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* taxname = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(taxname,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TAXLIN:{
            status = carquet_schema_add_column(schema, "taxlin", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* taxlin = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(taxlin,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_EMPTY:{
            status = carquet_schema_add_column(schema, "empty", CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* empty = typToPointer(CARQUET_PHYSICAL_BYTE_ARRAY);
            record_batch[record_counter] = Record(empty,CARQUET_PHYSICAL_BYTE_ARRAY);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QORFSTART:{
            status = carquet_schema_add_column(schema, "qorfstart", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qorfstart = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(qorfstart,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QORFEND:{
            status = carquet_schema_add_column(schema, "qorfend", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qorfend = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(qorfend,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TORFSTART:{
            status = carquet_schema_add_column(schema, "torfstart", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* torfstart = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(torfstart,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TORFEND:{
            status = carquet_schema_add_column(schema, "torfend", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* torfend = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(torfend,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_PPOS:{
            status = carquet_schema_add_column(schema, "ppos", CARQUET_PHYSICAL_FLOAT, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* ppos = typToPointer(CARQUET_PHYSICAL_FLOAT);
            record_batch[record_counter] = Record(ppos,CARQUET_PHYSICAL_FLOAT);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_QFRAME:{
            status = carquet_schema_add_column(schema, "qframe", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* qframe = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(qframe,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        case Parameters::OUTFMT_TFRAME:{
            status = carquet_schema_add_column(schema, "tframe", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
            void* tframe = typToPointer(CARQUET_PHYSICAL_INT32);
            record_batch[record_counter] = Record(tframe,CARQUET_PHYSICAL_INT32);
            record_counter++;
            break;
        }
        default:
            status = CARQUET_ERROR_INTERNAL;
            checkStatus("Tried to create schema for unkown field type, outcode number:" + std::to_string(outcode));
            break;
    }
}

void ParquetDBWriter::addSeqBasedOnAln(const int record_num, const char *seq, unsigned int offset,
                        const std::string &bt, bool reverse, bool isReverseStrand,
                        bool translateSequence, const TranslateNucl &translateNucl) {
    void* record_pointer = record_batch[record_num].first;
    StringColumn* carbyte_pointer = reinterpret_cast<StringColumn*>(record_pointer);
    // write directly into the column arena instead of a temporary string
    const size_t startSize = carbyte_pointer->arena.size();
    carbyte_pointer->startValue();
    unsigned int seqPos = 0;
    char codon[3];
    for (uint32_t i = 0; i < bt.size(); ++i) {
        char seqChar = (isReverseStrand == true) ? Orf::complement(seq[offset - seqPos]) : seq[offset + seqPos];
        if (translateSequence) {
            codon[0] = (isReverseStrand == true) ? Orf::complement(seq[offset - seqPos])     : seq[offset + seqPos];
            codon[1] = (isReverseStrand == true) ? Orf::complement(seq[offset - (seqPos+1)]) : seq[offset + (seqPos+1)];
            codon[2] = (isReverseStrand == true) ? Orf::complement(seq[offset - (seqPos+2)]) : seq[offset + (seqPos+2)];
            seqChar = translateNucl.translateSingleCodon(codon);
        }
        switch (bt[i]) {
            case 'M':
                carbyte_pointer->append(seqChar);
                seqPos += (translateSequence) ?  3 : 1;
                break;
            case 'I':
                if (reverse == true) {
                    carbyte_pointer->append('-');
                } else {
                    carbyte_pointer->append(seqChar);
                    seqPos += (translateSequence) ?  3 : 1;
                }
                break;
            case 'D':
                if (reverse == true) {
                    carbyte_pointer->append(seqChar);
                    seqPos += (translateSequence) ?  3 : 1;
                } else {
                    carbyte_pointer->append('-');
                }
                break;
        }
    }
    pending_bytes += carbyte_pointer->arena.size() - startSize;
}


void ParquetDBWriter::clearRecordBatch(){
    for (const auto& Record : record_batch) {
        void* record_pointer = Record.second.first;
        carquet_physical_type type = Record.second.second;
        switch (type){
            case CARQUET_PHYSICAL_BYTE_ARRAY:{
                StringColumn* carbyte_pointer_a = reinterpret_cast<StringColumn*>(record_pointer);
                carbyte_pointer_a->clear(2 * flush_bytes);
                break;
            }
            case CARQUET_PHYSICAL_DOUBLE:{
                std::vector<double>* carbyte_pointer_b = reinterpret_cast<std::vector<double>*>(record_pointer);
                carbyte_pointer_b->clear();
                break;
            }
            case CARQUET_PHYSICAL_INT32:{
                std::vector<int32_t>* carbyte_pointer_c = reinterpret_cast<std::vector<int32_t>*>(record_pointer);
                carbyte_pointer_c->clear();
                break;
            }
            case CARQUET_PHYSICAL_INT64:{
                std::vector<int64_t>* carbyte_pointer_c = reinterpret_cast<std::vector<int64_t>*>(record_pointer);
                carbyte_pointer_c->clear();
                break;
            }
            case CARQUET_PHYSICAL_FLOAT:{
                std::vector<float>* carbyte_pointer_d = reinterpret_cast<std::vector<float>*>(record_pointer);
                carbyte_pointer_d->clear();
                break;
            }
            default:
                status = CARQUET_ERROR_INTERNAL;
                checkStatus("Clear record of unkown field type");
                break;
        }
    }
}

void ParquetDBWriter::freeRecordBatch(){
    for (const auto& record : record_batch) {
        void* record_pointer = record.second.first;
        carquet_physical_type type = record.second.second;
        switch (type){
            case CARQUET_PHYSICAL_BYTE_ARRAY:{
                StringColumn* carbyte_pointer_a = reinterpret_cast<StringColumn*>(record_pointer);
                delete carbyte_pointer_a; 
                break;
            }
            case CARQUET_PHYSICAL_DOUBLE:{
                std::vector<double>* carbyte_pointer_b = reinterpret_cast<std::vector<double>*>(record_pointer);
                delete carbyte_pointer_b;
                break;
            }
            case CARQUET_PHYSICAL_INT32:{
                std::vector<int32_t>* carbyte_pointer_c = reinterpret_cast<std::vector<int32_t>*>(record_pointer);
                delete carbyte_pointer_c;
                break;
            }
            case CARQUET_PHYSICAL_INT64:{
                std::vector<int64_t>* carbyte_pointer_c = reinterpret_cast<std::vector<int64_t>*>(record_pointer);
                delete carbyte_pointer_c;
                break;
            }
            case CARQUET_PHYSICAL_FLOAT:{
                std::vector<float>* carbyte_pointer_d = reinterpret_cast<std::vector<float>*>(record_pointer);
                delete carbyte_pointer_d;
                break;
            }
            default:
                status = CARQUET_ERROR_INTERNAL;
                checkStatus("Free record of unkown field type");
                break;
        }
    }
}


void ParquetDBWriter::writeCell(const std::string& record_value, const int record_num){
    void* record_pointer = record_batch[record_num].first;
    carquet_physical_type type = record_batch[record_num].second;
    StringColumn* carbyte_pointer_a = reinterpret_cast<StringColumn*>(record_pointer);
    carbyte_pointer_a->add(record_value);
    pending_bytes += record_value.size();
}

void ParquetDBWriter::writeCell(const double& record_value, const int record_num){
    void* record_pointer = record_batch[record_num].first;
    carquet_physical_type type = record_batch[record_num].second;
    std::vector<double>* carbyte_pointer_b = reinterpret_cast<std::vector<double>*>(record_pointer);
    carbyte_pointer_b->emplace_back(record_value);
    pending_bytes += sizeof(double);
}

void ParquetDBWriter::writeCell(const int& record_value, const int record_num){
    void* record_pointer = record_batch[record_num].first;
    carquet_physical_type type = record_batch[record_num].second;
    std::vector<int32_t>* carbyte_pointer_c = reinterpret_cast<std::vector<int32_t>*>(record_pointer);
    carbyte_pointer_c->emplace_back(record_value);
    pending_bytes += sizeof(int32_t);
}

void ParquetDBWriter::writeCell(const int64_t& record_value, const int record_num){
    void* record_pointer = record_batch[record_num].first;
    carquet_physical_type type = record_batch[record_num].second;
    std::vector<int64_t>* carbyte_pointer_c = reinterpret_cast<std::vector<int64_t>*>(record_pointer);
    carbyte_pointer_c->emplace_back(record_value);
    pending_bytes += sizeof(int64_t);
}

void ParquetDBWriter::writeCell(const float& record_value, const int record_num){
    void* record_pointer = record_batch[record_num].first;
    carquet_physical_type type = record_batch[record_num].second;
    std::vector<float>* carbyte_pointer_d = reinterpret_cast<std::vector<float>*>(record_pointer);
    carbyte_pointer_d->emplace_back(record_value);
    pending_bytes += sizeof(float);
}



void ParquetDBWriter::writeColumn(int column_number){
    
    void* record_pointer = record_batch[column_number].first;
    carquet_physical_type type = record_batch[column_number].second;
    
    switch (type){
        case CARQUET_PHYSICAL_BYTE_ARRAY: {
            StringColumn* carbyte_pointer_a = reinterpret_cast<StringColumn*>(record_pointer);
            carbyte_pointer_a->finalize(); 
            status = carquet_writer_write_batch(writer, column_number, carbyte_pointer_a->values.data(),static_cast<int64_t>(carbyte_pointer_a->values.size()), NULL, NULL);
            checkStatus("Column write" + std::to_string(column_number) + "error");
            break;
        }
        case CARQUET_PHYSICAL_DOUBLE:{
            std::vector<double>* carbyte_pointer_b = reinterpret_cast<std::vector<double>*>(record_pointer);
            status = carquet_writer_write_batch(writer, column_number, carbyte_pointer_b->data(),static_cast<int64_t>(carbyte_pointer_b->size()), NULL, NULL);
            checkStatus("Column write" + std::to_string(column_number) + "error");
            break;
        }
        case CARQUET_PHYSICAL_INT32:{
            std::vector<int32_t>* carbyte_pointer_c = reinterpret_cast<std::vector<int32_t>*>(record_pointer);
            status = carquet_writer_write_batch(writer, column_number, carbyte_pointer_c->data(),static_cast<int64_t>(carbyte_pointer_c->size()), NULL, NULL);
            checkStatus("Column write" + std::to_string(column_number) + "error");
            break;
        }
        case CARQUET_PHYSICAL_INT64:{
            std::vector<int64_t>* carbyte_pointer_c = reinterpret_cast<std::vector<int64_t>*>(record_pointer);
            status = carquet_writer_write_batch(writer, column_number, carbyte_pointer_c->data(),static_cast<int64_t>(carbyte_pointer_c->size()), NULL, NULL);
            checkStatus("Column write" + std::to_string(column_number) + "error");
            break;
        }
        case CARQUET_PHYSICAL_FLOAT:{
            std::vector<float>* carbyte_pointer_d = reinterpret_cast<std::vector<float>*>(record_pointer);
            status = carquet_writer_write_batch(writer, column_number, carbyte_pointer_d->data(),static_cast<int64_t>(carbyte_pointer_d->size()), NULL, NULL);
            checkStatus("Column write" + std::to_string(column_number) + "error");
            break;
        }
        default:
            status = carquet_status::CARQUET_ERROR_INTERNAL;;
            checkStatus("Column write" + std::to_string(column_number) + "error"); 
            break;
    }
}

void ParquetDBWriter::writeBatchToFile(const std::vector<int>& outcodes){
    for(size_t i = 0; i < outcodes.size(); i++){
        writeColumn(i);
    }
    clearRecordBatch();
    pending_bytes = 0;
    pending_rows = 0;
}

void ParquetDBWriter::endRow(){
    pending_rows++;
    if (pending_bytes >= flush_bytes) {
        flush();
    }
}

void ParquetDBWriter::flush(){
    if (pending_rows == 0 && pending_bytes == 0) {
        return;
    }
    for(size_t i = 0; i < record_batch.size(); i++){
        writeColumn(i);
    }
    clearRecordBatch();
    pending_bytes = 0;
    pending_rows = 0;
}

void ParquetDBWriter::setFlushBytes(const size_t bytes){
    flush_bytes = bytes;
}

void ParquetDBWriter::close(){
    if (writer_close_status == CARQUET_ERROR_ALREADY_CLOSED){
        Debug(Debug::WARNING) << "Writer already closed" << "\n";
        return;
    }
    flush();
    writer_close_status = carquet_writer_close(writer);
    if (writer_close_status != carquet_status_t::CARQUET_OK){
        carquet_writer_abort(writer);
        writer_close_status = CARQUET_ERROR_ALREADY_CLOSED;
        checkStatus("Closing parquet writer failed");
        
    } else {
        writer_close_status = CARQUET_ERROR_ALREADY_CLOSED;
    }
}

// ---------------------------------------------------------------------------------------------
// DB -> parquet conversion (convert2parquet)
// ---------------------------------------------------------------------------------------------

ParquetDBWriter::WriterOptions ParquetDBWriter::writerOptionsFromParameters(const Parameters& par){
    // order of --parquet-compression; LZ4_RAW, the legacy LZ4 codec is not readable by most readers (e.g. DuckDB)
    // 5 (brotli) is only offered when carquet is built with -DCARQUET_ENABLE_BROTLI=ON: extend the option's
    // description and regex ("^[0-5]{1}$") in Parameters.cpp to enable it again
    static const int codecs[] = {CARQUET_COMPRESSION_UNCOMPRESSED, CARQUET_COMPRESSION_ZSTD, CARQUET_COMPRESSION_SNAPPY,
                                 CARQUET_COMPRESSION_GZIP, CARQUET_COMPRESSION_LZ4_RAW, CARQUET_COMPRESSION_BROTLI};
    WriterOptions options;
    options.compression = codecs[par.parquetCompression];
    options.compressionLevel = par.parquetCompressionLevel;
    options.rowGroupSize = par.parquetRowGroupSize;
    return options;
}

void ParquetDBWriter::setWriterOptions(const WriterOptions& options){
    writerOptions = options;
    // carquet starts a new row group only between batches, so a batch must not exceed the row group size
    flush_bytes = std::min<size_t>(DEFAULT_FLUSH_BYTES, options.rowGroupSize);
}

bool ParquetDBWriter::isSupportedDbType(const int dbtype){
    switch (dbtype & Parameters::DBTYPE_MASK) {
        case Parameters::DBTYPE_AMINO_ACIDS:
        case Parameters::DBTYPE_NUCLEOTIDES:
        case Parameters::DBTYPE_HMM_PROFILE:
        case Parameters::DBTYPE_ALIGNMENT_RES:
        case Parameters::DBTYPE_CLUSTER_RES:
        case Parameters::DBTYPE_PREFILTER_RES:
        case Parameters::DBTYPE_PREFILTER_REV_RES:
        case Parameters::DBTYPE_TAXONOMICAL_RESULT:
        case Parameters::DBTYPE_MSA_DB:
        case Parameters::DBTYPE_GENERIC_DB:
            return true;
        default:
            // CA3M (binary msa column still missing), INDEX_DB (binary k-mer tables, single values of
            // hundreds of MB that carquet does not compress), OMIT_FILE, OFFSETDB and the validation-only types
            return false;
    }
}

static bool isLineageField(const char* s, const char* end){
    // lineage is ';' separated and either all taxids (--tax-lineage 2) or all "<rank char>_<name>" (--tax-lineage 1)
    if (s == end) {
        return false;
    }
    bool allNumeric = true;
    bool allRankPrefixed = true;
    const char* e = s;
    while (e <= end) {
        if (e == end || *e == ';') {
            bool numeric = (e > s);
            for (const char* c = s; c < e; c++) {
                numeric = numeric && (*c >= '0' && *c <= '9');
            }
            allNumeric = allNumeric && numeric;
            allRankPrefixed = allRankPrefixed && (e - s >= 2 && s[1] == '_');
            s = e + 1;
        }
        e++;
    }
    return allNumeric || allRankPrefixed;
}

// taxonomy layout of one entry: taxid rank name [num_seqs assigned agree support] [ranks] [lineage]
unsigned int ParquetDBWriter::detectTaxonomyColumns(const char* entry){
    std::vector<std::pair<const char*, const char*>> c;
    const char* start = entry;
    for (const char* p = entry; ; p++) {
        if (*p == '\t' || *p == '\n' || *p == '\0') {
            c.emplace_back(start, p);
            if (*p != '\t') break;
            start = p + 1;
        }
    }
    auto isInteger = [](const std::pair<const char*, const char*>& f) {
        if (f.first == f.second) return false;
        for (const char* p = f.first; p < f.second; p++) if (*p < '0' || *p > '9') return false;
        return true;
    };
    unsigned int optional = DBOPT_NONE;
    size_t base = 3;
    if (c.size() >= 7 && isInteger(c[3]) && isInteger(c[4]) && isInteger(c[5])) {
        optional |= DBOPT_TAX_AGGREGATE;
        base = 7;
    }
    const size_t extra = c.size() > base ? c.size() - base : 0;
    if (extra >= 2) {
        optional |= DBOPT_TAX_RANKS | DBOPT_TAX_LINEAGE;
    } else if (extra == 1) {
        optional |= isLineageField(c[base].first, c[base].second) ? DBOPT_TAX_LINEAGE : DBOPT_TAX_RANKS;
    }
    return optional;
}

// strict integer parse: optional '-', digits only, within [lo, hi]
static bool parseInteger(const char* s, const size_t n, const int64_t lo, const int64_t hi, int64_t& value){
    size_t i = 0;
    bool negative = false;
    if (i < n && s[i] == '-') {
        negative = true;
        i++;
    }
    if (i == n) {
        return false;
    }
    uint64_t v = 0;
    for (; i < n; i++) {
        if (s[i] < '0' || s[i] > '9' || v > (UINT64_MAX - 9) / 10) {
            return false;
        }
        v = v * 10 + static_cast<uint64_t>(s[i] - '0');
    }
    if (negative) {
        if (v > static_cast<uint64_t>(INT64_MAX) + 1) return false;
        value = (v == static_cast<uint64_t>(INT64_MAX) + 1) ? INT64_MIN : -static_cast<int64_t>(v);
    } else {
        if (v > static_cast<uint64_t>(INT64_MAX)) return false;
        value = static_cast<int64_t>(v);
    }
    return value >= lo && value <= hi;
}

void ParquetDBWriter::writeField(const int column, const char* s, const size_t n){
    const ColumnType type = columnTypes[column];
    bool ok = true;
    int64_t i = 0;
    switch (type) {
        case COL_INT8:   ok = parseInteger(s, n, INT8_MIN, INT8_MAX, i);   if (ok) writeCell(static_cast<int>(i), column); break;
        case COL_UINT8:  ok = parseInteger(s, n, 0, UINT8_MAX, i);         if (ok) writeCell(static_cast<int>(i), column); break;
        case COL_INT16:  ok = parseInteger(s, n, INT16_MIN, INT16_MAX, i); if (ok) writeCell(static_cast<int>(i), column); break;
        case COL_INT32:  ok = parseInteger(s, n, INT32_MIN, INT32_MAX, i); if (ok) writeCell(static_cast<int>(i), column); break;
        case COL_INT64:  ok = parseInteger(s, n, INT64_MIN, INT64_MAX, i); if (ok) writeCell(i, column); break;
        case COL_FLOAT: {
            float v = 0;
            const fast_float::from_chars_result r = fast_float::from_chars(s, s + n, v);
            ok = (n > 0 && r.ec == std::errc() && r.ptr == s + n);
            if (ok) writeCell(v, column);
            break;
        }
        case COL_DOUBLE: {
            double v = 0;
            const fast_float::from_chars_result r = fast_float::from_chars(s, s + n, v);
            ok = (n > 0 && r.ec == std::errc() && r.ptr == s + n);
            if (ok) writeCell(v, column);
            break;
        }
        case COL_STRING:
        case COL_BINARY:
            writeCell(std::string(s, n), column);
            break;
    }
    if (ok == false) {
        Debug(Debug::ERROR) << "Entry " << currentKey << ": value '" << std::string(s, n) << "' of column "
                            << columnNames[column] << " is not a valid " << columnTypeName(type) << "\n";
        EXIT(EXIT_FAILURE);
    }
}

void ParquetDBWriter::writeEntry(const DBKeyType key, const char* data, const size_t len){
    const int64_t k = static_cast<int64_t>(key);
    currentKey = key;
    const int type = parseLines ? -1 : (dbtype & Parameters::DBTYPE_MASK);
    switch (type) {
        case Parameters::DBTYPE_AMINO_ACIDS:
        case Parameters::DBTYPE_NUCLEOTIDES: {
            size_t l = len;
            while (l > 0 && data[l - 1] == '\n') {
                l--;
            }
            put(0, k);
            put(1, std::string(data, l));
            endRow();
            break;
        }
        case Parameters::DBTYPE_GENERIC_DB:
        case Parameters::DBTYPE_MSA_DB:
        case Parameters::DBTYPE_INDEX_DB:
            put(0, k);
            put(1, std::string(data, len));
            endRow();
            break;
        case Parameters::DBTYPE_HMM_PROFILE: {
            // 25 bytes per position: 20 scores (int8), query residue, consensus, neff, gap del, gap ins (uint8)
            const size_t positions = len / Sequence::PROFILE_READIN_SIZE;
            for (size_t pos = 0; pos < positions; pos++) {
                const char* p = data + pos * Sequence::PROFILE_READIN_SIZE;
                int col = 0;
                put(col++, k);
                put(col++, static_cast<int>(pos));
                for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
                    put(col++, static_cast<int>(static_cast<signed char>(p[aa])));
                }
                for (size_t b = Sequence::PROFILE_AA_SIZE; b < Sequence::PROFILE_READIN_SIZE; b++) {
                    put(col++, static_cast<int>(static_cast<unsigned char>(p[b])));
                }
                endRow();
            }
            break;
        }
        default: {
            // tab separated lines: alignment, prefilter, cluster, taxonomy, and generic DBs with a schema;
            // native column 0 is the entry key, column c the (c-1)-th field of the line
            const char* p = data;
            const char* end = data + len;
            while (p < end && *p != '\0') {
                if (*p == '\n') {
                    p++;    // empty line
                    continue;
                }
                put(0, k);
                int col = 1;
                const char* start = p;
                while (true) {
                    if (p == end || *p == '\t' || *p == '\n' || *p == '\0') {
                        if (col < nativeColumns) {
                            const int out = projection.empty() ? col : projection[col];
                            if (out >= 0) {
                                writeField(out, start, p - start);
                            }
                        }
                        col++;
                        if (p == end || *p != '\t') {
                            break;
                        }
                        start = p + 1;
                    }
                    p++;
                }
                if (col != nativeColumns) {
                    Debug(Debug::ERROR) << "Entry " << key << " has a line with " << (col - 1) << " columns, expected "
                                        << (nativeColumns - 1) << "\n";
                    EXIT(EXIT_FAILURE);
                }
                endRow();
                if (p < end && *p == '\n') {
                    p++;
                }
            }
            break;
        }
    }
}

const char* ParquetDBWriter::columnTypeName(const ColumnType type){
    switch (type) {
        case COL_INT8:   return "int8";
        case COL_UINT8:  return "uint8";
        case COL_INT16:  return "int16";
        case COL_INT32:  return "int32";
        case COL_INT64:  return "int64";
        case COL_FLOAT:  return "float";
        case COL_DOUBLE: return "double";
        case COL_STRING: return "string";
        case COL_BINARY: return "binary";
    }
    return "unknown";
}

bool ParquetDBWriter::parseSchema(const std::string& text, std::vector<ColumnSpec>& columns, std::string& error){
    struct TypeName { const char* name; ColumnType type; };
    static const TypeName typeNames[] = {
        {"int8", COL_INT8}, {"tinyint", COL_INT8}, {"uint8", COL_UINT8}, {"utinyint", COL_UINT8},
        {"int16", COL_INT16}, {"smallint", COL_INT16}, {"int32", COL_INT32}, {"int", COL_INT32},
        {"integer", COL_INT32}, {"int64", COL_INT64}, {"bigint", COL_INT64}, {"float", COL_FLOAT},
        {"real", COL_FLOAT}, {"double", COL_DOUBLE}, {"string", COL_STRING}, {"varchar", COL_STRING},
        {"text", COL_STRING}
    };
    columns.clear();
    std::string t = text;
    // trim whitespace and the surrounding braces
    t.erase(0, t.find_first_not_of(" \t\n"));
    t.erase(t.find_last_not_of(" \t\n") + 1);
    if (t.size() < 2 || t.front() != '{' || t.back() != '}') {
        error = "schema must look like '{name:type, name:type, ...}'";
        return false;
    }
    t = t.substr(1, t.size() - 2);
    std::vector<std::string> entries = Util::split(t, ",");
    for (size_t i = 0; i < entries.size(); i++) {
        std::string e = entries[i];
        e.erase(0, e.find_first_not_of(" \t\n"));
        e.erase(e.find_last_not_of(" \t\n") + 1);
        const size_t colon = e.find(':');
        if (colon == std::string::npos) {
            error = "column '" + e + "' has no ':type'";
            return false;
        }
        std::string name = e.substr(0, colon);
        std::string type = e.substr(colon + 1);
        name.erase(name.find_last_not_of(" \t") + 1);
        type.erase(0, type.find_first_not_of(" \t"));
        if (name.empty() || (isalpha(name[0]) == 0 && name[0] != '_')) {
            error = "invalid column name '" + name + "'";
            return false;
        }
        for (size_t c = 0; c < name.size(); c++) {
            if (isalnum(name[c]) == 0 && name[c] != '_') {
                error = "invalid column name '" + name + "'";
                return false;
            }
        }
        std::transform(type.begin(), type.end(), type.begin(), ::tolower);
        bool found = false;
        for (size_t n = 0; n < sizeof(typeNames) / sizeof(typeNames[0]); n++) {
            if (type == typeNames[n].name) {
                columns.push_back(ColumnSpec{name, typeNames[n].type});
                found = true;
                break;
            }
        }
        if (found == false) {
            error = "unknown type '" + type + "' of column '" + name + "' (int8, uint8, int16, int32, int64, float, double, string)";
            return false;
        }
        for (size_t c = 0; c + 1 < columns.size(); c++) {
            if (columns[c].name == name) {
                error = "column '" + name + "' is given twice";
                return false;
            }
        }
    }
    if (columns.empty()) {
        error = "schema has no columns";
        return false;
    }
    return true;
}

void ParquetDBWriter::initFromSchema(const char* resultFile, const int dbtype, const unsigned int optional,
                                     const std::vector<ColumnSpec>& schemaColumns){
    this->dbtype = dbtype;
    fields.clear();
    std::vector<ColumnSpec> out;
    if (Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_GENERIC_DB)) {
        // the entry key is always the first column, the schema describes the fields of each line
        for (size_t i = 0; i < schemaColumns.size(); i++) {
            if (schemaColumns[i].name == "key") {
                Debug(Debug::ERROR) << "Column name 'key' is reserved for the entry key of generic DBs\n";
                EXIT(EXIT_FAILURE);
            }
        }
        out.push_back(ColumnSpec{"key", COL_INT64});
        out.insert(out.end(), schemaColumns.begin(), schemaColumns.end());
        nativeColumns = static_cast<int>(out.size());
        projection.clear();
        parseLines = true;
    } else {
        // select columns of the layout by name; the output keeps the schema order
        const std::vector<int> layout = dbTypeToFields(dbtype, optional);
        nativeColumns = static_cast<int>(layout.size());
        projection.assign(layout.size(), -1);
        for (size_t i = 0; i < schemaColumns.size(); i++) {
            int native = -1;
            for (size_t j = 0; j < layout.size(); j++) {
                if (fieldSpec(layout[j]).name == schemaColumns[i].name) {
                    native = static_cast<int>(j);
                    break;
                }
            }
            if (native < 0) {
                std::string names;
                for (size_t j = 0; j < layout.size(); j++) {
                    names += (j > 0 ? ", " : "") + fieldSpec(layout[j]).name;
                }
                Debug(Debug::ERROR) << "Column '" << schemaColumns[i].name << "' is not part of the "
                                    << Parameters::getDbTypeName(dbtype) << " layout: " << names << "\n";
                EXIT(EXIT_FAILURE);
            }
            const ColumnSpec spec = fieldSpec(layout[native]);
            if (spec.type != schemaColumns[i].type) {
                Debug(Debug::ERROR) << "Column '" << spec.name << "' is " << columnTypeName(spec.type)
                                    << ", conversion to " << columnTypeName(schemaColumns[i].type)
                                    << " is not supported for non-generic DBs\n";
                EXIT(EXIT_FAILURE);
            }
            projection[native] = static_cast<int>(i);
            out.push_back(spec);
        }
        parseLines = false;
    }
    initFromSpecs(resultFile, out);
}

static size_t availableMemory(){
    // MemAvailable counts reclaimable page cache, unlike free memory
    FILE* fp = fopen("/proc/meminfo", "r");
    if (fp != NULL) {
        char line[256];
        while (fgets(line, sizeof(line), fp) != NULL) {
            unsigned long long kb = 0;
            if (sscanf(line, "MemAvailable: %llu kB", &kb) == 1) {
                fclose(fp);
                return static_cast<size_t>(kb) * 1024;
            }
        }
        fclose(fp);
    }
    return Util::getTotalSystemMemory() / 2;
}

static inline size_t entryLength(DBReader<DBKeyType>& reader, size_t id, const char* data){
    // compressed DBs: the index length is the compressed size, the decompressed entry ends with \0
    return (reader.isCompressed() == DBReader<DBKeyType>::COMPRESSED) ? strlen(data) : reader.getEntryLen(id) - 1;
}

// starts of `parts` ranges of roughly equal data size over the given entry order
static std::vector<size_t> splitByBytes(DBReader<DBKeyType>& reader, const std::vector<size_t>* order, size_t parts){
    std::vector<size_t> bounds(1, 0);
    const size_t n = reader.getSize();
    const size_t target = reader.getTotalDataSize() / parts + 1;
    size_t acc = 0;
    for (size_t i = 0; i < n && bounds.size() < parts; i++) {
        acc += reader.getEntryLen(order != NULL ? (*order)[i] : i);
        if (acc >= target * bounds.size()) {
            bounds.push_back(i + 1);
        }
    }
    while (bounds.size() < parts + 1) {
        bounds.push_back(n);
    }
    bounds.back() = n;
    return bounds;
}

std::string ParquetDBWriter::lookupParquetName(const std::string& output, const bool multiFile){
    if (multiFile) {
        return output + "/_lookup.pq";
    }
    const size_t slash = output.find_last_of('/');
    const std::string dir = (slash == std::string::npos) ? "" : output.substr(0, slash + 1);
    const std::string base = (slash == std::string::npos) ? output : output.substr(slash + 1);
    return dir + "_" + base + ".lookup.pq";
}

// remove an earlier result of the same name: <output> file or <output>/part-*.parquet, and its lookup
static void removeStaleParquet(const std::string& output){
    const std::string singleLookup = ParquetDBWriter::lookupParquetName(output, false);
    if (FileUtil::fileExists(singleLookup.c_str())) {
        FileUtil::remove(singleLookup.c_str());
    }
    if (FileUtil::directoryExists(output.c_str())) {
        DIR* dir = opendir(output.c_str());
        if (dir != NULL) {
            struct dirent* ent;
            while ((ent = readdir(dir)) != NULL) {
                const std::string name = ent->d_name;
                if ((name.compare(0, 5, "part-") == 0 && name.find(".parquet") != std::string::npos)
                    || name.compare(0, 10, "_lookup.pq") == 0) {
                    FileUtil::remove((output + "/" + name).c_str());
                }
            }
            closedir(dir);
        }
        rmdir(output.c_str());   // only succeeds when nothing else is in it
    } else if (FileUtil::fileExists(output.c_str())) {
        FileUtil::remove(output.c_str());
    }
}

void ParquetDBWriter::removeNumberedParquet(const std::string& dir, const std::string& prefix){
    DIR* d = opendir(dir.c_str());
    if (d == NULL) {
        return;
    }
    struct dirent* ent;
    while ((ent = readdir(d)) != NULL) {
        const std::string name = ent->d_name;
        if (name.size() <= prefix.size() + 1 || name.compare(0, prefix.size() + 1, prefix + "_") != 0) {
            continue;
        }
        // <prefix>_<digits>.parquet or <prefix>_<digits>.parquet.tmp
        size_t p = prefix.size() + 1;
        const size_t digits = p;
        while (p < name.size() && isdigit(name[p])) {
            p++;
        }
        const std::string rest = name.substr(p);
        if (p > digits && (rest == ".parquet" || rest == ".parquet.tmp")) {
            FileUtil::remove((dir + "/" + name).c_str());
        }
    }
    closedir(d);
}

// key, accession, file_number of every lookup entry (reader opened with USE_LOOKUP)
static void writeLookupParquet(DBReader<DBKeyType>& reader, const std::string& name,
                               const ParquetDBWriter::WriterOptions& options){
    const std::string tmp = name + ".tmp";
    ParquetDBWriter w;
    w.setWriterOptions(options);
    w.initFromFields(tmp.c_str(), {ParquetDBWriter::DBFIELD_KEY, ParquetDBWriter::DBFIELD_ACCESSION,
                                   ParquetDBWriter::DBFIELD_FILE_NUMBER});
    DBReader<DBKeyType>::LookupEntry* lookup = reader.getLookup();
    for (size_t i = 0; i < reader.getLookupSize(); i++) {
        w.writeCell(static_cast<int64_t>(lookup[i].id), 0);
        w.writeCell(lookup[i].entryName, 1);
        w.writeCell(static_cast<int64_t>(lookup[i].fileNumber), 2);
        w.endRow();
    }
    w.close();
    FileUtil::move(tmp.c_str(), name.c_str());
}

void ParquetDBWriter::writeDbAsParquet(DBReader<DBKeyType>& reader, const std::string& output, const int sortMode,
                                       const int threads, const WriterOptions& options, const bool withLookup,
                                       const std::vector<ColumnSpec>& schemaColumns, const ConvertHooks& hooks){
    const size_t TARGET_FILE_BYTES = 2ULL * 1024 * 1024 * 1024;   // raw DB bytes per output file
    const size_t MIN_FILE_BYTES = 256ULL * 1024 * 1024;           // smallest file worth a separate worker
    const size_t SUB_BYTES = 256ULL * 1024 * 1024;                // spillsort key sub-range
    const size_t WORKER_BYTES = 512ULL * 1024 * 1024;             // memory estimate per writer

    const size_t rawBytes = reader.getTotalDataSize();
    const size_t entries = reader.getSize();
    const unsigned int optional = detectOptionalColumns(reader, 0);

    int mode = sortMode;
    if (mode == SORT_MODE_AUTO) {
        mode = (rawBytes < availableMemory() / 2) ? SORT_MODE_KEYORDER : SORT_MODE_SPILLSORT;
    }
    size_t files = std::max((rawBytes + TARGET_FILE_BYTES - 1) / TARGET_FILE_BYTES,
                            std::min<size_t>(threads, (rawBytes + MIN_FILE_BYTES - 1) / MIN_FILE_BYTES));
    files = std::max<size_t>(1, std::min(files, entries));
    files = std::max(files, std::min(hooks.minFiles, entries));
    const size_t perWorker = (mode == SORT_MODE_SPILLSORT) ? 3 * SUB_BYTES : WORKER_BYTES;
    const int workers = static_cast<int>(std::max<size_t>(1, std::min<size_t>(std::min<size_t>(threads, files), availableMemory() / perWorker)));

    std::vector<std::string> names;
    if (hooks.partNames) {
        names = hooks.partNames(files);
    } else if (files == 1) {
        removeStaleParquet(output);
        names.push_back(output);
    } else {
        removeStaleParquet(output);
        FileUtil::makeDir(output.c_str());
        for (size_t f = 0; f < files; f++) {
            char buf[32];
            snprintf(buf, sizeof(buf), "/part-%05zu.parquet", f);
            names.push_back(output + buf);
        }
    }
    const char* modeName[] = {"auto", "linear", "keyorder", "spillsort"};
    Debug(Debug::INFO) << "Writing " << files << " parquet file(s) with " << workers << " worker(s), sort mode "
                       << modeName[mode] << "\n";
    if (withLookup) {
        const std::string lookupName = lookupParquetName(output, files > 1);
        Debug(Debug::INFO) << "Writing lookup to " << lookupName << "\n";
        writeLookupParquet(reader, lookupName, options);
    }

    // each output file is written by one worker; its writer only lives while the file is written
    auto openPart = [&](ParquetDBWriter& w, size_t f) {
        w.setWriterOptions(options);
        if (hooks.openPart) {
            hooks.openPart(w, names[f] + ".tmp");
        } else if (schemaColumns.empty()) {
            w.initFromDbType((names[f] + ".tmp").c_str(), reader.getDbtype(), optional);
        } else {
            w.initFromSchema((names[f] + ".tmp").c_str(), reader.getDbtype(), optional, schemaColumns);
        }
    };
    auto writeEntryRows = [&](ParquetDBWriter& w, DBKeyType key, char* data, size_t len, unsigned int thread_idx) {
        if (hooks.writeEntry) {
            hooks.writeEntry(w, key, data, len, thread_idx);
        } else {
            w.writeEntry(key, data, len);
        }
    };
    auto closePart = [&](ParquetDBWriter& w, size_t f) {
        w.close();
        FileUtil::move((names[f] + ".tmp").c_str(), names[f].c_str());
    };

    Debug::Progress progress(entries);
    if (mode == SORT_MODE_KEYORDER) {
        // reader is sorted by key: contiguous id ranges are key ranges
        const std::vector<size_t> bounds = splitByBytes(reader, NULL, files);
#pragma omp parallel for schedule(dynamic, 1) num_threads(workers)
        for (size_t f = 0; f < files; f++) {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            ParquetDBWriter w;
            openPart(w, f);
            for (size_t id = bounds[f]; id < bounds[f + 1]; id++) {
                progress.updateProgress();
                char* data = reader.getData(id, thread_idx);
                writeEntryRows(w, reader.getDbKey(id), data, entryLength(reader, id, data), thread_idx);
            }
            closePart(w, f);
        }
        return;
    }

    // linear and spillsort read the data in file order
    std::vector<size_t> byOffset(entries);
    for (size_t i = 0; i < entries; i++) {
        byOffset[i] = i;
    }
    std::sort(byOffset.begin(), byOffset.end(),
              [&](size_t a, size_t b) { return reader.getOffset(a) < reader.getOffset(b); });
    reader.setSequentialAdvice();

    if (mode == SORT_MODE_LINEAR) {
        const std::vector<size_t> bounds = splitByBytes(reader, &byOffset, files);
#pragma omp parallel for schedule(dynamic, 1) num_threads(workers)
        for (size_t f = 0; f < files; f++) {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            ParquetDBWriter w;
            openPart(w, f);
            for (size_t i = bounds[f]; i < bounds[f + 1]; i++) {
                progress.updateProgress();
                const size_t id = byOffset[i];
                char* data = reader.getData(id, thread_idx);
                writeEntryRows(w, reader.getDbKey(id), data, entryLength(reader, id, data), thread_idx);
            }
            closePart(w, f);
        }
        return;
    }

    // spillsort phase 1: scatter entries in file order into one temporary file per key sub-range
    // record: uint64 key, uint64 length, data
    const size_t subs = std::max(files, (rawBytes + SUB_BYTES - 1) / SUB_BYTES);
    const std::vector<size_t> subBounds = splitByBytes(reader, NULL, subs);
    std::vector<DBKeyType> subFirstKey(subs);
    for (size_t s = 0; s < subs; s++) {
        subFirstKey[s] = reader.getDbKey(std::min(subBounds[s], entries - 1));
    }
    const std::string spillDir = output + ".spill";
    FileUtil::makeDir(spillDir.c_str());
    struct Spill {
        std::mutex lock;
        FILE* fp;
        std::vector<char> buf;
    };
    std::vector<Spill> spills(subs);
    std::vector<std::string> spillNames(subs);
    for (size_t s = 0; s < subs; s++) {
        spillNames[s] = spillDir + "/sub-" + SSTR(s) + ".bin";
        spills[s].fp = FileUtil::openFileOrDie(spillNames[s].c_str(), "wb", false);
    }
    const size_t SPILL_BUFFER = 4 * 1024 * 1024;
#pragma omp parallel for schedule(dynamic, 1000) num_threads(threads)
    for (size_t i = 0; i < entries; i++) {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        progress.updateProgress();
        const size_t id = byOffset[i];
        const uint64_t key = reader.getDbKey(id);
        const char* data = reader.getData(id, thread_idx);
        const uint64_t len = entryLength(reader, id, data);
        const size_t s = std::upper_bound(subFirstKey.begin(), subFirstKey.end(), key) - subFirstKey.begin() - 1;
        std::lock_guard<std::mutex> guard(spills[s].lock);
        std::vector<char>& buf = spills[s].buf;
        buf.insert(buf.end(), reinterpret_cast<const char*>(&key), reinterpret_cast<const char*>(&key) + sizeof(key));
        buf.insert(buf.end(), reinterpret_cast<const char*>(&len), reinterpret_cast<const char*>(&len) + sizeof(len));
        buf.insert(buf.end(), data, data + len);
        buf.push_back('\0');   // entries handed to the hooks are \0 terminated
        if (buf.size() >= SPILL_BUFFER) {
            fwrite(buf.data(), 1, buf.size(), spills[s].fp);
            buf.clear();
        }
    }
    for (size_t s = 0; s < subs; s++) {
        fwrite(spills[s].buf.data(), 1, spills[s].buf.size(), spills[s].fp);
        std::vector<char>().swap(spills[s].buf);
        fclose(spills[s].fp);
    }

    // phase 2: each file gets consecutive sub-ranges; load (next one prefetched), sort by key, write
    auto loadSpill = [&](size_t s, std::vector<char>& out) {
        FILE* fp = FileUtil::openFileOrDie(spillNames[s].c_str(), "rb", true);
        fseek(fp, 0, SEEK_END);
        const size_t size = ftell(fp);
        fseek(fp, 0, SEEK_SET);
        out.resize(size);
        if (size > 0 && fread(out.data(), 1, size, fp) != size) {
            Debug(Debug::ERROR) << "Could not read " << spillNames[s] << "\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(fp);
    };
    Debug(Debug::INFO) << "Sorting and writing " << subs << " key ranges\n";
#pragma omp parallel for schedule(dynamic, 1) num_threads(workers)
    for (size_t f = 0; f < files; f++) {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        struct Ent { uint64_t key; size_t off; uint64_t len; };
        std::vector<char> current, next;
        std::vector<Ent> ents;
        const size_t first = f * subs / files;
        const size_t last = (f + 1) * subs / files;
        ParquetDBWriter w;
        openPart(w, f);
        loadSpill(first, current);
        for (size_t s = first; s < last; s++) {
            std::thread prefetch;
            if (s + 1 < last) {
                prefetch = std::thread(loadSpill, s + 1, std::ref(next));
            }
            ents.clear();
            for (size_t o = 0; o < current.size();) {
                Ent e;
                memcpy(&e.key, current.data() + o, sizeof(e.key));
                memcpy(&e.len, current.data() + o + sizeof(e.key), sizeof(e.len));
                e.off = o + sizeof(e.key) + sizeof(e.len);
                ents.push_back(e);
                o = e.off + e.len + 1;
            }
            std::sort(ents.begin(), ents.end(), [](const Ent& a, const Ent& b) { return a.key < b.key; });
            for (size_t i = 0; i < ents.size(); i++) {
                writeEntryRows(w, static_cast<DBKeyType>(ents[i].key), current.data() + ents[i].off, ents[i].len, thread_idx);
            }
            if (prefetch.joinable()) {
                prefetch.join();
            }
            current.swap(next);
            FileUtil::remove(spillNames[s].c_str());
        }
        closePart(w, f);
    }
    rmdir(spillDir.c_str());
}

