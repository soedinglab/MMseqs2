/**
 * @file fuzz_writer.c
 * @brief Fuzz target for carquet Parquet writer
 *
 * Tests the writer with random flat schemas, the full primitive type surface,
 * logical annotations, opt-in encodings, compression codecs, nullable columns,
 * metadata knobs, buffer/file outputs, and readback through the newer reader
 * APIs.
 */

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <carquet/carquet.h>

#define MAX_COLUMNS 16
#define MAX_ROWS 500
#define MAX_STRING_LEN 128
#define MAX_TYPE_LENGTH 16

typedef struct {
    const uint8_t* data;
    size_t size;
    size_t pos;
} fuzz_input_t;

static uint8_t consume_byte(fuzz_input_t* in) {
    if (in->pos >= in->size) return 0;
    return in->data[in->pos++];
}
static uint16_t consume_u16(fuzz_input_t* in) {
    uint16_t lo = consume_byte(in), hi = consume_byte(in);
    return (hi << 8) | lo;
}
static uint32_t consume_u32(fuzz_input_t* in) {
    uint32_t a = consume_byte(in), b = consume_byte(in);
    uint32_t c = consume_byte(in), d = consume_byte(in);
    return a | (b << 8) | (c << 16) | (d << 24);
}
static uint64_t consume_u64(fuzz_input_t* in) {
    uint64_t lo = consume_u32(in), hi = consume_u32(in);
    return lo | (hi << 32);
}
static float consume_float(fuzz_input_t* in) {
    union { uint32_t u; float f; } v; v.u = consume_u32(in); return v.f;
}
static double consume_double(fuzz_input_t* in) {
    union { uint64_t u; double d; } v; v.u = consume_u64(in); return v.d;
}

static const carquet_physical_type_t FUZZ_TYPES[] = {
    CARQUET_PHYSICAL_BOOLEAN,
    CARQUET_PHYSICAL_INT32,
    CARQUET_PHYSICAL_INT64,
    CARQUET_PHYSICAL_INT96,
    CARQUET_PHYSICAL_FLOAT,
    CARQUET_PHYSICAL_DOUBLE,
    CARQUET_PHYSICAL_BYTE_ARRAY,
    CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY,
};
#define NUM_FUZZ_TYPES (sizeof(FUZZ_TYPES) / sizeof(FUZZ_TYPES[0]))

static const carquet_compression_t FUZZ_CODECS[] = {
    CARQUET_COMPRESSION_UNCOMPRESSED,
    CARQUET_COMPRESSION_SNAPPY,
    CARQUET_COMPRESSION_GZIP,
    CARQUET_COMPRESSION_LZ4,
    CARQUET_COMPRESSION_ZSTD,
};
#define NUM_FUZZ_CODECS (sizeof(FUZZ_CODECS) / sizeof(FUZZ_CODECS[0]))

typedef struct {
    carquet_physical_type_t type;
    carquet_encoding_t encoding;
    carquet_logical_type_t logical_type;
    bool has_logical_type;
    bool nullable;
    int32_t type_length;
} fuzz_column_plan_t;

static int32_t choose_type_length(fuzz_input_t* in, carquet_physical_type_t type) {
    static const int32_t lengths[] = { 2, 6, 8, 12, 16 };
    if (type != CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY) return 0;
    return lengths[consume_byte(in) % (sizeof(lengths) / sizeof(lengths[0]))];
}

static const carquet_logical_type_t* choose_logical_type(
    fuzz_input_t* in,
    carquet_physical_type_t type,
    int32_t type_length,
    carquet_logical_type_t* out) {

    memset(out, 0, sizeof(*out));
    uint8_t selector = consume_byte(in);

    switch (type) {
        case CARQUET_PHYSICAL_INT32:
            switch (selector % 5) {
                case 0: return NULL;
                case 1: out->id = CARQUET_LOGICAL_DATE; break;
                case 2:
                    out->id = CARQUET_LOGICAL_TIME;
                    out->params.time.unit = CARQUET_TIME_UNIT_MILLIS;
                    out->params.time.is_adjusted_to_utc = (selector & 0x80) != 0;
                    break;
                case 3:
                    out->id = CARQUET_LOGICAL_INTEGER;
                    out->params.integer.bit_width = (selector & 1) ? 16 : 32;
                    out->params.integer.is_signed = (selector & 2) != 0;
                    break;
                default:
                    out->id = CARQUET_LOGICAL_DECIMAL;
                    out->params.decimal.precision = 9;
                    out->params.decimal.scale = selector % 5;
                    break;
            }
            return out;
        case CARQUET_PHYSICAL_INT64:
            switch (selector % 5) {
                case 0: return NULL;
                case 1:
                    out->id = CARQUET_LOGICAL_TIME;
                    out->params.time.unit = (selector & 1)
                        ? CARQUET_TIME_UNIT_MICROS : CARQUET_TIME_UNIT_NANOS;
                    out->params.time.is_adjusted_to_utc = (selector & 0x80) != 0;
                    break;
                case 2:
                    out->id = CARQUET_LOGICAL_TIMESTAMP;
                    out->params.timestamp.unit = (selector & 1)
                        ? CARQUET_TIME_UNIT_MICROS : CARQUET_TIME_UNIT_NANOS;
                    out->params.timestamp.is_adjusted_to_utc = (selector & 0x80) != 0;
                    break;
                case 3:
                    out->id = CARQUET_LOGICAL_INTEGER;
                    out->params.integer.bit_width = 64;
                    out->params.integer.is_signed = (selector & 2) != 0;
                    break;
                default:
                    out->id = CARQUET_LOGICAL_DECIMAL;
                    out->params.decimal.precision = 18;
                    out->params.decimal.scale = selector % 9;
                    break;
            }
            return out;
        case CARQUET_PHYSICAL_BYTE_ARRAY:
            switch (selector % 7) {
                case 0: return NULL;
                case 1: out->id = CARQUET_LOGICAL_STRING; break;
                case 2: out->id = CARQUET_LOGICAL_ENUM; break;
                case 3: out->id = CARQUET_LOGICAL_JSON; break;
                case 4: out->id = CARQUET_LOGICAL_BSON; break;
                case 5:
                    out->id = CARQUET_LOGICAL_GEOMETRY;
                    memcpy(out->params.geometry.crs, "OGC:CRS84", 10);
                    break;
                default:
                    out->id = CARQUET_LOGICAL_GEOGRAPHY;
                    memcpy(out->params.geography.crs, "OGC:CRS84", 10);
                    out->params.geography.has_algorithm = true;
                    out->params.geography.algorithm =
                        (carquet_geospatial_edge_algorithm_t)(selector % 5);
                    break;
            }
            return out;
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
            if (type_length == 2 && (selector & 1)) {
                out->id = CARQUET_LOGICAL_FLOAT16;
                return out;
            }
            if (type_length == 12 && (selector & 1)) {
                out->id = CARQUET_LOGICAL_INTERVAL;
                return out;
            }
            if (type_length == 16 && (selector & 1)) {
                out->id = CARQUET_LOGICAL_UUID;
                return out;
            }
            return NULL;
        default:
            return NULL;
    }
}

static carquet_encoding_t choose_encoding(
    fuzz_input_t* in,
    carquet_physical_type_t type) {

    switch (type) {
        case CARQUET_PHYSICAL_INT32:
        case CARQUET_PHYSICAL_INT64: {
            static const carquet_encoding_t choices[] = {
                CARQUET_ENCODING_PLAIN,
                CARQUET_ENCODING_DELTA_BINARY_PACKED,
                CARQUET_ENCODING_BYTE_STREAM_SPLIT,
                CARQUET_ENCODING_RLE_DICTIONARY,
            };
            return choices[consume_byte(in) % (sizeof(choices) / sizeof(choices[0]))];
        }
        case CARQUET_PHYSICAL_FLOAT:
        case CARQUET_PHYSICAL_DOUBLE: {
            static const carquet_encoding_t choices[] = {
                CARQUET_ENCODING_PLAIN,
                CARQUET_ENCODING_BYTE_STREAM_SPLIT,
                CARQUET_ENCODING_RLE_DICTIONARY,
            };
            return choices[consume_byte(in) % (sizeof(choices) / sizeof(choices[0]))];
        }
        case CARQUET_PHYSICAL_BYTE_ARRAY: {
            static const carquet_encoding_t choices[] = {
                CARQUET_ENCODING_PLAIN,
                CARQUET_ENCODING_DELTA_LENGTH_BYTE_ARRAY,
                CARQUET_ENCODING_DELTA_BYTE_ARRAY,
                CARQUET_ENCODING_RLE_DICTIONARY,
            };
            return choices[consume_byte(in) % (sizeof(choices) / sizeof(choices[0]))];
        }
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY: {
            static const carquet_encoding_t choices[] = {
                CARQUET_ENCODING_PLAIN,
                CARQUET_ENCODING_DELTA_BYTE_ARRAY,
                CARQUET_ENCODING_BYTE_STREAM_SPLIT,
            };
            return choices[consume_byte(in) % (sizeof(choices) / sizeof(choices[0]))];
        }
        default:
            return CARQUET_ENCODING_PLAIN;
    }
}

static bool fuzz_row_group_filter(
    const carquet_reader_t* reader,
    int32_t row_group_index,
    void* user_data) {

    (void)reader;
    uint8_t selector = user_data ? *(const uint8_t*)user_data : 0;
    return ((uint8_t)row_group_index & 1u) == (selector & 1u);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    if (size < 10) return 0;
    (void)carquet_init();

    fuzz_input_t input = { data, size, 0 };
    carquet_error_t err = CARQUET_ERROR_INIT;

    uint8_t num_columns = (consume_byte(&input) % MAX_COLUMNS) + 1;
    uint16_t num_rows = (consume_u16(&input) % MAX_ROWS) + 1;
    uint8_t codec_idx = consume_byte(&input) % NUM_FUZZ_CODECS;
    uint8_t nullable_mask = consume_byte(&input);
    uint8_t config_byte = consume_byte(&input);
    uint8_t feature_byte = consume_byte(&input);

    carquet_compression_t codec = FUZZ_CODECS[codec_idx];

    /* Create schema */
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) return 0;

    fuzz_column_plan_t* columns = calloc(num_columns, sizeof(fuzz_column_plan_t));
    if (!columns) { carquet_schema_free(schema); return 0; }

    for (int i = 0; i < num_columns; i++) {
        char col_name[32];
        snprintf(col_name, sizeof(col_name), "col_%d", i);

        uint8_t type_idx = consume_byte(&input) % NUM_FUZZ_TYPES;
        columns[i].type = FUZZ_TYPES[type_idx];
        columns[i].type_length = choose_type_length(&input, columns[i].type);
        columns[i].nullable = (nullable_mask & (1 << (i % 8))) != 0;
        columns[i].encoding = choose_encoding(&input, columns[i].type);
        const carquet_logical_type_t* logical_type = choose_logical_type(
            &input, columns[i].type, columns[i].type_length, &columns[i].logical_type);
        columns[i].has_logical_type = logical_type != NULL;

        carquet_field_repetition_t rep = columns[i].nullable
            ? CARQUET_REPETITION_OPTIONAL : CARQUET_REPETITION_REQUIRED;

        if (carquet_schema_add_column(
                schema, col_name, columns[i].type, logical_type, rep,
                columns[i].type_length, 0) != CARQUET_OK) {
            free(columns);
            carquet_schema_free(schema);
            return 0;
        }
    }

    /* Temp file */
    char tmp_path[] = "/tmp/fuzz_writer_XXXXXX";
    int fd = mkstemp(tmp_path);
    if (fd < 0) { free(columns); carquet_schema_free(schema); return 0; }
    close(fd);

    /* Writer options — vary based on fuzz input */
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = codec;

    /* Vary page and row group sizes */
    uint32_t page_sizes[] = { 1024, 4096, 8192, 16384 };
    uint32_t rg_sizes[] = { 16384, 65536, 131072, 262144 };
    opts.page_size = page_sizes[config_byte & 0x03];
    opts.row_group_size = rg_sizes[(config_byte >> 2) & 0x03];
    opts.write_page_index = (config_byte & 0x10) != 0;
    opts.write_bloom_filters = (config_byte & 0x20) != 0;
    opts.write_arrow_schema = (config_byte & 0x40) != 0;
    opts.data_page_version = (config_byte & 0x80) ? 2 : 1;
    static const int64_t max_rows_choices[] = { 0, 1, 8, 64 };
    opts.max_rows_per_page =
        max_rows_choices[feature_byte & 0x03];
    static const int64_t dict_page_choices[] = { 64, 256, 4096, 1024 * 1024 };
    opts.dictionary_page_size =
        dict_page_choices[(feature_byte >> 2) & 0x03];

    bool use_buffer_writer = (feature_byte & 0x10) != 0;
    carquet_writer_t* writer = use_buffer_writer
        ? carquet_writer_create_buffer(schema, &opts, &err)
        : carquet_writer_create(tmp_path, schema, &opts, &err);
    if (!writer) {
        remove(tmp_path); free(columns); carquet_schema_free(schema);
        return 0;
    }

    if (feature_byte & 0x20) {
        int32_t sorting_count = num_columns < 3 ? num_columns : 3;
        carquet_sorting_column_t sorting[3];
        for (int32_t i = 0; i < sorting_count; i++) {
            sorting[i].column_index = i;
            sorting[i].descending = (consume_byte(&input) & 1) != 0;
            sorting[i].nulls_first = (consume_byte(&input) & 1) != 0;
        }
        (void)carquet_writer_set_sorting_columns(writer, sorting, sorting_count);
    }
    if (feature_byte & 0x40) {
        (void)carquet_writer_add_metadata(writer, "fuzz", "writer");
    }
    if (feature_byte & 0x80) {
        (void)carquet_writer_add_metadata(writer, "ARROW:schema", "already-present");
    }
    for (int i = 0; i < num_columns; i++) {
        (void)carquet_writer_set_column_encoding(writer, i, columns[i].encoding);
        if (consume_byte(&input) & 1) {
            carquet_compression_t col_codec =
                FUZZ_CODECS[consume_byte(&input) % NUM_FUZZ_CODECS];
            (void)carquet_writer_set_column_compression(writer, i, col_codec, 0);
        }
        (void)carquet_writer_set_column_statistics(
            writer, i, (consume_byte(&input) & 1) != 0);
        if (consume_byte(&input) & 1) {
            int64_t ndv = (int64_t)((consume_u16(&input) % 2048) + 1);
            double fpp = ((double)((consume_byte(&input) % 90) + 1)) / 100.0;
            (void)carquet_writer_set_column_bloom_filter_options(
                writer, i, true, ndv, fpp);
        } else {
            (void)carquet_writer_set_column_bloom_filter(
                writer, i, (consume_byte(&input) & 1) != 0);
        }
    }

    /* Allocate value buffers */
    void* values = malloc((size_t)num_rows * MAX_TYPE_LENGTH);
    int16_t* def_levels = malloc(num_rows * sizeof(int16_t));
    carquet_byte_array_t* ba_values = NULL;
    char* string_pool = NULL;

    if (!values || !def_levels) goto cleanup_write;

    /* Write each column */
    for (int col = 0; col < num_columns; col++) {
        bool is_nullable = columns[col].nullable;
        void* write_ptr = values;
        int64_t non_null_count = 0;

        switch (columns[col].type) {
            case CARQUET_PHYSICAL_BOOLEAN: {
                uint8_t* bools = (uint8_t*)values;
                for (int i = 0; i < num_rows; i++) {
                    bool present = !is_nullable || ((consume_byte(&input) & 1) != 0);
                    def_levels[i] = present ? 1 : 0;
                    if (present) {
                        bools[non_null_count++] = consume_byte(&input) & 1;
                    }
                }
                break;
            }
            case CARQUET_PHYSICAL_INT32: {
                int32_t* ints = (int32_t*)values;
                for (int i = 0; i < num_rows; i++) {
                    bool present = !is_nullable || ((consume_byte(&input) & 1) != 0);
                    def_levels[i] = present ? 1 : 0;
                    if (present) {
                        int32_t value = (int32_t)consume_u32(&input);
                        if (columns[col].encoding == CARQUET_ENCODING_RLE_DICTIONARY)
                            value %= 32;
                        ints[non_null_count++] = value;
                    }
                }
                break;
            }
            case CARQUET_PHYSICAL_INT64: {
                int64_t* longs = (int64_t*)values;
                for (int i = 0; i < num_rows; i++) {
                    bool present = !is_nullable || ((consume_byte(&input) & 1) != 0);
                    def_levels[i] = present ? 1 : 0;
                    if (present) {
                        int64_t value = (int64_t)consume_u64(&input);
                        if (columns[col].encoding == CARQUET_ENCODING_RLE_DICTIONARY)
                            value %= 32;
                        longs[non_null_count++] = value;
                    }
                }
                break;
            }
            case CARQUET_PHYSICAL_INT96: {
                carquet_int96_t* ints96 = (carquet_int96_t*)values;
                for (int i = 0; i < num_rows; i++) {
                    bool present = !is_nullable || ((consume_byte(&input) & 1) != 0);
                    def_levels[i] = present ? 1 : 0;
                    if (present) {
                        ints96[non_null_count].value[0] = consume_u32(&input);
                        ints96[non_null_count].value[1] = consume_u32(&input);
                        ints96[non_null_count].value[2] = consume_u32(&input);
                        non_null_count++;
                    }
                }
                break;
            }
            case CARQUET_PHYSICAL_FLOAT: {
                float* floats = (float*)values;
                for (int i = 0; i < num_rows; i++) {
                    bool present = !is_nullable || ((consume_byte(&input) & 1) != 0);
                    def_levels[i] = present ? 1 : 0;
                    if (present) {
                        float value = consume_float(&input);
                        if (columns[col].encoding == CARQUET_ENCODING_RLE_DICTIONARY)
                            value = (float)((int32_t)consume_byte(&input) % 16);
                        floats[non_null_count++] = value;
                    }
                }
                break;
            }
            case CARQUET_PHYSICAL_DOUBLE: {
                double* doubles = (double*)values;
                for (int i = 0; i < num_rows; i++) {
                    bool present = !is_nullable || ((consume_byte(&input) & 1) != 0);
                    def_levels[i] = present ? 1 : 0;
                    if (present) {
                        double value = consume_double(&input);
                        if (columns[col].encoding == CARQUET_ENCODING_RLE_DICTIONARY)
                            value = (double)((int32_t)consume_byte(&input) % 16);
                        doubles[non_null_count++] = value;
                    }
                }
                break;
            }
            case CARQUET_PHYSICAL_BYTE_ARRAY: {
                /* Real BYTE_ARRAY testing */
                if (!ba_values) {
                    ba_values = calloc(MAX_ROWS, sizeof(carquet_byte_array_t));
                    string_pool = malloc(MAX_ROWS * MAX_STRING_LEN);
                }
                if (!ba_values || !string_pool) goto cleanup_write;

                for (int i = 0; i < num_rows; i++) {
                    bool present = !is_nullable || ((consume_byte(&input) & 1) != 0);
                    def_levels[i] = present ? 1 : 0;
                    if (present) {
                        char* s = string_pool + ((size_t)non_null_count * MAX_STRING_LEN);
                        int slen = 0;
                        if (columns[col].encoding == CARQUET_ENCODING_DELTA_BYTE_ARRAY) {
                            slen = snprintf(s, MAX_STRING_LEN, "shared/prefix/%u/%lld",
                                (unsigned)(consume_byte(&input) % 8),
                                (long long)non_null_count);
                        } else if (columns[col].encoding == CARQUET_ENCODING_RLE_DICTIONARY) {
                            slen = snprintf(s, MAX_STRING_LEN, "dict_%u",
                                (unsigned)(consume_byte(&input) % 16));
                        } else {
                            slen = consume_byte(&input) % MAX_STRING_LEN;
                            for (int j = 0; j < slen; j++)
                                s[j] = (char)(consume_byte(&input) % 95 + 32);
                        }
                        ba_values[non_null_count].data = (uint8_t*)s;
                        ba_values[non_null_count].length = slen;
                        non_null_count++;
                    }
                }
                write_ptr = ba_values;
                break;
            }
            case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY: {
                uint8_t* fixed = (uint8_t*)values;
                int32_t len = columns[col].type_length;
                for (int i = 0; i < num_rows; i++) {
                    bool present = !is_nullable || ((consume_byte(&input) & 1) != 0);
                    def_levels[i] = present ? 1 : 0;
                    if (!present) continue;
                    uint8_t* dst = fixed + ((size_t)non_null_count * (size_t)len);
                    for (int32_t j = 0; j < len; j++) {
                        dst[j] = (uint8_t)(
                            columns[col].encoding == CARQUET_ENCODING_DELTA_BYTE_ARRAY && j < len / 2
                                ? j
                                : consume_byte(&input));
                    }
                    non_null_count++;
                }
                break;
            }
            default:
                break;
        }

        carquet_status_t status = carquet_writer_write_batch(
            writer, col, write_ptr, num_rows,
            is_nullable ? def_levels : NULL, NULL);

        if (status != CARQUET_OK) goto cleanup_write;
    }

    {
        carquet_status_t close_status = carquet_writer_close(writer);
        if (!use_buffer_writer) {
            writer = NULL; /* file writers are consumed by close */
        }

        free(values); values = NULL;
        free(def_levels); def_levels = NULL;
        free(ba_values); ba_values = NULL;
        free(string_pool); string_pool = NULL;

        if (close_status != CARQUET_OK) {
            if (use_buffer_writer && writer) {
                carquet_writer_abort(writer);
                writer = NULL;
            }
            free(columns);
            carquet_schema_free(schema);
            remove(tmp_path);
            return 0;
        }
    }

    /* Roundtrip: read back and verify metadata */
    {
        void* buffer_data = NULL;
        size_t buffer_size = 0;
        carquet_reader_t* reader = NULL;
        if (use_buffer_writer) {
            if (carquet_writer_get_buffer(writer, &buffer_data, &buffer_size) == CARQUET_OK) {
                writer = NULL; /* get_buffer consumes the closed buffer writer */
                reader = carquet_reader_open_buffer(buffer_data, buffer_size, NULL, &err);
            }
        } else {
            reader = carquet_reader_open(tmp_path, NULL, &err);
        }
        if (reader) {
            int64_t read_rows = carquet_reader_num_rows(reader);
            int32_t read_cols = carquet_reader_num_columns(reader);

            if (read_rows != num_rows || read_cols != num_columns) {
                carquet_reader_close(reader);
                free(buffer_data);
                free(columns);
                carquet_schema_free(schema);
                remove(tmp_path);
                __builtin_trap();
            }

            if (!use_buffer_writer && read_cols > 0) {
                int32_t cols_to_prebuffer[4];
                int32_t npre = read_cols < 4 ? read_cols : 4;
                for (int32_t i = 0; i < npre; i++) cols_to_prebuffer[i] = i;
                (void)carquet_reader_prebuffer(reader, 0, cols_to_prebuffer, npre, &err);
            }

            /* Read all data via batch reader */
            carquet_batch_reader_config_t config;
            carquet_batch_reader_config_init(&config);
            config.batch_size = 100;
            config.preserve_dictionaries = true;

            carquet_batch_reader_t* batch_reader = carquet_batch_reader_create(reader, &config, &err);
            if (batch_reader) {
                carquet_row_batch_t* batch = NULL;
                int64_t total_read = 0;
                while (carquet_batch_reader_next(batch_reader, &batch) == CARQUET_OK && batch) {
                    int64_t n = carquet_row_batch_num_rows(batch);
                    /* Access column data to exercise decoders */
                    for (int32_t c = 0; c < read_cols; c++) {
                        const void* col_data;
                        const uint8_t* nulls;
                        int64_t count;
                        (void)carquet_row_batch_column(batch, c, &col_data, &nulls, &count);
                        const uint32_t* indices;
                        const uint8_t* dict_data;
                        int32_t dict_count;
                        const uint32_t* dict_offsets;
                        (void)carquet_row_batch_column_dictionary(
                            batch, c, &indices, &nulls, &count,
                            &dict_data, &dict_count, &dict_offsets);
                    }
                    total_read += n;
                    carquet_row_batch_free(batch);
                    batch = NULL;
                }
                carquet_batch_reader_free(batch_reader);

                if (total_read != num_rows) {
                    carquet_reader_close(reader);
                    free(buffer_data);
                    free(columns);
                    carquet_schema_free(schema);
                    remove(tmp_path);
                    __builtin_trap();
                }
            }

            /* Exercise row-group filter handling even when it skips the only RG. */
            carquet_batch_reader_config_init(&config);
            config.batch_size = 64;
            config.row_group_filter = fuzz_row_group_filter;
            config.row_group_filter_ctx = &feature_byte;
            batch_reader = carquet_batch_reader_create(reader, &config, &err);
            if (batch_reader) {
                carquet_row_batch_t* batch = NULL;
                int batches = 0;
                while (batches < 4 &&
                       carquet_batch_reader_next(batch_reader, &batch) == CARQUET_OK && batch) {
                    carquet_row_batch_free(batch);
                    batch = NULL;
                    batches++;
                }
                carquet_batch_reader_free(batch_reader);
            }

            /* Also exercise low-level column reader */
            int32_t num_rg = carquet_reader_num_row_groups(reader);
            for (int32_t rg = 0; rg < num_rg && rg < 5; rg++) {
                carquet_row_group_metadata_t rg_meta;
                if (carquet_reader_row_group_metadata(reader, rg, &rg_meta) != CARQUET_OK)
                    continue;
                for (int32_t c = 0; c < read_cols && c < 8; c++) {
                    carquet_column_chunk_metadata_t meta;
                    (void)carquet_reader_column_chunk_metadata(reader, rg, c, &meta);
                    carquet_column_reader_t* col_reader =
                        carquet_reader_get_column(reader, rg, c, &err);
                    if (col_reader) {
                        void* col_vals = malloc((size_t)rg_meta.num_rows * MAX_TYPE_LENGTH);
                        int16_t* dl = malloc((size_t)rg_meta.num_rows * sizeof(int16_t));
                        if (col_vals && dl) {
                            (void)carquet_column_read_batch(col_reader, col_vals,
                                                           rg_meta.num_rows, dl, NULL);
                        }
                        free(col_vals);
                        free(dl);
                        carquet_column_reader_free(col_reader);
                    }
                }
            }

            carquet_reader_release_prebuffer(reader);
            carquet_reader_close(reader);
        }
        free(buffer_data);
    }

    free(columns);
    carquet_schema_free(schema);
    remove(tmp_path);
    return 0;

cleanup_write:
    free(values);
    free(def_levels);
    free(ba_values);
    free(string_pool);
    if (writer) carquet_writer_abort(writer);
    free(columns);
    carquet_schema_free(schema);
    remove(tmp_path);
    return 0;
}

#ifdef AFL_MAIN
#include <sys/stat.h>
int main(int argc, char** argv) {
    if (argc != 2) { fprintf(stderr, "Usage: %s <input_file>\n", argv[0]); return 1; }
    FILE* f = fopen(argv[1], "rb");
    if (!f) { perror("fopen"); return 1; }
    struct stat st; fstat(fileno(f), &st);
    uint8_t* d = malloc((size_t)st.st_size);
    if (!d) { fclose(f); return 1; }
    fread(d, 1, (size_t)st.st_size, f); fclose(f);
    int r = LLVMFuzzerTestOneInput(d, (size_t)st.st_size);
    free(d); return r;
}
#endif
