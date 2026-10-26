/**
 * @file error.c
 * @brief Error handling implementation
 */

#include <carquet/error.h>
#include <carquet/types.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

/* ============================================================================
 * Error Functions
 * ============================================================================
 */

void carquet_error_init(carquet_error_t* error) {
    if (!error) return;

    error->code = CARQUET_OK;
    error->message[0] = '\0';
    error->file = NULL;
    error->line = 0;
    error->function = NULL;
    error->offset = -1;
    error->column_index = -1;
    error->row_group_index = -1;
}

void carquet_error_clear(carquet_error_t* error) {
    carquet_error_init(error);
}

void carquet_error_set(carquet_error_t* error,
                       carquet_status_t code,
                       const char* file,
                       int line,
                       const char* function,
                       const char* format, ...) {
    if (!error) return;

    error->code = code;
    error->file = file;
    error->line = line;
    error->function = function;

    if (format) {
        va_list args;
        va_start(args, format);
        vsnprintf(error->message, CARQUET_ERROR_MESSAGE_MAX, format, args);
        va_end(args);
    } else {
        error->message[0] = '\0';
    }
}

void carquet_error_copy(carquet_error_t* dest, const carquet_error_t* src) {
    if (!dest || !src) return;
    *dest = *src;
}

const char* carquet_status_string(carquet_status_t status) {
    switch (status) {
        case CARQUET_OK:
            return "Success";
        case CARQUET_ERROR_INVALID_ARGUMENT:
            return "Invalid argument";
        case CARQUET_ERROR_OUT_OF_MEMORY:
            return "Out of memory";
        case CARQUET_ERROR_NOT_IMPLEMENTED:
            return "Not implemented";
        case CARQUET_ERROR_INTERNAL:
            return "Internal error";
        case CARQUET_ERROR_FILE_NOT_FOUND:
            return "File not found";
        case CARQUET_ERROR_FILE_OPEN:
            return "Failed to open file";
        case CARQUET_ERROR_FILE_READ:
            return "Failed to read file";
        case CARQUET_ERROR_FILE_WRITE:
            return "Failed to write file";
        case CARQUET_ERROR_FILE_SEEK:
            return "Failed to seek in file";
        case CARQUET_ERROR_FILE_TRUNCATED:
            return "File truncated or incomplete";
        case CARQUET_ERROR_INVALID_MAGIC:
            return "Invalid magic bytes";
        case CARQUET_ERROR_INVALID_FOOTER:
            return "Invalid file footer";
        case CARQUET_ERROR_INVALID_SCHEMA:
            return "Invalid schema";
        case CARQUET_ERROR_INVALID_METADATA:
            return "Invalid metadata";
        case CARQUET_ERROR_INVALID_PAGE:
            return "Invalid page";
        case CARQUET_ERROR_INVALID_ENCODING:
            return "Invalid or unsupported encoding";
        case CARQUET_ERROR_VERSION_NOT_SUPPORTED:
            return "Version not supported";
        case CARQUET_ERROR_THRIFT_DECODE:
            return "Thrift decode error";
        case CARQUET_ERROR_THRIFT_ENCODE:
            return "Thrift encode error";
        case CARQUET_ERROR_THRIFT_INVALID_TYPE:
            return "Invalid Thrift type";
        case CARQUET_ERROR_THRIFT_TRUNCATED:
            return "Truncated Thrift data";
        case CARQUET_ERROR_DECODE:
            return "Decode error";
        case CARQUET_ERROR_ENCODE:
            return "Encode error";
        case CARQUET_ERROR_DICTIONARY_NOT_FOUND:
            return "Dictionary not found";
        case CARQUET_ERROR_INVALID_RLE:
            return "Invalid RLE data";
        case CARQUET_ERROR_INVALID_DELTA:
            return "Invalid delta encoding data";
        case CARQUET_ERROR_COMPRESSION:
            return "Compression error";
        case CARQUET_ERROR_DECOMPRESSION:
            return "Decompression error";
        case CARQUET_ERROR_UNSUPPORTED_CODEC:
            return "Unsupported compression codec";
        case CARQUET_ERROR_INVALID_COMPRESSED_DATA:
            return "Invalid compressed data";
        case CARQUET_ERROR_TYPE_MISMATCH:
            return "Type mismatch";
        case CARQUET_ERROR_COLUMN_NOT_FOUND:
            return "Column not found";
        case CARQUET_ERROR_ROW_GROUP_NOT_FOUND:
            return "Row group not found";
        case CARQUET_ERROR_END_OF_DATA:
            return "End of data";
        case CARQUET_ERROR_CHECKSUM:
            return "Checksum error";
        case CARQUET_ERROR_CRC_MISMATCH:
            return "CRC mismatch";
        case CARQUET_ERROR_INVALID_STATE:
            return "Invalid state";
        case CARQUET_ERROR_ALREADY_CLOSED:
            return "Already closed";
        case CARQUET_ERROR_NOT_OPEN:
            return "Not open";
        case CARQUET_ERROR_PAGE_INDEX_REQUIRED:
            return "Page index required but absent for filtered column";
        default:
            return "Unknown error";
    }
}

/* ============================================================================
 * Type Name Functions
 * ============================================================================
 */

const char* carquet_physical_type_name(carquet_physical_type_t type) {
    switch (type) {
        case CARQUET_PHYSICAL_BOOLEAN:
            return "BOOLEAN";
        case CARQUET_PHYSICAL_INT32:
            return "INT32";
        case CARQUET_PHYSICAL_INT64:
            return "INT64";
        case CARQUET_PHYSICAL_INT96:
            return "INT96";
        case CARQUET_PHYSICAL_FLOAT:
            return "FLOAT";
        case CARQUET_PHYSICAL_DOUBLE:
            return "DOUBLE";
        case CARQUET_PHYSICAL_BYTE_ARRAY:
            return "BYTE_ARRAY";
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
            return "FIXED_LEN_BYTE_ARRAY";
        default:
            return "UNKNOWN";
    }
}

const char* carquet_compression_name(carquet_compression_t codec) {
    switch (codec) {
        case CARQUET_COMPRESSION_UNCOMPRESSED:
            return "UNCOMPRESSED";
        case CARQUET_COMPRESSION_SNAPPY:
            return "SNAPPY";
        case CARQUET_COMPRESSION_GZIP:
            return "GZIP";
        case CARQUET_COMPRESSION_LZO:
            return "LZO";
        case CARQUET_COMPRESSION_BROTLI:
            return "BROTLI";
        case CARQUET_COMPRESSION_LZ4:
            return "LZ4";
        case CARQUET_COMPRESSION_ZSTD:
            return "ZSTD";
        case CARQUET_COMPRESSION_LZ4_RAW:
            return "LZ4_RAW";
        default:
            return "UNKNOWN";
    }
}

const char* carquet_encoding_name(carquet_encoding_t encoding) {
    switch (encoding) {
        case CARQUET_ENCODING_PLAIN:
            return "PLAIN";
        case CARQUET_ENCODING_PLAIN_DICTIONARY:
            return "PLAIN_DICTIONARY";
        case CARQUET_ENCODING_RLE:
            return "RLE";
        case CARQUET_ENCODING_BIT_PACKED:
            return "BIT_PACKED";
        case CARQUET_ENCODING_DELTA_BINARY_PACKED:
            return "DELTA_BINARY_PACKED";
        case CARQUET_ENCODING_DELTA_LENGTH_BYTE_ARRAY:
            return "DELTA_LENGTH_BYTE_ARRAY";
        case CARQUET_ENCODING_DELTA_BYTE_ARRAY:
            return "DELTA_BYTE_ARRAY";
        case CARQUET_ENCODING_RLE_DICTIONARY:
            return "RLE_DICTIONARY";
        case CARQUET_ENCODING_BYTE_STREAM_SPLIT:
            return "BYTE_STREAM_SPLIT";
        default:
            return "UNKNOWN";
    }
}

/* ============================================================================
 * Enhanced Error Reporting
 * ============================================================================
 */

const char* carquet_error_recovery_hint(carquet_status_t status) {
    switch (status) {
        case CARQUET_OK:
            return NULL;

        case CARQUET_ERROR_INVALID_MAGIC:
            return "Ensure the file is a valid Parquet file (should start with 'PAR1')";

        case CARQUET_ERROR_INVALID_FOOTER:
            return "The file may be corrupted or incomplete. Try re-downloading or regenerating it";

        case CARQUET_ERROR_FILE_TRUNCATED:
            return "The file appears incomplete. Check if the write operation completed successfully";

        case CARQUET_ERROR_CRC_MISMATCH:
            return "Data integrity check failed. The file may be corrupted during transfer or storage";

        case CARQUET_ERROR_UNSUPPORTED_CODEC:
            return "This compression codec is not supported. Supported: UNCOMPRESSED, SNAPPY, GZIP, LZ4, ZSTD";

        case CARQUET_ERROR_INVALID_ENCODING:
            return "Encoding not supported. Supported: PLAIN, RLE, DICTIONARY, DELTA_*, BYTE_STREAM_SPLIT";

        case CARQUET_ERROR_OUT_OF_MEMORY:
            return "Not enough memory. Try processing data in smaller batches or free system memory";

        case CARQUET_ERROR_DICTIONARY_NOT_FOUND:
            return "Dictionary page missing for dictionary-encoded column. File may be malformed";

        case CARQUET_ERROR_VERSION_NOT_SUPPORTED:
            return "Parquet file uses unsupported features. Try with a different Parquet writer";

        case CARQUET_ERROR_COLUMN_NOT_FOUND:
            return "Verify column name or index is correct for this file's schema";

        case CARQUET_ERROR_ROW_GROUP_NOT_FOUND:
            return "Row group index is out of range. Check carquet_reader_num_row_groups()";

        case CARQUET_ERROR_DECOMPRESSION:
            return "Failed to decompress data. The file may be corrupted or use an unsupported variant";

        case CARQUET_ERROR_TYPE_MISMATCH:
            return "Requested type doesn't match column physical type. Check schema before reading";

        default:
            return NULL;
    }
}

int carquet_error_format(const carquet_error_t* error, char* buffer, size_t buffer_size) {
    if (!error || !buffer || buffer_size == 0) return 0;

    int written = 0;

    /* Basic error info */
    written = snprintf(buffer, buffer_size, "[%s] %s",
                       carquet_status_string(error->code),
                       error->message[0] ? error->message : "(no details)");

    if (written < 0 || (size_t)written >= buffer_size) {
        return written < 0 ? -1 : (int)buffer_size - 1;
    }

    /* Add location context if available */
    if (error->offset >= 0) {
        int len = snprintf(buffer + written, buffer_size - written,
                           " (file offset: %lld)", (long long)error->offset);
        if (len > 0 && (size_t)(written + len) < buffer_size) {
            written += len;
        }
    }

    if (error->row_group_index >= 0) {
        int len = snprintf(buffer + written, buffer_size - written,
                           " (row group: %d)", error->row_group_index);
        if (len > 0 && (size_t)(written + len) < buffer_size) {
            written += len;
        }
    }

    if (error->column_index >= 0) {
        int len = snprintf(buffer + written, buffer_size - written,
                           " (column: %d)", error->column_index);
        if (len > 0 && (size_t)(written + len) < buffer_size) {
            written += len;
        }
    }

    /* Add recovery hint */
    const char* hint = carquet_error_recovery_hint(error->code);
    if (hint) {
        int len = snprintf(buffer + written, buffer_size - written,
                           "\n  Hint: %s", hint);
        if (len > 0 && (size_t)(written + len) < buffer_size) {
            written += len;
        }
    }

    return written;
}

void carquet_error_set_context(carquet_error_t* error,
                                int64_t offset,
                                int32_t row_group_index,
                                int32_t column_index) {
    if (!error) return;

    if (offset >= 0) error->offset = offset;
    if (row_group_index >= 0) error->row_group_index = row_group_index;
    if (column_index >= 0) error->column_index = column_index;
}

bool carquet_error_is_recoverable(carquet_status_t status) {
    switch (status) {
        /* These are generally not recoverable without user intervention */
        case CARQUET_ERROR_INVALID_MAGIC:
        case CARQUET_ERROR_INVALID_FOOTER:
        case CARQUET_ERROR_FILE_TRUNCATED:
        case CARQUET_ERROR_CRC_MISMATCH:
        case CARQUET_ERROR_VERSION_NOT_SUPPORTED:
            return false;

        /* These might be recoverable by skipping or retrying */
        case CARQUET_ERROR_DECOMPRESSION:
        case CARQUET_ERROR_DECODE:
        case CARQUET_ERROR_INVALID_PAGE:
            return true;

        /* Resource errors - might resolve with retry */
        case CARQUET_ERROR_OUT_OF_MEMORY:
        case CARQUET_ERROR_FILE_READ:
        case CARQUET_ERROR_FILE_SEEK:
            return true;

        /* Generally not recoverable */
        default:
            return false;
    }
}
