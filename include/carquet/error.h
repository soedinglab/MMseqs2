/**
 * @file error.h
 * @brief Error handling for Carquet library
 *
 * This header provides error codes and error handling utilities.
 * All Carquet functions that can fail return an error code or use
 * the carquet_error_t structure for detailed error information.
 */

#ifndef CARQUET_ERROR_H
#define CARQUET_ERROR_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Error Codes
 * ============================================================================
 */

typedef enum carquet_status {
    /* Success */
    CARQUET_OK = 0,

    /* General errors */
    CARQUET_ERROR_INVALID_ARGUMENT = 1,
    CARQUET_ERROR_OUT_OF_MEMORY = 2,
    CARQUET_ERROR_NOT_IMPLEMENTED = 3,
    CARQUET_ERROR_INTERNAL = 4,

    /* File I/O errors */
    CARQUET_ERROR_FILE_NOT_FOUND = 10,
    CARQUET_ERROR_FILE_OPEN = 11,
    CARQUET_ERROR_FILE_READ = 12,
    CARQUET_ERROR_FILE_WRITE = 13,
    CARQUET_ERROR_FILE_SEEK = 14,
    CARQUET_ERROR_FILE_TRUNCATED = 15,

    /* Format errors */
    CARQUET_ERROR_INVALID_MAGIC = 20,
    CARQUET_ERROR_INVALID_FOOTER = 21,
    CARQUET_ERROR_INVALID_SCHEMA = 22,
    CARQUET_ERROR_INVALID_METADATA = 23,
    CARQUET_ERROR_INVALID_PAGE = 24,
    CARQUET_ERROR_INVALID_ENCODING = 25,
    CARQUET_ERROR_VERSION_NOT_SUPPORTED = 26,

    /* Thrift parsing errors */
    CARQUET_ERROR_THRIFT_DECODE = 30,
    CARQUET_ERROR_THRIFT_ENCODE = 31,
    CARQUET_ERROR_THRIFT_INVALID_TYPE = 32,
    CARQUET_ERROR_THRIFT_TRUNCATED = 33,

    /* Encoding/decoding errors */
    CARQUET_ERROR_DECODE = 40,
    CARQUET_ERROR_ENCODE = 41,
    CARQUET_ERROR_DICTIONARY_NOT_FOUND = 42,
    CARQUET_ERROR_INVALID_RLE = 43,
    CARQUET_ERROR_INVALID_DELTA = 44,

    /* Compression errors */
    CARQUET_ERROR_COMPRESSION = 50,
    CARQUET_ERROR_DECOMPRESSION = 51,
    CARQUET_ERROR_UNSUPPORTED_CODEC = 52,
    CARQUET_ERROR_INVALID_COMPRESSED_DATA = 53,

    /* Data errors */
    CARQUET_ERROR_TYPE_MISMATCH = 60,
    CARQUET_ERROR_COLUMN_NOT_FOUND = 61,
    CARQUET_ERROR_ROW_GROUP_NOT_FOUND = 62,
    CARQUET_ERROR_END_OF_DATA = 63,

    /* Checksum errors */
    CARQUET_ERROR_CHECKSUM = 70,
    CARQUET_ERROR_CRC_MISMATCH = 71,

    /* State errors */
    CARQUET_ERROR_INVALID_STATE = 80,
    CARQUET_ERROR_ALREADY_CLOSED = 81,
    CARQUET_ERROR_NOT_OPEN = 82,

    /* Filter / page-index errors */
    CARQUET_ERROR_PAGE_INDEX_REQUIRED = 90,

} carquet_status_t;

/* ============================================================================
 * Error Context
 * ============================================================================
 * Detailed error information for debugging.
 */

#define CARQUET_ERROR_MESSAGE_MAX 256

typedef struct carquet_error {
    carquet_status_t code;
    char message[CARQUET_ERROR_MESSAGE_MAX];

    /* Location information (optional) */
    const char* file;
    int line;
    const char* function;

    /* Additional context */
    int64_t offset;           /* File offset where error occurred */
    int32_t column_index;     /* Column index if applicable */
    int32_t row_group_index;  /* Row group index if applicable */
} carquet_error_t;

/* ============================================================================
 * Error Handling Macros
 * ============================================================================
 */

/**
 * Initialize an error structure to success state. Matches carquet_error_init():
 * the context fields use -1 for "not applicable", since 0 is a valid file
 * offset, column and row group.
 */
#define CARQUET_ERROR_INIT { .code = CARQUET_OK, .message = {0}, .offset = -1, \
                             .column_index = -1, .row_group_index = -1 }

/**
 * Check if status indicates success.
 */
#define CARQUET_SUCCEEDED(status) ((status) == CARQUET_OK)

/**
 * Check if status indicates failure.
 */
#define CARQUET_FAILED(status) ((status) != CARQUET_OK)

/**
 * Return early if status is not OK.
 */
#define CARQUET_RETURN_IF_ERROR(status) \
    do { \
        carquet_status_t _status = (status); \
        if (CARQUET_FAILED(_status)) return _status; \
    } while (0)

/**
 * Set error with location information.
 * Format string is included in variadic args to avoid C23 extension warnings.
 */
#define CARQUET_SET_ERROR(err, status_code, ...) \
    carquet_error_set((err), (status_code), __FILE__, __LINE__, __func__, __VA_ARGS__)

/**
 * Set error if condition is false, return status.
 */
#define CARQUET_CHECK(cond, err, status_code, ...) \
    do { \
        if (!(cond)) { \
            CARQUET_SET_ERROR((err), (status_code), __VA_ARGS__); \
            return (status_code); \
        } \
    } while (0)

/* ============================================================================
 * Error Functions
 * ============================================================================
 */

/**
 * Initialize an error structure.
 */
void carquet_error_init(carquet_error_t* error);

/**
 * Clear an error structure (reset to success state).
 */
void carquet_error_clear(carquet_error_t* error);

/**
 * Set error information.
 */
#if defined(__GNUC__) || defined(__clang__)
__attribute__((format(printf, 6, 7)))
#endif
void carquet_error_set(carquet_error_t* error,
                       carquet_status_t code,
                       const char* file,
                       int line,
                       const char* function,
                       const char* format, ...);

/**
 * Copy error from source to destination.
 */
void carquet_error_copy(carquet_error_t* dest, const carquet_error_t* src);

/**
 * Get a human-readable description of a status code.
 */
const char* carquet_status_string(carquet_status_t status);

/**
 * Check if error is set (not OK).
 */
static inline bool carquet_error_is_set(const carquet_error_t* error) {
    return error && error->code != CARQUET_OK;
}

/**
 * Get error code from error structure.
 */
static inline carquet_status_t carquet_error_code(const carquet_error_t* error) {
    return error ? error->code : CARQUET_OK;
}

/**
 * Get error message from error structure.
 */
static inline const char* carquet_error_message(const carquet_error_t* error) {
    return error ? error->message : "";
}

/**
 * Get a recovery hint for a status code.
 * Returns NULL if no hint is available.
 */
const char* carquet_error_recovery_hint(carquet_status_t status);

/**
 * Format an error into a human-readable string.
 *
 * The output includes:
 * - Status code name and message
 * - File offset, row group, and column context (if set)
 * - Recovery hint (if available)
 *
 * @param error The error to format
 * @param buffer Output buffer
 * @param buffer_size Size of output buffer
 * @return Number of characters written (excluding null terminator)
 */
int carquet_error_format(const carquet_error_t* error, char* buffer, size_t buffer_size);

/**
 * Set additional context on an error.
 *
 * @param error The error to modify
 * @param offset File offset where error occurred (-1 to skip)
 * @param row_group_index Row group index (-1 to skip)
 * @param column_index Column index (-1 to skip)
 */
void carquet_error_set_context(carquet_error_t* error,
                                int64_t offset,
                                int32_t row_group_index,
                                int32_t column_index);

/**
 * Check if an error might be recoverable.
 *
 * Some errors (like file corruption) are not recoverable, while
 * others (like temporary I/O errors) might succeed on retry.
 *
 * @param status The status code to check
 * @return true if the error might be recoverable
 */
bool carquet_error_is_recoverable(carquet_status_t status);

/* ============================================================================
 * Result Type Pattern
 * ============================================================================
 * For functions that return a value or an error.
 */

#define CARQUET_RESULT(type) \
    struct { \
        carquet_status_t status; \
        type value; \
    }

/* Common result types */
typedef struct carquet_result_i32 {
    carquet_status_t status;
    int32_t value;
} carquet_result_i32_t;

typedef struct carquet_result_i64 {
    carquet_status_t status;
    int64_t value;
} carquet_result_i64_t;

typedef struct carquet_result_size {
    carquet_status_t status;
    size_t value;
} carquet_result_size_t;

typedef struct carquet_result_ptr {
    carquet_status_t status;
    void* value;
} carquet_result_ptr_t;

#ifdef __cplusplus
}
#endif

#endif /* CARQUET_ERROR_H */
