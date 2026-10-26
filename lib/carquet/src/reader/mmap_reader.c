/**
 * @file mmap_reader.c
 * @brief Memory-mapped I/O support for zero-copy reads
 *
 * Provides memory-mapped file access for improved performance when reading
 * large Parquet files. Memory mapping allows the OS to handle paging and
 * caching efficiently.
 */

#include "core/allocator.h"
#include <carquet/carquet.h>
#include "reader_internal.h"
#include "../core/endian.h"
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#endif

/* ============================================================================
 * Platform-specific Implementation
 * ============================================================================
 */

#ifdef _WIN32

carquet_mmap_info_t* carquet_mmap_open(const char* path, carquet_error_t* error) {
    carquet_mmap_info_t* mmap_info = carquet_mem_calloc(1, sizeof(carquet_mmap_info_t));
    if (!mmap_info) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate mmap info");
        return NULL;
    }

    /* Open file */
    mmap_info->file_handle = CreateFileA(
        path,
        GENERIC_READ,
        FILE_SHARE_READ,
        NULL,
        OPEN_EXISTING,
        FILE_ATTRIBUTE_NORMAL,
        NULL);

    if (mmap_info->file_handle == INVALID_HANDLE_VALUE) {
        carquet_mem_free(mmap_info);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_OPEN, "Failed to open file for mmap");
        return NULL;
    }

    /* Get file size */
    LARGE_INTEGER file_size;
    if (!GetFileSizeEx(mmap_info->file_handle, &file_size)) {
        CloseHandle(mmap_info->file_handle);
        carquet_mem_free(mmap_info);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_READ, "Failed to get file size");
        return NULL;
    }
    mmap_info->size = (size_t)file_size.QuadPart;

    /* Create file mapping */
    mmap_info->mapping_handle = CreateFileMappingA(
        mmap_info->file_handle,
        NULL,
        PAGE_READONLY,
        0, 0,
        NULL);

    if (!mmap_info->mapping_handle) {
        CloseHandle(mmap_info->file_handle);
        carquet_mem_free(mmap_info);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_READ, "Failed to create file mapping");
        return NULL;
    }

    /* Map view */
    mmap_info->data = (uint8_t*)MapViewOfFile(
        mmap_info->mapping_handle,
        FILE_MAP_READ,
        0, 0, 0);

    if (!mmap_info->data) {
        CloseHandle(mmap_info->mapping_handle);
        CloseHandle(mmap_info->file_handle);
        carquet_mem_free(mmap_info);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_READ, "Failed to map file view");
        return NULL;
    }

    mmap_info->is_valid = true;
    return mmap_info;
}

void carquet_mmap_close(carquet_mmap_info_t* mmap_info) {
    if (!mmap_info) return;

    if (mmap_info->data) {
        UnmapViewOfFile(mmap_info->data);
    }
    if (mmap_info->mapping_handle) {
        CloseHandle(mmap_info->mapping_handle);
    }
    if (mmap_info->file_handle != INVALID_HANDLE_VALUE) {
        CloseHandle(mmap_info->file_handle);
    }
    mmap_info->is_valid = false;
    carquet_mem_free(mmap_info);
}

#else /* POSIX */

carquet_mmap_info_t* carquet_mmap_open(const char* path, carquet_error_t* error) {
    carquet_mmap_info_t* mmap_info = carquet_mem_calloc(1, sizeof(carquet_mmap_info_t));
    if (!mmap_info) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate mmap info");
        return NULL;
    }

    /* Open file */
    mmap_info->fd = open(path, O_RDONLY);
    if (mmap_info->fd < 0) {
        carquet_mem_free(mmap_info);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_OPEN, "Failed to open file for mmap: %s", path);
        return NULL;
    }

    /* Get file size */
    struct stat st;
    if (fstat(mmap_info->fd, &st) < 0) {
        close(mmap_info->fd);
        carquet_mem_free(mmap_info);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_READ, "Failed to stat file");
        return NULL;
    }
    mmap_info->size = (size_t)st.st_size;

    /* Memory map the file */
    /* Read-only, so MAP_SHARED and MAP_PRIVATE see the same bytes; a private
     * mapping additionally sets up copy-on-write shadow state that costs
     * noticeable time at mmap/munmap for a large file. */
    mmap_info->data = mmap(NULL, mmap_info->size, PROT_READ, MAP_SHARED, mmap_info->fd, 0);
    if (mmap_info->data == MAP_FAILED) {
        close(mmap_info->fd);
        carquet_mem_free(mmap_info);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_READ, "Failed to mmap file");
        return NULL;
    }

    /* Analytics scans walk column chunks in file-offset order and benefit from
     * sequential readahead more than random-page heuristics. */
    madvise(mmap_info->data, mmap_info->size, MADV_SEQUENTIAL);

    mmap_info->is_valid = true;
    return mmap_info;
}

void carquet_mmap_close(carquet_mmap_info_t* mmap_info) {
    if (!mmap_info) return;

    if (mmap_info->data && mmap_info->data != MAP_FAILED) {
        munmap(mmap_info->data, mmap_info->size);
    }
    if (mmap_info->fd >= 0) {
        close(mmap_info->fd);
    }
    mmap_info->is_valid = false;
    carquet_mem_free(mmap_info);
}

#endif

/* ============================================================================
 * Public API for Memory-Mapped Reading
 * ============================================================================
 */

/**
 * Internal function to open a file with memory mapping.
 * This is called from file_reader.c when use_mmap is true.
 */
carquet_status_t carquet_reader_open_mmap_internal(
    carquet_reader_t* reader,
    const char* path,
    carquet_error_t* error) {

    carquet_mmap_info_t* mmap_info = carquet_mmap_open(path, error);
    if (!mmap_info) {
        return error ? error->code : CARQUET_ERROR_FILE_OPEN;
    }

    reader->mmap_data = mmap_info->data;
    reader->file_size = mmap_info->size;
    reader->mmap_info = mmap_info;  /* Store for cleanup in close() */

    return CARQUET_OK;
}

/* ============================================================================
 * Zero-Copy Eligibility Check
 * ============================================================================
 */

/**
 * Check if a page is eligible for zero-copy reading.
 * Zero-copy requires:
 * - Little-endian system (Parquet stores values in little-endian)
 * - Uncompressed data (no decompression needed)
 * - PLAIN encoding (no decoding needed)
 * - Fixed-size type (predictable layout)
 */
bool carquet_page_is_zero_copy_eligible(
    carquet_compression_t codec,
    carquet_encoding_t encoding,
    carquet_physical_type_t type) {

#if !CARQUET_LITTLE_ENDIAN
    /* Big-endian systems cannot use zero-copy for numeric types
     * because Parquet stores values in little-endian format */
    (void)codec;
    (void)encoding;
    (void)type;
    return false;
#else
    /* Must be uncompressed */
    if (codec != CARQUET_COMPRESSION_UNCOMPRESSED) {
        return false;
    }

    /* Must be PLAIN encoding */
    if (encoding != CARQUET_ENCODING_PLAIN) {
        return false;
    }

    /* Must be fixed-size type */
    switch (type) {
        case CARQUET_PHYSICAL_INT32:
        case CARQUET_PHYSICAL_INT64:
        case CARQUET_PHYSICAL_FLOAT:
        case CARQUET_PHYSICAL_DOUBLE:
        case CARQUET_PHYSICAL_INT96:
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
            return true;

        case CARQUET_PHYSICAL_BOOLEAN:
            /* Boolean is bit-packed, not directly mappable */
            return false;

        case CARQUET_PHYSICAL_BYTE_ARRAY:
            /* Variable length, requires length parsing */
            return false;

        default:
            return false;
    }
#endif
}

/**
 * Open a Parquet file from a memory buffer.
 */
carquet_reader_t* carquet_reader_open_buffer(
    const void* buffer,
    size_t size,
    const carquet_reader_options_t* options,
    carquet_error_t* error) {

    /* buffer is nonnull per API contract */
    if (size == 0) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT, "Invalid buffer size");
        return NULL;
    }

    carquet_reader_t* reader = carquet_mem_calloc(1, sizeof(carquet_reader_t));
    if (!reader) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate reader");
        return NULL;
    }

    reader->mmap_data = (const uint8_t*)buffer;
    reader->file_size = size;
    reader->owns_file = false;  /* We don't own the buffer */

    if (options) {
        reader->options = *options;
    } else {
        carquet_reader_options_init(&reader->options);
    }

    /* Initialize arena */
    if (carquet_arena_init(&reader->arena) != CARQUET_OK) {
        carquet_mem_free(reader);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to initialize arena");
        return NULL;
    }

    /* Parse footer from buffer */
    /* Minimum size check */
    if (size < 12) {  /* 4 (magic) + 4 (footer size) + 4 (magic) */
        carquet_arena_destroy(&reader->arena);
        carquet_mem_free(reader);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_FOOTER, "Buffer too small");
        return NULL;
    }

    /* Check magic bytes */
    if (memcmp(buffer, "PAR1", 4) != 0) {
        carquet_arena_destroy(&reader->arena);
        carquet_mem_free(reader);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_MAGIC, "Invalid header magic");
        return NULL;
    }

    const uint8_t* end = (const uint8_t*)buffer + size;
    if (memcmp(end - 4, "PAR1", 4) != 0) {
        bool encrypted = memcmp(end - 4, "PARE", 4) == 0;
        carquet_arena_destroy(&reader->arena);
        carquet_mem_free(reader);
        if (encrypted) {
            /* Parquet modular encryption: say so instead of "corrupt". */
            CARQUET_SET_ERROR(error, CARQUET_ERROR_NOT_IMPLEMENTED,
                "File uses Parquet modular encryption (encrypted footer, magic "
                "\"PARE\"), which carquet does not support");
            return NULL;
        }
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_MAGIC, "Invalid footer magic");
        return NULL;
    }

    /* Get footer size */
    uint32_t footer_size = carquet_read_u32_le(end - 8);
    if (footer_size > size - 8) {
        carquet_arena_destroy(&reader->arena);
        carquet_mem_free(reader);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_FOOTER, "Footer size too large");
        return NULL;
    }

    /* Parse footer */
    const uint8_t* footer_data = end - 8 - footer_size;
    carquet_status_t status = parquet_parse_file_metadata(
        footer_data, footer_size, &reader->arena, &reader->metadata, error);

    if (status != CARQUET_OK) {
        carquet_arena_destroy(&reader->arena);
        carquet_mem_free(reader);
        return NULL;
    }

    /* Build schema - declared in reader_internal.h */
    reader->schema = build_schema(&reader->arena, &reader->metadata, error);
    if (!reader->schema) {
        carquet_arena_destroy(&reader->arena);
        carquet_mem_free(reader);
        return NULL;
    }

    reader->is_open = true;
    return reader;
}
