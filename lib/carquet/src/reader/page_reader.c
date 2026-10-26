/**
 * @file page_reader.c
 * @brief Page reading implementation
 *
 * Handles reading and decoding of Parquet data pages.
 */

#include "core/allocator.h"
#include "core/compat.h"
#include <carquet/carquet.h>
#include "reader_internal.h"
#include "thrift/parquet_types.h"
#include "encoding/plain.h"
#include "encoding/rle.h"
#include "compression/custom.h"
#include "compression/snappy.h"
#include "core/endian.h"
#include "core/bitpack.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>

#if defined(CARQUET_ARCH_ARM) && defined(CARQUET_ENABLE_NEON) && \
    (defined(__ARM_NEON) || defined(__ARM_NEON__))
#include <arm_neon.h>
#endif

/* CRC32 verification */
extern uint32_t carquet_crc32(const uint8_t* data, size_t length);

/* SIMD dispatch functions for dictionary gather */
extern void carquet_dispatch_gather_i32(const int32_t* dict, const uint32_t* indices,
                                         int64_t count, int32_t* output);
extern void carquet_dispatch_gather_i64(const int64_t* dict, const uint32_t* indices,
                                         int64_t count, int64_t* output);
extern void carquet_dispatch_gather_float(const float* dict, const uint32_t* indices,
                                           int64_t count, float* output);
extern void carquet_dispatch_gather_double(const double* dict, const uint32_t* indices,
                                            int64_t count, double* output);
extern bool carquet_dispatch_checked_gather_i32(const int32_t* dict, int32_t dict_count,
                                                 const uint32_t* indices, int64_t count,
                                                 int32_t* output);
extern bool carquet_dispatch_checked_gather_i64(const int64_t* dict, int32_t dict_count,
                                                 const uint32_t* indices, int64_t count,
                                                 int64_t* output);
extern bool carquet_dispatch_checked_gather_float(const float* dict, int32_t dict_count,
                                                   const uint32_t* indices, int64_t count,
                                                   float* output);
extern bool carquet_dispatch_checked_gather_double(const double* dict, int32_t dict_count,
                                                    const uint32_t* indices, int64_t count,
                                                    double* output);

/* SIMD dispatch functions for definition level processing */
extern int64_t carquet_dispatch_count_non_nulls(const int16_t* def_levels, int64_t count,
                                                  int16_t max_def_level);
extern void carquet_dispatch_fill_def_levels(int16_t* def_levels, int64_t count, int16_t value);

/* Forward declarations for compression functions */
extern carquet_status_t carquet_lz4_decompress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);
extern carquet_status_t carquet_lz4_hadoop_decompress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);
/* carquet_snappy_decompress / _diag come from compression/snappy.h */
extern int carquet_gzip_decompress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);
extern int carquet_zstd_decompress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);
extern int carquet_brotli_decompress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);
extern carquet_status_t carquet_byte_stream_split_decode_float(
    const uint8_t* data,
    size_t data_size,
    float* values,
    int64_t count);
extern carquet_status_t carquet_byte_stream_split_decode_double(
    const uint8_t* data,
    size_t data_size,
    double* values,
    int64_t count);
extern carquet_status_t carquet_byte_stream_split_decode(
    const uint8_t* data,
    size_t data_size,
    int32_t type_length,
    uint8_t* values,
    int64_t count);
extern carquet_status_t carquet_delta_decode_int32(
    const uint8_t* data,
    size_t data_size,
    int32_t* values,
    int32_t num_values,
    size_t* bytes_consumed);
extern carquet_status_t carquet_delta_decode_int64(
    const uint8_t* data,
    size_t data_size,
    int64_t* values,
    int32_t num_values,
    size_t* bytes_consumed);
extern carquet_status_t carquet_delta_length_decode(
    const uint8_t* data,
    size_t data_size,
    carquet_byte_array_t* values,
    int32_t num_values,
    size_t* bytes_consumed);
extern carquet_status_t carquet_delta_strings_decode(
    const uint8_t* data,
    size_t data_size,
    carquet_byte_array_t* values,
    int32_t num_values,
    uint8_t* work_buffer,
    size_t work_buffer_size,
    size_t* bytes_consumed);
extern carquet_status_t carquet_delta_strings_decode_lengths(
    const uint8_t* data, size_t data_size, int32_t num_values,
    int32_t** prefix_out, int32_t** suffix_out,
    size_t* header_bytes, size_t* total_size);
extern carquet_status_t carquet_delta_strings_reconstruct(
    const uint8_t* data, size_t data_size, size_t header_bytes,
    const int32_t* prefix_lengths, const int32_t* suffix_lengths,
    carquet_byte_array_t* values, int32_t num_values,
    uint8_t* work_buffer, size_t work_buffer_size, size_t* bytes_consumed);
extern carquet_status_t carquet_delta_strings_decoded_size(
    const uint8_t* data,
    size_t data_size,
    int32_t num_values,
    size_t* required_size);

static bool checked_add_size(size_t a, size_t b, size_t* out);
static bool checked_mul_size(size_t a, size_t b, size_t* out);
static bool checked_add_i64(int64_t a, int64_t b, int64_t* out);

/* ============================================================================
 * Pre-buffered I/O Helper
 * ============================================================================
 */

/**
 * Read data from a file offset, using the prebuffer cache if available.
 * Returns bytes read (0 on failure).
 */
static size_t prebuf_read_at(carquet_reader_t* file_reader,
                             int64_t offset, void* buf, size_t size) {
    if (offset < 0 || size > (size_t)INT64_MAX) {
        return 0;
    }

    size_t start = (size_t)offset;
    size_t end = 0;
    if (!checked_add_size(start, size, &end)) {
        return 0;
    }
    size_t prebuf_start = file_reader->prebuffer.file_offset >= 0
                        ? (size_t)file_reader->prebuffer.file_offset : 0;
    size_t prebuf_end = 0;
    bool prebuf_end_ok = checked_add_size(prebuf_start, file_reader->prebuffer.size,
                                          &prebuf_end);

    /* Check prebuffer cache first */
    if (file_reader->prebuffer.data &&
        offset >= file_reader->prebuffer.file_offset &&
        prebuf_end_ok &&
        start >= prebuf_start &&
        end <= prebuf_end) {
        memcpy(buf, file_reader->prebuffer.data +
               (offset - file_reader->prebuffer.file_offset), size);
        return size;
    }

    /* Fall back to 64-bit aware seek + fread */
    if (carquet_fseek64(file_reader->file, (int64_t)offset, SEEK_SET) != 0) return 0;
    return fread(buf, 1, size, file_reader->file);
}

/* ============================================================================
 * Decompression
 * ============================================================================
 */

carquet_status_t carquet_decompress_page_diag(
    carquet_compression_t codec,
    const uint8_t* compressed,
    size_t compressed_size,
    uint8_t* decompressed,
    size_t decompressed_capacity,
    size_t* decompressed_size,
    char* diag,
    size_t diag_size) {

    if (diag && diag_size > 0) diag[0] = '\0';

    /* Only the built-in Snappy decoder reports a detailed reason; a
     * user-registered replacement keeps the plain path. */
    carquet_custom_codec_t custom_override;
    if (codec == CARQUET_COMPRESSION_SNAPPY && diag && diag_size > 0 &&
        !carquet_custom_codec_lookup(codec, &custom_override)) {
        carquet_snappy_diag_t sd;
        carquet_status_t st = carquet_snappy_decompress_diag(
            compressed, compressed_size, decompressed, decompressed_capacity,
            decompressed_size, &sd);
        if (st != CARQUET_OK && sd.reason != CARQUET_SNAPPY_OK) {
            carquet_snappy_diag_format(&sd, diag, diag_size);
        }
        return st;
    }

    return carquet_decompress_page(codec, compressed, compressed_size,
                                   decompressed, decompressed_capacity,
                                   decompressed_size);
}

carquet_status_t carquet_decompress_page(
    carquet_compression_t codec,
    const uint8_t* compressed,
    size_t compressed_size,
    uint8_t* decompressed,
    size_t decompressed_capacity,
    size_t* decompressed_size) {

    /* Honor user-registered codec implementations before falling through to
     * the built-ins, so callers can replace e.g. GZIP with a HW-accelerated
     * variant or fill the LZO slot (and BROTLI when built without libbrotli). */
    carquet_custom_codec_t custom;
    if (carquet_custom_codec_lookup(codec, &custom)) {
        return custom.decompress(compressed, compressed_size,
                                 decompressed, decompressed_capacity,
                                 decompressed_size, custom.user_data);
    }

    switch (codec) {
        case CARQUET_COMPRESSION_UNCOMPRESSED:
            if (compressed_size > decompressed_capacity) {
                return CARQUET_ERROR_DECOMPRESSION;
            }
            memcpy(decompressed, compressed, compressed_size);
            *decompressed_size = compressed_size;
            return CARQUET_OK;

        case CARQUET_COMPRESSION_SNAPPY:
            return carquet_snappy_decompress(
                compressed, compressed_size,
                decompressed, decompressed_capacity, decompressed_size);

        case CARQUET_COMPRESSION_LZ4:
            /* Codec 5 is the deprecated Hadoop-framed LZ4. */
            return carquet_lz4_hadoop_decompress(
                compressed, compressed_size,
                decompressed, decompressed_capacity, decompressed_size);

        case CARQUET_COMPRESSION_LZ4_RAW:
            return carquet_lz4_decompress(
                compressed, compressed_size,
                decompressed, decompressed_capacity, decompressed_size);

        case CARQUET_COMPRESSION_GZIP:
            return carquet_gzip_decompress(
                compressed, compressed_size,
                decompressed, decompressed_capacity, decompressed_size);

        case CARQUET_COMPRESSION_ZSTD:
            return carquet_zstd_decompress(
                compressed, compressed_size,
                decompressed, decompressed_capacity, decompressed_size);

        case CARQUET_COMPRESSION_BROTLI:
            /* UNSUPPORTED_CODEC when built without libbrotli. */
            return carquet_brotli_decompress(
                compressed, compressed_size,
                decompressed, decompressed_capacity, decompressed_size);

        default:
            return CARQUET_ERROR_UNSUPPORTED_CODEC;
    }
}

/* ============================================================================
 * Utility Functions
 * ============================================================================
 */

#define CARQUET_MAX_PAGE_PAYLOAD_SIZE (256ULL * 1024 * 1024)

static bool checked_add_size(size_t a, size_t b, size_t* out) {
    if (a > SIZE_MAX - b) {
        return false;
    }
    *out = a + b;
    return true;
}

static bool checked_mul_size(size_t a, size_t b, size_t* out) {
    if (a != 0 && b > SIZE_MAX / a) {
        return false;
    }
    *out = a * b;
    return true;
}

static bool checked_add_i64(int64_t a, int64_t b, int64_t* out) {
    if ((b > 0 && a > INT64_MAX - b) ||
        (b < 0 && a < INT64_MIN - b)) {
        return false;
    }
    *out = a + b;
    return true;
}

static carquet_status_t validate_page_payload_size(
    const parquet_page_header_t* header,
    bool allow_empty,
    carquet_error_t* error) {

    if (header->compressed_page_size < 0 ||
        (!allow_empty && header->compressed_page_size == 0) ||
        (size_t)header->compressed_page_size > CARQUET_MAX_PAGE_PAYLOAD_SIZE ||
        header->uncompressed_page_size < 0 ||
        (size_t)header->uncompressed_page_size > CARQUET_MAX_PAGE_PAYLOAD_SIZE) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE, "Page size out of range");
        return CARQUET_ERROR_INVALID_PAGE;
    }

    if (header->type == CARQUET_PAGE_DATA_V2) {
        const parquet_data_page_header_v2_t* v2h = &header->data_page_header_v2;
        if (v2h->num_values <= 0 ||
            v2h->repetition_levels_byte_length < 0 ||
            v2h->definition_levels_byte_length < 0) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE, "Invalid V2 page header");
            return CARQUET_ERROR_INVALID_PAGE;
        }

        size_t levels_size;
        if (!checked_add_size((size_t)v2h->repetition_levels_byte_length,
                              (size_t)v2h->definition_levels_byte_length,
                              &levels_size) ||
            levels_size > (size_t)header->compressed_page_size ||
            levels_size > (size_t)header->uncompressed_page_size) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE, "V2 level sizes exceed page size");
            return CARQUET_ERROR_DECODE;
        }
    } else if (header->type == CARQUET_PAGE_DATA) {
        if (header->data_page_header.num_values <= 0) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE, "Invalid data page value count");
            return CARQUET_ERROR_INVALID_PAGE;
        }
    } else if (header->type == CARQUET_PAGE_DICTIONARY) {
        if (header->dictionary_page_header.num_values < 0) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE, "Invalid dictionary value count");
            return CARQUET_ERROR_INVALID_PAGE;
        }
    }

    return CARQUET_OK;
}

static carquet_status_t validate_page_payload_span(
    const carquet_reader_t* file_reader,
    int64_t page_offset,
    size_t header_size,
    int32_t compressed_size,
    carquet_error_t* error) {

    if (page_offset < 0 || (size_t)page_offset > file_reader->file_size) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE, "Page offset out of range");
        return CARQUET_ERROR_INVALID_PAGE;
    }

    size_t offset = (size_t)page_offset;
    size_t payload_start;
    size_t payload_end;
    if (!checked_add_size(offset, header_size, &payload_start) ||
        !checked_add_size(payload_start, (size_t)compressed_size, &payload_end) ||
        payload_end > file_reader->file_size) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE, "Page payload exceeds file size");
        return CARQUET_ERROR_INVALID_PAGE;
    }

    return CARQUET_OK;
}

/**
 * Resolve the offset of the first data page of a dictionary-encoded chunk.
 *
 * `data_page_offset` is the spec-mandated answer and is what every conforming
 * writer records. Some writers (e.g. DuckDB) instead point it at the dictionary
 * page, so carquet used to *always* recompute the offset as
 * `dictionary_page_offset + dict_header + dict_payload`.
 *
 * That override is wrong whenever anything sits between the dictionary page and
 * the first data page — an index page (ColumnMetaData.index_page_offset is a
 * legal, if rare, field), or writer-inserted alignment padding. In those files
 * the recomputed offset lands inside or short of the first data page, so the
 * page header parses out of a payload and the reader hands garbage bytes to the
 * codec, which then reports the block as corrupt even though the file is fine.
 *
 * Only recompute when `data_page_offset` is *provably* wrong, i.e. when it
 * points at or before the end of the dictionary page. Conforming files keep
 * their own offset.
 */
int64_t carquet_resolve_data_start_offset(int64_t data_page_offset,
                                          int64_t dict_end) {
    if (data_page_offset >= dict_end) {
        return data_page_offset;
    }
    return dict_end;
}

static carquet_status_t ensure_decompress_capacity(
    carquet_column_reader_t* reader,
    size_t needed,
    const char* message,
    carquet_error_t* error) {

    if (needed > CARQUET_MAX_PAGE_PAYLOAD_SIZE) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE, "Page size out of range");
        return CARQUET_ERROR_INVALID_PAGE;
    }

    if (needed > reader->decompress_capacity) {
        uint8_t* new_buf = carquet_mem_realloc(reader->decompress_buffer, needed);
        if (!new_buf) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "%s", message);
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        reader->decompress_buffer = new_buf;
        reader->decompress_capacity = needed;
    }

    return CARQUET_OK;
}

static inline int bit_width_for_max(int max_val) {
    if (max_val == 0) return 0;
    int width = 0;
    while (max_val > 0) {
        width++;
        max_val >>= 1;
    }
    return width;
}

#if defined(CARQUET_ARCH_ARM) && defined(CARQUET_ENABLE_NEON) && \
    (defined(__ARM_NEON) || defined(__ARM_NEON__))
static bool gather_fixed_dictionary_values_neon(const uint8_t* dict_data,
                                                 int32_t dict_count,
                                                 const uint32_t* indices,
                                                 int32_t count,
                                                 size_t value_size,
                                                 uint8_t* output) {
    for (int32_t i = 0; i < count; i++) {
        uint32_t idx = indices[i];
        if (idx >= (uint32_t)dict_count) {
            return false;
        }

        const uint8_t* src = dict_data + (size_t)idx * value_size;
        uint8_t* dst = output + (size_t)i * value_size;

        switch (value_size) {
            case 12:
                vst1_u8(dst, vld1_u8(src));
                memcpy(dst + 8, src + 8, 4);
                break;
            case 16:
                vst1q_u8(dst, vld1q_u8(src));
                break;
            case 32:
                vst1q_u8(dst, vld1q_u8(src));
                vst1q_u8(dst + 16, vld1q_u8(src + 16));
                break;
            default:
                if ((value_size & 15U) == 0 && value_size <= 64) {
                    for (size_t off = 0; off < value_size; off += 16) {
                        vst1q_u8(dst + off, vld1q_u8(src + off));
                    }
                } else {
                    memcpy(dst, src, value_size);
                }
                break;
        }
    }

    return true;
}
#endif

static bool page_values_can_be_viewed_directly(
    const carquet_column_reader_t* reader,
    carquet_encoding_t encoding) {

#if !CARQUET_LITTLE_ENDIAN
    (void)reader;
    (void)encoding;
    return false;
#else
    if (encoding != CARQUET_ENCODING_PLAIN ||
        reader->max_def_level > 0 ||
        reader->max_rep_level > 0) {
        return false;
    }

    switch (reader->type) {
        case CARQUET_PHYSICAL_INT32:
        case CARQUET_PHYSICAL_INT64:
        case CARQUET_PHYSICAL_INT96:
        case CARQUET_PHYSICAL_FLOAT:
        case CARQUET_PHYSICAL_DOUBLE:
            return true;
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
            return reader->type_length > 0;
        case CARQUET_PHYSICAL_BOOLEAN:
        case CARQUET_PHYSICAL_BYTE_ARRAY:
        default:
            return false;
    }
#endif
}

void carquet_column_clear_retained_pages(carquet_column_reader_t* reader) {
    carquet_retained_page_t* p = reader->retained_pages;
    while (p) {
        carquet_retained_page_t* next = p->next;
        carquet_mem_free(p);
        p = next;
    }
    reader->retained_pages = NULL;
}

uint8_t* carquet_column_retain_page(
    carquet_column_reader_t* reader,
    const uint8_t* src,
    size_t size) {
    carquet_retained_page_t* node =
        (carquet_retained_page_t*)carquet_mem_malloc(sizeof(carquet_retained_page_t) + size);
    if (!node) return NULL;
    node->next = reader->retained_pages;
    node->size = size;
    if (size > 0 && src != NULL) {
        memcpy(node->data, src, size);
    }
    reader->retained_pages = node;
    return node->data;
}

static void release_decoded_level_buffers(carquet_column_reader_t* reader) {
    carquet_mem_free(reader->decoded_def_levels);
    carquet_mem_free(reader->decoded_rep_levels);
    reader->decoded_def_levels = NULL;
    reader->decoded_rep_levels = NULL;
}

static carquet_status_t ensure_decoded_page_buffers(
    carquet_column_reader_t* reader,
    int32_t num_values,
    size_t value_size,
    carquet_error_t* error) {

    bool need_def_levels = reader->max_def_level > 0;
    bool need_rep_levels = reader->max_rep_level > 0;
    size_t values_buffer_size = 0;

    if (num_values < 0 ||
        !checked_mul_size(value_size, (size_t)num_values, &values_buffer_size)) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Decode buffer size overflow");
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    /* A single data page can never legitimately contain more values than the
     * column chunk declares. The page header's num_values is attacker-
     * controlled, so without this bound a crafted header drives a multi-GB
     * decode-buffer allocation from a tiny file. */
    if (reader->col_meta &&
        reader->col_meta->num_values >= 0 &&
        (int64_t)num_values > reader->col_meta->num_values) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE,
            "Data page value count (%d) exceeds column chunk total (%lld)",
            num_values, (long long)reader->col_meta->num_values);
        return CARQUET_ERROR_INVALID_PAGE;
    }

    if (reader->decoded_ownership == CARQUET_DATA_VIEW) {
        reader->decoded_values = NULL;
        reader->decoded_capacity = 0;
    }
    reader->decoded_ownership = CARQUET_DATA_OWNED;

    if (!need_def_levels && reader->decoded_def_levels) {
        carquet_mem_free(reader->decoded_def_levels);
        reader->decoded_def_levels = NULL;
    }
    if (!need_rep_levels && reader->decoded_rep_levels) {
        carquet_mem_free(reader->decoded_rep_levels);
        reader->decoded_rep_levels = NULL;
    }

    /* Realloc when the value COUNT grows, or when the per-value byte width
     * changed. The width can change for the same count when a column reader's
     * preserve_dictionary flips on (values decoded as 2/8-byte physical values
     * vs. 4-byte uint32 dictionary indices) — has_dictionary, and therefore
     * preserve mode, only becomes known after the first page loads. Without
     * the width check, the count-keyed capacity would reuse a too-small buffer
     * and the wider decode/copy would overflow it. */
    if ((size_t)num_values > reader->decoded_capacity ||
        value_size != reader->decoded_value_size) {
        carquet_mem_free(reader->decoded_values);
        reader->decoded_values = NULL;
        release_decoded_level_buffers(reader);

        if (values_buffer_size > 0) {
            reader->decoded_values = carquet_mem_calloc(1, values_buffer_size);
        }
        if (need_def_levels && num_values > 0) {
            reader->decoded_def_levels = carquet_mem_malloc(sizeof(int16_t) * (size_t)num_values);
        }
        if (need_rep_levels && num_values > 0) {
            reader->decoded_rep_levels = carquet_mem_malloc(sizeof(int16_t) * (size_t)num_values);
        }
        reader->decoded_capacity = (size_t)num_values;
        reader->decoded_value_size = value_size;
    } else {
        if (need_def_levels && !reader->decoded_def_levels && num_values > 0) {
            reader->decoded_def_levels = carquet_mem_malloc(sizeof(int16_t) * (size_t)num_values);
        }
        if (need_rep_levels && !reader->decoded_rep_levels && num_values > 0) {
            reader->decoded_rep_levels = carquet_mem_malloc(sizeof(int16_t) * (size_t)num_values);
        }
    }

    if ((values_buffer_size > 0 && !reader->decoded_values) ||
        (need_def_levels && num_values > 0 && !reader->decoded_def_levels) ||
        (need_rep_levels && num_values > 0 && !reader->decoded_rep_levels)) {
        carquet_mem_free(reader->decoded_values);
        reader->decoded_values = NULL;
        release_decoded_level_buffers(reader);
        reader->decoded_capacity = 0;
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate decode buffers");
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    return CARQUET_OK;
}

/* ============================================================================
 * Level Decoding
 * ============================================================================
 */

static carquet_status_t decode_levels_rle(
    const uint8_t* data,
    size_t data_size,
    int bit_width,
    int32_t num_values,
    int16_t* levels,
    size_t* bytes_consumed) {

    if (bit_width == 0) {
        /* All zeros */
        memset(levels, 0, num_values * sizeof(int16_t));
        *bytes_consumed = 0;
        return CARQUET_OK;
    }

    /* Use the convenience function for decoding levels */
    int64_t decoded = carquet_rle_decode_levels(
        data, data_size, bit_width, levels, num_values);

    if (decoded < 0) {
        return CARQUET_ERROR_DECODE;
    }

    if (decoded != num_values) {
        return CARQUET_ERROR_DECODE;
    }

    /* Estimate bytes consumed (not perfect, but good enough) */
    *bytes_consumed = data_size;
    return CARQUET_OK;
}

/* Decode the deprecated BIT_PACKED level encoding (Encoding=4). Unlike the
 * RLE/bit-packing hybrid, legacy BIT_PACKED packs values MSB-first with no
 * run headers and no length prefix; the byte length is implied by
 * ceil(num_values * bit_width / 8). Only valid for V1 page levels. */
static carquet_status_t decode_levels_bitpacked(
    const uint8_t* data,
    size_t data_size,
    int bit_width,
    int32_t num_values,
    int16_t* levels,
    size_t* bytes_consumed) {
    return carquet_decode_bitpacked_levels(data, data_size, bit_width,
               num_values, levels, bytes_consumed) == 0
           ? CARQUET_OK : CARQUET_ERROR_DECODE;
}

/* Decode a V1 page level section, dispatching on the (deprecated) BIT_PACKED
 * vs RLE encoding declared in the page header. RLE is length-prefixed (4-byte
 * LE) in V1; BIT_PACKED is not. Advances *ptr / *remaining past the section. */
static carquet_status_t decode_v1_level_section(
    const uint8_t** ptr,
    size_t* remaining,
    carquet_encoding_t level_encoding,
    int bit_width,
    int32_t num_values,
    int16_t* levels,
    carquet_error_t* error) {

    size_t bytes_consumed = 0;
    carquet_status_t status;

    if (level_encoding == CARQUET_ENCODING_BIT_PACKED) {
        status = decode_levels_bitpacked(*ptr, *remaining, bit_width,
                                         num_values, levels, &bytes_consumed);
        if (status != CARQUET_OK) {
            CARQUET_SET_ERROR(error, status, "Failed to decode BIT_PACKED levels");
            return status;
        }
        *ptr += bytes_consumed;
        *remaining -= bytes_consumed;
        return CARQUET_OK;
    }

    /* RLE (default): 4-byte LE length prefix, then the hybrid stream. */
    if (*remaining < 4) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE, "Truncated levels");
        return CARQUET_ERROR_DECODE;
    }
    uint32_t sz = carquet_read_u32_le(*ptr);
    *ptr += 4;
    *remaining -= 4;
    if (sz > *remaining) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE, "Invalid level size");
        return CARQUET_ERROR_DECODE;
    }
    status = decode_levels_rle(*ptr, sz, bit_width, num_values, levels,
                               &bytes_consumed);
    if (status != CARQUET_OK) {
        CARQUET_SET_ERROR(error, status, "Failed to decode RLE levels");
        return status;
    }
    *ptr += sz;
    *remaining -= sz;
    return CARQUET_OK;
}

/* ============================================================================
 * Dictionary Page Reading
 * ============================================================================
 */

carquet_status_t carquet_read_dictionary_page(
    carquet_column_reader_t* reader,
    uint8_t* page_data,
    size_t page_size,
    const parquet_dictionary_page_header_t* header,
    carquet_data_ownership_t ownership,
    carquet_error_t* error) {

    /* Ownership contract: page_data is owned by the CALLER. On error this
     * function never frees page_data (it only frees its own allocations and
     * NULLs reader->dictionary_data); the caller frees page_data on a non-OK
     * return. On success the reader adopts page_data (dictionary_data) and
     * frees it later via reset/close. Freeing page_data here previously caused
     * a double free with the caller's error-path free. */
    if (!reader || !page_data || !header) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT, "NULL argument");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    if (header->num_values < 0) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE, "Invalid dictionary value count");
        return CARQUET_ERROR_INVALID_PAGE;
    }

    /* Allocate dictionary storage */
    size_t value_size = 0;
    switch (reader->type) {
        case CARQUET_PHYSICAL_INT32:
        case CARQUET_PHYSICAL_FLOAT:
            value_size = 4;
            break;
        case CARQUET_PHYSICAL_INT64:
        case CARQUET_PHYSICAL_DOUBLE:
            value_size = 8;
            break;
        case CARQUET_PHYSICAL_INT96:
            value_size = 12;
            break;
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
            value_size = reader->type_length;
            break;
        case CARQUET_PHYSICAL_BYTE_ARRAY:
            /* Variable length - will be handled differently */
            break;
        default:
            break;
    }

    reader->dictionary_count = header->num_values;
    reader->dictionary_ownership = ownership;

    if (reader->type == CARQUET_PHYSICAL_BYTE_ARRAY) {
        /* For variable length, keep raw dictionary bytes and build offsets once */
        reader->dictionary_data = page_data;
        reader->dictionary_size = page_size;

        /* Build offset table for O(1) BYTE_ARRAY lookup. Every entry needs
         * at least its 4-byte length prefix, so an (untrusted) num_values
         * above page_size / 4 is malformed; reject it before it sizes an
         * allocation (fuzzer: num_values = 1.7e9 asked for 6.9 GB). */
        if ((size_t)header->num_values > page_size / 4) {
            reader->dictionary_data = NULL;  /* caller frees page_data */
            CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE,
                "Dictionary declares %d entries but holds %llu bytes",
                header->num_values, (unsigned long long)page_size);
            return CARQUET_ERROR_DECODE;
        }
        reader->dictionary_offsets = carquet_mem_malloc((size_t)header->num_values * sizeof(uint32_t));
        if (!reader->dictionary_offsets) {
            reader->dictionary_data = NULL;  /* caller frees page_data */
            CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate offset table");
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }

        /* Scan dictionary once to build offset table */
        const uint8_t* dict_ptr = page_data;
        size_t dict_remaining = page_size;
        for (int32_t i = 0; i < header->num_values; i++) {
            if (dict_remaining < 4) {
                carquet_mem_free(reader->dictionary_offsets);
                reader->dictionary_data = NULL;  /* caller frees page_data */
                reader->dictionary_offsets = NULL;
                CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE, "Truncated dictionary");
                return CARQUET_ERROR_DECODE;
            }
            reader->dictionary_offsets[i] = (uint32_t)(dict_ptr - page_data);
            uint32_t len = carquet_read_u32_le(dict_ptr);
            /* Widen to size_t before adding: `4 + len` in 32-bit unsigned
             * arithmetic wraps for len >= 0xFFFFFFFC, which would defeat the
             * `dict_remaining < entry_size` bounds check below. */
            size_t entry_size = (size_t)4 + (size_t)len;
            if (dict_remaining < entry_size) {
                carquet_mem_free(reader->dictionary_offsets);
                reader->dictionary_data = NULL;  /* caller frees page_data */
                reader->dictionary_offsets = NULL;
                CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE, "Invalid dictionary entry");
                return CARQUET_ERROR_DECODE;
            }
            dict_ptr += entry_size;
            dict_remaining -= entry_size;
        }
    } else {
        /* Fixed size values */
        size_t dict_size = 0;
        if (!checked_mul_size(value_size, (size_t)header->num_values, &dict_size)) {
            /* page_data not yet adopted; caller frees it on this error. */
            CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE, "Dictionary size overflow");
            return CARQUET_ERROR_DECODE;
        }
        if (dict_size > page_size) {
            /* page_data not yet adopted; caller frees it on this error. */
            CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE, "Truncated dictionary");
            return CARQUET_ERROR_DECODE;
        }
        reader->dictionary_data = page_data;
        reader->dictionary_size = dict_size;
    }

    reader->has_dictionary = true;
    return CARQUET_OK;
}

/* ============================================================================
 * Phase 3 encoding decode dispatch (DELTA_*, BYTE_STREAM_SPLIT int/FLBA)
 * ============================================================================
 *
 * Returns true via *handled if the encoding was one of the Phase 3 encodings
 * this helper owns; in that case *status holds the decode result.
 *
 * Byte-array lifetime:
 *  - DELTA_LENGTH_BYTE_ARRAY: value .data pointers reference into `ptr` (the
 *    page payload), exactly like PLAIN BYTE_ARRAY. *needs_page_retain is set so
 *    the caller routes it through the same page-retain+fixup path.
 *  - DELTA_BYTE_ARRAY: strings are reconstructed into a scratch buffer that is
 *    allocated through carquet_column_retain_page() (so it lives until the
 *    batch is consumed / row-group reset / reader close). Values are decoded
 *    directly into that retained buffer, so no pointer fixup is required and
 *    there is no use-after-free.
 *  - DELTA_BINARY_PACKED / BSS int/FLBA: decoded into reader->decoded_values
 *    (fixed-size, owned) — no extra lifetime concern.
 */
static carquet_status_t decode_phase3_values(
    carquet_column_reader_t* reader,
    int32_t encoding,
    const uint8_t* ptr,
    size_t remaining,
    void* values,
    int32_t non_null_count,
    bool* handled,
    bool* needs_page_retain,
    carquet_error_t* error) {

    *handled = true;
    *needs_page_retain = false;
    size_t consumed = 0;

    switch (encoding) {
        case CARQUET_ENCODING_DELTA_BINARY_PACKED:
            if (reader->type == CARQUET_PHYSICAL_INT32) {
                return carquet_delta_decode_int32(
                    ptr, remaining, (int32_t*)values, non_null_count, &consumed);
            } else if (reader->type == CARQUET_PHYSICAL_INT64) {
                return carquet_delta_decode_int64(
                    ptr, remaining, (int64_t*)values, non_null_count, &consumed);
            }
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ENCODING,
                "DELTA_BINARY_PACKED requires INT32/INT64");
            return CARQUET_ERROR_INVALID_ENCODING;

        case CARQUET_ENCODING_DELTA_LENGTH_BYTE_ARRAY:
            if (reader->type != CARQUET_PHYSICAL_BYTE_ARRAY) {
                CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ENCODING,
                    "DELTA_LENGTH_BYTE_ARRAY requires BYTE_ARRAY");
                return CARQUET_ERROR_INVALID_ENCODING;
            }
            /* Values point into the page payload — same lifetime as PLAIN
             * BYTE_ARRAY, so the caller must retain the page buffer. */
            *needs_page_retain = true;
            return carquet_delta_length_decode(
                ptr, remaining, (carquet_byte_array_t*)values,
                non_null_count, &consumed);

        case CARQUET_ENCODING_DELTA_BYTE_ARRAY:
            if (reader->type != CARQUET_PHYSICAL_BYTE_ARRAY &&
                reader->type != CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY) {
                CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ENCODING,
                    "DELTA_BYTE_ARRAY requires BYTE_ARRAY/FIXED_LEN_BYTE_ARRAY");
                return CARQUET_ERROR_INVALID_ENCODING;
            }
            {
                if (non_null_count <= 0) {
                    return CARQUET_OK;
                }
                /* Decode the two length streams once: they give the exact
                 * reconstruction size (no guessing) and are reused below. */
                int32_t* pfx = NULL;
                int32_t* sfx = NULL;
                size_t hdr_bytes = 0, work_size = 0;
                carquet_status_t st = carquet_delta_strings_decode_lengths(
                    ptr, remaining, non_null_count, &pfx, &sfx, &hdr_bytes, &work_size);
                if (st != CARQUET_OK) {
                    return st;
                }
                /* Allocate the scratch through the retain list so it outlives
                 * the batch (freed on row-group reset / reader close). */
                uint8_t* work = carquet_column_retain_page(
                    reader, NULL, work_size == 0 ? 1 : work_size);
                if (!work) {
                    carquet_mem_free(pfx);
                    carquet_mem_free(sfx);
                    CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY,
                        "Failed to allocate DELTA_BYTE_ARRAY scratch");
                    return CARQUET_ERROR_OUT_OF_MEMORY;
                }
                /* For FIXED_LEN_BYTE_ARRAY the column buffer holds raw fixed
                 * width values; decode into a temporary byte-array view, then
                 * copy the reconstructed bytes out (values are length
                 * type_length). For BYTE_ARRAY decode straight into the
                 * carquet_byte_array_t output. */
                if (reader->type == CARQUET_PHYSICAL_BYTE_ARRAY) {
                    st = carquet_delta_strings_reconstruct(
                        ptr, remaining, hdr_bytes, pfx, sfx,
                        (carquet_byte_array_t*)values, non_null_count,
                        work, work_size, &consumed);
                    carquet_mem_free(pfx);
                    carquet_mem_free(sfx);
                    return st;
                } else {
                    carquet_byte_array_t* tmp = (carquet_byte_array_t*)carquet_mem_malloc(
                        (size_t)non_null_count * sizeof(carquet_byte_array_t));
                    if (!tmp) {
                        carquet_mem_free(pfx);
                        carquet_mem_free(sfx);
                        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY,
                            "Failed to allocate FLBA decode view");
                        return CARQUET_ERROR_OUT_OF_MEMORY;
                    }
                    st = carquet_delta_strings_reconstruct(
                        ptr, remaining, hdr_bytes, pfx, sfx, tmp, non_null_count,
                        work, work_size, &consumed);
                    carquet_mem_free(pfx);
                    carquet_mem_free(sfx);
                    if (st != CARQUET_OK) {
                        carquet_mem_free(tmp);
                        return st;
                    }
                    uint8_t* out = (uint8_t*)values;
                    size_t len = (size_t)reader->type_length;
                    for (int32_t i = 0; i < non_null_count; i++) {
                        if ((size_t)tmp[i].length != len) {
                            carquet_mem_free(tmp);
                            CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE,
                                "DELTA_BYTE_ARRAY FLBA value length mismatch");
                            return CARQUET_ERROR_DECODE;
                        }
                        memcpy(out + (size_t)i * len, tmp[i].data, len);
                    }
                    carquet_mem_free(tmp);
                    return CARQUET_OK;
                }
            }

        case CARQUET_ENCODING_BYTE_STREAM_SPLIT:
            switch (reader->type) {
                case CARQUET_PHYSICAL_FLOAT:
                    return carquet_byte_stream_split_decode_float(
                        ptr, remaining, (float*)values, non_null_count);
                case CARQUET_PHYSICAL_DOUBLE:
                    return carquet_byte_stream_split_decode_double(
                        ptr, remaining, (double*)values, non_null_count);
                case CARQUET_PHYSICAL_INT32:
                    return carquet_byte_stream_split_decode(
                        ptr, remaining, 4, (uint8_t*)values, non_null_count);
                case CARQUET_PHYSICAL_INT64:
                    return carquet_byte_stream_split_decode(
                        ptr, remaining, 8, (uint8_t*)values, non_null_count);
                case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
                    return carquet_byte_stream_split_decode(
                        ptr, remaining, reader->type_length,
                        (uint8_t*)values, non_null_count);
                default:
                    CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ENCODING,
                        "BYTE_STREAM_SPLIT unsupported type");
                    return CARQUET_ERROR_INVALID_ENCODING;
            }

        case CARQUET_ENCODING_RLE:
            /* RLE as a value encoding is defined only for BOOLEAN. Layout: a
             * 4-byte little-endian length prefix followed by the RLE/bit-packed
             * hybrid at bit width 1. */
            if (reader->type != CARQUET_PHYSICAL_BOOLEAN) {
                CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ENCODING,
                    "RLE value encoding requires BOOLEAN");
                return CARQUET_ERROR_INVALID_ENCODING;
            }
            if (non_null_count <= 0) {
                return CARQUET_OK;
            }
            if (remaining < 4) {
                CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE,
                    "RLE boolean: truncated length prefix");
                return CARQUET_ERROR_DECODE;
            }
            {
                uint32_t rle_len = (uint32_t)ptr[0] | ((uint32_t)ptr[1] << 8) |
                                   ((uint32_t)ptr[2] << 16) | ((uint32_t)ptr[3] << 24);
                if ((size_t)rle_len > remaining - 4) {
                    CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE,
                        "RLE boolean: data length exceeds page");
                    return CARQUET_ERROR_DECODE;
                }
                uint32_t* tmp = carquet_mem_malloc(
                    (size_t)non_null_count * sizeof(uint32_t));
                if (!tmp) {
                    CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY,
                        "RLE boolean scratch allocation failed");
                    return CARQUET_ERROR_OUT_OF_MEMORY;
                }
                int64_t decoded = carquet_rle_decode_all(
                    ptr + 4, rle_len, 1, tmp, non_null_count);
                if (decoded != non_null_count) {
                    carquet_mem_free(tmp);
                    CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE,
                        "RLE boolean: decoded count mismatch");
                    return CARQUET_ERROR_DECODE;
                }
                uint8_t* out = (uint8_t*)values;
                for (int64_t i = 0; i < non_null_count; i++) {
                    out[i] = (uint8_t)(tmp[i] & 1u);
                }
                carquet_mem_free(tmp);
                return CARQUET_OK;
            }

        default:
            *handled = false;
            return CARQUET_OK;
    }
}

/* ============================================================================
 * Data Page Reading
 * ============================================================================
 */

carquet_status_t carquet_read_data_page_v1(
    carquet_column_reader_t* reader,
    const uint8_t* page_data,
    size_t page_size,
    const parquet_data_page_header_t* header,
    void* values,
    int64_t max_values,
    int16_t* def_levels,
    int16_t* rep_levels,
    int64_t* values_read,
    carquet_error_t* error) {

    const uint8_t* ptr = page_data;
    size_t remaining = page_size;

    int32_t num_values = header->num_values;
    if (num_values > max_values) {
        num_values = (int32_t)max_values;
    }

    /* Decode repetition levels if needed */
    if (reader->max_rep_level > 0 && rep_levels) {
        int bit_width = bit_width_for_max(reader->max_rep_level);
        carquet_status_t status = decode_v1_level_section(
            &ptr, &remaining, header->repetition_level_encoding,
            bit_width, num_values, rep_levels, error);
        if (status != CARQUET_OK) {
            return status;
        }
    } else if (rep_levels) {
        memset(rep_levels, 0, num_values * sizeof(int16_t));
    }

    /* Decode definition levels if needed */
    if (reader->max_def_level > 0 && def_levels) {
        int bit_width = bit_width_for_max(reader->max_def_level);
        carquet_status_t status = decode_v1_level_section(
            &ptr, &remaining, header->definition_level_encoding,
            bit_width, num_values, def_levels, error);
        if (status != CARQUET_OK) {
            return status;
        }
    } else if (def_levels) {
        /* Set all to max level (all values present) - use SIMD dispatch */
        carquet_dispatch_fill_def_levels(def_levels, num_values, reader->max_def_level);
    }

    /* Count non-null values */
    int32_t non_null_count = num_values;
    if (def_levels && reader->max_def_level > 0) {
        non_null_count = (int32_t)carquet_dispatch_count_non_nulls(
            def_levels, num_values, reader->max_def_level);
    }

    /* Decode values based on encoding */
    carquet_status_t status = CARQUET_OK;

    switch (header->encoding) {
        case CARQUET_ENCODING_PLAIN:
            /* In dictionary-preserving mode the destination buffer is sized for
             * uint32_t indices (sizeof(uint32_t) per value). A column chunk that
             * starts dictionary-encoded but falls back to a PLAIN data page
             * mid-chunk would have carquet_decode_plain() write full physical
             * values (e.g. 16-byte carquet_byte_array_t for BYTE_ARRAY) into that
             * narrow buffer, overrunning the heap. Reject rather than corrupt. */
            if (reader->preserve_dictionary) {
                CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ENCODING,
                    "Cannot preserve dictionary: column chunk falls back to PLAIN "
                    "encoding mid-chunk (mixed encodings)");
                return CARQUET_ERROR_INVALID_ENCODING;
            }
            {
                int64_t bytes = carquet_decode_plain(
                    ptr, remaining, reader->type, reader->type_length,
                    values, non_null_count);
                if (bytes < 0) {
                    status = CARQUET_ERROR_DECODE;
                }
            }
            break;

        case CARQUET_ENCODING_RLE_DICTIONARY:
        case CARQUET_ENCODING_PLAIN_DICTIONARY:
            if (!reader->has_dictionary) {
                CARQUET_SET_ERROR(error, CARQUET_ERROR_DICTIONARY_NOT_FOUND,
                    "Dictionary encoding without dictionary");
                return CARQUET_ERROR_DICTIONARY_NOT_FOUND;
            }
            /* Decode dictionary indices using RLE */
            {
                /* Read bit width byte */
                if (remaining < 1) {
                    CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE, "Missing bit width");
                    return CARQUET_ERROR_DECODE;
                }
                int bit_width = ptr[0];
                ptr++;
                remaining--;
                if (bit_width > 32) {
                    /* The RLE decoder clamps a width above 32 to 0 and would
                     * then hand back `count` copies of dictionary entry 0 with
                     * no error. */
                    CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE,
                        "Dictionary index bit width %d exceeds 32", bit_width);
                    return CARQUET_ERROR_DECODE;
                }

                int32_t encoded_count = non_null_count;

                /* Dictionary preservation: decode indices directly into output */
                if (reader->preserve_dictionary) {
                    uint32_t* out_indices = (uint32_t*)values;
                    int64_t decoded = carquet_rle_decode_all(
                        ptr, remaining, bit_width, out_indices, encoded_count);
                    /* A short count means a truncated index stream; the tail
                     * of the output would be whatever the buffer held before. */
                    if (decoded != encoded_count) {
                        CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE,
                            "Dictionary indices truncated: decoded %lld of %lld",
                            (long long)decoded, (long long)encoded_count);
                        return CARQUET_ERROR_DECODE;
                    }
                    break;
                }

                /* Use reusable indices buffer to avoid per-page allocation */
                uint32_t* indices;
                if ((size_t)encoded_count <= reader->indices_capacity) {
                    indices = reader->indices_buffer;
                } else {
                    /* Need larger buffer - reallocate */
                    carquet_mem_free(reader->indices_buffer);
                    reader->indices_buffer = carquet_mem_malloc((size_t)encoded_count * sizeof(uint32_t));
                    if (!reader->indices_buffer) {
                        reader->indices_capacity = 0;
                        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate indices");
                        return CARQUET_ERROR_OUT_OF_MEMORY;
                    }
                    reader->indices_capacity = encoded_count;
                    indices = reader->indices_buffer;
                }

                int64_t decoded = carquet_rle_decode_all(
                    ptr, remaining, bit_width, indices, encoded_count);

                if (decoded != encoded_count) {
                    CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE,
                        "Dictionary indices truncated: decoded %lld of %lld",
                        (long long)decoded, (long long)encoded_count);
                    return CARQUET_ERROR_DECODE;
                }

                /* Look up dictionary values for the dense non-null stream */
                if (reader->type == CARQUET_PHYSICAL_BYTE_ARRAY) {
                    /* BYTE_ARRAY: dictionary is stored as length-prefixed values */
                    carquet_byte_array_t* out = (carquet_byte_array_t*)values;

                    /* Use O(1) offset table lookup (built when dictionary was read) */
                    if (reader->dictionary_offsets) {
                        for (int32_t i = 0; i < encoded_count; i++) {
                            int32_t idx = (int32_t)indices[i];
                            if (idx < 0 || idx >= reader->dictionary_count) {
                                status = CARQUET_ERROR_DECODE;
                                break;
                            }

                            /* Direct O(1) lookup using offset table */
                            uint32_t offset = reader->dictionary_offsets[idx];
                            const uint8_t* dict_ptr = reader->dictionary_data + offset;
                            uint32_t len = carquet_read_u32_le(dict_ptr);
                            out[i].data = (uint8_t*)(dict_ptr + 4);
                            out[i].length = (int32_t)len;
                        }
                    } else {
                        /* Fallback: scan each time (shouldn't happen for new readers).
                         * Bounds-checked and size_t-widened as defense-in-depth: a
                         * 32-bit `4 + len` could otherwise wrap and walk past the
                         * dictionary buffer. */
                        const uint8_t* dict_end =
                            reader->dictionary_data + reader->dictionary_size;
                        for (int32_t i = 0; i < encoded_count; i++) {
                            int32_t idx = (int32_t)indices[i];
                            if (idx < 0 || idx >= reader->dictionary_count) {
                                status = CARQUET_ERROR_DECODE;
                                break;
                            }

                            const uint8_t* dict_ptr = reader->dictionary_data;
                            bool oob = false;
                            for (int32_t j = 0; j <= idx; j++) {
                                if ((size_t)(dict_end - dict_ptr) < 4) {
                                    oob = true;
                                    break;
                                }
                                uint32_t len = carquet_read_u32_le(dict_ptr);
                                if ((size_t)(dict_end - dict_ptr) - 4 < (size_t)len) {
                                    oob = true;
                                    break;
                                }
                                if (j == idx) {
                                    out[i].data = (uint8_t*)(dict_ptr + 4);
                                    out[i].length = (int32_t)len;
                                    break;
                                }
                                dict_ptr += (size_t)4 + (size_t)len;
                            }
                            if (oob) {
                                status = CARQUET_ERROR_DECODE;
                                break;
                            }
                        }
                    }
                } else {
                    /* Use SIMD-optimized gather for common types */
                    switch (reader->type) {
                        case CARQUET_PHYSICAL_INT32:
                            if (!carquet_dispatch_checked_gather_i32(
                                    (const int32_t*)reader->dictionary_data,
                                    reader->dictionary_count,
                                    indices, encoded_count, (int32_t*)values)) {
                                CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE,
                                                  "Dictionary index out of bounds");
                                return CARQUET_ERROR_DECODE;
                            }
                            break;
                        case CARQUET_PHYSICAL_INT64:
                            if (!carquet_dispatch_checked_gather_i64(
                                    (const int64_t*)reader->dictionary_data,
                                    reader->dictionary_count,
                                    indices, encoded_count, (int64_t*)values)) {
                                CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE,
                                                  "Dictionary index out of bounds");
                                return CARQUET_ERROR_DECODE;
                            }
                            break;
                        case CARQUET_PHYSICAL_FLOAT:
                            if (!carquet_dispatch_checked_gather_float(
                                    (const float*)reader->dictionary_data,
                                    reader->dictionary_count,
                                    indices, encoded_count, (float*)values)) {
                                CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE,
                                                  "Dictionary index out of bounds");
                                return CARQUET_ERROR_DECODE;
                            }
                            break;
                        case CARQUET_PHYSICAL_DOUBLE:
                            if (!carquet_dispatch_checked_gather_double(
                                    (const double*)reader->dictionary_data,
                                    reader->dictionary_count,
                                    indices, encoded_count, (double*)values)) {
                                CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE,
                                                  "Dictionary index out of bounds");
                                return CARQUET_ERROR_DECODE;
                            }
                            break;
                        case CARQUET_PHYSICAL_INT96:
                        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
                            {
                                size_t value_size = (reader->type == CARQUET_PHYSICAL_INT96)
                                    ? 12 : (size_t)reader->type_length;
                                uint8_t* out = (uint8_t*)values;
                                bool ok;
#if defined(CARQUET_ARCH_ARM) && defined(CARQUET_ENABLE_NEON) && \
    (defined(__ARM_NEON) || defined(__ARM_NEON__))
                                ok = gather_fixed_dictionary_values_neon(
                                    reader->dictionary_data,
                                    reader->dictionary_count,
                                    indices,
                                    encoded_count,
                                    value_size,
                                    out);
#else
                                ok = true;
                                for (int32_t i = 0; i < encoded_count; i++) {
                                    uint32_t idx = indices[i];
                                    if (idx >= (uint32_t)reader->dictionary_count) {
                                        ok = false;
                                        break;
                                    }
                                    memcpy(out + (size_t)i * value_size,
                                           reader->dictionary_data + (size_t)idx * value_size,
                                           value_size);
                                }
#endif
                                if (!ok) {
                                    CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE,
                                                      "Dictionary index out of bounds");
                                    return CARQUET_ERROR_DECODE;
                                }
                            }
                            break;
                        default:
                            break;
                    }
                }
                /* indices buffer is reused, don't free */
            }
            break;

        default:
            {
                /* Phase 3 encodings: DELTA_BINARY_PACKED,
                 * DELTA_LENGTH_BYTE_ARRAY, DELTA_BYTE_ARRAY, and
                 * BYTE_STREAM_SPLIT for FLOAT/DOUBLE/INT32/INT64/FLBA. */
                bool handled = false;
                bool needs_page_retain = false;
                status = decode_phase3_values(
                    reader, header->encoding, ptr, remaining,
                    values, non_null_count,
                    &handled, &needs_page_retain, error);
                if (!handled) {
                    CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ENCODING,
                        "Unsupported encoding: %d", header->encoding);
                    return CARQUET_ERROR_INVALID_ENCODING;
                }
            }
            break;
    }

    if (status != CARQUET_OK) {
        CARQUET_SET_ERROR(error, status, "Failed to decode values");
        return status;
    }

    *values_read = num_values;
    return CARQUET_OK;
}

/* ============================================================================
 * Data Page V2 Reading
 * ============================================================================
 *
 * V2 page layout: [rep_levels_bytes | def_levels_bytes | data_bytes]
 * - Rep/def levels are NOT compressed (stored before compressed data)
 * - No 4-byte length prefixes for levels (byte lengths come from header)
 * - Levels are RLE-encoded (same as V1, but without length prefix)
 * - Data portion may or may not be compressed (header.is_compressed)
 */

carquet_status_t carquet_read_data_page_v2(
    carquet_column_reader_t* reader,
    const uint8_t* page_data,
    size_t page_size,
    const parquet_data_page_header_v2_t* header,
    void* values,
    int64_t max_values,
    int16_t* def_levels,
    int16_t* rep_levels,
    int64_t* values_read,
    carquet_error_t* error) {

    const uint8_t* ptr = page_data;
    size_t remaining = page_size;
    size_t bytes_consumed;

    int32_t num_values = header->num_values;
    if (num_values > max_values) {
        num_values = (int32_t)max_values;
    }

    /* V2: Repetition levels come first, with known byte length (no length prefix) */
    if (reader->max_rep_level > 0 && rep_levels) {
        int32_t rep_bytes = header->repetition_levels_byte_length;
        if (rep_bytes < 0 || (size_t)rep_bytes > remaining) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE, "Invalid V2 rep level size");
            return CARQUET_ERROR_DECODE;
        }

        if (rep_bytes > 0) {
            int bit_width = bit_width_for_max(reader->max_rep_level);
            carquet_status_t status = decode_levels_rle(
                ptr, (size_t)rep_bytes, bit_width, num_values, rep_levels, &bytes_consumed);
            if (status != CARQUET_OK) {
                CARQUET_SET_ERROR(error, status, "Failed to decode V2 rep levels");
                return status;
            }
        } else {
            memset(rep_levels, 0, num_values * sizeof(int16_t));
        }
        ptr += rep_bytes;
        remaining -= (size_t)rep_bytes;
    } else {
        /* Skip rep level bytes even if we don't need them */
        int32_t rep_bytes = header->repetition_levels_byte_length;
        if (rep_bytes < 0 || (size_t)rep_bytes > remaining) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE, "Invalid V2 rep level size");
            return CARQUET_ERROR_DECODE;
        }
        ptr += rep_bytes;
        remaining -= (size_t)rep_bytes;
        if (rep_levels) {
            memset(rep_levels, 0, num_values * sizeof(int16_t));
        }
    }

    /* V2: Definition levels come second, with known byte length (no length prefix) */
    if (reader->max_def_level > 0 && def_levels) {
        int32_t def_bytes = header->definition_levels_byte_length;
        if (def_bytes < 0 || (size_t)def_bytes > remaining) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE, "Invalid V2 def level size");
            return CARQUET_ERROR_DECODE;
        }

        if (def_bytes > 0) {
            int bit_width = bit_width_for_max(reader->max_def_level);
            carquet_status_t status = decode_levels_rle(
                ptr, (size_t)def_bytes, bit_width, num_values, def_levels, &bytes_consumed);
            if (status != CARQUET_OK) {
                CARQUET_SET_ERROR(error, status, "Failed to decode V2 def levels");
                return status;
            }
        } else {
            memset(def_levels, 0, num_values * sizeof(int16_t));
        }
        ptr += def_bytes;
        remaining -= (size_t)def_bytes;
    } else {
        /* Skip def level bytes even if we don't need them */
        int32_t def_bytes = header->definition_levels_byte_length;
        if (def_bytes < 0 || (size_t)def_bytes > remaining) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE, "Invalid V2 def level size");
            return CARQUET_ERROR_DECODE;
        }
        ptr += def_bytes;
        remaining -= (size_t)def_bytes;
        if (def_levels) {
            /* Set all to max level (all values present) - use SIMD dispatch */
            carquet_dispatch_fill_def_levels(def_levels, num_values, reader->max_def_level);
        }
    }

    /* V2: Remaining bytes are the data payload.
     * Note: For V2, decompression of the data portion is handled by the caller
     * (load_next_page_mmap/fread) BEFORE calling this function, since the caller
     * must decompress only the data portion while leaving levels uncompressed.
     * By the time we get here, ptr points to uncompressed data. */

    /* Count non-null values */
    int32_t non_null_count = num_values;
    if (def_levels && reader->max_def_level > 0) {
        non_null_count = (int32_t)carquet_dispatch_count_non_nulls(
            def_levels, num_values, reader->max_def_level);
    }

    /* Decode values based on encoding - reuse V1 value decoding logic */
    carquet_status_t status = CARQUET_OK;

    switch (header->encoding) {
        case CARQUET_ENCODING_PLAIN:
            /* In dictionary-preserving mode the destination buffer is sized for
             * uint32_t indices (sizeof(uint32_t) per value). A column chunk that
             * starts dictionary-encoded but falls back to a PLAIN data page
             * mid-chunk would have carquet_decode_plain() write full physical
             * values (e.g. 16-byte carquet_byte_array_t for BYTE_ARRAY) into that
             * narrow buffer, overrunning the heap. Reject rather than corrupt. */
            if (reader->preserve_dictionary) {
                CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ENCODING,
                    "Cannot preserve dictionary: column chunk falls back to PLAIN "
                    "encoding mid-chunk (mixed encodings)");
                return CARQUET_ERROR_INVALID_ENCODING;
            }
            {
                int64_t bytes = carquet_decode_plain(
                    ptr, remaining, reader->type, reader->type_length,
                    values, non_null_count);
                if (bytes < 0) {
                    status = CARQUET_ERROR_DECODE;
                }
            }
            break;

        case CARQUET_ENCODING_RLE_DICTIONARY:
        case CARQUET_ENCODING_PLAIN_DICTIONARY:
            if (!reader->has_dictionary) {
                CARQUET_SET_ERROR(error, CARQUET_ERROR_DICTIONARY_NOT_FOUND,
                    "Dictionary encoding without dictionary");
                return CARQUET_ERROR_DICTIONARY_NOT_FOUND;
            }
            /* Decode dictionary indices using RLE */
            {
                /* Read bit width byte */
                if (remaining < 1) {
                    CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE, "Missing bit width");
                    return CARQUET_ERROR_DECODE;
                }
                int bit_width = ptr[0];
                ptr++;
                remaining--;
                if (bit_width > 32) {
                    /* The RLE decoder clamps a width above 32 to 0 and would
                     * then hand back `count` copies of dictionary entry 0 with
                     * no error. */
                    CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE,
                        "Dictionary index bit width %d exceeds 32", bit_width);
                    return CARQUET_ERROR_DECODE;
                }
                int32_t encoded_count = non_null_count;

                /* Dictionary preservation: decode indices directly into output */
                if (reader->preserve_dictionary) {
                    uint32_t* out_indices = (uint32_t*)values;
                    int64_t decoded = carquet_rle_decode_all(
                        ptr, remaining, bit_width, out_indices, encoded_count);
                    /* A short count means a truncated index stream; the tail
                     * of the output would be whatever the buffer held before. */
                    if (decoded != encoded_count) {
                        CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE,
                            "Dictionary indices truncated: decoded %lld of %lld",
                            (long long)decoded, (long long)encoded_count);
                        return CARQUET_ERROR_DECODE;
                    }
                    break;
                }

                /* Use reusable indices buffer */
                uint32_t* indices;
                if ((size_t)encoded_count <= reader->indices_capacity) {
                    indices = reader->indices_buffer;
                } else {
                    carquet_mem_free(reader->indices_buffer);
                    reader->indices_buffer = carquet_mem_malloc((size_t)encoded_count * sizeof(uint32_t));
                    if (!reader->indices_buffer) {
                        reader->indices_capacity = 0;
                        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate indices");
                        return CARQUET_ERROR_OUT_OF_MEMORY;
                    }
                    reader->indices_capacity = encoded_count;
                    indices = reader->indices_buffer;
                }

                int64_t decoded = carquet_rle_decode_all(
                    ptr, remaining, bit_width, indices, encoded_count);

                if (decoded != encoded_count) {
                    CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE,
                        "Dictionary indices truncated: decoded %lld of %lld",
                        (long long)decoded, (long long)encoded_count);
                    return CARQUET_ERROR_DECODE;
                }

                /* Look up values from dictionary — same logic as V1 */
                if (reader->type == CARQUET_PHYSICAL_BYTE_ARRAY) {
                    carquet_byte_array_t* out = (carquet_byte_array_t*)values;

                    if (reader->dictionary_offsets) {
                        for (int32_t i = 0; i < encoded_count; i++) {
                            int32_t idx = (int32_t)indices[i];
                            if (idx < 0 || idx >= reader->dictionary_count) {
                                status = CARQUET_ERROR_DECODE;
                                break;
                            }
                            uint32_t offset = reader->dictionary_offsets[idx];
                            const uint8_t* dict_ptr = reader->dictionary_data + offset;
                            uint32_t len = carquet_read_u32_le(dict_ptr);
                            out[i].data = (uint8_t*)(dict_ptr + 4);
                            out[i].length = (int32_t)len;
                        }
                    } else {
                        /* Fallback scan (unreachable for current readers); bounds-checked
                         * and size_t-widened as defense-in-depth against a 32-bit
                         * `4 + len` wrap walking past the dictionary buffer. */
                        const uint8_t* dict_end =
                            reader->dictionary_data + reader->dictionary_size;
                        for (int32_t i = 0; i < encoded_count; i++) {
                            int32_t idx = (int32_t)indices[i];
                            if (idx < 0 || idx >= reader->dictionary_count) {
                                status = CARQUET_ERROR_DECODE;
                                break;
                            }
                            const uint8_t* dict_ptr = reader->dictionary_data;
                            bool oob = false;
                            for (int32_t j = 0; j <= idx; j++) {
                                if ((size_t)(dict_end - dict_ptr) < 4) {
                                    oob = true;
                                    break;
                                }
                                uint32_t len = carquet_read_u32_le(dict_ptr);
                                if ((size_t)(dict_end - dict_ptr) - 4 < (size_t)len) {
                                    oob = true;
                                    break;
                                }
                                if (j == idx) {
                                    out[i].data = (uint8_t*)(dict_ptr + 4);
                                    out[i].length = (int32_t)len;
                                    break;
                                }
                                dict_ptr += (size_t)4 + (size_t)len;
                            }
                            if (oob) {
                                status = CARQUET_ERROR_DECODE;
                                break;
                            }
                        }
                    }
                } else {
                    switch (reader->type) {
                        case CARQUET_PHYSICAL_INT32:
                            if (!carquet_dispatch_checked_gather_i32(
                                    (const int32_t*)reader->dictionary_data,
                                    reader->dictionary_count,
                                    indices, encoded_count, (int32_t*)values)) {
                                CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE,
                                                  "Dictionary index out of bounds");
                                return CARQUET_ERROR_DECODE;
                            }
                            break;
                        case CARQUET_PHYSICAL_INT64:
                            if (!carquet_dispatch_checked_gather_i64(
                                    (const int64_t*)reader->dictionary_data,
                                    reader->dictionary_count,
                                    indices, encoded_count, (int64_t*)values)) {
                                CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE,
                                                  "Dictionary index out of bounds");
                                return CARQUET_ERROR_DECODE;
                            }
                            break;
                        case CARQUET_PHYSICAL_FLOAT:
                            if (!carquet_dispatch_checked_gather_float(
                                    (const float*)reader->dictionary_data,
                                    reader->dictionary_count,
                                    indices, encoded_count, (float*)values)) {
                                CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE,
                                                  "Dictionary index out of bounds");
                                return CARQUET_ERROR_DECODE;
                            }
                            break;
                        case CARQUET_PHYSICAL_DOUBLE:
                            if (!carquet_dispatch_checked_gather_double(
                                    (const double*)reader->dictionary_data,
                                    reader->dictionary_count,
                                    indices, encoded_count, (double*)values)) {
                                CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE,
                                                  "Dictionary index out of bounds");
                                return CARQUET_ERROR_DECODE;
                            }
                            break;
                        case CARQUET_PHYSICAL_INT96:
                        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
                            {
                                size_t value_size = (reader->type == CARQUET_PHYSICAL_INT96)
                                    ? 12 : (size_t)reader->type_length;
                                uint8_t* out = (uint8_t*)values;
                                bool ok;
#if defined(CARQUET_ARCH_ARM) && defined(CARQUET_ENABLE_NEON) && \
    (defined(__ARM_NEON) || defined(__ARM_NEON__))
                                ok = gather_fixed_dictionary_values_neon(
                                    reader->dictionary_data,
                                    reader->dictionary_count,
                                    indices, encoded_count,
                                    value_size, out);
#else
                                ok = true;
                                for (int32_t i = 0; i < encoded_count; i++) {
                                    uint32_t idx = indices[i];
                                    if (idx >= (uint32_t)reader->dictionary_count) {
                                        ok = false;
                                        break;
                                    }
                                    memcpy(out + (size_t)i * value_size,
                                           reader->dictionary_data + (size_t)idx * value_size,
                                           value_size);
                                }
#endif
                                if (!ok) {
                                    CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE,
                                                      "Dictionary index out of bounds");
                                    return CARQUET_ERROR_DECODE;
                                }
                            }
                            break;
                        default:
                            break;
                    }
                }
            }
            break;

        default:
            {
                /* Phase 3 encodings (shared with V1 path). */
                bool handled = false;
                bool needs_page_retain = false;
                status = decode_phase3_values(
                    reader, header->encoding, ptr, remaining,
                    values, non_null_count,
                    &handled, &needs_page_retain, error);
                if (!handled) {
                    CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ENCODING,
                        "Unsupported encoding: %d", header->encoding);
                    return CARQUET_ERROR_INVALID_ENCODING;
                }
            }
            break;
    }

    if (status != CARQUET_OK) {
        CARQUET_SET_ERROR(error, status, "Failed to decode values");
        return status;
    }

    *values_read = num_values;
    return CARQUET_OK;
}

/* ============================================================================
 * Helper: Get value size for a physical type
 * ============================================================================
 */

static size_t get_value_size(carquet_physical_type_t type, int32_t type_length) {
    switch (type) {
        case CARQUET_PHYSICAL_BOOLEAN:
            return 1;
        case CARQUET_PHYSICAL_INT32:
        case CARQUET_PHYSICAL_FLOAT:
            return 4;
        case CARQUET_PHYSICAL_INT64:
        case CARQUET_PHYSICAL_DOUBLE:
            return 8;
        case CARQUET_PHYSICAL_INT96:
            return 12;
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
            return type_length;
        case CARQUET_PHYSICAL_BYTE_ARRAY:
            return sizeof(carquet_byte_array_t);
        default:
            return 0;
    }
}

static int32_t count_present_levels(
    const int16_t* def_levels,
    int32_t count,
    int16_t max_def_level) {
    return (int32_t)carquet_dispatch_count_non_nulls(def_levels, count, max_def_level);
}

/* Absolute file offset of the payload currently being decompressed, used only
 * to make decompression failures locatable in the source file. */
static int64_t current_page_payload_offset(const carquet_column_reader_t* reader) {
    return reader->data_start_offset + reader->current_page;
}

static carquet_status_t prepare_data_page_payload(
    carquet_column_reader_t* reader,
    const parquet_column_metadata_t* col_meta,
    const parquet_page_header_t* page_header,
    const uint8_t* compressed,
    const uint8_t** page_data,
    size_t* page_size,
    bool* used_decompress_buffer,
    carquet_error_t* error) {

    bool is_v2 = (page_header->type == CARQUET_PAGE_DATA_V2);
    carquet_status_t status;

    *page_data = NULL;
    *page_size = 0;
    *used_decompress_buffer = false;

    status = validate_page_payload_size(page_header, false, error);
    if (status != CARQUET_OK) {
        return status;
    }

    if (is_v2) {
        const parquet_data_page_header_v2_t* v2h = &page_header->data_page_header_v2;
        size_t levels_size = (size_t)v2h->repetition_levels_byte_length +
                             (size_t)v2h->definition_levels_byte_length;
        size_t compressed_data_size = (size_t)page_header->compressed_page_size - levels_size;
        bool data_is_compressed = v2h->is_compressed &&
                                  col_meta->codec != CARQUET_COMPRESSION_UNCOMPRESSED;

        if (data_is_compressed) {
            size_t uncompressed_data_size =
                (size_t)page_header->uncompressed_page_size - levels_size;
            size_t total_needed = levels_size + uncompressed_data_size;

            status = ensure_decompress_capacity(
                reader, total_needed, "Failed to allocate V2 decompress buffer", error);
            if (status != CARQUET_OK) {
                return status;
            }

            if (levels_size > 0) {
                memcpy(reader->decompress_buffer, compressed, levels_size);
            }

            size_t decompressed_data_size = 0;
            char diag[CARQUET_DECOMPRESS_DIAG_SIZE];
            status = carquet_decompress_page_diag(col_meta->codec,
                compressed + levels_size, compressed_data_size,
                reader->decompress_buffer + levels_size, uncompressed_data_size,
                &decompressed_data_size, diag, sizeof(diag));
            if (status != CARQUET_OK) {
                /* The diagnostic goes first: carquet_error_t messages are
                 * capped at CARQUET_ERROR_MESSAGE_MAX, so the detail that
                 * identifies the failure must not be what gets truncated.
                 * %llu, not %zu — the MinGW/msvcrt printf has no z modifier. */
                CARQUET_SET_ERROR(error, status,
                    "Failed to decompress V2 page at file offset %lld%s%s "
                    "(codec=%d, levels=%llu, compressed=%llu, uncompressed=%llu)",
                    (long long)current_page_payload_offset(reader),
                    diag[0] ? ": " : "", diag,
                    (int)col_meta->codec,
                    (unsigned long long)levels_size,
                    (unsigned long long)compressed_data_size,
                    (unsigned long long)uncompressed_data_size);
                return status;
            }

            *page_data = reader->decompress_buffer;
            *page_size = levels_size + decompressed_data_size;
            *used_decompress_buffer = true;
            return CARQUET_OK;
        }

        *page_data = compressed;
        *page_size = (size_t)page_header->compressed_page_size;
        return CARQUET_OK;
    }

    if (col_meta->codec == CARQUET_COMPRESSION_UNCOMPRESSED) {
        *page_data = compressed;
        *page_size = (size_t)page_header->compressed_page_size;
        return CARQUET_OK;
    }

    status = ensure_decompress_capacity(
        reader, (size_t)page_header->uncompressed_page_size,
        "Failed to allocate decompress buffer", error);
    if (status != CARQUET_OK) {
        return status;
    }

    char diag[CARQUET_DECOMPRESS_DIAG_SIZE];
    status = carquet_decompress_page_diag(col_meta->codec,
        compressed, (size_t)page_header->compressed_page_size,
        reader->decompress_buffer, (size_t)page_header->uncompressed_page_size,
        page_size, diag, sizeof(diag));
    if (status != CARQUET_OK) {
        CARQUET_SET_ERROR(error, status,
            "Failed to decompress data page at file offset %lld%s%s "
            "(codec=%d, compressed=%d, uncompressed=%d)",
            (long long)current_page_payload_offset(reader),
            diag[0] ? ": " : "", diag,
            (int)col_meta->codec,
            page_header->compressed_page_size, page_header->uncompressed_page_size);
        return status;
    }

    *page_data = reader->decompress_buffer;
    *used_decompress_buffer = true;
    return CARQUET_OK;
}

/* ============================================================================
 * Helper: Load dictionary page (mmap path)
 * ============================================================================
 */

static carquet_status_t load_dictionary_page_mmap(
    carquet_column_reader_t* reader,
    carquet_error_t* error) {

    carquet_reader_t* file_reader = reader->file_reader;
    const uint8_t* mmap_data = file_reader->mmap_data;
    const parquet_column_metadata_t* col_meta = reader->col_meta;

    /* Parse page header directly from mmap */
    int64_t dict_offset = col_meta->dictionary_page_offset;
    if (dict_offset < 0 || (size_t)dict_offset >= file_reader->file_size) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE, "Dictionary page offset out of range");
        return CARQUET_ERROR_INVALID_PAGE;
    }
    const uint8_t* header_ptr = mmap_data + dict_offset;
    /* Page headers have no spec size limit (large statistics can exceed any
     * fixed guess). The whole file is mapped, so let the thrift parser read
     * the full remaining span; it stops at the struct end. */
    size_t max_header = file_reader->file_size - (size_t)dict_offset;

    parquet_page_header_t page_header;
    size_t header_size;
    carquet_status_t status = parquet_parse_page_header(
        header_ptr, max_header, &page_header, &header_size, error);
    if (status != CARQUET_OK) {
        return status;
    }

    if (page_header.type != CARQUET_PAGE_DICTIONARY) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE, "Expected dictionary page");
        return CARQUET_ERROR_INVALID_PAGE;
    }

    status = validate_page_payload_size(&page_header, true, error);
    if (status != CARQUET_OK) {
        return status;
    }
    status = validate_page_payload_span(file_reader, dict_offset, header_size,
                                        page_header.compressed_page_size, error);
    if (status != CARQUET_OK) {
        return status;
    }

    /* Get pointer to compressed data */
    const uint8_t* compressed = header_ptr + header_size;

    /* Verify CRC32 if present */
    if (page_header.has_crc && file_reader->options.verify_checksums) {
        uint32_t computed_crc = carquet_crc32(compressed, page_header.compressed_page_size);
        uint32_t expected_crc = (uint32_t)page_header.crc;
        if (computed_crc != expected_crc) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_CRC_MISMATCH,
                "Dictionary page CRC mismatch: expected 0x%08X, got 0x%08X",
                expected_crc, computed_crc);
            return CARQUET_ERROR_CRC_MISMATCH;
        }
    }

    /* Process dictionary data */
    const uint8_t* page_data;
    size_t page_size;
    uint8_t* decompressed = NULL;

    if (col_meta->codec == CARQUET_COMPRESSION_UNCOMPRESSED) {
        /* Zero-copy: point directly to mmap data */
        page_data = compressed;
        page_size = page_header.compressed_page_size;
    } else {
        /* Must decompress */
        decompressed = carquet_mem_malloc(page_header.uncompressed_page_size);
        if (!decompressed) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate decompress buffer");
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }

        char diag[CARQUET_DECOMPRESS_DIAG_SIZE];
        status = carquet_decompress_page_diag(col_meta->codec,
            compressed, page_header.compressed_page_size,
            decompressed, page_header.uncompressed_page_size, &page_size,
            diag, sizeof(diag));

        if (status != CARQUET_OK) {
            carquet_mem_free(decompressed);
            CARQUET_SET_ERROR(error, status,
                "Failed to decompress dictionary page at file offset %lld%s%s "
                "(codec=%d, compressed=%d, uncompressed=%d)",
                (long long)dict_offset,
                diag[0] ? ": " : "", diag,
                (int)col_meta->codec,
                page_header.compressed_page_size, page_header.uncompressed_page_size);
            return status;
        }
        page_data = decompressed;
    }

    /* Parse dictionary */
    status = carquet_read_dictionary_page(
        reader, (uint8_t*)page_data, page_size,
        &page_header.dictionary_page_header,
        col_meta->codec == CARQUET_COMPRESSION_UNCOMPRESSED
            ? CARQUET_DATA_VIEW
            : CARQUET_DATA_OWNED,
        error);

    if (status == CARQUET_OK) {
        int64_t dict_end;
        if (!checked_add_i64(dict_offset, (int64_t)header_size, &dict_end) ||
            !checked_add_i64(dict_end, (int64_t)page_header.compressed_page_size,
                             &dict_end)) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE, "Dictionary page offset overflow");
            return CARQUET_ERROR_INVALID_PAGE;
        }
        reader->data_start_offset =
            carquet_resolve_data_start_offset(col_meta->data_page_offset, dict_end);
    }

    if (col_meta->codec != CARQUET_COMPRESSION_UNCOMPRESSED && status != CARQUET_OK) {
        carquet_mem_free(decompressed);
    }
    return status;
}

/* Parquet page headers have no spec size limit: large column statistics
 * (BYTE_ARRAY min/max) can push a header past any fixed guess. Read from the
 * prebuffer/file growing the window whenever the thrift parser reports
 * truncation, until it parses or a hard ceiling is reached. */
static carquet_status_t read_and_parse_page_header_fread(
    carquet_reader_t* file_reader,
    int64_t offset,
    parquet_page_header_t* page_header,
    size_t* header_size,
    carquet_error_t* error) {

    uint8_t stackbuf[256];
    uint8_t* buf = stackbuf;
    size_t cap = sizeof(stackbuf);

    for (;;) {
        size_t n = prebuf_read_at(file_reader, offset, buf, cap);
        if (n < 8) {
            if (buf != stackbuf) carquet_mem_free(buf);
            CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_READ,
                              "Failed to read page header");
            return CARQUET_ERROR_FILE_READ;
        }

        carquet_status_t status = parquet_parse_page_header(
            buf, n, page_header, header_size, error);
        if (status == CARQUET_OK) {
            if (buf != stackbuf) carquet_mem_free(buf);
            return CARQUET_OK;
        }

        /* Only a truncated parse on a completely filled buffer means the
         * header may extend further; anything else is a real error or we
         * already have all available bytes. */
        if (status != CARQUET_ERROR_THRIFT_TRUNCATED || n < cap ||
            cap >= CARQUET_MAX_PAGE_PAYLOAD_SIZE) {
            if (buf != stackbuf) carquet_mem_free(buf);
            return status;
        }

        size_t new_cap = cap * 4;
        if (new_cap > CARQUET_MAX_PAGE_PAYLOAD_SIZE) {
            new_cap = CARQUET_MAX_PAGE_PAYLOAD_SIZE;
        }
        uint8_t* nb = carquet_mem_malloc(new_cap);
        if (!nb) {
            if (buf != stackbuf) carquet_mem_free(buf);
            CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY,
                              "Failed to allocate page header buffer");
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        if (buf != stackbuf) carquet_mem_free(buf);
        buf = nb;
        cap = new_cap;
    }
}

/* ============================================================================
 * Helper: Load dictionary page (fread path)
 * ============================================================================
 */

static carquet_status_t load_dictionary_page_fread(
    carquet_column_reader_t* reader,
    carquet_error_t* error) {

    carquet_reader_t* file_reader = reader->file_reader;
    const parquet_column_metadata_t* col_meta = reader->col_meta;
    int64_t dict_offset = col_meta->dictionary_page_offset;

    /* Read page header (from prebuffer cache or file) */
    parquet_page_header_t page_header;
    size_t header_size;
    carquet_status_t status = read_and_parse_page_header_fread(
        file_reader, dict_offset, &page_header, &header_size, error);
    if (status != CARQUET_OK) {
        return status;
    }

    if (page_header.type != CARQUET_PAGE_DICTIONARY) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE, "Expected dictionary page");
        return CARQUET_ERROR_INVALID_PAGE;
    }

    status = validate_page_payload_size(&page_header, true, error);
    if (status != CARQUET_OK) {
        return status;
    }
    status = validate_page_payload_span(file_reader, dict_offset, header_size,
                                        page_header.compressed_page_size, error);
    if (status != CARQUET_OK) {
        return status;
    }

    /* Read compressed data (from prebuffer cache or file) */
    uint8_t* compressed = carquet_mem_malloc(page_header.compressed_page_size);
    if (!compressed) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate compressed buffer");
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    int64_t dict_data_offset;
    if (!checked_add_i64(dict_offset, (int64_t)header_size, &dict_data_offset)) {
        carquet_mem_free(compressed);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE, "Dictionary page offset overflow");
        return CARQUET_ERROR_INVALID_PAGE;
    }
    if (prebuf_read_at(file_reader, dict_data_offset,
                       compressed, page_header.compressed_page_size) !=
        (size_t)page_header.compressed_page_size) {
        carquet_mem_free(compressed);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_READ, "Failed to read dictionary data");
        return CARQUET_ERROR_FILE_READ;
    }

    /* Verify CRC32 if present */
    if (page_header.has_crc && file_reader->options.verify_checksums) {
        uint32_t computed_crc = carquet_crc32(compressed, page_header.compressed_page_size);
        uint32_t expected_crc = (uint32_t)page_header.crc;
        if (computed_crc != expected_crc) {
            carquet_mem_free(compressed);
            CARQUET_SET_ERROR(error, CARQUET_ERROR_CRC_MISMATCH,
                "Dictionary page CRC mismatch: expected 0x%08X, got 0x%08X",
                expected_crc, computed_crc);
            return CARQUET_ERROR_CRC_MISMATCH;
        }
    }

    /* Decompress if needed */
    uint8_t* page_data;
    size_t page_size;

    if (col_meta->codec == CARQUET_COMPRESSION_UNCOMPRESSED) {
        page_data = compressed;
        page_size = page_header.compressed_page_size;
    } else {
        page_data = carquet_mem_malloc(page_header.uncompressed_page_size);
        if (!page_data) {
            carquet_mem_free(compressed);
            CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate decompress buffer");
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }

        char diag[CARQUET_DECOMPRESS_DIAG_SIZE];
        status = carquet_decompress_page_diag(col_meta->codec,
            compressed, page_header.compressed_page_size,
            page_data, page_header.uncompressed_page_size, &page_size,
            diag, sizeof(diag));
        carquet_mem_free(compressed);

        if (status != CARQUET_OK) {
            carquet_mem_free(page_data);
            CARQUET_SET_ERROR(error, status,
                "Failed to decompress dictionary page at file offset %lld%s%s "
                "(codec=%d, compressed=%d, uncompressed=%d)",
                (long long)col_meta->dictionary_page_offset,
                diag[0] ? ": " : "", diag,
                (int)col_meta->codec,
                page_header.compressed_page_size, page_header.uncompressed_page_size);
            return status;
        }
    }

    /* Parse dictionary */
    status = carquet_read_dictionary_page(
        reader, page_data, page_size,
        &page_header.dictionary_page_header,
        CARQUET_DATA_OWNED, error);

    if (status == CARQUET_OK) {
        int64_t dict_end;
        if (!checked_add_i64(col_meta->dictionary_page_offset,
                             (int64_t)header_size, &dict_end) ||
            !checked_add_i64(dict_end, (int64_t)page_header.compressed_page_size,
                             &dict_end)) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE, "Dictionary page offset overflow");
            return CARQUET_ERROR_INVALID_PAGE;
        }
        reader->data_start_offset =
            carquet_resolve_data_start_offset(col_meta->data_page_offset, dict_end);
    }

    if (status != CARQUET_OK) {
        if (page_data != compressed) {
            carquet_mem_free(page_data);
        } else {
            carquet_mem_free(compressed);
        }
    }

    return status;
}

/* ============================================================================
 * Helper: Load and decode a new page (mmap path with zero-copy support)
 * ============================================================================
 */

static carquet_status_t load_next_page_mmap(
    carquet_column_reader_t* reader,
    carquet_error_t* error) {

    carquet_reader_t* file_reader = reader->file_reader;
    const uint8_t* mmap_data = file_reader->mmap_data;
    const parquet_column_metadata_t* col_meta = reader->col_meta;

    /* Load dictionary if needed (may update data_start_offset) */
    if (col_meta->has_dictionary_page_offset && !reader->has_dictionary) {
        carquet_status_t status = load_dictionary_page_mmap(reader, error);
        if (status != CARQUET_OK) {
            return status;
        }
    }

    /* Parse page header directly from mmap */
    int64_t page_offset;
    if (!checked_add_i64(reader->data_start_offset, reader->current_page, &page_offset)) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE, "Data page offset overflow");
        return CARQUET_ERROR_INVALID_PAGE;
    }
    if (page_offset < 0 || (size_t)page_offset >= file_reader->file_size) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE, "Data page offset out of range");
        return CARQUET_ERROR_INVALID_PAGE;
    }
    const uint8_t* header_ptr = mmap_data + page_offset;
    /* Page headers have no spec size limit (large statistics can exceed any
     * fixed guess). The whole file is mapped, so let the thrift parser read
     * the full remaining span; it stops at the struct end. */
    size_t max_hdr = file_reader->file_size - (size_t)page_offset;

    parquet_page_header_t page_header;
    size_t header_size;
    carquet_status_t status = parquet_parse_page_header(
        header_ptr, max_hdr, &page_header, &header_size, error);
    if (status != CARQUET_OK) {
        return status;
    }

    if (page_header.type != CARQUET_PAGE_DATA && page_header.type != CARQUET_PAGE_DATA_V2) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE, "Expected data page");
        return CARQUET_ERROR_INVALID_PAGE;
    }

    status = validate_page_payload_size(&page_header, false, error);
    if (status != CARQUET_OK) {
        return status;
    }
    status = validate_page_payload_span(file_reader, page_offset, header_size,
                                        page_header.compressed_page_size, error);
    if (status != CARQUET_OK) {
        return status;
    }

    /* Get pointer to page data in mmap */
    const uint8_t* page_data_ptr = header_ptr + header_size;

    /* Verify CRC32 if present */
    if (page_header.has_crc && file_reader->options.verify_checksums) {
        uint32_t computed_crc = carquet_crc32(page_data_ptr, page_header.compressed_page_size);
        uint32_t expected_crc = (uint32_t)page_header.crc;
        if (computed_crc != expected_crc) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_CRC_MISMATCH,
                "Page CRC mismatch: expected 0x%08X, got 0x%08X at offset %lld",
                expected_crc, computed_crc, (long long)page_offset);
            return CARQUET_ERROR_CRC_MISMATCH;
        }
    }

    /* Extract num_values and encoding from the correct header union member */
    bool is_v2 = (page_header.type == CARQUET_PAGE_DATA_V2);
    int32_t num_values = is_v2 ? page_header.data_page_header_v2.num_values
                               : page_header.data_page_header.num_values;
    carquet_encoding_t page_encoding = is_v2 ? page_header.data_page_header_v2.encoding
                                             : page_header.data_page_header.encoding;
    /* In dictionary-preserving mode, decoded_values contains uint32_t indices
     * rather than materialized physical values. The batch reader sizes its
     * destination buffer accordingly; keep the page-copy path symmetric or a
     * preserved INT64/DOUBLE page will overrun a uint32_t output buffer. */
    size_t value_size = reader->preserve_dictionary
        ? sizeof(uint32_t)
        : get_value_size(reader->type, reader->type_length);

    /* Check if zero-copy is possible (V1 only — V2 has levels interleaved) */
    bool zero_copy_eligible = !is_v2 && carquet_page_is_zero_copy_eligible(
        col_meta->codec, page_encoding, reader->type);

    /* Additional constraint: no definition/repetition levels for zero-copy
     * (levels require RLE decoding which modifies data layout) */
    bool has_levels = (reader->max_def_level > 0 || reader->max_rep_level > 0);

    if (zero_copy_eligible && !has_levels) {
        /* ====== ZERO-COPY PATH ====== */

        /* Validate that the page payload actually holds num_values fixed-width
         * values before viewing it directly. Without this a crafted header that
         * declares more values than the payload contains causes an
         * out-of-bounds read when the batch reader memcpy's num_values *
         * value_size bytes out of decoded_values (which points straight into
         * the mapped file). Mirrors the same check on the buffered PLAIN
         * view-directly paths below. The zero-copy codec is always
         * uncompressed, so compressed_page_size is the exact payload length. */
        size_t required_bytes = 0;
        if (!checked_mul_size(value_size, (size_t)num_values, &required_bytes)) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE, "Page payload size overflow");
            return CARQUET_ERROR_DECODE;
        }
        if ((size_t)page_header.compressed_page_size < required_bytes) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE, "Truncated PLAIN page payload");
            return CARQUET_ERROR_DECODE;
        }

        /* Free previous owned buffer if any */
        if (reader->decoded_ownership == CARQUET_DATA_OWNED) {
            carquet_mem_free(reader->decoded_values);
        }

        /* Point directly to mmap data - no copy! */
        reader->decoded_values = (uint8_t*)page_data_ptr;
        reader->decoded_ownership = CARQUET_DATA_VIEW;

        /* Zero-copy path only triggers when max_def/rep == 0 (REQUIRED columns).
         * Level buffers are unused by callers for REQUIRED columns, so set to NULL
         * to avoid unnecessary allocation and memset overhead. */
        if (reader->decoded_def_levels) {
            carquet_mem_free(reader->decoded_def_levels);
            reader->decoded_def_levels = NULL;
        }
        if (reader->decoded_rep_levels) {
            carquet_mem_free(reader->decoded_rep_levels);
            reader->decoded_rep_levels = NULL;
        }
        reader->decoded_capacity = 0;

        reader->page_loaded = true;
        reader->page_num_values = num_values;
        reader->page_values_read = 0;
        reader->page_header_size = (int32_t)header_size;
        reader->page_compressed_size = page_header.compressed_page_size;

        return CARQUET_OK;
    }

    /* ====== STANDARD PATH (with decompression/decoding) ====== */

    const uint8_t* page_data;
    size_t page_size;
    bool used_decompress_buffer = false;

    status = prepare_data_page_payload(reader, col_meta, &page_header, page_data_ptr,
                                       &page_data, &page_size, &used_decompress_buffer,
                                       error);
    if (status != CARQUET_OK) {
        return status;
    }

    if (!is_v2 && page_values_can_be_viewed_directly(reader, page_encoding)) {
        size_t required_bytes = 0;
        if (!checked_mul_size(value_size, (size_t)num_values, &required_bytes)) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE, "Page payload size overflow");
            return CARQUET_ERROR_DECODE;
        }
        if (page_size < required_bytes) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE, "Truncated PLAIN page payload");
            return CARQUET_ERROR_DECODE;
        }

        if (reader->decoded_ownership == CARQUET_DATA_OWNED) {
            carquet_mem_free(reader->decoded_values);
        }
        reader->decoded_values = (uint8_t*)page_data;
        reader->decoded_ownership = CARQUET_DATA_VIEW;
        reader->decoded_capacity = 0;
        release_decoded_level_buffers(reader);

        reader->page_loaded = true;
        reader->page_num_values = num_values;
        reader->page_values_read = 0;
        reader->page_header_size = (int32_t)header_size;
        reader->page_compressed_size = page_header.compressed_page_size;

        return CARQUET_OK;
    }

    status = ensure_decoded_page_buffers(reader, num_values, value_size, error);
    if (status != CARQUET_OK) {
        return status;
    }

    /* Decode the page */
    int64_t decoded_count;
    if (is_v2) {
        status = carquet_read_data_page_v2(
            reader, page_data, page_size,
            &page_header.data_page_header_v2,
            reader->decoded_values, num_values,
            reader->decoded_def_levels, reader->decoded_rep_levels,
            &decoded_count, error);
    } else {
        status = carquet_read_data_page_v1(
            reader, page_data, page_size,
            &page_header.data_page_header,
            reader->decoded_values, num_values,
            reader->decoded_def_levels, reader->decoded_rep_levels,
            &decoded_count, error);
    }

    if (status != CARQUET_OK) {
        return status;
    }

    /* For BYTE_ARRAY PLAIN columns with compressed data, retain a copy of the
     * decompressed buffer since carquet_byte_array_t.data pointers reference it.
     * The decompression buffer is reused across pages, AND a single batch read
     * may span multiple pages, so every page must be retained until the batch
     * is consumed (list is flushed on row-group reset / reader close). */
    if (used_decompress_buffer && reader->type == CARQUET_PHYSICAL_BYTE_ARRAY &&
        (page_encoding == CARQUET_ENCODING_PLAIN ||
         page_encoding == CARQUET_ENCODING_DELTA_LENGTH_BYTE_ARRAY)) {
        uint8_t* retained = carquet_column_retain_page(reader, page_data, page_size);
        if (!retained) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to retain page data");
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        /* Fixup BYTE_ARRAY pointers to reference the retained copy */
        ptrdiff_t offset = retained - page_data;
        carquet_byte_array_t* ba = (carquet_byte_array_t*)reader->decoded_values;
        for (int64_t i = 0; i < decoded_count; i++) {
            if (ba[i].data) {
                ba[i].data = ba[i].data + offset;
            }
        }
    }

    reader->page_loaded = true;
    reader->page_num_values = (int32_t)decoded_count;
    reader->page_values_read = 0;
    reader->page_header_size = (int32_t)header_size;
    reader->page_compressed_size = page_header.compressed_page_size;

    return CARQUET_OK;
}

/* ============================================================================
 * Helper: Load and decode a new page (fread path)
 * ============================================================================
 */

static carquet_status_t load_next_page_fread(
    carquet_column_reader_t* reader,
    carquet_error_t* error) {

    carquet_reader_t* file_reader = reader->file_reader;
    const parquet_column_metadata_t* col_meta = reader->col_meta;

    /* Load dictionary if needed (may update data_start_offset) */
    if (col_meta->has_dictionary_page_offset && !reader->has_dictionary) {
        carquet_status_t status = load_dictionary_page_fread(reader, error);
        if (status != CARQUET_OK) {
            return status;
        }
    }

    /* Read page header (from prebuffer cache or file) */
    int64_t data_offset = reader->data_start_offset;
    int64_t page_file_offset;
    if (!checked_add_i64(data_offset, reader->current_page, &page_file_offset)) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE, "Data page offset overflow");
        return CARQUET_ERROR_INVALID_PAGE;
    }

    parquet_page_header_t page_header;
    size_t header_size;
    carquet_status_t status = read_and_parse_page_header_fread(
        file_reader, page_file_offset, &page_header, &header_size, error);
    if (status != CARQUET_OK) {
        return status;
    }

    if (page_header.type != CARQUET_PAGE_DATA && page_header.type != CARQUET_PAGE_DATA_V2) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE, "Expected data page");
        return CARQUET_ERROR_INVALID_PAGE;
    }

    status = validate_page_payload_size(&page_header, false, error);
    if (status != CARQUET_OK) {
        return status;
    }

    status = validate_page_payload_span(file_reader, page_file_offset, header_size,
                                        page_header.compressed_page_size, error);
    if (status != CARQUET_OK) {
        return status;
    }

    /* Read compressed page data into a reusable buffer (from prebuffer or file) */
    if ((size_t)page_header.compressed_page_size > reader->page_buffer_capacity) {
        uint8_t* new_buffer = carquet_mem_realloc(reader->page_buffer, (size_t)page_header.compressed_page_size);
        if (!new_buffer) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate page buffer");
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        reader->page_buffer = new_buffer;
        reader->page_buffer_capacity = (size_t)page_header.compressed_page_size;
    }
    reader->page_buffer_size = (size_t)page_header.compressed_page_size;
    uint8_t* compressed = reader->page_buffer;

    int64_t page_data_offset;
    if (!checked_add_i64(page_file_offset, (int64_t)header_size, &page_data_offset)) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE, "Data page offset overflow");
        return CARQUET_ERROR_INVALID_PAGE;
    }
    if (prebuf_read_at(file_reader, page_data_offset,
                       compressed, page_header.compressed_page_size) !=
        (size_t)page_header.compressed_page_size) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_READ, "Failed to read page data");
        return CARQUET_ERROR_FILE_READ;
    }

    /* Verify CRC32 if present */
    if (page_header.has_crc && file_reader->options.verify_checksums) {
        uint32_t computed_crc = carquet_crc32(compressed, page_header.compressed_page_size);
        uint32_t expected_crc = (uint32_t)page_header.crc;
        if (computed_crc != expected_crc) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_CRC_MISMATCH,
                "Page CRC mismatch: expected 0x%08X, got 0x%08X at offset %lld",
                expected_crc, computed_crc, (long long)(data_offset + reader->current_page));
            return CARQUET_ERROR_CRC_MISMATCH;
        }
    }

    /* Extract num_values and encoding from the correct header union member */
    bool is_v2 = (page_header.type == CARQUET_PAGE_DATA_V2);
    int32_t num_values = is_v2 ? page_header.data_page_header_v2.num_values
                               : page_header.data_page_header.num_values;
    carquet_encoding_t page_encoding = is_v2 ? page_header.data_page_header_v2.encoding
                                             : page_header.data_page_header.encoding;
    size_t value_size = reader->preserve_dictionary
        ? sizeof(uint32_t)
        : get_value_size(reader->type, reader->type_length);

    /* Decompress if needed. V2 pages keep level bytes uncompressed and only
     * decompress the data tail, so both fread and mmap paths share this helper. */
    const uint8_t* page_data;
    size_t page_size;
    bool used_decompress_buffer = false;
    status = prepare_data_page_payload(reader, col_meta, &page_header, compressed,
                                       &page_data, &page_size, &used_decompress_buffer,
                                       error);
    if (status != CARQUET_OK) {
        return status;
    }

    if (!is_v2 && page_values_can_be_viewed_directly(reader, page_encoding)) {
        size_t required_bytes = 0;
        if (!checked_mul_size(value_size, (size_t)num_values, &required_bytes)) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE, "Page payload size overflow");
            return CARQUET_ERROR_DECODE;
        }
        if (page_size < required_bytes) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_DECODE, "Truncated PLAIN page payload");
            return CARQUET_ERROR_DECODE;
        }

        if (reader->decoded_ownership == CARQUET_DATA_OWNED) {
            carquet_mem_free(reader->decoded_values);
        }
        reader->decoded_values = (uint8_t*)page_data;
        reader->decoded_ownership = CARQUET_DATA_VIEW;
        reader->decoded_capacity = 0;
        release_decoded_level_buffers(reader);

        reader->page_loaded = true;
        reader->page_num_values = num_values;
        reader->page_values_read = 0;
        reader->page_header_size = (int32_t)header_size;
        reader->page_compressed_size = page_header.compressed_page_size;
        return CARQUET_OK;
    }

    status = ensure_decoded_page_buffers(reader, num_values, value_size, error);
    if (status != CARQUET_OK) {
        return status;
    }

    /* Decode the entire page into our buffers */
    int64_t decoded_count;
    if (is_v2) {
        status = carquet_read_data_page_v2(
            reader, page_data, page_size,
            &page_header.data_page_header_v2,
            reader->decoded_values, num_values,
            reader->decoded_def_levels, reader->decoded_rep_levels,
            &decoded_count, error);
    } else {
        status = carquet_read_data_page_v1(
            reader, page_data, page_size,
            &page_header.data_page_header,
            reader->decoded_values, num_values,
            reader->decoded_def_levels, reader->decoded_rep_levels,
            &decoded_count, error);
    }

    if (status != CARQUET_OK) {
        return status;
    }

    /* For BYTE_ARRAY PLAIN columns, the decoded carquet_byte_array_t structs
     * have .data pointers into the page data buffer. Retain every page buffer
     * so those pointers stay valid across page boundaries within a batch;
     * the retention list is flushed on row-group reset / reader close. */
    bool retain = (reader->type == CARQUET_PHYSICAL_BYTE_ARRAY &&
                   (page_encoding == CARQUET_ENCODING_PLAIN ||
                    page_encoding == CARQUET_ENCODING_DELTA_LENGTH_BYTE_ARRAY));

    if (retain) {
        uint8_t* retained = carquet_column_retain_page(reader, page_data, page_size);
        if (!retained) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to retain page data");
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }

        ptrdiff_t offset = retained - page_data;
        carquet_byte_array_t* ba = (carquet_byte_array_t*)reader->decoded_values;
        for (int64_t i = 0; i < decoded_count; i++) {
            if (ba[i].data) {
                ba[i].data = ba[i].data + offset;
            }
        }
    }

    /* Update page tracking state */
    reader->page_loaded = true;
    reader->page_num_values = (int32_t)decoded_count;
    reader->page_values_read = 0;
    reader->page_header_size = (int32_t)header_size;
    reader->page_compressed_size = page_header.compressed_page_size;

    return CARQUET_OK;
}

/* ============================================================================
 * Helper: Load and decode a new page (dispatcher)
 * ============================================================================
 */

static carquet_status_t load_next_page(
    carquet_column_reader_t* reader,
    carquet_error_t* error) {

    carquet_reader_t* file_reader = reader->file_reader;

    /* Use mmap/buffer path if memory-mapped or buffer-based reader */
    if (file_reader->mmap_data != NULL) {
        return load_next_page_mmap(reader, error);
    }

    /* Fall back to fread path (requires valid file handle) */
    if (file_reader->file == NULL) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_STATE, "No data source available");
        return CARQUET_ERROR_INVALID_STATE;
    }
    return load_next_page_fread(reader, error);
}

/* ============================================================================
 * Page Loading Helper
 * ============================================================================
 */

carquet_status_t carquet_column_ensure_page_loaded(
    carquet_column_reader_t* reader,
    carquet_error_t* error) {

    if (!reader) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT, "NULL reader");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    if (!reader->page_loaded || reader->page_values_read >= reader->page_num_values) {
        if (reader->page_loaded) {
            reader->current_page += reader->page_header_size + reader->page_compressed_size;
            reader->page_loaded = false;
        }

        return load_next_page(reader, error);
    }

    return CARQUET_OK;
}

/* ============================================================================
 * Dictionary Pre-Load (shared between data-page seek and standard read)
 * ============================================================================
 */

carquet_status_t carquet_column_ensure_dictionary_loaded(
    carquet_column_reader_t* reader,
    carquet_error_t* error) {

    if (!reader || !reader->col_meta) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT, "NULL reader");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    if (!reader->col_meta->has_dictionary_page_offset || reader->has_dictionary) {
        return CARQUET_OK;
    }

    if (reader->file_reader->mmap_data != NULL) {
        return load_dictionary_page_mmap(reader, error);
    }
    if (reader->file_reader->file == NULL) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_STATE, "No data source");
        return CARQUET_ERROR_INVALID_STATE;
    }
    return load_dictionary_page_fread(reader, error);
}

/* ============================================================================
 * Page Seek (used by page-filter row-range iteration)
 * ============================================================================
 *
 * Repositions the reader so that the next page load decodes the data page
 * at absolute file offset `page_file_offset`. We walk the page headers
 * between the current position and the target, accumulating num_values so
 * values_remaining stays consistent — this is robust to OPTIONAL columns
 * where row count and value count diverge.
 *
 * Dictionary state is preserved across the seek; the dictionary is shared
 * across all data pages in the column chunk.
 */

static carquet_status_t parse_header_at_offset(
    carquet_reader_t* file_reader, int64_t offset,
    parquet_page_header_t* hdr, size_t* hdr_size,
    carquet_error_t* error) {

    if (offset < 0 || (size_t)offset >= file_reader->file_size) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE,
            "Page header offset out of file");
        return CARQUET_ERROR_INVALID_PAGE;
    }

    if (file_reader->mmap_data != NULL) {
        size_t max_hdr = file_reader->file_size - (size_t)offset;
        return parquet_parse_page_header(
            file_reader->mmap_data + offset, max_hdr, hdr, hdr_size, error);
    }
    if (file_reader->file == NULL) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_STATE, "No data source");
        return CARQUET_ERROR_INVALID_STATE;
    }
    return read_and_parse_page_header_fread(
        file_reader, offset, hdr, hdr_size, error);
}

carquet_status_t carquet_column_reader_seek_to_data_page(
    carquet_column_reader_t* reader,
    int64_t page_file_offset,
    int64_t values_before_page,
    carquet_error_t* error) {

    (void)values_before_page;  /* Computed internally via header walk. */

    if (!reader || !reader->col_meta) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT, "NULL reader");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    /* Ensure dictionary is loaded so data_start_offset is final. */
    carquet_status_t st = carquet_column_ensure_dictionary_loaded(reader, error);
    if (st != CARQUET_OK) return st;

    if (page_file_offset < reader->data_start_offset) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
            "Target page offset %lld before chunk data start %lld",
            (long long)page_file_offset,
            (long long)reader->data_start_offset);
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    int64_t target_offset = page_file_offset - reader->data_start_offset;
    int64_t chunk_total_values = reader->col_meta->num_values;

    /* By API contract, the caller has already consumed any prior batch
     * before seeking. Flush retained BYTE_ARRAY page buffers. */
    carquet_column_clear_retained_pages(reader);

    /* If a VIEW-owned decoded buffer is held (mmap), release the
     * reference; owned buffers are kept for reuse. */
    if (reader->decoded_ownership == CARQUET_DATA_VIEW) {
        reader->decoded_values = NULL;
        reader->decoded_capacity = 0;
        reader->decoded_ownership = CARQUET_DATA_OWNED;
    }

    /* Always walk from the chunk start so accumulated_values counts every
     * value in pages before target_offset. Walking only from the current
     * position would omit pages [0, current_page) on a forward seek (e.g. the
     * second of two ranges in a page filter), leaving values_remaining too
     * high; a later unbounded read would then run current_page past the chunk
     * end and parse the next column's bytes as pages (load_next_page bounds
     * the offset only by file size, not by the chunk extent). */
    /* Always walk from the chunk start so accumulated_values counts every
     * value in pages before target_offset. Walking only from the current
     * position would omit pages [0, current_page) on a forward seek (e.g. the
     * second of two ranges in a page filter), leaving values_remaining too
     * high; a later unbounded read would then run current_page past the chunk
     * end and parse the next column's bytes as pages (load_next_page bounds
     * the offset only by file size, not by the chunk extent). */
    int64_t walk_pos = 0;
    int64_t accumulated_values = 0;

    /* Clear any mid-page state; values_remaining is recomputed from the walk
     * below, so the previously-decremented per-page counts no longer apply. */
    if (reader->page_loaded) {
        reader->page_loaded = false;
        reader->page_num_values = 0;
        reader->page_values_read = 0;
        reader->page_header_size = 0;
        reader->page_compressed_size = 0;
    }

    /* Walk page headers from walk_pos up to (but not including)
     * target_offset, summing each page's num_values into the
     * accumulator. */
    while (walk_pos < target_offset) {
        int64_t abs_off = reader->data_start_offset + walk_pos;
        parquet_page_header_t hdr;
        size_t hdr_size = 0;
        st = parse_header_at_offset(reader->file_reader, abs_off,
                                    &hdr, &hdr_size, error);
        if (st != CARQUET_OK) return st;

        int32_t num_values = 0;
        if (hdr.type == CARQUET_PAGE_DATA) {
            num_values = hdr.data_page_header.num_values;
        } else if (hdr.type == CARQUET_PAGE_DATA_V2) {
            num_values = hdr.data_page_header_v2.num_values;
        }
        /* Dictionary pages between data pages are not expected — we
         * arrive here only for in-chunk data pages. */
        if (num_values < 0) num_values = 0;
        accumulated_values += num_values;

        int64_t advance = (int64_t)hdr_size + (int64_t)hdr.compressed_page_size;
        if (advance <= 0) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_PAGE,
                "Non-positive page advance during seek walk");
            return CARQUET_ERROR_INVALID_PAGE;
        }
        walk_pos += advance;
    }

    if (walk_pos != target_offset) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
            "Target page offset %lld does not align with a page boundary "
            "(walked past to %lld)",
            (long long)target_offset, (long long)walk_pos);
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    reader->current_page = target_offset;
    reader->values_remaining = chunk_total_values - accumulated_values;
    if (reader->values_remaining < 0) reader->values_remaining = 0;
    return CARQUET_OK;
}

/* ============================================================================
 * Page Reading Entry Point
 * ============================================================================
 */

carquet_status_t carquet_read_next_page(
    carquet_column_reader_t* reader,
    void* values,
    int64_t max_values,
    int16_t* def_levels,
    int16_t* rep_levels,
    int64_t* values_read,
    carquet_error_t* error) {

    if (!reader || !values || !values_read) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT, "NULL argument");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    {
        carquet_status_t status = carquet_column_ensure_page_loaded(reader, error);
        if (status != CARQUET_OK) {
            return status;
        }
    }

    if (max_values == 0) {
        *values_read = 0;
        return CARQUET_OK;
    }

    /* Calculate how many values to return from the current page */
    int32_t available = reader->page_num_values - reader->page_values_read;
    int32_t to_copy = (int32_t)max_values;
    if (to_copy > available) {
        to_copy = available;
    }
    if (to_copy <= 0) {
        *values_read = 0;
        return CARQUET_OK;
    }

    /* Copy values from decoded buffers. Optional columns store a dense stream
     * of present values; definition levels preserve the logical row shape. */
    size_t value_size = reader->preserve_dictionary
        ? sizeof(uint32_t)
        : get_value_size(reader->type, reader->type_length);
    size_t offset = (size_t)reader->page_values_read * value_size;
    int32_t values_to_copy = to_copy;

    if (reader->decoded_def_levels && reader->max_def_level > 0) {
        /* Dense (present-value) offset of the first unread level. Recounting
         * [0, page_values_read) on every call made a page read in small
         * batches quadratic; keep a cursor and only recount when a skip
         * moved page_values_read behind our back. */
        int32_t dense_start;
        if (reader->page_values_read == 0) {
            dense_start = 0;
        } else if (reader->dense_cursor_pos == reader->page_values_read) {
            dense_start = reader->dense_cursor_val;
        } else {
            dense_start = count_present_levels(
                reader->decoded_def_levels,
                reader->page_values_read,
                reader->max_def_level);
        }
        values_to_copy = count_present_levels(
            reader->decoded_def_levels + reader->page_values_read,
            to_copy,
            reader->max_def_level);
        offset = (size_t)dense_start * value_size;
        reader->dense_cursor_pos = reader->page_values_read + to_copy;
        reader->dense_cursor_val = dense_start + values_to_copy;
    }

    memcpy(values, (uint8_t*)reader->decoded_values + offset,
           (size_t)values_to_copy * value_size);

    if (def_levels) {
        if (reader->decoded_def_levels) {
            memcpy(def_levels, reader->decoded_def_levels + reader->page_values_read,
                   (size_t)to_copy * sizeof(int16_t));
        } else {
            memset(def_levels, 0, (size_t)to_copy * sizeof(int16_t));
        }
    }
    if (rep_levels) {
        if (reader->decoded_rep_levels) {
            memcpy(rep_levels, reader->decoded_rep_levels + reader->page_values_read,
                   (size_t)to_copy * sizeof(int16_t));
        } else {
            memset(rep_levels, 0, (size_t)to_copy * sizeof(int16_t));
        }
    }

    /* Update state */
    reader->page_values_read += to_copy;
    reader->values_remaining -= to_copy;
    *values_read = to_copy;

    return CARQUET_OK;
}

/* ============================================================================
 * Skip Values (non-materializing where possible)
 * ============================================================================
 *
 * Advances the reader past up to num_values logical values. Any whole page that
 * fits entirely within the remaining skip count is advanced by parsing only its
 * page header — the compressed payload is never read, decompressed, or decoded.
 * Only a final partial page is decoded (into the reader's reusable decoded
 * buffer; no caller buffer is needed), and values already decoded in a loaded
 * page are dropped for free. This keeps offset-index-driven seeks cheap: a skip
 * of a million rows no longer pays to decode a million values.
 *
 * Returns the number of values actually skipped (may be less than requested if
 * the column chunk is exhausted first).
 */
int64_t carquet_column_skip(
    carquet_column_reader_t* reader,
    int64_t num_values) {

    /* reader is nonnull per API contract */
    if (num_values <= 0 || reader->values_remaining <= 0) {
        return 0;
    }

    carquet_error_t error = CARQUET_ERROR_INIT;

    /* Finalize data_start_offset (past any dictionary page) so that
     * current_page == 0 refers to the first data page. */
    if (carquet_column_ensure_dictionary_loaded(reader, &error) != CARQUET_OK) {
        return 0;
    }

    int64_t total_skipped = 0;
    while (total_skipped < num_values && reader->values_remaining > 0) {
        /* 1. Values already decoded in the loaded page: dropping them is free. */
        if (reader->page_loaded &&
            reader->page_values_read < reader->page_num_values) {
            int64_t avail = (int64_t)reader->page_num_values -
                            (int64_t)reader->page_values_read;
            int64_t want = num_values - total_skipped;
            int64_t n = avail < want ? avail : want;
            if (n > reader->values_remaining) n = reader->values_remaining;
            reader->page_values_read += (int32_t)n;
            reader->values_remaining -= n;
            total_skipped += n;
            continue;
        }

        /* 2. Loaded page fully consumed: advance past it (no I/O). */
        if (reader->page_loaded) {
            reader->current_page += (int64_t)reader->page_header_size +
                                    (int64_t)reader->page_compressed_size;
            reader->page_loaded = false;
            reader->page_num_values = 0;
            reader->page_values_read = 0;
            reader->page_header_size = 0;
            reader->page_compressed_size = 0;
            continue;
        }

        /* 3. No page loaded: peek the next page header (no payload read). */
        int64_t abs_off;
        if (!checked_add_i64(reader->data_start_offset, reader->current_page,
                             &abs_off)) {
            break;
        }
        parquet_page_header_t hdr;
        size_t hdr_size = 0;
        if (parse_header_at_offset(reader->file_reader, abs_off,
                                   &hdr, &hdr_size, &error) != CARQUET_OK) {
            break;  /* End of chunk or malformed: stop, return what we have. */
        }
        if (hdr.type != CARQUET_PAGE_DATA && hdr.type != CARQUET_PAGE_DATA_V2) {
            break;  /* Walked past the chunk's data pages. */
        }
        int32_t page_vals = (hdr.type == CARQUET_PAGE_DATA_V2)
            ? hdr.data_page_header_v2.num_values
            : hdr.data_page_header.num_values;
        if (page_vals < 0) page_vals = 0;

        int64_t advance = (int64_t)hdr_size + (int64_t)hdr.compressed_page_size;
        if (advance <= 0) {
            break;  /* Defensive: never spin on a degenerate header. */
        }

        int64_t want = num_values - total_skipped;
        if ((int64_t)page_vals <= want &&
            (int64_t)page_vals <= reader->values_remaining) {
            /* Whole page fits in the skip: advance by header size only. */
            reader->current_page += advance;
            reader->values_remaining -= page_vals;
            total_skipped += page_vals;
            continue;
        }

        /* 4. Partial page: decode just this page; iteration 1 drains it. */
        if (carquet_column_ensure_page_loaded(reader, &error) != CARQUET_OK ||
            !reader->page_loaded) {
            break;
        }
    }

    return total_skipped;
}
