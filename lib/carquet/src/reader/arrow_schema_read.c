/**
 * @file arrow_schema_read.c
 * @brief Minimal bounds-checked FlatBuffer reader for the "ARROW:schema" blob.
 *
 * Navigates Message -> Schema -> [Field] -> Field.custom_metadata -> [KeyValue]
 * and attaches each field's metadata to the matching Parquet schema element.
 * The input is untrusted (it comes from the file), so every offset, length and
 * vector span is validated before use; anything inconsistent aborts the parse
 * of that sub-tree without touching the rest.
 */

#include "arrow_schema_read.h"
#include "core/allocator.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Arrow FlatBuffer field ids (format/Message.fbs, Schema.fbs). A union in the
 * schema occupies two vtable slots (type, then value), which is why the Schema
 * table sits at Message slot 2. */
enum { MSG_HEADER_SLOT = 2 };      /* Message.header (union value) */
enum { SCHEMA_FIELDS_SLOT = 1 };   /* Schema.fields */
enum { FIELD_NAME_SLOT = 0 };      /* Field.name */
enum { FIELD_TYPETYPE_SLOT = 2 };  /* Field.type_type (union tag, u8) */
enum { FIELD_META_SLOT = 6 };      /* Field.custom_metadata */
enum { KV_KEY_SLOT = 0, KV_VALUE_SLOT = 1 };

/* Arrow Type union tags for the 64-bit-offset variants Parquet can't express
 * (format/Type.fbs). Kept in sync with carquet_arrow_type_refinement_t. */
enum { AT_LARGEBINARY = 19, AT_LARGEUTF8 = 20, AT_LARGELIST = 21 };

#include "fb_read.h"

/* Map a Field's inline type_type (union tag, u8) to a carquet refinement.
 * Returns 0 (CARQUET_ARROW_REFINE_NONE) when absent or not a 64-bit variant. */
static int32_t field_type_refinement(const fbr* r, size_t field_tbl) {
    size_t f;
    if (!fbr_field(r, field_tbl, FIELD_TYPETYPE_SLOT, &f) || f == 0) return 0;
    uint8_t tag;
    if (!fbr_u8(r, f, &tag)) return 0;
    switch (tag) {
        case AT_LARGEUTF8:   return 1;  /* CARQUET_ARROW_REFINE_LARGE_UTF8   */
        case AT_LARGEBINARY: return 2;  /* CARQUET_ARROW_REFINE_LARGE_BINARY */
        case AT_LARGELIST:   return 3;  /* CARQUET_ARROW_REFINE_LARGE_LIST   */
        default:             return 0;
    }
}

/* ---- attach metadata to the matching top-level schema element ---- */

static int32_t match_element(const char* name, size_t name_len,
                             parquet_schema_element_t* elements,
                             int32_t num_elements,
                             const int32_t* parent_indices,
                             const uint8_t* used) {
    if (!name) return -1;
    /* Match by name among direct children of the root that have no metadata
     * yet (used[] guards against duplicate names being over-written). */
    for (int32_t i = 1; i < num_elements; i++) {
        int32_t parent = parent_indices ? parent_indices[i] : 0;
        if (parent != 0) continue;
        if (used[i]) continue;
        if (elements[i].name && strlen(elements[i].name) == name_len &&
            memcmp(elements[i].name, name, name_len) == 0) return i;
    }
    return -1;
}

void carquet_apply_arrow_field_metadata(
    const char* b64_value,
    parquet_schema_element_t* elements,
    int32_t num_elements,
    const int32_t* parent_indices,
    carquet_arena_t* arena) {

    if (!b64_value || !elements || num_elements < 2 || !arena) return;

    size_t raw_len = 0;
    uint8_t* raw = base64_decode(b64_value, &raw_len);
    if (!raw) return;

    /* Skip the Arrow IPC encapsulation prefix: either the modern
     * 0xFFFFFFFF continuation marker + u32 length, or the legacy bare u32. */
    size_t fb_start;
    if (raw_len >= 8 && raw[0] == 0xFF && raw[1] == 0xFF &&
        raw[2] == 0xFF && raw[3] == 0xFF) {
        fb_start = 8;
    } else if (raw_len >= 4) {
        fb_start = 4;
    } else {
        carquet_mem_free(raw);
        return;
    }

    fbr r = { raw + fb_start, raw_len - fb_start };

    uint8_t* used = (uint8_t*)carquet_mem_calloc((size_t)num_elements, 1);
    if (!used) { carquet_mem_free(raw); return; }

    size_t msg, hdr_field, schema_tbl, fields_data;
    uint32_t nfields;
    if (!fbr_indirect(&r, 0, &msg)) goto done;                         /* root Message */
    if (!fbr_field(&r, msg, MSG_HEADER_SLOT, &hdr_field) || hdr_field == 0) goto done;
    if (!fbr_indirect(&r, hdr_field, &schema_tbl)) goto done;          /* Schema */
    if (!fbr_field_vector(&r, schema_tbl, SCHEMA_FIELDS_SLOT, 4,
                          &fields_data, &nfields)) goto done;

    /* Only root children can match, so a legitimate schema has fewer
     * top-level fields than Parquet elements; ignore any excess rather than
     * scanning (num_elements x 2^20) candidates. */
    if (nfields > (uint32_t)num_elements) nfields = (uint32_t)num_elements;

    /* Every string of a well-formed blob is stored (and so copied) once;
     * allow generous headroom for shared strings but bound amplification. */
    size_t copy_budget = (r.len <= (SIZE_MAX - 4096) / 4)
                             ? r.len * 4 + 4096 : SIZE_MAX;

    for (uint32_t j = 0; j < nfields; j++) {
        size_t field_tbl;
        if (!fbr_indirect(&r, fields_data + (size_t)j * 4, &field_tbl)) continue;

        const char* name = NULL;
        size_t name_len = 0;
        if (!fbr_field_string_view(&r, field_tbl, FIELD_NAME_SLOT,
                                   &name, &name_len)) continue;

        /* Match the Arrow field to its Parquet element up front so both the
         * type refinement and the custom_metadata can be attached, even when
         * the field carries no metadata. */
        int32_t idx = match_element(name, name_len, elements, num_elements,
                                    parent_indices, used);
        if (idx < 0) continue;
        used[idx] = 1;

        /* Type refinement: a 64-bit-offset Arrow type Parquet can't express. */
        int32_t refine = field_type_refinement(&r, field_tbl);
        if (refine != 0) elements[idx].arrow_type_refinement = refine;

        /* custom_metadata (variable labels/descriptions). */
        size_t meta_data;
        uint32_t nmeta;
        if (!fbr_field_vector(&r, field_tbl, FIELD_META_SLOT, 4,
                              &meta_data, &nmeta) || nmeta == 0) continue;
        /* Fields may share one metadata vector: charge the array too. */
        if ((size_t)nmeta * sizeof(parquet_key_value_t) > copy_budget) continue;
        copy_budget -= (size_t)nmeta * sizeof(parquet_key_value_t);

        parquet_key_value_t* kvs = (parquet_key_value_t*)carquet_arena_calloc(
            arena, nmeta, sizeof(parquet_key_value_t));
        if (!kvs) continue;

        int32_t got = 0;
        for (uint32_t k = 0; k < nmeta; k++) {
            size_t kv_tbl;
            if (!fbr_indirect(&r, meta_data + (size_t)k * 4, &kv_tbl)) continue;
            char* key = fbr_field_string(&r, kv_tbl, KV_KEY_SLOT, arena,
                                         &copy_budget);
            if (!key) continue;  /* KeyValue.key is required to be useful */
            kvs[got].key = key;
            kvs[got].value = fbr_field_string(&r, kv_tbl, KV_VALUE_SLOT, arena,
                                              &copy_budget);
            got++;
        }
        if (got == 0) continue;

        elements[idx].field_metadata = kvs;
        elements[idx].num_field_metadata = got;
    }

done:
    carquet_mem_free(used);
    carquet_mem_free(raw);
}
