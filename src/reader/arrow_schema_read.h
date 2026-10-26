/**
 * @file arrow_schema_read.h
 * @brief Parse the "ARROW:schema" footer blob and recover per-field metadata.
 *
 * PyArrow / Arrow C++ (and carquet's own writer) store the original Arrow
 * schema in the Parquet footer under "ARROW:schema" as a base64-encoded,
 * encapsulated Arrow IPC Schema message. That blob is the only place Arrow's
 * per-field `custom_metadata` (variable labels/descriptions) lives — the
 * Parquet SchemaElement wire format cannot express it.
 *
 * This module contains a minimal, bounds-checked FlatBuffer *reader* (the
 * counterpart to the writer in src/writer/arrow_schema.c) that extracts each
 * field's custom_metadata and attaches it to the matching Parquet schema
 * element. It is best-effort: malformed or unexpected input is ignored rather
 * than failing the file open, and only flat top-level fields are matched.
 */
#ifndef CARQUET_ARROW_SCHEMA_READ_H
#define CARQUET_ARROW_SCHEMA_READ_H

#include "core/arena.h"
#include "thrift/parquet_types.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Parse @p b64_value (the "ARROW:schema" metadata value) and populate
 * `field_metadata` / `num_field_metadata` on the top-level schema elements
 * whose names match the Arrow fields. Copies are made in @p arena.
 *
 * @param b64_value      base64 "ARROW:schema" value (may be NULL → no-op).
 * @param elements       Parsed schema elements (element 0 is the root group).
 * @param num_elements   Number of schema elements.
 * @param parent_indices Parent element index per element (-1/0 for root).
 * @param arena          Arena for the copied key/value strings and arrays.
 */
void carquet_apply_arrow_field_metadata(
    const char* b64_value,
    parquet_schema_element_t* elements,
    int32_t num_elements,
    const int32_t* parent_indices,
    carquet_arena_t* arena);

#ifdef __cplusplus
}
#endif

#endif /* CARQUET_ARROW_SCHEMA_READ_H */
