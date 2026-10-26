# API map

Every public function of `include/carquet/carquet.h` (v0.7.x), grouped by task. Signatures are abbreviated to argument types; `err` is a `carquet_error_t*` out-parameter that may be NULL. Functions returning a pointer return NULL on failure and fill `err`; functions returning `carquet_status_t` return `CARQUET_OK` (0) on success. Error helpers live in `include/carquet/error.h` (listed at the end).

## Library

| Function | Returns | Arguments | Purpose |
|---|---|---|---|
| `carquet_version` | `const char*` | `void` | Library version string ("0.7.x"). |
| `carquet_version_components` | `void` | `int* major, int* minor, int* patch` | Get individual version components. |
| `carquet_init` | `carquet_status_t` | `void` | Initialize the Carquet library. |
| `carquet_cleanup` | `void` | `void` | Release library-level resources. |
| `carquet_get_cpu_info` | `const carquet_cpu_info_t*` | `void` | CPU feature information detected at runtime. |
| `carquet_set_allocator` | `void` | `const carquet_allocator_t* allocator` | Custom memory allocator interface. |
| `carquet_get_allocator` | `const carquet_allocator_t*` | `void` | Get the current memory allocator. |
| `carquet_register_codec` | `carquet_status_t` | `carquet_compression_t codec, const carquet_custom_codec_t* impl` | Pluggable compress/decompress implementation for one codec slot. |

## Schema building

| Function | Returns | Arguments | Purpose |
|---|---|---|---|
| `carquet_schema_create` | `carquet_schema_t*` | `err` | Create an empty schema (root group at element index 0). |
| `carquet_schema_free` | `void` | `carquet_schema_t* schema` | Free a schema and all associated resources. |
| `carquet_schema_add_column` | `carquet_status_t` | `carquet_schema_t* schema, const char* name, carquet_physical_type_t physical_type, const carquet_logical_type_t* logical_type, carquet_field_repetition_t repetition, int32_t type_length, int32_t parent_index` | Add a primitive (leaf) column to the schema. |
| `carquet_schema_add_group` | `int32_t` | `carquet_schema_t* schema, const char* name, carquet_field_repetition_t repetition, int32_t parent_index` | Add a group (struct) to the schema for nested structures. |
| `carquet_schema_add_list` | `int32_t` | `carquet_schema_t* schema, const char* name, carquet_physical_type_t element_type, const carquet_logical_type_t* element_logical_type, carquet_field_repetition_t list_repetition, int32_t type_length, int32_t parent_index` | Add a LIST column to the schema using the standard 3-level encoding. |
| `carquet_schema_add_map` | `int32_t` | `carquet_schema_t* schema, const char* name, carquet_physical_type_t key_type, const carquet_logical_type_t* key_logical_type, int32_t key_type_length, carquet_physical_type_t value_type, const carquet_logical_type_t* value_logical_type, int32_t value_type_length, carquet_field_repetition_t map_repetition, int32_t parent_index` | Add a MAP column to the schema using the standard encoding. |
| `carquet_schema_add_list_group` | `int32_t` | `carquet_schema_t* schema, const char* name, carquet_field_repetition_t list_repetition, int32_t parent_index` | Add a LIST container whose element is an arbitrary nested subtree. |
| `carquet_schema_add_map_group` | `int32_t` | `carquet_schema_t* schema, const char* name, carquet_field_repetition_t map_repetition, int32_t parent_index` | Add a MAP container whose key/value are arbitrary nested subtrees. |
| `carquet_schema_add_variant` | `int32_t` | `carquet_schema_t* schema, const char* name, carquet_field_repetition_t variant_repetition, int32_t parent_index` | Add an unshredded VARIANT group to the schema. |
| `carquet_schema_set_field_metadata` | `carquet_status_t` | `carquet_schema_t* schema, int32_t element_index, const char* key, const char* value` | Attach Arrow-style per-field metadata (label, description) to a schema element; emitted with write_arrow_schema. |

## Schema introspection

| Function | Returns | Arguments | Purpose |
|---|---|---|---|
| `carquet_schema_num_columns` | `int32_t` | `const carquet_schema_t* schema` | Get the number of leaf columns in the schema. |
| `carquet_schema_num_elements` | `int32_t` | `const carquet_schema_t* schema` | Get the total number of schema elements (columns + groups). |
| `carquet_schema_get_element` | `const carquet_schema_node_t*` | `const carquet_schema_t* schema, int32_t index` | Get a schema element by index. |
| `carquet_schema_find_column` | `int32_t` | `const carquet_schema_t* schema, const char* name` | Find a column by name. |
| `carquet_schema_column_name` | `const char*` | `const carquet_schema_t* schema, int32_t leaf_index` | Get the name of a leaf column by index. |
| `carquet_schema_column_type` | `carquet_physical_type_t` | `const carquet_schema_t* schema, int32_t leaf_index` | Get the physical type of a leaf column by index. |
| `carquet_schema_column_path` | `int32_t` | `const carquet_schema_t* schema, int32_t leaf_index, const char** path_out, int32_t max_depth` | Get the full schema path for a leaf column. |
| `carquet_schema_max_def_level` | `int16_t` | `const carquet_schema_t* schema, int32_t leaf_index` | Get the accumulated maximum definition level for a leaf column. |
| `carquet_schema_max_rep_level` | `int16_t` | `const carquet_schema_t* schema, int32_t leaf_index` | Get the accumulated maximum repetition level for a leaf column. |
| `carquet_schema_node_name` | `const char*` | `const carquet_schema_node_t* node` | Get the name of a schema node. |
| `carquet_schema_node_is_leaf` | `bool` | `const carquet_schema_node_t* node` | Check if a schema node is a leaf (column) or group. |
| `carquet_schema_node_physical_type` | `carquet_physical_type_t` | `const carquet_schema_node_t* node` | Get the physical type of a leaf node. |
| `carquet_schema_node_logical_type` | `const carquet_logical_type_t*` | `const carquet_schema_node_t* node` | Get the logical type annotation of a node. |
| `carquet_schema_node_repetition` | `carquet_field_repetition_t` | `const carquet_schema_node_t* node` | Get the repetition level of a node. |
| `carquet_schema_node_max_def_level` | `int16_t` | `const carquet_schema_node_t* node` | Get the maximum definition level for a column. |
| `carquet_schema_node_max_rep_level` | `int16_t` | `const carquet_schema_node_t* node` | Get the maximum repetition level for a column. |
| `carquet_schema_node_type_length` | `int32_t` | `const carquet_schema_node_t* node` | Get the type length for a FIXED_LEN_BYTE_ARRAY column. |

## Writer

| Function | Returns | Arguments | Purpose |
|---|---|---|---|
| `carquet_writer_options_init` | `void` | `carquet_writer_options_t* options` | Fill writer options with defaults (UNCOMPRESSED, 128 MB row groups, 1 MB pages, stats + CRC on, page index / bloom off, async_io on). |
| `carquet_writer_create` | `carquet_writer_t*` | `const char* path, const carquet_schema_t* schema, const carquet_writer_options_t* options, err` | Create a new Parquet file for writing. |
| `carquet_writer_create_file` | `carquet_writer_t*` | `FILE* file, const carquet_schema_t* schema, const carquet_writer_options_t* options, err` | Create a writer to a FILE handle. |
| `carquet_writer_create_buffer` | `carquet_writer_t*` | `const carquet_schema_t* schema, const carquet_writer_options_t* options, err` | Create a writer that writes to an internal memory buffer. |
| `carquet_writer_open_append` | `carquet_writer_t*` | `const char* path, const carquet_schema_t* schema, const carquet_writer_options_t* options, err` | Open an existing file and append row groups; schema must match leaf names/types/repetitions. |
| `carquet_writer_write_batch` | `carquet_status_t` | `carquet_writer_t* writer, int32_t column_index, const void* values, int64_t num_values, const int16_t* def_levels, const int16_t* rep_levels` | Write one leaf column: num_values logical rows, dense values, optional def/rep levels. |
| `carquet_writer_write_list_column` | `carquet_status_t` | `carquet_writer_t* writer, int32_t column_index, int64_t num_lists, const int32_t* offsets, const uint8_t* list_validity, const void* values, const uint8_t* value_validity, err` | Write a single-level LIST/MAP leaf from Arrow layout (int32 offsets + list/element validity bitmaps). |
| `carquet_writer_write_arrow` | `carquet_status_t` | `carquet_writer_t* writer, ArrowArray* array, ArrowSchema* schema, err` | Shred and write a top-level Arrow struct array of any nesting; consumes array and schema. |
| `carquet_writer_new_row_group` | `carquet_status_t` | `carquet_writer_t* writer` | Start a new row group. |
| `carquet_writer_num_columns` | `int32_t` | `const carquet_writer_t* writer` | Get the number of leaf columns the writer expects. |
| `carquet_writer_add_metadata` | `carquet_status_t` | `carquet_writer_t* writer, const char* key, const char* value` | Add key-value metadata to the file being written. |
| `carquet_writer_close` | `carquet_status_t` | `carquet_writer_t* writer` | Close the writer and finalize the file. |
| `carquet_writer_abort` | `void` | `carquet_writer_t* writer` | Abort writing and clean up without finalizing the file. |
| `carquet_writer_get_buffer` | `carquet_status_t` | `carquet_writer_t* writer, void** buffer, size_t* size` | Get the buffer contents after writing. |

## Writer per-column overrides (after create, before the first write)

| Function | Returns | Arguments | Purpose |
|---|---|---|---|
| `carquet_writer_set_column_encoding` | `carquet_status_t` | `carquet_writer_t* writer, int32_t column_index, carquet_encoding_t encoding` | Set encoding for a specific column. |
| `carquet_writer_set_column_compression` | `carquet_status_t` | `carquet_writer_t* writer, int32_t column_index, carquet_compression_t codec, int32_t level` | Set compression for a specific column. |
| `carquet_writer_set_column_page_size` | `carquet_status_t` | `carquet_writer_t* writer, int32_t column_index, int64_t bytes` | Override the target byte-based page-flush size for one column. |
| `carquet_writer_set_column_statistics` | `carquet_status_t` | `carquet_writer_t* writer, int32_t column_index, bool enabled` | Enable or disable statistics for a specific column. |
| `carquet_writer_set_max_statistics_size` | `carquet_status_t` | `carquet_writer_t* writer, int64_t bytes` | Maximum stored size of variable-length min/max statistics. |
| `carquet_writer_set_column_bloom_filter` | `carquet_status_t` | `carquet_writer_t* writer, int32_t column_index, bool enabled` | Enable or disable bloom filter for a specific column. |
| `carquet_writer_set_column_bloom_filter_options` | `carquet_status_t` | `carquet_writer_t* writer, int32_t column_index, bool enabled, int64_t ndv, double fpp` | Enable a bloom filter for a column with explicit NDV and false-positive rate. |
| `carquet_writer_set_sorting_columns` | `carquet_status_t` | `carquet_writer_t* writer, const carquet_sorting_column_t* columns, int32_t count` | Declare the row-group sort order in metadata (declaration only; data is not sorted). |

## Reader: open and inspect

| Function | Returns | Arguments | Purpose |
|---|---|---|---|
| `carquet_reader_options_init` | `void` | `carquet_reader_options_t* options` | Fill reader options with defaults (use_mmap off, verify_checksums on, threads auto). |
| `carquet_reader_open` | `carquet_reader_t*` | `const char* path, const carquet_reader_options_t* options, err` | Open a Parquet file for reading. |
| `carquet_reader_open_file` | `carquet_reader_t*` | `FILE* file, const carquet_reader_options_t* options, err` | Open a Parquet file from a FILE handle. |
| `carquet_reader_open_buffer` | `carquet_reader_t*` | `const void* buffer, size_t size, const carquet_reader_options_t* options, err` | Open a Parquet file from a memory buffer. |
| `carquet_reader_close` | `void` | `carquet_reader_t* reader` | Close a reader and release all resources. |
| `carquet_reader_schema` | `const carquet_schema_t*` | `const carquet_reader_t* reader` | Get the file schema. |
| `carquet_reader_num_rows` | `int64_t` | `const carquet_reader_t* reader` | Get the total number of rows in the file. |
| `carquet_reader_num_row_groups` | `int32_t` | `const carquet_reader_t* reader` | Get the number of row groups in the file. |
| `carquet_reader_num_columns` | `int32_t` | `const carquet_reader_t* reader` | Get the number of columns in the file. |
| `carquet_reader_is_mmap` | `bool` | `const carquet_reader_t* reader` | Check if reader is using memory-mapped I/O. |
| `carquet_reader_can_zero_copy` | `bool` | `const carquet_reader_t* reader, int32_t row_group_index, int32_t column_index` | True when a (row group, column) can be served straight from the mmap (uncompressed PLAIN fixed-width REQUIRED). |
| `carquet_reader_row_group_metadata` | `carquet_status_t` | `const carquet_reader_t* reader, int32_t row_group_index, carquet_row_group_metadata_t* metadata` | Row-group metadata: num_rows, total byte size, file offset. |
| `carquet_reader_column_chunk_metadata` | `carquet_status_t` | `const carquet_reader_t* reader, int32_t row_group_index, int32_t column_index, carquet_column_chunk_metadata_t* metadata` | Detailed metadata for a column chunk. |
| `carquet_get_file_info` | `carquet_status_t` | `const char* path, carquet_file_info_t* info, err` | Footer-only summary (size, rows, row groups, columns, version, created_by) without a reader. |
| `carquet_validate_file` | `carquet_status_t` | `const char* path, err` | Validate a Parquet file structure. |
| `carquet_reader_prebuffer` | `carquet_status_t` | `carquet_reader_t* reader, int32_t row_group_index, const int32_t* column_indices, int32_t num_columns, err` | Pre-buffer column data for I/O coalescing. |
| `carquet_reader_release_prebuffer` | `void` | `carquet_reader_t* reader` | Release pre-buffered data. |

## Batch reader (default scan API)

| Function | Returns | Arguments | Purpose |
|---|---|---|---|
| `carquet_batch_reader_config_init` | `void` | `carquet_batch_reader_config_t* config` | Fill batch-reader config with defaults (batch_size 65536, all columns, threads auto). |
| `carquet_batch_reader_create` | `carquet_batch_reader_t*` | `carquet_reader_t* reader, const carquet_batch_reader_config_t* config, err` | Create a batch reader for efficient columnar reading. |
| `carquet_batch_reader_next` | `carquet_status_t` | `carquet_batch_reader_t* batch_reader, carquet_row_batch_t** batch` | Read the next batch of rows. |
| `carquet_batch_reader_last_error` | `const carquet_error_t*` | `const carquet_batch_reader_t* batch_reader` | Detail for the most recent carquet_batch_reader_next() failure. |
| `carquet_batch_reader_set_page_filter` | `carquet_status_t` | `carquet_batch_reader_t* reader, const carquet_filter_clause_t* clauses, int32_t count` | Install AND-ed predicate clauses: row groups pruned by stats/bloom, pages by the column index; NULL/0 clears. |
| `carquet_batch_reader_rows_skipped` | `int64_t` | `const carquet_batch_reader_t* reader` | Number of rows skipped by the active page filter so far. |
| `carquet_batch_reader_set_row_range` | `carquet_status_t` | `carquet_batch_reader_t* reader, int64_t offset, int64_t limit` | Restrict the scan to logical rows [offset, offset+limit); row groups outside are skipped without I/O. |
| `carquet_batch_reader_free` | `void` | `carquet_batch_reader_t* batch_reader` | Free a batch reader. |
| `carquet_row_batch_num_rows` | `int64_t` | `const carquet_row_batch_t* batch` | Get the number of rows in a batch. |
| `carquet_row_batch_num_columns` | `int32_t` | `const carquet_row_batch_t* batch` | Get the number of columns in a batch. |
| `carquet_row_batch_column` | `carquet_status_t` | `const carquet_row_batch_t* batch, int32_t column_index, const void** data, const uint8_t** null_bitmap, int64_t* num_values` | Get column data from a batch. |
| `carquet_row_batch_column_list` | `carquet_status_t` | `const carquet_row_batch_t* batch, int32_t column_index, const int32_t** offsets, int64_t* num_lists, const void** values, const uint8_t** value_validity, int64_t* num_values, const uint8_t** list_validity` | Single-level LIST/MAP column as Arrow layout: offsets, flattened values, element and list validity. |
| `carquet_row_batch_column_dictionary` | `carquet_status_t` | `const carquet_row_batch_t* batch, int32_t column_index, const uint32_t** indices, const uint8_t** null_bitmap, int64_t* num_values, const uint8_t** dictionary_data, int32_t* dictionary_count, const uint32_t** dictionary_offsets` | Dictionary-preserved column: uint32 indices + dictionary bytes/offsets (needs preserve_dictionaries). |
| `carquet_row_batch_free` | `void` | `carquet_row_batch_t* batch` | Free a row batch. |
| `carquet_thread_pool_create` | `carquet_thread_pool_t*` | `int32_t num_threads` | Create a reusable thread pool for parallel reading. |
| `carquet_thread_pool_destroy` | `void` | `carquet_thread_pool_t* pool` | Destroy a thread pool. |

## Column reader (levels, one row group)

| Function | Returns | Arguments | Purpose |
|---|---|---|---|
| `carquet_reader_get_column` | `carquet_column_reader_t*` | `carquet_reader_t* reader, int32_t row_group_index, int32_t column_index, err` | Get a column reader for a specific row group and column. |
| `carquet_column_read_batch` | `int64_t` | `carquet_column_reader_t* reader, void* values, int64_t max_values, int16_t* def_levels, int16_t* rep_levels` | Read a batch of values from a column. |
| `carquet_column_read_batch_ex` | `int64_t` | `carquet_column_reader_t* reader, void* values, int64_t max_values, int16_t* def_levels, int16_t* rep_levels, err` | Like read_batch but reports partial reads and failures through a carquet_error_t. |
| `carquet_column_skip` | `int64_t` | `carquet_column_reader_t* reader, int64_t num_values` | Skip values in a column without reading them. |
| `carquet_column_has_next` | `bool` | `const carquet_column_reader_t* reader` | Check if there are more values to read. |
| `carquet_column_remaining` | `int64_t` | `const carquet_column_reader_t* reader` | Get the number of remaining values in the column. |
| `carquet_column_reader_free` | `void` | `carquet_column_reader_t* reader` | Free a column reader. |
| `carquet_count_rows` | `int64_t` | `const int16_t* rep_levels, int64_t num_values` | Count logical rows from repetition levels. |
| `carquet_list_offsets` | `int64_t` | `const int16_t* rep_levels, int64_t num_values, int16_t list_rep_level, int64_t* offsets_out, int64_t max_offsets` | Compute list offsets from repetition levels. |

## Statistics and pruning

| Function | Returns | Arguments | Purpose |
|---|---|---|---|
| `carquet_reader_column_statistics` | `carquet_status_t` | `const carquet_reader_t* reader, int32_t row_group_index, int32_t column_index, carquet_column_statistics_t* stats` | Chunk statistics for one (row group, column): has_min_max, min/max bytes, null_count, distinct_count, num_values. |
| `carquet_reader_row_group_matches` | `carquet_status_t` | `const carquet_reader_t* reader, int32_t row_group_index, int32_t column_index, carquet_compare_op_t op, const void* value, int32_t value_size, bool* might_match` | Might this row group satisfy `column OP value`? Conservative (true when unknown). |
| `carquet_reader_filter_row_groups` | `int32_t` | `const carquet_reader_t* reader, int32_t column_index, carquet_compare_op_t op, const void* value, int32_t value_size, int32_t* matching_indices, int32_t max_indices` | Collect the indices of row groups that might satisfy `column OP value`. |
| `carquet_reader_geospatial_statistics` | `carquet_status_t` | `const carquet_reader_t* reader, int32_t row_group_index, int32_t column_index, carquet_geospatial_statistics_t* stats` | Bounding box and WKB geometry types for a GEOMETRY/GEOGRAPHY chunk (INVALID_METADATA when absent). |

## Bloom filters and page indexes

| Function | Returns | Arguments | Purpose |
|---|---|---|---|
| `carquet_reader_get_bloom_filter` | `carquet_bloom_filter_t*` | `carquet_reader_t* reader, int32_t row_group_index, int32_t column_index, err` | Read a bloom filter for a column in a row group. |
| `carquet_bloom_filter_check_i32` | `bool` | `const carquet_bloom_filter_t* filter, int32_t value` | Check if a bloom filter might contain an int32 value. |
| `carquet_bloom_filter_check_i64` | `bool` | `const carquet_bloom_filter_t* filter, int64_t value` | Check if a bloom filter might contain an int64 value. |
| `carquet_bloom_filter_check_float` | `bool` | `const carquet_bloom_filter_t* filter, float value` | Check if a bloom filter might contain a float value. |
| `carquet_bloom_filter_check_double` | `bool` | `const carquet_bloom_filter_t* filter, double value` | Check if a bloom filter might contain a double value. |
| `carquet_bloom_filter_check_bytes` | `bool` | `const carquet_bloom_filter_t* filter, const uint8_t* data, size_t len` | Check if a bloom filter might contain a byte sequence. |
| `carquet_bloom_filter_size` | `size_t` | `const carquet_bloom_filter_t* filter` | Get bloom filter size in bytes. |
| `carquet_bloom_filter_destroy` | `void` | `carquet_bloom_filter_t* filter` | Free a bloom filter. |
| `carquet_reader_get_column_index` | `carquet_column_index_t*` | `carquet_reader_t* reader, int32_t row_group_index, int32_t column_index, err` | Load a column chunk ColumnIndex (per-page min/max/null counts); free with carquet_column_index_free. |
| `carquet_column_index_num_pages` | `int32_t` | `const carquet_column_index_t* index` | Get the number of pages in a column index. |
| `carquet_column_index_get_page_stats` | `carquet_status_t` | `const carquet_column_index_t* index, int32_t page_index, carquet_page_stats_t* stats` | Get per-page statistics from a column index. |
| `carquet_column_index_boundary_order` | `int32_t` | `const carquet_column_index_t* index` | Get boundary order of a column index. |
| `carquet_column_index_free` | `void` | `carquet_column_index_t* index` | Free a column index. |
| `carquet_reader_get_offset_index` | `carquet_offset_index_t*` | `carquet_reader_t* reader, int32_t row_group_index, int32_t column_index, err` | Read offset index (per-page locations) for a column chunk. |
| `carquet_offset_index_num_pages` | `int32_t` | `const carquet_offset_index_t* index` | Get the number of pages in an offset index. |
| `carquet_offset_index_get_page_location` | `carquet_status_t` | `const carquet_offset_index_t* index, int32_t page_index, carquet_page_location_t* location` | Get page location from an offset index. |
| `carquet_offset_index_free` | `void` | `carquet_offset_index_t* index` | Free an offset index. |

## Footer and per-field metadata

| Function | Returns | Arguments | Purpose |
|---|---|---|---|
| `carquet_reader_num_metadata` | `int32_t` | `const carquet_reader_t* reader` | Get the number of key-value metadata entries in the file. |
| `carquet_reader_get_metadata` | `carquet_status_t` | `const carquet_reader_t* reader, int32_t index, const char** key, const char** value` | Get a key-value metadata entry by index. |
| `carquet_reader_find_metadata` | `const char*` | `const carquet_reader_t* reader, const char* key` | Find a metadata value by key. |
| `carquet_reader_column_num_metadata` | `int32_t` | `const carquet_reader_t* reader, int32_t column_index` | Number of Arrow per-field metadata entries for a leaf column. |
| `carquet_reader_column_get_metadata` | `carquet_status_t` | `const carquet_reader_t* reader, int32_t column_index, int32_t index, const char** key, const char** value` | Get an Arrow per-field metadata entry for a leaf column by index. |
| `carquet_reader_column_find_metadata` | `const char*` | `const carquet_reader_t* reader, int32_t column_index, const char* key` | Find an Arrow per-field metadata value for a leaf column by key. |
| `carquet_reader_column_arrow_type_refinement` | `carquet_arrow_type_refinement_t` | `const carquet_reader_t* reader, int32_t column_index` | Arrow refinement recovered from ARROW:schema (LargeUtf8 / LargeBinary), or NONE. |

## Arrow C Data Interface bridge

| Function | Returns | Arguments | Purpose |
|---|---|---|---|
| `carquet_arrow_export_schema` | `carquet_status_t` | `const carquet_schema_t* schema, ArrowSchema* out, err` | Export a flat Carquet schema as an Arrow C Data Interface schema. |
| `carquet_arrow_export_batch` | `carquet_status_t` | `const carquet_row_batch_t* batch, const carquet_schema_t* schema, ArrowSchema* out_schema, ArrowArray* out_array, err` | Export a flat row batch (plus single-level LIST/MAP) as ArrowSchema + ArrowArray (owned copies; release once). |
| `carquet_arrow_import_schema` | `carquet_status_t` | `ArrowSchema* schema, carquet_schema_t** out, err` | Build a carquet schema from an Arrow "+s" schema (nested types, decimal/temporal coercions); consumes it. |
| `carquet_writer_write_arrow` | `carquet_status_t` | `carquet_writer_t* writer, ArrowArray* array, ArrowSchema* schema, err` | Shred and write a top-level Arrow struct array of any nesting; consumes array and schema. |
| `carquet_reader_read_arrow` | `carquet_status_t` | `carquet_reader_t* reader, int32_t row_group_index, ArrowSchema* out_schema, ArrowArray* out_array, err` | Read one row group into a nested ArrowArray tree (any depth of struct/list/map); release once. |

## Error helpers (`carquet/error.h`)

| Function / macro | Purpose |
|---|---|
| `CARQUET_ERROR_INIT` | Initialiser for a `carquet_error_t` on the stack (`carquet_error_t err = CARQUET_ERROR_INIT;`; in C++ use `carquet_error_t err{};`) |
| `carquet_status_string(code)` | Enumerator name for a status code |
| `carquet_error_format(&err, buf, n)` | One-line rendering: code, message, source location, file offset, column / row group |
| `carquet_error_recovery_hint(code)` | Human hint for a code, or NULL |
| `carquet_error_is_set(&err)` / `carquet_error_code(&err)` / `carquet_error_message(&err)` | Inline accessors |
| `carquet_error_init` / `carquet_error_clear` / `carquet_error_copy` | Lifecycle helpers |
| `CARQUET_RETURN_IF_ERROR(status)` | Early return on failure (for code extending carquet) |

## Key structs (fields you set)

- `carquet_writer_options_t`: `compression`, `compression_level`, `row_group_size` (bytes, default 128 MB), `page_size` (default 1 MB), `write_statistics`, `write_crc`, `write_page_index`, `write_bloom_filters`, `dictionary_encoding`, `dictionary_page_size`, `max_rows_per_page`, `write_arrow_schema`, `data_page_version` (1 or 2), `coerce_timestamps` + `coerce_timestamp_unit` + `allow_timestamp_truncation`, `write_batch_size`, `file_format_version`, `async_io`.
- `carquet_reader_options_t`: `use_mmap`, `verify_checksums`, `buffer_size`, `num_threads`.
- `carquet_batch_reader_config_t`: `batch_size`, `num_threads`, `use_mmap`, `column_indices` + `num_columns` (take precedence) or `column_names` + `num_column_names`, `row_group_filter` + `row_group_filter_ctx`, `preserve_dictionaries`, `thread_pool`.
- `carquet_filter_clause_t`: `column_index`, `op` (`CARQUET_FILTER_EQ/NE/LT/LE/GT/GE/RANGE/IN/IS_NULL/IS_NOT_NULL`), `value` + `value_size`; for RANGE `has_lo`/`lo`/`lo_size` and `has_hi`/`hi`/`hi_size`; for IN `values` + `value_count`.
- `carquet_logical_type_t`: `id` plus `params.decimal {precision, scale}`, `params.integer {bit_width, is_signed}`, `params.time` / `params.timestamp {unit (CARQUET_TIME_UNIT_MILLIS/MICROS/NANOS), is_adjusted_to_utc}`, `params.geography {crs, algorithm, has_algorithm}`, `params.geometry {crs}`, `params.variant {specification_version}`.
- `carquet_byte_array_t`: `uint8_t* data; int32_t length;` (used for BYTE_ARRAY values in both directions).
- `carquet_column_statistics_t`: `has_min_max`, `has_null_count`, `has_distinct_count`, `null_count`, `distinct_count`, `num_values`, `min_value`/`max_value` (raw bytes) + sizes.
- `carquet_column_chunk_metadata_t`: `type`, `codec`, `num_values`, `total_compressed_size`, `total_uncompressed_size`, `data_page_offset`, `has_dictionary_page`, `dictionary_page_offset`, `encodings[]`, `has_bloom_filter`, `has_column_index`, `has_offset_index`.
- `carquet_file_info_t`: `file_size`, `num_rows`, `num_row_groups`, `num_columns`, `version`, `created_by`.
- `carquet_custom_codec_t`: `compress(src, src_size, dst, dst_capacity, &out_size, level, user_data)`, `decompress(...)`, `compress_bound(src_size, user_data)`, `user_data`.
