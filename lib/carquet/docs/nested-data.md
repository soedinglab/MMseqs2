# Nested and Nullable Data

## The Mental Model

Carquet gives you direct Parquet semantics:

- schema tree
- leaf columns
- definition levels
- repetition levels

There is no row-object layer in the C API. Even for nested data, you still write and read leaf columns plus level streams.

## Prefer the Schema Helpers

Use the helpers when your data matches standard Parquet layouts:

- `carquet_schema_add_group()`
- `carquet_schema_add_list()`
- `carquet_schema_add_map()`

Example:

```c
carquet_schema_t* schema = carquet_schema_create(NULL);
carquet_schema_add_column(schema, "id",
    CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

int32_t scores_group = carquet_schema_add_list(schema, "scores",
    CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_OPTIONAL, 0, 0);

int32_t tags_group = carquet_schema_add_map(schema, "tags",
    CARQUET_PHYSICAL_BYTE_ARRAY, NULL, 0,
    CARQUET_PHYSICAL_BYTE_ARRAY, NULL, 0,
    CARQUET_REPETITION_OPTIONAL, 0);
```

These helpers build the standard 3-level Parquet encodings for LIST and MAP. That keeps the schema interoperable with Arrow, DuckDB, Spark, and other readers.

## Definition Levels

Definition levels answer "how much of this path exists?"

Common cases:

- required scalar: no `def_levels`
- optional scalar: level at max definition means present, lower means null
- nested optional structures: higher levels mean more of the path is defined

Useful schema queries:

- `carquet_schema_max_def_level()`
- `carquet_schema_node_max_def_level()`

For flat optional columns, the pattern is simple:

```c
int16_t def_levels[] = {1, 0, 1, 1, 0};
```

For lists and maps, compute the exact max levels from the schema rather than hard-coding them.

## Repetition Levels

Repetition levels answer "did this value continue the current repeated container, or start a new one?"

Useful schema queries:

- `carquet_schema_max_rep_level()`
- `carquet_schema_node_max_rep_level()`

When reading repeated data, ask for `rep_levels` from `carquet_column_read_batch()` and reconstruct list boundaries from them.

## Reconstruct Lists While Reading

Top-level `LIST` columns created with `carquet_schema_add_list()` usually use `list_rep_level = 1`.

```c
int32_t values[256];
int16_t rep_levels[256];
int64_t offsets[257];

int64_t n = carquet_column_read_batch(col, values, 256, NULL, rep_levels);
int64_t num_lists = carquet_list_offsets(rep_levels, n, 1, offsets, 257);

for (int64_t i = 0; i < num_lists; i++) {
    int64_t begin = offsets[i];
    int64_t end = offsets[i + 1];
    /* values[begin..end-1] belongs to list i */
}
```

Helpers:

- `carquet_count_rows()`: count top-level logical rows from repetition levels
- `carquet_list_offsets()`: build Arrow-style offsets for one repeated level

## Writing Repeated Data

For single-level `LIST<T>` and `MAP<K,V>` (schemas built with `add_list()` / `add_map()`), use the high-level helper `carquet_writer_write_list_column()` — it takes Arrow-style `int32` offsets plus optional list/element validity bitmaps and computes the levels for you:

```c
/* [[10, 20], [], NULL, [30]] on the leaf element column. */
int32_t offsets[]     = {0, 2, 2, 2, 3};
int32_t values[]      = {10, 20, 30};
uint8_t list_validity = 0x0B;            /* row 2 (bit 2) is a null list */
carquet_writer_write_list_column(writer, elem_col, 4, offsets, &list_validity,
                                 values, /*value_validity=*/NULL, &err);
```

For a `MAP<K,V>`, call it twice with the same `offsets` and `list_validity` — once for the key leaf and once for the value leaf. It is the write-side inverse of `carquet_row_batch_column_list()`; see [`writing.md`](./writing.md#write-nested-columns-without-computing-levels).

If you need full control (deeper nesting, or a non-standard encoding), write the leaf values plus matching `def_levels` and `rep_levels` yourself:

1. Build the schema with `add_list()` / `add_map()` or explicit groups.
2. Query max def/rep levels from the schema.
3. Emit one leaf stream per column with Parquet-correct levels via `carquet_writer_write_batch()`.

## Schema Introspection

Use schema introspection when you need to validate generated level streams or map projected columns back to paths:

- `carquet_schema_find_column()`
- `carquet_schema_column_name()`
- `carquet_schema_column_path()`
- `carquet_schema_get_element()`
- `carquet_schema_node_name()`
- `carquet_schema_node_is_leaf()`
- `carquet_schema_node_physical_type()`
- `carquet_schema_node_logical_type()`
- `carquet_schema_node_repetition()`
- `carquet_schema_node_type_length()`

This is especially useful for generic readers and schema-driven code generation.

## Arbitrary-Depth Nesting via the Arrow C Data Interface

For deeply nested data — `LIST<LIST<T>>`, `LIST<STRUCT<...>>`, `MAP<K, LIST<V>>`, structs of lists, and any composition — the simplest path is the Arrow C Data Interface bridge, which shreds and reassembles the full tree for you:

- **Write**: `carquet_writer_write_arrow(writer, array, schema, &err)` takes a top-level Arrow struct array (`format = "+s"`) and runs a generic Dremel record-shredder, computing each leaf column's repetition/definition levels and dense values internally. It handles `struct` (`+s`), `list` (`+l`), `large_list` (`+L`) and `map` (`+m`) composed to any depth. Both structs are consumed (Arrow move semantics).
- **Read**: `carquet_reader_read_arrow(reader, row_group_index, &out_schema, &out_array, &err)` reassembles one row group into the matching nested Arrow struct array. Every buffer is an owned copy released via the standard `release` callback, so the result outlives the reader.
- **Schema mapping** (`carquet_arrow_import_schema()` / `carquet_arrow_export_schema()`) is recursive: a Parquet 3-level `LIST` ↔ Arrow `List<element>`, a Parquet `MAP` ↔ Arrow `Map<entries: struct<key, value>>`, and a plain group ↔ `Struct`.

To build a nested schema by hand (without going through Arrow), use the composable builders, which return the inner `REPEATED` group so you can attach an arbitrary element/key/value subtree:

- `carquet_schema_add_list_group()` — returns the inner `list` group; add exactly one child (the element), which may itself be a leaf, a struct (`add_group()`), or another list/map.
- `carquet_schema_add_map_group()` — returns the inner `key_value` group; add a required `key` and a `value`, each of which may be nested.

```c
/* list<list<int32>> */
int32_t inner = carquet_schema_add_list_group(schema, "matrix",
    CARQUET_REPETITION_OPTIONAL, 0);
int32_t inner2 = carquet_schema_add_list_group(schema, "element",
    CARQUET_REPETITION_OPTIONAL, inner);
carquet_schema_add_column(schema, "element", CARQUET_PHYSICAL_INT32, NULL,
    CARQUET_REPETITION_OPTIONAL, 0, inner2);
```

The single-level `carquet_writer_write_list_column()` / `carquet_row_batch_column_list()` helpers remain the lightweight path when your data is exactly one level of `LIST`/`MAP`.
