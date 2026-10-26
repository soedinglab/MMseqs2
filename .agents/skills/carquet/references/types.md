# Type reference

## Physical type ↔ C buffer (write_batch input, column/batch reader output)

| Physical type | C element | Bytes | Notes |
|---|---|---|---|
| `CARQUET_PHYSICAL_BOOLEAN` | `uint8_t` (0 / 1) | 1 | one byte per value in the API, bit-packed on disk |
| `CARQUET_PHYSICAL_INT32` | `int32_t` | 4 | also carries INTEGER(8/16/32, signed or unsigned), DATE, TIME(MILLIS), DECIMAL ≤ 9 digits |
| `CARQUET_PHYSICAL_INT64` | `int64_t` | 8 | also INTEGER(64, signed or unsigned), TIME(MICROS/NANOS), TIMESTAMP, DECIMAL ≤ 18 digits |
| `CARQUET_PHYSICAL_INT96` | `carquet_int96_t` | 12 | deprecated Impala/Hive timestamps; read and written PLAIN only, no statistics, cannot be filtered |
| `CARQUET_PHYSICAL_FLOAT` | `float` | 4 | BYTE_STREAM_SPLIT applied automatically under a codec |
| `CARQUET_PHYSICAL_DOUBLE` | `double` | 8 | idem |
| `CARQUET_PHYSICAL_BYTE_ARRAY` | `carquet_byte_array_t { uint8_t* data; int32_t length; }` | ptr + len | STRING, JSON, BSON, ENUM, DECIMAL (big-endian two's complement, variable length), GEOMETRY/GEOGRAPHY (WKB) |
| `CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY` | `uint8_t[type_length]` packed | N | UUID (16), FLOAT16 (2), INTERVAL (12), DECIMAL (big-endian two's complement, N bytes), arbitrary fixed binary |

Read-side layout differs by API: `carquet_row_batch_column()` is row-aligned (one slot per row, null slots zero-filled, bitmap bit `i` = row `i` present, LSB first, NULL bitmap means all present); `carquet_column_read_batch()` is dense (present values only) with one definition level per row.

## Logical annotations (`carquet_logical_type_t`)

| `id` | Physical | `params` | Meaning of the stored integer / bytes |
|---|---|---|---|
| `CARQUET_LOGICAL_STRING` | BYTE_ARRAY | none | UTF-8 |
| `CARQUET_LOGICAL_ENUM`, `JSON`, `BSON` | BYTE_ARRAY | none | bytes as named |
| `CARQUET_LOGICAL_DATE` | INT32 | none | days since 1970-01-01 |
| `CARQUET_LOGICAL_TIME` | INT32 (MILLIS) or INT64 (MICROS, NANOS) | `.time.unit`, `.time.is_adjusted_to_utc` | time of day in the unit |
| `CARQUET_LOGICAL_TIMESTAMP` | INT64 | `.timestamp.unit`, `.timestamp.is_adjusted_to_utc` | since epoch in the unit (`CARQUET_TIME_UNIT_MILLIS` / `MICROS` / `NANOS`) |
| `CARQUET_LOGICAL_INTEGER` | INT32 (8/16/32 bits) or INT64 (64) | `.integer.bit_width`, `.integer.is_signed` | value in the full-width C type; unsigned ordering applied to statistics |
| `CARQUET_LOGICAL_DECIMAL` | INT32, INT64, FLBA, BYTE_ARRAY | `.decimal.precision`, `.decimal.scale` | unscaled integer (`12.34` with scale 2 is `1234`); byte-backed forms are big-endian two's complement |
| `CARQUET_LOGICAL_UUID` | FLBA(16) | none | 16 raw bytes |
| `CARQUET_LOGICAL_FLOAT16` | FLBA(2) | none | IEEE half, little-endian; statistics ordered numerically |
| `CARQUET_LOGICAL_INTERVAL` | FLBA(12) | none | months, days, millis (3 × uint32 LE) |
| `CARQUET_LOGICAL_GEOMETRY` / `GEOGRAPHY` | BYTE_ARRAY | `.geometry.crs` / `.geography.crs`, `.algorithm`, `.has_algorithm` | ISO WKB; writer computes bounding-box statistics |
| `CARQUET_LOGICAL_VARIANT` | group | via `carquet_schema_add_variant()` | unshredded `{metadata, value}` group |
| `CARQUET_LOGICAL_LIST` / `MAP` | group | via `carquet_schema_add_list()` / `add_map()` (3-level standard layout) | |
| `CARQUET_LOGICAL_NULL` | any | none | all-null column |

Rules enforced by `carquet_schema_add_column()` (error `CARQUET_ERROR_INVALID_ARGUMENT` otherwise): STRING/ENUM/JSON/BSON/GEOMETRY/GEOGRAPHY need BYTE_ARRAY; DATE needs INT32; TIME(MILLIS) INT32, TIME(MICROS/NANOS) INT64; TIMESTAMP INT64; DECIMAL precision ≤ 9 on INT32, ≤ 18 on INT64, ≤ floor((8N−1)·log10 2) on FLBA(N) (16 bytes = 38 digits), `0 ≤ scale ≤ precision`; UUID FLBA(16); FLOAT16 FLBA(2); INTERVAL FLBA(12). The writer also emits the legacy `ConvertedType` for every annotation that has one, so old readers (Hive, Spark 2) understand the files.

## Encodings and codecs: what is allowed where

| Encoding (`carquet_writer_set_column_encoding`) | Physical types | Notes |
|---|---|---|
| `CARQUET_ENCODING_PLAIN` | all | default; fastest reads, zero-copy eligible |
| `CARQUET_ENCODING_RLE_DICTIONARY` | INT32, INT64, FLOAT, DOUBLE, BYTE_ARRAY, FLBA | PLAIN dictionary page + RLE indices; automatic fallback to PLAIN when the dictionary exceeds `dictionary_page_size` or the data is all-unique (per row group) |
| `CARQUET_ENCODING_BYTE_STREAM_SPLIT` | FLOAT, DOUBLE, INT32, INT64, FLBA | automatic for FLOAT/DOUBLE under a codec |
| `CARQUET_ENCODING_DELTA_BINARY_PACKED` | INT32, INT64 | sorted or slowly varying integers |
| `CARQUET_ENCODING_DELTA_LENGTH_BYTE_ARRAY`, `DELTA_BYTE_ARRAY` | BYTE_ARRAY (DELTA_BYTE_ARRAY also FLBA) | strings with shared prefixes |
| `CARQUET_ENCODING_RLE` | BOOLEAN | long boolean runs |

| Codec (`opts.compression` / `set_column_compression`) | Built in | Level |
|---|---|---|
| `CARQUET_COMPRESSION_UNCOMPRESSED` | yes | zero-copy reads possible |
| `CARQUET_COMPRESSION_SNAPPY` | yes (in-house) | none |
| `CARQUET_COMPRESSION_GZIP` | zlib | 1–9 |
| `CARQUET_COMPRESSION_LZ4_RAW` (codec 7) | liblz4 | prefer over LZ4 for new files |
| `CARQUET_COMPRESSION_LZ4` (codec 5, Hadoop-framed) | liblz4 | legacy readers |
| `CARQUET_COMPRESSION_ZSTD` | libzstd | 1 default (best speed/size); 3–9 smaller, slower |
| `CARQUET_COMPRESSION_BROTLI` | libbrotli (optional, `CARQUET_ENABLE_BROTLI`) | quality 0–11, default 8 (Arrow's); pyarrow / pandas `compression="brotli"` files read natively |
| `CARQUET_COMPRESSION_LZO` | no (GPL reference implementation; Arrow has none either) | plug in with `carquet_register_codec()` |

## Arrow C Data Interface format mapping (import and export)

| Arrow format | Parquet (import) | Export from Parquet |
|---|---|---|
| `b` | BOOLEAN | BOOLEAN → `b` |
| `c` `C` `s` `S` | INT32 + INTEGER(8/16, signed/unsigned) | INTEGER(8/16) → same codes |
| `i` `I` `l` `L` | INT32 / INT64 (+ INTEGER unsigned for `I`/`L`) | |
| `f` `g` | FLOAT / DOUBLE | |
| `e` | FLBA(2) FLOAT16 | |
| `u` `U` `z` `Z` | BYTE_ARRAY (STRING for `u`/`U`) | 32-bit offsets on export; `U`/`Z` origin recoverable via `carquet_reader_column_arrow_type_refinement()` |
| `w:N` | FLBA(N) | FLBA → `w:N` (UUID, DECIMAL on FLBA export as `w:N`) |
| `tdD` / `tdm` | DATE (date64 converted to days) | DATE → `tdD` |
| `ttm` `ttu` `ttn` / `tts` | TIME(MILLIS/MICROS/NANOS) (`tts` → MILLIS) | |
| `tsm:tz` `tsu:tz` `tsn:tz` / `tss:tz` | TIMESTAMP (seconds → MILLIS); non-empty tz sets `is_adjusted_to_utc` | |
| `d:p,s[,bits]` | FLBA DECIMAL(p, s), parquet-cpp byte length, LE→BE conversion | |
| `+s` | group (struct); top level must be `+s` | |
| `+l` `+L` `+w:n` | 3-level LIST | LIST → `+l` |
| `+m` | MAP (`key_value` struct child) | MAP → `+m` |
| dictionary-encoded | resolved to the value type | |
| `+ud` `+us` run-end | rejected (`NOT_IMPLEMENTED`) | |

Sliced arrays (`offset != 0`) are supported at every level on import. On export every buffer is an owned copy; release once.
