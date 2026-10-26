/**
 * @file nested_data.c
 * @brief Write and read nested LIST and MAP columns
 *
 * There is no row-object layer in the C API: even for nested data you write
 * and read leaf columns plus their definition/repetition level streams.
 *
 *   - definition level: "how much of this path exists?" (null vs empty vs present)
 *   - repetition level: "did this value continue the current list, or start a new row?"
 *
 * Values passed to write_batch are dense (only the present leaf values); the
 * def/rep arrays carry one entry per logical slot. This example writes a
 * list<int32> with empty/null cases and a map<string,string>, then reads both
 * back and reconstructs the original rows.
 */

#include <carquet/carquet.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(expr) do { if ((expr) != CARQUET_OK) { \
    fprintf(stderr, "FAIL at %s:%d\n", __FILE__, __LINE__); return 1; } } while(0)

/* Leaf column indices in the schema built below. */
enum { COL_SCORES = 0, COL_TAGS_KEY = 1, COL_TAGS_VALUE = 2 };

static carquet_schema_t* make_schema(void) {
    carquet_schema_t* s = carquet_schema_create(NULL);
    /* list<int32>, list and element both OPTIONAL (so we can express null
     * lists, empty lists, and null elements). */
    carquet_schema_add_list(s, "scores", CARQUET_PHYSICAL_INT32, NULL,
                            CARQUET_REPETITION_OPTIONAL, 0, 0);
    /* map<string,string>, map OPTIONAL. */
    carquet_schema_add_map(s, "tags",
                           CARQUET_PHYSICAL_BYTE_ARRAY, NULL, 0,   /* key   */
                           CARQUET_PHYSICAL_BYTE_ARRAY, NULL, 0,   /* value */
                           CARQUET_REPETITION_OPTIONAL, 0);
    return s;
}

static carquet_byte_array_t ba(const char* s) {
    carquet_byte_array_t b;
    b.data = (uint8_t*)s;
    b.length = (int64_t)strlen(s);
    return b;
}

static int write_file(const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = make_schema();

    /* Report the max levels carquet derived for each leaf — these are what the
     * hand-written def/rep arrays below are built against. */
    printf("  scores  leaf: max_def=%d max_rep=%d\n",
           carquet_schema_max_def_level(schema, COL_SCORES),
           carquet_schema_max_rep_level(schema, COL_SCORES));
    printf("  tags.key   : max_def=%d max_rep=%d\n",
           carquet_schema_max_def_level(schema, COL_TAGS_KEY),
           carquet_schema_max_rep_level(schema, COL_TAGS_KEY));
    printf("  tags.value : max_def=%d max_rep=%d\n",
           carquet_schema_max_def_level(schema, COL_TAGS_VALUE),
           carquet_schema_max_rep_level(schema, COL_TAGS_VALUE));

    carquet_writer_t* w = carquet_writer_create(path, schema, NULL, &err);
    if (!w) { fprintf(stderr, "create: %s\n", err.message); carquet_schema_free(schema); return 1; }

    /* scores, 4 logical rows:
     *   row0 = [10, 20]
     *   row1 = []            (empty list)
     *   row2 = NULL          (null list)
     *   row3 = [30, NULL, 50]
     * def: 3=element present, 2=null element, 1=empty list, 0=null list.
     * rep: 0=new row, 1=continue current list. */
    int32_t scores[]    = {10, 20, 30, 50};            /* dense: only def==3 slots */
    int16_t scores_def[] = {3, 3, 1, 0, 3, 2, 3};
    int16_t scores_rep[] = {0, 1, 0, 0, 0, 1, 1};
    CHECK(carquet_writer_write_batch(w, COL_SCORES, scores,
                                     (int64_t)(sizeof(scores_def) / sizeof(scores_def[0])),
                                     scores_def, scores_rep));

    /* tags, same 4 rows, all maps present and non-empty:
     *   row0 = {a:1}
     *   row1 = {b:2}
     *   row2 = {c:3}
     *   row3 = {d:4, e:5}   */
    carquet_byte_array_t keys[]   = { ba("a"), ba("b"), ba("c"), ba("d"), ba("e") };
    carquet_byte_array_t values[] = { ba("1"), ba("2"), ba("3"), ba("4"), ba("5") };
    /* key leaf: max_def=2 (all keys present); value leaf: max_def=3. */
    int16_t key_def[]   = {2, 2, 2, 2, 2};
    int16_t key_rep[]   = {0, 0, 0, 0, 1};
    int16_t value_def[] = {3, 3, 3, 3, 3};
    int16_t value_rep[] = {0, 0, 0, 0, 1};
    CHECK(carquet_writer_write_batch(w, COL_TAGS_KEY,   keys,   5, key_def,   key_rep));
    CHECK(carquet_writer_write_batch(w, COL_TAGS_VALUE, values, 5, value_def, value_rep));

    CHECK(carquet_writer_close(w));
    carquet_schema_free(schema);
    return 0;
}

static int read_scores(carquet_reader_t* r) {
    printf("\n=== scores: list<int32> ===\n");
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_column_reader_t* col = carquet_reader_get_column(r, 0, COL_SCORES, &err);
    if (!col) { fprintf(stderr, "get_column: %s\n", err.message); return 1; }

    int32_t values[64];
    int16_t def[64], rep[64];
    int64_t n = carquet_column_read_batch(col, values, 64, def, rep);
    if (n < 0) { fprintf(stderr, "read_batch failed\n"); carquet_column_reader_free(col); return 1; }

    /* Walk slots: rep==0 starts a new row; def selects null/empty/null-elem/value.
     * `open` tracks whether the current row has an unterminated "[". */
    int vi = 0;       /* index into dense values (advances only on def==3) */
    int row = -1;
    bool open = false;
    for (int64_t i = 0; i < n; i++) {
        if (rep[i] == 0) {
            if (open) printf("]");          /* close previous non-empty list */
            if (row >= 0) printf("\n");
            open = false;
            row++;
            printf("  row%d = ", row);
            if (def[i] == 0) { printf("NULL"); continue; }   /* null list  */
            if (def[i] == 1) { printf("[]");   continue; }   /* empty list */
            printf("[");
            open = true;
        } else {
            printf(", ");
        }
        if (def[i] == 2)      printf("NULL");        /* null element */
        else if (def[i] == 3) printf("%d", values[vi++]);
    }
    if (open) printf("]");
    printf("\n");

    carquet_column_reader_free(col);
    return 0;
}

static int read_tags(carquet_reader_t* r) {
    printf("\n=== tags: map<string,string> ===\n");
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_column_reader_t* kcol = carquet_reader_get_column(r, 0, COL_TAGS_KEY,   &err);
    carquet_column_reader_t* vcol = carquet_reader_get_column(r, 0, COL_TAGS_VALUE, &err);
    if (!kcol || !vcol) { fprintf(stderr, "get_column tags failed\n"); return 1; }

    carquet_byte_array_t keys[64], vals[64];
    int16_t kdef[64], krep[64], vdef[64], vrep[64];
    int64_t nk = carquet_column_read_batch(kcol, keys, 64, kdef, krep);
    int64_t nv = carquet_column_read_batch(vcol, vals, 64, vdef, vrep);
    if (nk < 0 || nv < 0) { fprintf(stderr, "read tags failed\n"); goto done; }

    /* Keys and values share the same key_value repetition, so the rep stream
     * lines up entry-for-entry; rep==0 opens a new row. */
    int row = -1;
    for (int64_t i = 0; i < nk; i++) {
        if (krep[i] == 0) {
            if (row >= 0) printf(" }\n");
            row++;
            printf("  row%d = {", row);
        } else {
            printf(",");
        }
        printf(" %.*s:%.*s",
               (int)keys[i].length, (const char*)keys[i].data,
               (int)vals[i].length, (const char*)vals[i].data);
    }
    printf(" }\n");

done:
    carquet_column_reader_free(kcol);
    carquet_column_reader_free(vcol);
    return 0;
}

int main(void) {
    carquet_init();
    const char* path = "/tmp/carquet_nested_example.parquet";

    printf("Writing nested file...\n");
    if (write_file(path)) return 1;

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) { fprintf(stderr, "open: %s\n", err.message); return 1; }
    printf("\n%lld logical rows, %d leaf columns\n",
           (long long)carquet_reader_num_rows(r), carquet_reader_num_columns(r));

    if (read_scores(r)) { carquet_reader_close(r); return 1; }
    if (read_tags(r))   { carquet_reader_close(r); return 1; }

    carquet_reader_close(r);
    remove(path);
    printf("\nDone.\n");
    return 0;
}
