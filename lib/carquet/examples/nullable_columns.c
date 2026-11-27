/**
 * @file nullable_columns.c
 * @brief Example demonstrating optional/nullable columns with definition levels
 *
 * This example demonstrates:
 * - OPTIONAL (nullable) columns
 * - Definition levels to track null values
 * - Reading back null values correctly
 *
 * Parquet Definition Levels:
 * - For REQUIRED columns: No definition levels needed
 * - For OPTIONAL columns: def_level = 0 means NULL, def_level = 1 means value present
 * - For nested structures: Higher levels indicate deeper nesting
 *
 * Build:
 *   gcc -o nullable_columns nullable_columns.c -lcarquet -I../include
 *
 * Run:
 *   ./nullable_columns
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <carquet/carquet.h>
#include "example_paths.h"

#define NUM_ROWS 20

/**
 * Create schema with both required and optional columns
 */
static carquet_schema_t* create_schema_with_nullables(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) {
        fprintf(stderr, "Failed to create schema: %s\n", err.message);
        return NULL;
    }

    carquet_logical_type_t string_type = { .id = CARQUET_LOGICAL_STRING };

    /* Required column (always has a value) */
(void)carquet_schema_add_column(schema, "id",
        CARQUET_PHYSICAL_INT64, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    /* Required string column */
(void)carquet_schema_add_column(schema, "name",
        CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);

    /* Optional integer column (can be NULL) */
(void)carquet_schema_add_column(schema, "age",
        CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_OPTIONAL, 0, 0);

    /* Optional double column (can be NULL) */
(void)carquet_schema_add_column(schema, "score",
        CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_OPTIONAL, 0, 0);

    /* Optional string column (can be NULL) */
(void)carquet_schema_add_column(schema, "email",
        CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_OPTIONAL, 0, 0);

    printf("Schema created:\n");
    printf("  - id: INT64 (REQUIRED)\n");
    printf("  - name: STRING (REQUIRED)\n");
    printf("  - age: INT32 (OPTIONAL - can be NULL)\n");
    printf("  - score: DOUBLE (OPTIONAL - can be NULL)\n");
    printf("  - email: STRING (OPTIONAL - can be NULL)\n\n");

    return schema;
}

/**
 * Write data with some NULL values
 */
static int write_nullable_data(const char* filename, carquet_schema_t* schema) {
    printf("Writing data with NULL values to: %s\n\n", filename);

    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_SNAPPY;

    carquet_writer_t* writer = carquet_writer_create(filename, schema, &opts, &err);
    if (!writer) {
        fprintf(stderr, "Failed to create writer: %s\n", err.message);
        return -1;
    }

    /* Allocate data arrays */
    int64_t ids[NUM_ROWS];
    carquet_byte_array_t names[NUM_ROWS];
    char name_buffers[NUM_ROWS][32];

    int32_t ages[NUM_ROWS];
    int16_t age_def_levels[NUM_ROWS];

    double scores[NUM_ROWS];
    int16_t score_def_levels[NUM_ROWS];

    carquet_byte_array_t emails[NUM_ROWS];
    char email_buffers[NUM_ROWS][64];
    int16_t email_def_levels[NUM_ROWS];

    printf("Data being written:\n");
    printf("%-4s %-12s %-6s %-8s %-25s\n", "ID", "Name", "Age", "Score", "Email");
    printf("%-4s %-12s %-6s %-8s %-25s\n", "--", "----", "---", "-----", "-----");

    int age_value_count = 0;
    int score_value_count = 0;
    int email_value_count = 0;

    for (int i = 0; i < NUM_ROWS; i++) {
        /* Required fields - always have values */
        ids[i] = i + 1;

        snprintf(name_buffers[i], sizeof(name_buffers[i]), "User_%02d", i + 1);
        names[i].data = (uint8_t*)name_buffers[i];
        names[i].length = (int32_t)strlen(name_buffers[i]);

        /* Age: NULL for every 3rd row */
        if (i % 3 == 2) {
            age_def_levels[i] = 0;  /* NULL - no value stored */
        } else {
            ages[age_value_count] = 20 + (i % 50);
            age_def_levels[i] = 1;  /* Value present */
            age_value_count++;
        }

        /* Score: NULL for every 4th row */
        if (i % 4 == 3) {
            score_def_levels[i] = 0;  /* NULL */
        } else {
            scores[score_value_count] = 75.0 + (i % 25) * 1.0;
            score_def_levels[i] = 1;  /* Value present */
            score_value_count++;
        }

        /* Email: NULL for every 5th row */
        if (i % 5 == 4) {
            email_def_levels[i] = 0;  /* NULL */
        } else {
            snprintf(email_buffers[email_value_count], sizeof(email_buffers[0]),
                     "user%02d@example.com", i + 1);
            emails[email_value_count].data = (uint8_t*)email_buffers[email_value_count];
            emails[email_value_count].length = (int32_t)strlen(email_buffers[email_value_count]);
            email_def_levels[i] = 1;
            email_value_count++;
        }

        /* Print the row for visualization */
        printf("%-4lld %-12s ", (long long)ids[i], name_buffers[i]);

        if (age_def_levels[i] == 0) {
            printf("%-6s ", "NULL");
        } else {
            printf("%-6d ", ages[age_value_count - 1]);
        }

        if (score_def_levels[i] == 0) {
            printf("%-8s ", "NULL");
        } else {
            printf("%-8.1f ", scores[score_value_count - 1]);
        }

        if (email_def_levels[i] == 0) {
            printf("%-25s", "NULL");
        } else {
            printf("%-25s", email_buffers[email_value_count - 1]);
        }
        printf("\n");
    }

    printf("\nValue counts:\n");
    printf("  ages: %d values, %d NULLs\n", age_value_count, NUM_ROWS - age_value_count);
    printf("  scores: %d values, %d NULLs\n", score_value_count, NUM_ROWS - score_value_count);
    printf("  emails: %d values, %d NULLs\n", email_value_count, NUM_ROWS - email_value_count);

    carquet_status_t status;

    /* Write required columns (no def_levels needed) */
    status = carquet_writer_write_batch(writer, 0, ids, NUM_ROWS, NULL, NULL);
    if (status != CARQUET_OK) goto error;

    status = carquet_writer_write_batch(writer, 1, names, NUM_ROWS, NULL, NULL);
    if (status != CARQUET_OK) goto error;

    /* Write optional columns WITH definition levels */
    /* Note: Only non-null values are in the values array, but num_values
     * is the total number of logical rows (including NULLs). The def_levels
     * array has one entry per logical row indicating null vs present. */
    status = carquet_writer_write_batch(writer, 2, ages, NUM_ROWS,
                                        age_def_levels, NULL);
    if (status != CARQUET_OK) goto error;

    status = carquet_writer_write_batch(writer, 3, scores, NUM_ROWS,
                                        score_def_levels, NULL);
    if (status != CARQUET_OK) goto error;

    status = carquet_writer_write_batch(writer, 4, emails, NUM_ROWS,
                                        email_def_levels, NULL);
    if (status != CARQUET_OK) goto error;

    status = carquet_writer_close(writer);
    if (status != CARQUET_OK) {
        fprintf(stderr, "Failed to close writer\n");
        return -1;
    }

    printf("\nSuccessfully wrote %d rows\n", NUM_ROWS);
    return 0;

error:
    fprintf(stderr, "Failed to write data\n");
    carquet_writer_abort(writer);
    return -1;
}

/**
 * Read data and interpret NULL values
 */
static int read_nullable_data(const char* filename) {
    printf("\nReading data with NULL values from: %s\n\n", filename);

    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_reader_t* reader = carquet_reader_open(filename, NULL, &err);
    if (!reader) {
        fprintf(stderr, "Failed to open reader: %s\n", err.message);
        return -1;
    }

    int64_t num_rows = carquet_reader_num_rows(reader);
    printf("File contains %lld rows\n\n", (long long)num_rows);

    /* Read from first row group */
    enum { BATCH_SIZE = 10 };

    printf("First %d rows:\n", BATCH_SIZE);
    printf("%-4s %-12s %-6s %-8s %-25s\n", "ID", "Name", "Age", "Score", "Email");
    printf("%-4s %-12s %-6s %-8s %-25s\n", "--", "----", "---", "-----", "-----");

    /* Read IDs (required) */
    int64_t ids[BATCH_SIZE];
    carquet_column_reader_t* id_col = carquet_reader_get_column(reader, 0, 0, &err);
    int64_t id_count = 0;
    if (id_col) {
        id_count = carquet_column_read_batch(id_col, ids, BATCH_SIZE, NULL, NULL);
        carquet_column_reader_free(id_col);
    }

    /* Read ages (optional) with definition levels */
    int32_t ages[BATCH_SIZE];
    int16_t age_def[BATCH_SIZE];
    carquet_column_reader_t* age_col = carquet_reader_get_column(reader, 0, 2, &err);
    int64_t age_count = 0;
    if (age_col) {
        age_count = carquet_column_read_batch(age_col, ages, BATCH_SIZE, age_def, NULL);
        carquet_column_reader_free(age_col);
    }

    /* Read scores (optional) with definition levels */
    double scores[BATCH_SIZE];
    int16_t score_def[BATCH_SIZE];
    carquet_column_reader_t* score_col = carquet_reader_get_column(reader, 0, 3, &err);
    int64_t score_count = 0;
    if (score_col) {
        score_count = carquet_column_read_batch(score_col, scores, BATCH_SIZE, score_def, NULL);
        carquet_column_reader_free(score_col);
    }

    /* Display with NULL handling */
    int age_idx = 0;
    int score_idx = 0;

    for (int i = 0; i < id_count; i++) {
        printf("%-4lld %-12s ", (long long)ids[i], "...");  /* Name omitted for brevity */

        /* Age: check definition level */
        if (i < age_count && age_def[i] > 0) {
            printf("%-6d ", ages[age_idx++]);
        } else {
            printf("%-6s ", "NULL");
        }

        /* Score: check definition level */
        if (i < score_count && score_def[i] > 0) {
            printf("%-8.1f ", scores[score_idx++]);
        } else {
            printf("%-8s ", "NULL");
        }

        printf("%-25s\n", "...");  /* Email omitted */
    }

    carquet_reader_close(reader);

    printf("\nSuccessfully read nullable data\n");
    return 0;
}

/**
 * Explain definition levels concept
 */
static void explain_definition_levels(void) {
    printf("\n=== Understanding Definition Levels ===\n\n");

    printf("Definition levels track the 'depth' at which a value is defined.\n");
    printf("For a simple OPTIONAL column:\n");
    printf("  - def_level = 0: Value is NULL\n");
    printf("  - def_level = 1: Value is present\n\n");

    printf("For nested structures (e.g., OPTIONAL struct with OPTIONAL field):\n");
    printf("  - def_level = 0: Outer struct is NULL\n");
    printf("  - def_level = 1: Outer struct present, inner field is NULL\n");
    printf("  - def_level = 2: Both outer struct and inner field are present\n\n");

    printf("When writing, only provide non-NULL values in the values array\n");
    printf("(sparse encoding). The definition levels array has one entry per\n");
    printf("logical row. num_values = number of logical rows.\n\n");

    printf("Example:\n");
    printf("  Logical rows:  [10, NULL, 20, NULL, 30]\n");
    printf("  Values array:  [10, 20, 30]       (3 non-null values, packed)\n");
    printf("  Def levels:    [1, 0, 1, 0, 1]    (5 entries, one per row)\n");
    printf("  num_values:    5                   (logical row count)\n\n");
}

int main(int argc, char* argv[]) {
    char default_path[1024];
    example_temp_path(default_path, sizeof(default_path), "example_nullable.parquet");
    const char* filename = default_path;

    printf("=== Carquet Nullable Columns Example ===\n");
    printf("Library version: %s\n\n", carquet_version());

    if (argc > 1) {
        filename = argv[1];
    }

    /* Explain the concepts */
    explain_definition_levels();

    /* Create schema */
    carquet_schema_t* schema = create_schema_with_nullables();
    if (!schema) {
        return 1;
    }

    /* Write data with NULLs */
    if (write_nullable_data(filename, schema) != 0) {
        carquet_schema_free(schema);
        return 1;
    }

    carquet_schema_free(schema);

    /* Read data back */
    if (read_nullable_data(filename) != 0) {
        return 1;
    }

    printf("\n=== Nullable columns example completed ===\n");

    /* Clean up */
    if (argc <= 1) {
        remove(filename);
        printf("(Removed temporary file)\n");
    }

    return 0;
}
