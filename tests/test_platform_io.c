/**
 * @file test_platform_io.c
 * @brief Platform file-API regressions (first exercised on Windows).
 *
 * - UTF-8 paths with non-ASCII characters: the CRT fopen()/CreateFileA read a
 *   narrow path in the ANSI code page, so such a file could not be written,
 *   read (fread or mmap), appended to, or cleaned up by writer abort.
 * - An empty file opened with use_mmap: CreateFileMapping rejects 0-byte
 *   files; the open must fail cleanly (no leaked handles, proper error).
 * - A memory-mapped file can be deleted while the reader is open
 *   (FILE_SHARE_DELETE on Windows, unlink-while-open on POSIX).
 *
 * Runs on every platform: on POSIX the same calls exercise the native UTF-8
 * path handling.
 */

#include <carquet/carquet.h>
#include "test_helpers.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#endif

/* "carquet_" + e-acute + CJK "数据" + ".parquet", spelled as UTF-8 bytes so the
 * test does not depend on the compiler's source/execution charset. */
#define UTF8_NAME "carquet_\xC3\xA9\xE6\x95\xB0\xE6\x8D\xAE"

static char g_dir[1024];

/* Temp directory as UTF-8. On Windows %TEMP% via getenv() is ANSI, so take
 * the wide temp path and convert it. */
static const char* temp_dir_utf8(void) {
#ifdef _WIN32
    wchar_t wtmp[MAX_PATH + 1];
    DWORD n = GetTempPathW(MAX_PATH + 1, wtmp);
    if (n > 0 && n <= MAX_PATH &&
        WideCharToMultiByte(CP_UTF8, 0, wtmp, -1, g_dir, (int)sizeof(g_dir), NULL, NULL) > 0) {
        return g_dir;
    }
#endif
    snprintf(g_dir, sizeof(g_dir), "%s", carquet_test_get_temp_dir());
    return g_dir;
}

/* Existence check that does not go through the library under test. */
static int utf8_file_exists(const char* path) {
#ifdef _WIN32
    wchar_t wpath[2048];
    if (MultiByteToWideChar(CP_UTF8, 0, path, -1, wpath, 2048) <= 0) return 0;
    return GetFileAttributesW(wpath) != INVALID_FILE_ATTRIBUTES;
#else
    FILE* f = fopen(path, "rb");
    if (f) fclose(f);
    return f != NULL;
#endif
}

static int utf8_remove(const char* path) {
#ifdef _WIN32
    wchar_t wpath[2048];
    if (MultiByteToWideChar(CP_UTF8, 0, path, -1, wpath, 2048) <= 0) return -1;
    return DeleteFileW(wpath) ? 0 : -1;
#else
    return remove(path);
#endif
}

static carquet_schema_t* make_schema(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) return NULL;
    if (carquet_schema_add_column(schema, "v", CARQUET_PHYSICAL_INT64, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK) {
        carquet_schema_free(schema);
        return NULL;
    }
    return schema;
}

static int write_rows(const char* path, int64_t base, int n) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = make_schema();
    if (!schema) return -1;
    carquet_writer_t* w = carquet_writer_create(path, schema, NULL, &err);
    carquet_schema_free(schema);
    if (!w) {
        printf("  writer_create: %s\n", err.message);
        return -1;
    }
    int64_t vals[16];
    for (int i = 0; i < n; i++) vals[i] = base + i;
    if (carquet_writer_write_batch(w, 0, vals, n, NULL, NULL) != CARQUET_OK) {
        carquet_writer_abort(w);
        return -1;
    }
    return carquet_writer_close(w) == CARQUET_OK ? 0 : -1;
}

static int64_t count_rows(const char* path, bool use_mmap) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_options_t opts;
    carquet_reader_options_init(&opts);
    opts.use_mmap = use_mmap;
    carquet_reader_t* r = carquet_reader_open(path, &opts, &err);
    if (!r) {
        printf("  reader_open(mmap=%d): %s\n", (int)use_mmap, err.message);
        return -1;
    }
    if (use_mmap && !carquet_reader_is_mmap(r)) {
        printf("  reader_open(mmap=1) fell back to fread\n");
        carquet_reader_close(r);
        return -1;
    }
    int64_t rows = carquet_reader_num_rows(r);
    carquet_reader_close(r);
    return rows;
}

static int test_utf8_path_roundtrip(void) {
    char path[1200];
    snprintf(path, sizeof(path), "%s/" UTF8_NAME "_%d.parquet",
             temp_dir_utf8(), carquet_test_getpid());
    utf8_remove(path);

    if (write_rows(path, 0, 10) != 0) TEST_FAIL("utf8_path_roundtrip", "write failed");
    if (!utf8_file_exists(path)) {
        TEST_FAIL("utf8_path_roundtrip", "file not created under its UTF-8 name");
    }
    if (count_rows(path, false) != 10) {
        utf8_remove(path);
        TEST_FAIL("utf8_path_roundtrip", "fread reader failed");
    }
    if (count_rows(path, true) != 10) {
        utf8_remove(path);
        TEST_FAIL("utf8_path_roundtrip", "mmap reader failed");
    }

    /* Append through the same path. */
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = make_schema();
    carquet_writer_t* w = carquet_writer_open_append(path, schema, NULL, &err);
    carquet_schema_free(schema);
    if (!w) {
        utf8_remove(path);
        printf("  open_append: %s\n", err.message);
        TEST_FAIL("utf8_path_roundtrip", "append open failed");
    }
    int64_t more[3] = {100, 101, 102};
    if (carquet_writer_write_batch(w, 0, more, 3, NULL, NULL) != CARQUET_OK ||
        carquet_writer_close(w) != CARQUET_OK) {
        utf8_remove(path);
        TEST_FAIL("utf8_path_roundtrip", "append write failed");
    }
    if (count_rows(path, false) != 13) {
        utf8_remove(path);
        TEST_FAIL("utf8_path_roundtrip", "wrong row count after append");
    }

    utf8_remove(path);
    TEST_PASS("utf8_path_roundtrip");
    return 0;
}

static int test_utf8_path_abort_removes_file(void) {
    char path[1200];
    snprintf(path, sizeof(path), "%s/" UTF8_NAME "_abort_%d.parquet",
             temp_dir_utf8(), carquet_test_getpid());
    utf8_remove(path);

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = make_schema();
    carquet_writer_t* w = carquet_writer_create(path, schema, NULL, &err);
    carquet_schema_free(schema);
    if (!w) TEST_FAIL("utf8_path_abort", "writer_create failed");
    if (!utf8_file_exists(path)) {
        carquet_writer_abort(w);
        TEST_FAIL("utf8_path_abort", "file not created under its UTF-8 name");
    }
    carquet_writer_abort(w);
    if (utf8_file_exists(path)) {
        utf8_remove(path);
        TEST_FAIL("utf8_path_abort", "abort left the partial file behind");
    }
    TEST_PASS("utf8_path_abort");
    return 0;
}

static int test_mmap_empty_file(void) {
    char path[1200];
    carquet_test_temp_path(path, sizeof(path), "platform_empty");
    FILE* f = fopen(path, "wb");
    if (!f) TEST_FAIL("mmap_empty_file", "cannot create empty file");
    fclose(f);

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_options_t opts;
    carquet_reader_options_init(&opts);
    opts.use_mmap = true;
    carquet_reader_t* r = carquet_reader_open(path, &opts, &err);
    if (r) {
        carquet_reader_close(r);
        remove(path);
        TEST_FAIL("mmap_empty_file", "empty file opened successfully");
    }
    if (err.code != CARQUET_ERROR_INVALID_FOOTER) {
        printf("  got %d (%s)\n", (int)err.code, err.message);
        remove(path);
        TEST_FAIL("mmap_empty_file", "expected INVALID_FOOTER");
    }
    /* No handle may be left open: the file must be removable right away. */
    if (remove(path) != 0) TEST_FAIL("mmap_empty_file", "file still locked after failed open");
    TEST_PASS("mmap_empty_file");
    return 0;
}

static int test_mmap_file_deletable_while_open(void) {
    char path[1200];
    carquet_test_temp_path(path, sizeof(path), "platform_unlink");
    remove(path);
    if (write_rows(path, 7, 5) != 0) TEST_FAIL("mmap_deletable", "write failed");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_options_t opts;
    carquet_reader_options_init(&opts);
    opts.use_mmap = true;
    carquet_reader_t* r = carquet_reader_open(path, &opts, &err);
    if (!r || !carquet_reader_is_mmap(r)) {
        if (r) carquet_reader_close(r);
        remove(path);
        TEST_FAIL("mmap_deletable", "mmap open failed");
    }
    int rc = remove(path);
    /* The mapping stays readable after the name is gone. */
    int64_t rows = carquet_reader_num_rows(r);
    carquet_reader_close(r);
    if (rc != 0) {
        remove(path);
        TEST_FAIL("mmap_deletable", "could not delete a memory-mapped file");
    }
    if (rows != 5) TEST_FAIL("mmap_deletable", "wrong row count");
    TEST_PASS("mmap_deletable");
    return 0;
}

static int test_buffer_writer(void) {
    /* carquet_writer_create_buffer() goes through tmpfile(), which on Windows
     * may need write access to the drive root; the fallback must work. */
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = make_schema();
    carquet_writer_t* w = carquet_writer_create_buffer(schema, NULL, &err);
    carquet_schema_free(schema);
    if (!w) TEST_FAIL("buffer_writer", err.message);
    int64_t vals[4] = {1, 2, 3, 4};
    if (carquet_writer_write_batch(w, 0, vals, 4, NULL, NULL) != CARQUET_OK ||
        carquet_writer_close(w) != CARQUET_OK) {
        carquet_writer_abort(w);
        TEST_FAIL("buffer_writer", "write/close failed");
    }
    void* buf = NULL;
    size_t size = 0;
    /* On success get_buffer() also releases the writer. */
    if (carquet_writer_get_buffer(w, &buf, &size) != CARQUET_OK || size < 12) {
        carquet_writer_abort(w);
        TEST_FAIL("buffer_writer", "get_buffer failed");
    }
    carquet_reader_t* r = carquet_reader_open_buffer(buf, size, NULL, &err);
    int64_t rows = r ? carquet_reader_num_rows(r) : -1;
    if (r) carquet_reader_close(r);
    const carquet_allocator_t* a = carquet_get_allocator();
    a->free(buf, a->ctx);
    if (rows != 4) TEST_FAIL("buffer_writer", "buffer did not round-trip");
    TEST_PASS("buffer_writer");
    return 0;
}

int main(void) {
    printf("=== Platform I/O Tests ===\n\n");
    carquet_init();

    int failures = 0;
    failures += test_utf8_path_roundtrip();
    failures += test_utf8_path_abort_removes_file();
    failures += test_mmap_empty_file();
    failures += test_mmap_file_deletable_while_open();
    failures += test_buffer_writer();

    printf("\n%s (%d failure(s))\n", failures ? "FAILED" : "All platform I/O tests passed", failures);
    return failures ? 1 : 0;
}
