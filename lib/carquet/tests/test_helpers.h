/**
 * @file test_helpers.h
 * @brief Portable test helper utilities for carquet tests
 *
 * Provides cross-platform utilities for test file paths and common test macros.
 */

#ifndef CARQUET_TEST_HELPERS_H
#define CARQUET_TEST_HELPERS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Platform-specific includes */
#ifdef _WIN32
#include <process.h>
#define carquet_test_getpid _getpid
#else
#include <unistd.h>
#define carquet_test_getpid getpid
#endif

/* Test macros */
#define TEST_PASS(name) printf("[PASS] %s\n", name)
#define TEST_FAIL(name, msg) do { printf("[FAIL] %s: %s\n", name, msg); return 1; } while(0)

/**
 * Get the platform-appropriate temp directory.
 * Returns /tmp on Unix, %TEMP% or %TMP% on Windows.
 */
static inline const char* carquet_test_get_temp_dir(void) {
#ifdef _WIN32
    static char temp_dir[512] = {0};
    if (temp_dir[0] == 0) {
        const char* tmp = getenv("TEMP");
        if (!tmp) tmp = getenv("TMP");
        if (!tmp) tmp = ".";
        snprintf(temp_dir, sizeof(temp_dir), "%s", tmp);
    }
    return temp_dir;
#else
    return "/tmp";
#endif
}

/**
 * Build a temp file path with the given base name.
 * Uses process ID to avoid conflicts between concurrent test runs.
 *
 * @param buffer Output buffer for the path
 * @param buffer_size Size of the output buffer
 * @param base_name Base name for the temp file (without .parquet extension)
 */
static inline void carquet_test_temp_path(char* buffer, size_t buffer_size, const char* base_name) {
    snprintf(buffer, buffer_size, "%s/carquet_%s_%d.parquet",
             carquet_test_get_temp_dir(), base_name, carquet_test_getpid());
}

/**
 * Build a temp file path with a custom extension.
 */
static inline void carquet_test_temp_path_ext(char* buffer, size_t buffer_size,
                                               const char* base_name, const char* ext) {
    snprintf(buffer, buffer_size, "%s/carquet_%s_%d.%s",
             carquet_test_get_temp_dir(), base_name, carquet_test_getpid(), ext);
}

/**
 * Remove a test file (cleanup helper).
 */
static inline void carquet_test_cleanup(const char* path) {
    remove(path);
}

#endif /* CARQUET_TEST_HELPERS_H */
