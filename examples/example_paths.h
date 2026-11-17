/**
 * @file example_paths.h
 * @brief Portable scratch-file paths for the example programs.
 *
 * The examples used to hard-code "/tmp/...", which does not exist on Windows
 * (it resolves to C:\tmp on the current drive). This picks the platform temp
 * directory instead: $TMPDIR (else /tmp) on POSIX, %TEMP% / %TMP% on Windows.
 */
#ifndef CARQUET_EXAMPLE_PATHS_H
#define CARQUET_EXAMPLE_PATHS_H

#include <stdio.h>
#include <stdlib.h>

static void example_temp_path(char* buf, size_t size, const char* name) {
#ifdef _WIN32
    const char* dir = getenv("TEMP");
    if (!dir || !*dir) dir = getenv("TMP");
    if (!dir || !*dir) dir = ".";
    snprintf(buf, size, "%s\\%s", dir, name);
#else
    const char* dir = getenv("TMPDIR");
    if (!dir || !*dir) dir = "/tmp";
    snprintf(buf, size, "%s/%s", dir, name);
#endif
}

#endif /* CARQUET_EXAMPLE_PATHS_H */
