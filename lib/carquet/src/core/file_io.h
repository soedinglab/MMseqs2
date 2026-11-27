/**
 * @file file_io.h
 * @brief Path-taking file helpers that accept UTF-8 paths on Windows.
 *
 * The CRT's fopen()/remove() and the Win32 "A" APIs interpret a narrow path in
 * the process ANSI code page, so a UTF-8 path with non-ASCII characters (the
 * encoding every other platform, and most callers, use) fails to open or
 * names the wrong file. These helpers convert a path that is valid UTF-8 to
 * UTF-16 and use the wide API; a path that is not valid UTF-8 (a legacy
 * ANSI-encoded path) is passed through to the narrow API unchanged, so
 * existing callers keep working. On non-Windows platforms they are thin
 * wrappers over fopen()/remove().
 */
#ifndef CARQUET_CORE_FILE_IO_H
#define CARQUET_CORE_FILE_IO_H

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>

/* Positional reads on a FILE's descriptor: pread() on POSIX, ReadFile() with
 * an explicit offset on Windows. They neither use nor move the stream
 * position, so several threads may read one file concurrently, which
 * fseek()+fread() on a shared FILE cannot do. carquet_file_supports_pread()
 * is true only for a regular disk file with a descriptor (not a pipe,
 * terminal or memory stream); carquet_file_pread() returns the bytes read,
 * short only at end of file or on error. */
bool carquet_file_supports_pread(FILE* file);
size_t carquet_file_pread(FILE* file, void* buf, size_t size, int64_t offset);

#ifdef _WIN32
#include <wchar.h>

FILE* carquet_fopen(const char* path, const char* mode);
int carquet_remove(const char* path);

/* UTF-16 copy of `path` when it contains non-ASCII bytes and is valid UTF-8;
 * NULL otherwise (use the narrow API). Free with carquet_mem_free(). */
wchar_t* carquet_path_to_wide(const char* path);

#else

static inline FILE* carquet_fopen(const char* path, const char* mode) {
    return fopen(path, mode);
}

static inline int carquet_remove(const char* path) {
    return remove(path);
}

#endif

#endif /* CARQUET_CORE_FILE_IO_H */
