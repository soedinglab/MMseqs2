/**
 * @file file_io.c
 * @brief UTF-8 path support for Windows file APIs and positional reads
 *        (see file_io.h).
 */

/* pread()/fileno() are POSIX; strict -std=c11 headers hide them otherwise. */
#if !defined(_WIN32) && !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200809L
#endif

#include "file_io.h"

#ifdef _WIN32

#include "allocator.h"
#include <windows.h>
#include <io.h>
#include <string.h>

wchar_t* carquet_path_to_wide(const char* path) {
    if (!path) return NULL;

    /* Pure ASCII reads the same in every ANSI code page: keep the narrow API. */
    const unsigned char* p = (const unsigned char*)path;
    while (*p && *p < 0x80) p++;
    if (!*p) return NULL;

    /* MB_ERR_INVALID_CHARS: a path that is not well-formed UTF-8 is most
     * likely ANSI-encoded; let the caller fall back to the narrow API. */
    int n = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path, -1, NULL, 0);
    if (n <= 0) return NULL;

    wchar_t* wide = (wchar_t*)carquet_mem_malloc((size_t)n * sizeof(wchar_t));
    if (!wide) return NULL;
    if (MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path, -1, wide, n) != n) {
        carquet_mem_free(wide);
        return NULL;
    }
    return wide;
}

FILE* carquet_fopen(const char* path, const char* mode) {
    wchar_t* wpath = carquet_path_to_wide(path);
    if (!wpath) return fopen(path, mode);

    wchar_t wmode[16];
    size_t i = 0;
    for (; mode[i] && i + 1 < sizeof(wmode) / sizeof(wmode[0]); i++) {
        wmode[i] = (wchar_t)(unsigned char)mode[i];
    }
    wmode[i] = L'\0';

    FILE* f = _wfopen(wpath, wmode);
    carquet_mem_free(wpath);
    return f;
}

int carquet_remove(const char* path) {
    wchar_t* wpath = carquet_path_to_wide(path);
    if (!wpath) return remove(path);
    int rc = _wremove(wpath);
    carquet_mem_free(wpath);
    return rc;
}

static HANDLE file_handle(FILE* file) {
    int fd = file ? _fileno(file) : -1;
    if (fd < 0) return INVALID_HANDLE_VALUE;
    return (HANDLE)_get_osfhandle(fd);
}

bool carquet_file_supports_pread(FILE* file) {
    HANDLE h = file_handle(file);
    return h != INVALID_HANDLE_VALUE && GetFileType(h) == FILE_TYPE_DISK;
}

size_t carquet_file_pread(FILE* file, void* buf, size_t size, int64_t offset) {
    HANDLE h = file_handle(file);
    if (h == INVALID_HANDLE_VALUE || offset < 0) return 0;
    size_t done = 0;
    while (done < size) {
        size_t want = size - done;
        DWORD chunk = want > 0x40000000u ? 0x40000000u : (DWORD)want;
        uint64_t pos = (uint64_t)offset + done;
        OVERLAPPED ov;
        memset(&ov, 0, sizeof(ov));
        ov.Offset = (DWORD)(pos & 0xFFFFFFFFu);
        ov.OffsetHigh = (DWORD)(pos >> 32);
        DWORD got = 0;
        if (!ReadFile(h, (char*)buf + done, chunk, &got, &ov) || got == 0) break;
        done += got;
    }
    return done;
}

#else

#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

bool carquet_file_supports_pread(FILE* file) {
    int fd = file ? fileno(file) : -1;
    struct stat st;
    return fd >= 0 && fstat(fd, &st) == 0 && S_ISREG(st.st_mode);
}

size_t carquet_file_pread(FILE* file, void* buf, size_t size, int64_t offset) {
    int fd = file ? fileno(file) : -1;
    if (fd < 0 || offset < 0) return 0;
    size_t done = 0;
    while (done < size) {
        int64_t pos = offset + (int64_t)done;
        if ((int64_t)(off_t)pos != pos) break;   /* 32-bit off_t */
        ssize_t n = pread(fd, (char*)buf + done, size - done, (off_t)pos);
        if (n < 0) {
            if (errno == EINTR) continue;
            break;
        }
        if (n == 0) break;
        done += (size_t)n;
    }
    return done;
}

#endif
