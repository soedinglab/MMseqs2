/**
 * @file bench_purge_win32.h
 * @brief Evict one file from the Windows file cache (benchmarks only).
 *
 * Windows has no per-file equivalent of posix_fadvise(DONTNEED). As on macOS,
 * the file is copied to a new path with unbuffered I/O (FILE_FLAG_NO_BUFFERING
 * on both handles, so neither the source read nor the copy enters the cache)
 * and the copy replaces the original. Deleting the original drops its cached
 * pages, and the copy was never cached, so the next open reads from storage.
 *
 * Shared by benchmark_carquet.c and benchmark_arrow_cpp.cpp so both libraries
 * read under the same conditions. Include after <windows.h>.
 */
#ifndef CARQUET_BENCH_PURGE_WIN32_H
#define CARQUET_BENCH_PURGE_WIN32_H

#ifdef _WIN32

#include <stdio.h>
#include <string.h>

static void bench_purge_file_cache_win32(const char* filename) {
    /* Unbuffered I/O needs sector-aligned offsets, sizes and buffers; 4 KiB
     * covers 512-byte and 4K-native disks, VirtualAlloc returns page-aligned
     * memory. */
    enum { BENCH_ALIGN = 4096, BENCH_CHUNK = 1 << 20 };
    char tmp[1024];
    HANDLE src, dst;
    LARGE_INTEGER size;
    void* buf;
    int ok = 1;

    if (snprintf(tmp, sizeof(tmp), "%s.nocache", filename) >= (int)sizeof(tmp)) return;

    src = CreateFileA(filename, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING,
                      FILE_FLAG_NO_BUFFERING | FILE_FLAG_SEQUENTIAL_SCAN, NULL);
    if (src == INVALID_HANDLE_VALUE) return;
    if (!GetFileSizeEx(src, &size)) { CloseHandle(src); return; }

    dst = CreateFileA(tmp, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS,
                      FILE_FLAG_NO_BUFFERING | FILE_FLAG_WRITE_THROUGH, NULL);
    if (dst == INVALID_HANDLE_VALUE) { CloseHandle(src); return; }

    buf = VirtualAlloc(NULL, BENCH_CHUNK, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    if (!buf) ok = 0;

    while (ok) {
        DWORD got = 0, put = 0, padded;
        if (!ReadFile(src, buf, BENCH_CHUNK, &got, NULL)) { ok = 0; break; }
        if (got == 0) break;
        /* The last read is short; write a whole sector and trim below. */
        padded = (got + (BENCH_ALIGN - 1)) & ~(DWORD)(BENCH_ALIGN - 1);
        memset((char*)buf + got, 0, padded - got);
        if (!WriteFile(dst, buf, padded, &put, NULL) || put != padded) ok = 0;
    }

    if (ok) {
        FILE_END_OF_FILE_INFO eof;
        eof.EndOfFile = size;
        ok = SetFileInformationByHandle(dst, FileEndOfFileInfo, &eof, sizeof(eof)) != 0;
    }

    if (buf) VirtualFree(buf, 0, MEM_RELEASE);
    CloseHandle(src);
    CloseHandle(dst);

    if (!ok || !MoveFileExA(tmp, filename, MOVEFILE_REPLACE_EXISTING)) {
        DeleteFileA(tmp);  /* keep the original; the read is then warm */
    }
}

#endif /* _WIN32 */

#endif /* CARQUET_BENCH_PURGE_WIN32_H */
