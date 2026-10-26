#ifndef CARQUET_CORE_COMPAT_H
#define CARQUET_CORE_COMPAT_H

#include "allocator.h"
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if !defined(_WIN32)
#include <sys/types.h>
#endif

static inline char* carquet_heap_strdup(const char* str) {
    if (!str) {
        return NULL;
    }

    size_t len = strlen(str) + 1;
    char* copy = (char*)carquet_mem_malloc(len);
    if (!copy) {
        return NULL;
    }

    memcpy(copy, str, len);
    return copy;
}

/* 64-bit file positioning wrappers.
 * `long` is 32-bit on 64-bit Windows, so fseek/ftell silently fail (or wrap)
 * for files larger than 2 GiB. Use platform-specific 64-bit variants. */
static inline int carquet_fseek64(FILE* file, int64_t offset, int whence) {
#if defined(_WIN32)
    return _fseeki64(file, (__int64)offset, whence);
#elif (defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE >= 200112L) || \
      defined(__linux__) || defined(__APPLE__) || defined(__FreeBSD__)
    return fseeko(file, (off_t)offset, whence);
#else
    if (offset > (int64_t)LONG_MAX || offset < (int64_t)LONG_MIN) {
        return -1;
    }
    return fseek(file, (long)offset, whence);
#endif
}

static inline int64_t carquet_ftell64(FILE* file) {
#if defined(_WIN32)
    return (int64_t)_ftelli64(file);
#elif (defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE >= 200112L) || \
      defined(__linux__) || defined(__APPLE__) || defined(__FreeBSD__)
    return (int64_t)ftello(file);
#else
    return (int64_t)ftell(file);
#endif
}

#endif /* CARQUET_CORE_COMPAT_H */
