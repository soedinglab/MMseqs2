#ifndef __FMEMOPEN_H__
#define __FMEMOPEN_H__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Include this file only for OSX / BSD compilations */
#ifdef OS_DARWIN
#define USE_FMEM_WRAPPER 1
#endif

#ifdef OS_FREEBSD
#define USE_FMEM_WRAPPER 1
#endif

#define USE_FMEM_WRAPPER 1
#ifdef USE_FMEM_WRAPPER
struct fmem {
    size_t pos;
    size_t size;
    char *buffer;
};
typedef struct fmem fmem_t;

FILE *fmemopen(void *, size_t, const char *);
#endif

#endif /* __FMEMOPEN_H__ */
