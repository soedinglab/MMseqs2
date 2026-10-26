/**
 * @file allocator.c
 * @brief Global memory allocator accessor.
 *
 * Stores the process-wide allocator configuration. The default is the C
 * standard library allocator. carquet_set_allocator() must be called before
 * any concurrent use (it is documented as not thread-safe).
 *
 * This provides the public allocator accessors declared in carquet.h and the
 * internal carquet_mem_* wrappers (see allocator.h) that the rest of the
 * library uses for every heap allocation, so a custom allocator is honored.
 */

#include <carquet/carquet.h>
#include "allocator.h"
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

static void* default_malloc(size_t size, void* ctx) {
    (void)ctx;
    return malloc(size);
}

static void* default_realloc(void* ptr, size_t size, void* ctx) {
    (void)ctx;
    return realloc(ptr, size);
}

static void default_free(void* ptr, void* ctx) {
    (void)ctx;
    free(ptr);
}

static const carquet_allocator_t g_default_allocator = {
    default_malloc,
    default_realloc,
    default_free,
    NULL
};

static carquet_allocator_t g_allocator = {
    default_malloc,
    default_realloc,
    default_free,
    NULL
};

void carquet_set_allocator(const carquet_allocator_t* allocator) {
    if (allocator == NULL ||
        allocator->malloc == NULL ||
        allocator->realloc == NULL ||
        allocator->free == NULL) {
        /* NULL or incomplete allocator resets to the libc default. */
        g_allocator = g_default_allocator;
        return;
    }
    g_allocator = *allocator;
}

const carquet_allocator_t* carquet_get_allocator(void) {
    return &g_allocator;
}

/* ============================================================================
 * Internal wrappers (see allocator.h)
 * ============================================================================
 */

void* carquet_mem_malloc(size_t size) {
    return g_allocator.malloc(size, g_allocator.ctx);
}

void* carquet_mem_calloc(size_t nmemb, size_t size) {
    size_t total;
    if (nmemb != 0 && size > SIZE_MAX / nmemb) {
        return NULL;  /* multiplication would overflow */
    }
    total = nmemb * size;
    void* p = g_allocator.malloc(total, g_allocator.ctx);
    if (p && total) {
        memset(p, 0, total);
    }
    return p;
}

void* carquet_mem_realloc(void* ptr, size_t size) {
    return g_allocator.realloc(ptr, size, g_allocator.ctx);
}

void carquet_mem_free(void* ptr) {
    if (ptr) {
        g_allocator.free(ptr, g_allocator.ctx);
    }
}
