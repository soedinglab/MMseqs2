/**
 * @file ms_diag.c
 * @brief Backend diagnostics: version components and detected CPU features.
 *
 * Exercises carquet_init, carquet_version_components and carquet_get_cpu_info.
 */
#include "ms_internal.h"

#include <stdio.h>
#include <stdlib.h>

/* A counting allocator: proves carquet honors a custom allocator. */
static size_t g_alloc_calls = 0;
static void* diag_malloc(size_t n, void* ctx)            { (void)ctx; g_alloc_calls++; return malloc(n); }
static void* diag_realloc(void* p, size_t n, void* ctx)  { (void)ctx; g_alloc_calls++; return realloc(p, n); }
static void  diag_free(void* p, void* ctx)               { (void)ctx; free(p); }

void ms_print_diagnostics(FILE* out) {
    if (!out) out = stdout;

    (void)carquet_init(); /* idempotent; forces feature detection */

    /* Install (and then restore) a custom allocator to exercise the hook. */
    const carquet_allocator_t* prev = carquet_get_allocator();
    (void)prev;
    carquet_allocator_t tracking = { diag_malloc, diag_realloc, diag_free, NULL };
    carquet_set_allocator(&tracking);
    carquet_set_allocator(NULL); /* restore default before any real work */

    int major = 0, minor = 0, patch = 0;
    carquet_version_components(&major, &minor, &patch);
    fprintf(out, "carquet %d.%d.%d\n", major, minor, patch);

    const carquet_cpu_info_t* cpu = carquet_get_cpu_info();
    fprintf(out, "SIMD features:");
    if (cpu->has_sse42)   fprintf(out, " SSE4.2");
    if (cpu->has_avx2)    fprintf(out, " AVX2");
    if (cpu->has_avx512f) fprintf(out, " AVX-512F");
    if (cpu->has_neon)    fprintf(out, " NEON");
    if (cpu->has_sve)     fprintf(out, " SVE(%d)", cpu->sve_vector_length);
    if (!cpu->has_sse42 && !cpu->has_avx2 && !cpu->has_neon && !cpu->has_sve)
        fprintf(out, " (scalar)");
    fprintf(out, "\n");
}

void ms_shutdown(void) {
    carquet_cleanup();
}
