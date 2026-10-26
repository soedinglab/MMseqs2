/**
 * @file allocator.h
 * @brief Internal allocation wrappers that route through the global allocator.
 *
 * Every heap allocation in the library goes through these wrappers so that a
 * custom allocator installed via carquet_set_allocator() is actually used.
 *
 * The public contract requires carquet_set_allocator() to be called before
 * any allocation and before concurrent use, so the active allocator is fixed
 * for the lifetime of all allocations: a block allocated through these
 * wrappers is always freed through them with the same allocator. Never mix
 * these with libc malloc/free for the same pointer.
 */
#ifndef CARQUET_CORE_ALLOCATOR_H
#define CARQUET_CORE_ALLOCATOR_H

#include <stddef.h>

/** Allocate @p size bytes (size 0 yields a unique freeable pointer or NULL). */
void* carquet_mem_malloc(size_t size);

/** Allocate @p nmemb * @p size zeroed bytes, with overflow check. */
void* carquet_mem_calloc(size_t nmemb, size_t size);

/** Resize @p ptr to @p size bytes (ptr may be NULL => malloc). */
void* carquet_mem_realloc(void* ptr, size_t size);

/** Free @p ptr (NULL is a no-op). */
void carquet_mem_free(void* ptr);

#endif /* CARQUET_CORE_ALLOCATOR_H */
