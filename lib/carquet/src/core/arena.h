/**
 * @file arena.h
 * @brief Arena (bump) memory allocator
 *
 * Arena allocators provide fast allocation by simply bumping a pointer.
 * Memory is freed all at once when the arena is reset or destroyed.
 * This is ideal for parsing where many small allocations are made
 * and then discarded together.
 */

#ifndef CARQUET_CORE_ARENA_H
#define CARQUET_CORE_ARENA_H

#include <carquet/error.h>
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Constants
 * ============================================================================
 */

#define CARQUET_ARENA_DEFAULT_BLOCK_SIZE (64 * 1024)  /* 64 KB */
#define CARQUET_ARENA_ALIGNMENT 16

/* ============================================================================
 * Types
 * ============================================================================
 */

/**
 * A single block in the arena.
 * Note: The union ensures data[] is properly aligned for all platforms,
 * including 32-bit systems where the struct fields alone would leave
 * data[] at an offset that's not 8-byte aligned.
 */
typedef struct carquet_arena_block {
    struct carquet_arena_block* next;
    size_t size;
    size_t used;
    union {
        uint8_t data[1];  /* Flexible array member (C89 compat) */
        /* Force alignment to match CARQUET_ARENA_ALIGNMENT (16) */
        double _align_double;
        void* _align_ptr;
        long long _align_ll;
    } u;
} carquet_arena_block_t;

/* Access data via u.data */
#define CARQUET_ARENA_BLOCK_DATA(block) ((block)->u.data)

/**
 * Arena allocator.
 */
typedef struct carquet_arena {
    carquet_arena_block_t* head;    /* First block */
    carquet_arena_block_t* current; /* Current block for allocation */
    size_t default_block_size;
    size_t total_allocated;         /* Total bytes allocated */
    size_t total_capacity;          /* Total capacity across all blocks */
} carquet_arena_t;

/* ============================================================================
 * Arena Operations
 * ============================================================================
 */

/**
 * Initialize an arena with default block size.
 */
carquet_status_t carquet_arena_init(carquet_arena_t* arena);

/**
 * Initialize an arena with custom block size.
 */
carquet_status_t carquet_arena_init_size(carquet_arena_t* arena, size_t block_size);

/**
 * Destroy an arena and free all memory.
 */
void carquet_arena_destroy(carquet_arena_t* arena);

/**
 * Reset an arena, freeing all allocations but keeping blocks.
 * This is more efficient than destroy + init for reuse.
 */
void carquet_arena_reset(carquet_arena_t* arena);

/**
 * Allocate memory from the arena.
 *
 * @param arena The arena
 * @param size Number of bytes to allocate
 * @return Pointer to allocated memory, or NULL on failure
 */
void* carquet_arena_alloc(carquet_arena_t* arena, size_t size);

/**
 * Allocate zeroed memory from the arena.
 */
void* carquet_arena_calloc(carquet_arena_t* arena, size_t count, size_t size);

/**
 * Allocate aligned memory from the arena.
 */
void* carquet_arena_alloc_aligned(carquet_arena_t* arena, size_t size, size_t alignment);

/**
 * Duplicate a string into the arena.
 */
char* carquet_arena_strdup(carquet_arena_t* arena, const char* str);

/**
 * Duplicate a string with maximum length into the arena.
 */
char* carquet_arena_strndup(carquet_arena_t* arena, const char* str, size_t max_len);

/**
 * Duplicate a memory region into the arena.
 */
void* carquet_arena_memdup(carquet_arena_t* arena, const void* src, size_t size);

/**
 * Get total bytes allocated from the arena.
 */
static inline size_t carquet_arena_allocated(const carquet_arena_t* arena) {
    return arena->total_allocated;
}

/**
 * Get total capacity of the arena.
 */
static inline size_t carquet_arena_capacity(const carquet_arena_t* arena) {
    return arena->total_capacity;
}

/* ============================================================================
 * Temporary Allocation (Save/Restore)
 * ============================================================================
 */

/**
 * Arena save point for temporary allocations.
 */
typedef struct carquet_arena_mark {
    carquet_arena_block_t* block;
    size_t used;
    size_t total_allocated;
} carquet_arena_mark_t;

/**
 * Save the current arena position.
 */
carquet_arena_mark_t carquet_arena_save(const carquet_arena_t* arena);

/**
 * Restore arena to a saved position, freeing newer allocations.
 */
void carquet_arena_restore(carquet_arena_t* arena, carquet_arena_mark_t mark);

#ifdef __cplusplus
}
#endif

#endif /* CARQUET_CORE_ARENA_H */
