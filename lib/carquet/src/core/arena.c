/**
 * @file arena.c
 * @brief Arena (bump) memory allocator implementation
 */

#include "allocator.h"
#include "arena.h"
#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * Internal Helpers
 * ============================================================================
 */

static inline size_t align_up(size_t value, size_t alignment) {
    if (alignment == 0 || value > SIZE_MAX - (alignment - 1)) {
        return SIZE_MAX;
    }
    return (value + alignment - 1) & ~(alignment - 1);
}

static int add_overflows_size(size_t a, size_t b, size_t* out) {
    if (a > SIZE_MAX - b) {
        return 1;
    }
    *out = a + b;
    return 0;
}

static carquet_arena_block_t* arena_new_block(size_t min_size) {
    size_t block_size = min_size < CARQUET_ARENA_DEFAULT_BLOCK_SIZE
                            ? CARQUET_ARENA_DEFAULT_BLOCK_SIZE
                            : align_up(min_size, CARQUET_ARENA_DEFAULT_BLOCK_SIZE);
    if (block_size == SIZE_MAX) {
        return NULL;
    }

    /* Allocate the header plus the data block.
     * Note: offsetof accounts for the union's alignment padding */
    size_t header_size = offsetof(carquet_arena_block_t, u);
    size_t alloc_size;
    if (add_overflows_size(header_size, block_size, &alloc_size)) {
        return NULL;
    }
    carquet_arena_block_t* block = (carquet_arena_block_t*)carquet_mem_malloc(alloc_size);

    if (!block) {
        return NULL;
    }

    block->next = NULL;
    block->size = block_size;
    block->used = 0;

    return block;
}

/* ============================================================================
 * Arena Operations
 * ============================================================================
 */

carquet_status_t carquet_arena_init(carquet_arena_t* arena) {
    return carquet_arena_init_size(arena, CARQUET_ARENA_DEFAULT_BLOCK_SIZE);
}

carquet_status_t carquet_arena_init_size(carquet_arena_t* arena, size_t block_size) {
    assert(arena != NULL);

    /* Zero-initialize the arena structure first */
    arena->head = NULL;
    arena->current = NULL;
    arena->default_block_size = block_size;
    arena->total_allocated = 0;
    arena->total_capacity = 0;

    arena->head = arena_new_block(block_size);
    if (!arena->head) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    arena->current = arena->head;
    arena->total_capacity = arena->head->size;

    return CARQUET_OK;
}

void carquet_arena_destroy(carquet_arena_t* arena) {
    assert(arena != NULL);

    carquet_arena_block_t* block = arena->head;
    while (block) {
        carquet_arena_block_t* next = block->next;
        carquet_mem_free(block);
        block = next;
    }

    arena->head = NULL;
    arena->current = NULL;
    arena->total_allocated = 0;
    arena->total_capacity = 0;
}

void carquet_arena_reset(carquet_arena_t* arena) {
    assert(arena != NULL);

    /* Reset all blocks to empty */
    carquet_arena_block_t* block = arena->head;
    while (block) {
        block->used = 0;
        block = block->next;
    }

    arena->current = arena->head;
    arena->total_allocated = 0;
}

void* carquet_arena_alloc(carquet_arena_t* arena, size_t size) {
    return carquet_arena_alloc_aligned(arena, size, CARQUET_ARENA_ALIGNMENT);
}

/**
 * Helper to calculate aligned offset within a block.
 * This calculates alignment based on absolute addresses, not just offsets,
 * which is necessary on 32-bit systems where malloc may not provide
 * sufficient alignment.
 */
static inline size_t arena_aligned_offset(carquet_arena_block_t* block,
                                           size_t current_used,
                                           size_t alignment) {
    uintptr_t base = (uintptr_t)CARQUET_ARENA_BLOCK_DATA(block);
    uintptr_t current_addr = base + current_used;
    uintptr_t aligned_addr = (current_addr + alignment - 1) & ~(alignment - 1);
    return (size_t)(aligned_addr - base);
}

void* carquet_arena_alloc_aligned(carquet_arena_t* arena, size_t size, size_t alignment) {
    assert(arena != NULL);
    if (size == 0) {
        return NULL;
    }

    /* Ensure alignment is power of 2 and at least 1 */
    if (alignment == 0) {
        alignment = 1;
    }

    carquet_arena_block_t* block = arena->current;
    assert(block != NULL);  /* Arena must be properly initialized */

    /* Calculate aligned offset based on absolute address */
    size_t aligned_offset = arena_aligned_offset(block, block->used, alignment);
    size_t new_used;
    if (add_overflows_size(aligned_offset, size, &new_used)) {
        return NULL;
    }

    /* Check if current block has space */
    if (new_used <= block->size) {
        void* ptr = CARQUET_ARENA_BLOCK_DATA(block) + aligned_offset;
        block->used = new_used;
        arena->total_allocated += size;
        return ptr;
    }

    /* Try next blocks */
    while (block->next) {
        block = block->next;
        aligned_offset = arena_aligned_offset(block, block->used, alignment);
        if (add_overflows_size(aligned_offset, size, &new_used)) {
            return NULL;
        }

        if (new_used <= block->size) {
            arena->current = block;
            void* ptr = CARQUET_ARENA_BLOCK_DATA(block) + aligned_offset;
            block->used = new_used;
            arena->total_allocated += size;
            return ptr;
        }
    }

    /* Need new block */
    size_t needed;  /* Worst case alignment overhead */
    if (add_overflows_size(size, alignment, &needed)) {
        return NULL;
    }
    size_t block_size = needed > arena->default_block_size
                            ? needed
                            : arena->default_block_size;

    carquet_arena_block_t* new_block = arena_new_block(block_size);
    if (!new_block) {
        return NULL;
    }

    /* Link new block */
    block->next = new_block;
    arena->current = new_block;
    arena->total_capacity += new_block->size;

    /* Allocate from new block */
    aligned_offset = arena_aligned_offset(new_block, new_block->used, alignment);
    if (add_overflows_size(aligned_offset, size, &new_used)) {
        return NULL;
    }
    new_block->used = new_used;
    arena->total_allocated += size;

    return CARQUET_ARENA_BLOCK_DATA(new_block) + aligned_offset;
}

void* carquet_arena_calloc(carquet_arena_t* arena, size_t count, size_t size) {
    size_t total = count * size;

    /* Check for overflow */
    if (count != 0 && total / count != size) {
        return NULL;
    }

    void* ptr = carquet_arena_alloc(arena, total);
    if (ptr) {
        memset(ptr, 0, total);
    }
    return ptr;
}

char* carquet_arena_strdup(carquet_arena_t* arena, const char* str) {
    if (!str) {
        return NULL;
    }
    return carquet_arena_strndup(arena, str, strlen(str));
}

char* carquet_arena_strndup(carquet_arena_t* arena, const char* str, size_t max_len) {
    if (!str) {
        return NULL;
    }

    size_t len = 0;
    while (len < max_len && str[len]) {
        len++;
    }

    char* copy = (char*)carquet_arena_alloc_aligned(arena, len + 1, 1);
    if (copy) {
        memcpy(copy, str, len);
        copy[len] = '\0';
    }
    return copy;
}

void* carquet_arena_memdup(carquet_arena_t* arena, const void* src, size_t size) {
    if (!src || size == 0) {
        return NULL;
    }

    void* copy = carquet_arena_alloc(arena, size);
    if (copy) {
        memcpy(copy, src, size);
    }
    return copy;
}

/* ============================================================================
 * Save/Restore
 * ============================================================================
 */

carquet_arena_mark_t carquet_arena_save(const carquet_arena_t* arena) {
    assert(arena != NULL);
    assert(arena->current != NULL);

    carquet_arena_mark_t mark = {
        .block = arena->current,
        .used = arena->current->used,
        .total_allocated = arena->total_allocated,
    };
    return mark;
}

void carquet_arena_restore(carquet_arena_t* arena, carquet_arena_mark_t mark) {
    assert(arena != NULL);
    assert(mark.block != NULL);

    /* Reset blocks after the marked block */
    carquet_arena_block_t* block = mark.block->next;
    while (block) {
        block->used = 0;
        block = block->next;
    }

    /* Restore marked block state */
    mark.block->used = mark.used;
    arena->current = mark.block;
    arena->total_allocated = mark.total_allocated;
}
