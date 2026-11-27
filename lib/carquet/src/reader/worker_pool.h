/**
 * @file worker_pool.h
 * @brief Persistent thread pool for batch reader parallelism
 *
 * Replaces per-batch OpenMP fork/join with a persistent pool that stays alive
 * across batch_reader_next() calls, eliminating barrier overhead.
 * Also supports row-group lookahead: while the current row group is being
 * consumed, workers pre-decompress pages for the next row group.
 */

#ifndef CARQUET_WORKER_POOL_H
#define CARQUET_WORKER_POOL_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <pthread.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Task and Pool Structures
 * ============================================================================ */

typedef void (*carquet_task_fn)(void* arg);

typedef struct carquet_task {
    carquet_task_fn fn;
    void* arg;
} carquet_task_t;

#define CARQUET_POOL_QUEUE_CAPACITY 512

typedef struct carquet_worker_pool {
#ifdef _WIN32
    HANDLE* threads;
#else
    pthread_t* threads;
#endif
    int32_t num_threads;

    /* Circular task queue protected by mutex */
    carquet_task_t queue[CARQUET_POOL_QUEUE_CAPACITY];
    int32_t queue_head;     /* Next slot to dequeue from */
    int32_t queue_tail;     /* Next slot to enqueue into */
    int32_t queue_count;    /* Number of tasks in queue */

    /* Synchronization */
#ifdef _WIN32
    CRITICAL_SECTION mutex;
    CONDITION_VARIABLE work_available;
    CONDITION_VARIABLE work_done;
    CONDITION_VARIABLE queue_not_full;
#else
    pthread_mutex_t mutex;
    pthread_cond_t work_available;
    pthread_cond_t work_done;
    pthread_cond_t queue_not_full;
#endif

    int32_t active_tasks;   /* Tasks currently being executed */
    bool shutdown;
} carquet_worker_pool_t;

/* ============================================================================
 * API
 * ============================================================================ */

/**
 * Create a worker pool with the given number of threads.
 * Returns NULL on failure.
 */
carquet_worker_pool_t* carquet_worker_pool_create(int32_t num_threads);

/**
 * Submit a task to the pool. The task function will be called with the
 * given argument on a worker thread. Non-blocking.
 */
void carquet_worker_pool_submit(carquet_worker_pool_t* pool,
                                 carquet_task_fn fn, void* arg);

/**
 * Submit N tasks with the same function but different arguments.
 * Acquires the lock once for the entire batch, reducing synchronization overhead.
 */
void carquet_worker_pool_submit_batch(carquet_worker_pool_t* pool,
                                       carquet_task_fn fn,
                                       void** args, int32_t count);

/**
 * Block until all submitted tasks have completed.
 */
void carquet_worker_pool_wait(carquet_worker_pool_t* pool);

/**
 * Submit N tasks and wait for all to complete.
 * Convenience wrapper for the common pattern of parallel-for.
 */
void carquet_worker_pool_parallel_for(carquet_worker_pool_t* pool,
                                       carquet_task_fn fn,
                                       void** args, int32_t count);

/**
 * Destroy the pool, joining all threads.
 */
void carquet_worker_pool_destroy(carquet_worker_pool_t* pool);

#ifdef __cplusplus
}
#endif

#endif /* CARQUET_WORKER_POOL_H */
