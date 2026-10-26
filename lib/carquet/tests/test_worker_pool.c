/**
 * @file test_worker_pool.c
 * @brief Concurrency tests for the persistent batch-reader thread pool.
 *
 * carquet_worker_pool_* (src/reader/worker_pool.h) is the parallelism engine
 * behind batch reading: it stays alive across batch_reader_next() calls and
 * fans page-decompression work out to N threads. Bugs here (lost tasks, races
 * on the shared counter, deadlock when the queue fills past capacity, unclean
 * shutdown) would show up as flaky, hard-to-reproduce reader failures. These
 * tests drive the pool directly with atomic accounting so a regression fails
 * deterministically instead of intermittently.
 */

#include <stdatomic.h>
#include <stdint.h>
#include <string.h>

#include "reader/worker_pool.h"
#include "test_helpers.h"

/* Each task bumps a shared atomic; the sum/count is checked after wait(). */
typedef struct {
    _Atomic long* counter;
    long add;
} add_task_t;

static void add_task_fn(void* arg) {
    add_task_t* t = (add_task_t*)arg;
    atomic_fetch_add_explicit(t->counter, t->add, memory_order_relaxed);
}

/* A per-slot writer task: proves every submitted arg was executed exactly once
 * (no dropped or duplicated tasks) without relying on ordering. */
typedef struct { _Atomic int* slot; } touch_task_t;
static void touch_task_fn(void* arg) {
    touch_task_t* t = (touch_task_t*)arg;
    atomic_fetch_add_explicit(t->slot, 1, memory_order_relaxed);
}

static int test_create_destroy(void) {
    carquet_worker_pool_t* p = carquet_worker_pool_create(4);
    if (!p) TEST_FAIL("create_destroy", "create returned NULL");
    /* Destroy with no work submitted must not hang or crash. */
    carquet_worker_pool_destroy(p);
    TEST_PASS("create_destroy");
    return 0;
}

static int test_single_thread_pool(void) {
    /* A one-thread pool still executes everything, just serially. */
    carquet_worker_pool_t* p = carquet_worker_pool_create(1);
    if (!p) TEST_FAIL("single_thread_pool", "create failed");
    _Atomic long counter = 0;
    add_task_t tasks[50];
    void* args[50];
    for (int i = 0; i < 50; i++) { tasks[i].counter = &counter; tasks[i].add = i; args[i] = &tasks[i]; }
    carquet_worker_pool_parallel_for(p, add_task_fn, args, 50);
    long expected = 50 * 49 / 2;
    carquet_worker_pool_destroy(p);
    if (atomic_load(&counter) != expected) TEST_FAIL("single_thread_pool", "sum wrong");
    TEST_PASS("single_thread_pool");
    return 0;
}

static int test_submit_wait(void) {
    carquet_worker_pool_t* p = carquet_worker_pool_create(4);
    if (!p) TEST_FAIL("submit_wait", "create failed");
    _Atomic long counter = 0;
    enum { N = 200 };
    add_task_t tasks[N];
    for (int i = 0; i < N; i++) {
        tasks[i].counter = &counter;
        tasks[i].add = 1;
        carquet_worker_pool_submit(p, add_task_fn, &tasks[i]);
    }
    carquet_worker_pool_wait(p);
    /* After wait() every submitted task must have run. */
    if (atomic_load(&counter) != N) TEST_FAIL("submit_wait", "not all tasks ran before wait returned");
    carquet_worker_pool_destroy(p);
    TEST_PASS("submit_wait");
    return 0;
}

static int test_parallel_for_exactly_once(void) {
    /* Each of N distinct slots must be touched exactly once — catches both
     * dropped tasks (slot==0) and duplicated dispatch (slot>1). */
    carquet_worker_pool_t* p = carquet_worker_pool_create(8);
    if (!p) TEST_FAIL("parallel_for_exactly_once", "create failed");
    enum { N = 1000 };
    static _Atomic int slots[N];
    static touch_task_t tasks[N];
    static void* args[N];
    for (int i = 0; i < N; i++) {
        atomic_store(&slots[i], 0);
        tasks[i].slot = &slots[i];
        args[i] = &tasks[i];
    }
    carquet_worker_pool_parallel_for(p, touch_task_fn, args, N);
    carquet_worker_pool_destroy(p);
    for (int i = 0; i < N; i++) {
        int v = atomic_load(&slots[i]);
        if (v != 1) { printf("  slot %d touched %d times\n", i, v);
                      TEST_FAIL("parallel_for_exactly_once", "task not run exactly once"); }
    }
    TEST_PASS("parallel_for_exactly_once");
    return 0;
}

static int test_submit_batch(void) {
    carquet_worker_pool_t* p = carquet_worker_pool_create(4);
    if (!p) TEST_FAIL("submit_batch", "create failed");
    _Atomic long counter = 0;
    enum { N = 300 };
    static add_task_t tasks[N];
    static void* args[N];
    for (int i = 0; i < N; i++) { tasks[i].counter = &counter; tasks[i].add = 2; args[i] = &tasks[i]; }
    carquet_worker_pool_submit_batch(p, add_task_fn, args, N);
    carquet_worker_pool_wait(p);
    if (atomic_load(&counter) != 2 * N) TEST_FAIL("submit_batch", "batch sum wrong");
    carquet_worker_pool_destroy(p);
    TEST_PASS("submit_batch");
    return 0;
}

static int test_exceeds_queue_capacity(void) {
    /* Submit far more tasks than CARQUET_POOL_QUEUE_CAPACITY (512) in one
     * batch. The queue must apply backpressure (queue_not_full) rather than
     * overflow or deadlock, and every task must still complete. */
    carquet_worker_pool_t* p = carquet_worker_pool_create(4);
    if (!p) TEST_FAIL("exceeds_queue_capacity", "create failed");
    _Atomic long counter = 0;
    enum { N = 4096 };  /* 8x queue capacity */
    static add_task_t tasks[N];
    static void* args[N];
    for (int i = 0; i < N; i++) { tasks[i].counter = &counter; tasks[i].add = 1; args[i] = &tasks[i]; }
    carquet_worker_pool_submit_batch(p, add_task_fn, args, N);
    carquet_worker_pool_wait(p);
    if (atomic_load(&counter) != N) TEST_FAIL("exceeds_queue_capacity", "lost tasks past queue capacity");
    carquet_worker_pool_destroy(p);
    TEST_PASS("exceeds_queue_capacity");
    return 0;
}

static int test_reuse_across_waves(void) {
    /* The pool is persistent: it must handle many submit/wait waves without
     * being recreated, mirroring one wave per row group. */
    carquet_worker_pool_t* p = carquet_worker_pool_create(6);
    if (!p) TEST_FAIL("reuse_across_waves", "create failed");
    _Atomic long counter = 0;
    enum { WAVES = 20, PER = 100 };
    static add_task_t tasks[PER];
    static void* args[PER];
    for (int i = 0; i < PER; i++) { tasks[i].counter = &counter; tasks[i].add = 1; args[i] = &tasks[i]; }
    for (int w = 0; w < WAVES; w++) {
        carquet_worker_pool_parallel_for(p, add_task_fn, args, PER);
    }
    carquet_worker_pool_destroy(p);
    if (atomic_load(&counter) != (long)WAVES * PER) TEST_FAIL("reuse_across_waves", "wave accounting wrong");
    TEST_PASS("reuse_across_waves");
    return 0;
}

static int test_zero_count(void) {
    /* parallel_for / submit_batch with count 0 must be a clean no-op. */
    carquet_worker_pool_t* p = carquet_worker_pool_create(4);
    if (!p) TEST_FAIL("zero_count", "create failed");
    carquet_worker_pool_parallel_for(p, add_task_fn, NULL, 0);
    carquet_worker_pool_submit_batch(p, add_task_fn, NULL, 0);
    carquet_worker_pool_wait(p);
    carquet_worker_pool_destroy(p);
    TEST_PASS("zero_count");
    return 0;
}

int main(void) {
    int failures = 0;
    failures += test_create_destroy();
    failures += test_single_thread_pool();
    failures += test_submit_wait();
    failures += test_parallel_for_exactly_once();
    failures += test_submit_batch();
    failures += test_exceeds_queue_capacity();
    failures += test_reuse_across_waves();
    failures += test_zero_count();
    if (failures) { printf("\n%d test(s) FAILED\n", failures); return 1; }
    printf("\nAll worker pool tests passed\n");
    return 0;
}
