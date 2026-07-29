// Tests for the shared-filesystem coordination layer.
//
// The interesting failure modes only appear across *processes*: fcntl locks are
// owned per-process, so a thread-only test would pass even if the locking were
// completely absent. The concurrency tests therefore fork real children, which is
// also how the layer is used in production (one worker process per node).

#include "ParallelCoordination.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <set>
#include <string>
#include <thread>
#include <vector>

#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>

const char* binary_name = "test_parallelcoordination";

static int failures = 0;

static void check(bool condition, const std::string &what) {
    if (condition == false) {
        fprintf(stderr, "FAIL: %s\n", what.c_str());
        failures++;
    } else {
        fprintf(stdout, "ok:   %s\n", what.c_str());
    }
}

static std::string makeTempDir() {
    char tmpl[] = "/tmp/mmseqs_coord_testXXXXXX";
    char *dir = mkdtemp(tmpl);
    if (dir == NULL) {
        perror("mkdtemp");
        exit(EXIT_FAILURE);
    }
    return std::string(dir);
}

static void removeTempDir(const std::string &dir) {
    std::string cmd = "rm -rf '" + dir + "'";
    if (system(cmd.c_str()) != 0) {
        fprintf(stderr, "warning: could not clean up %s\n", dir.c_str());
    }
}

// Children report the values they observed by appending them to their own file,
// which avoids needing any shared memory between parent and children.
static void writeValues(const std::string &path, const std::vector<int64_t> &values) {
    FILE *file = fopen(path.c_str(), "w");
    if (file == NULL) {
        perror(path.c_str());
        exit(EXIT_FAILURE);
    }
    for (size_t i = 0; i < values.size(); i++) {
        fprintf(file, "%lld\n", static_cast<long long>(values[i]));
    }
    fclose(file);
}

static std::vector<int64_t> readValues(const std::string &path) {
    std::vector<int64_t> values;
    FILE *file = fopen(path.c_str(), "r");
    if (file == NULL) {
        return values;
    }
    long long value;
    while (fscanf(file, "%lld", &value) == 1) {
        values.push_back(static_cast<int64_t>(value));
    }
    fclose(file);
    return values;
}

// Every fetchAdd must return a distinct value even when several processes race,
// because worker identity is assigned this way and duplicate ids would make two
// workers write to the same output shard.
static void testCounterAcrossProcesses(const std::string &dir) {
    const std::string counterPath = dir + "/worker_id";
    const int childCount = 8;
    const int perChild = 25;

    for (int child = 0; child < childCount; child++) {
        pid_t pid = fork();
        if (pid == 0) {
            SharedCounter counter(counterPath);
            std::vector<int64_t> mine;
            for (int i = 0; i < perChild; i++) {
                mine.push_back(counter.fetchAdd(1));
            }
            writeValues(dir + "/counter_child_" + std::to_string(child), mine);
            _exit(EXIT_SUCCESS);
        }
    }
    for (int child = 0; child < childCount; child++) {
        int status = 0;
        wait(&status);
    }

    std::set<int64_t> seen;
    for (int child = 0; child < childCount; child++) {
        std::vector<int64_t> values = readValues(dir + "/counter_child_" + std::to_string(child));
        for (size_t i = 0; i < values.size(); i++) {
            seen.insert(values[i]);
        }
    }

    SharedCounter counter(counterPath);
    check(counter.get() == childCount * perChild, "counter total is exact under process contention");
    check(seen.size() == static_cast<size_t>(childCount * perChild),
          "every fetchAdd returned a distinct value");
    check(*seen.begin() == 0 && *seen.rbegin() == childCount * perChild - 1,
          "fetchAdd values form a dense range");
}

// Threads inside one process share the fcntl lock (it is per-process), so this
// exercises the process-local mutex that FileLock adds on top.
static void testCounterAcrossThreads(const std::string &dir) {
    const std::string counterPath = dir + "/thread_counter";
    SharedCounter counter(counterPath);

    const int threadCount = 8;
    const int perThread = 25;
    std::vector<std::vector<int64_t> > perThreadValues(threadCount);
    std::vector<std::thread> threads;
    for (int t = 0; t < threadCount; t++) {
        threads.push_back(std::thread([&counter, &perThreadValues, t, perThread]() {
            for (int i = 0; i < perThread; i++) {
                perThreadValues[t].push_back(counter.fetchAdd(1));
            }
        }));
    }
    for (size_t t = 0; t < threads.size(); t++) {
        threads[t].join();
    }

    std::set<int64_t> seen;
    for (int t = 0; t < threadCount; t++) {
        for (size_t i = 0; i < perThreadValues[t].size(); i++) {
            seen.insert(perThreadValues[t][i]);
        }
    }
    check(seen.size() == static_cast<size_t>(threadCount * perThread),
          "every fetchAdd returned a distinct value across threads");
}

static void testQueueDrainsExactlyOnce(const std::string &dir) {
    const std::string queuePath = dir + "/single_queue";
    const int64_t itemCount = 64;
    WorkQueue queue(queuePath, itemCount);

    std::set<int64_t> claimed;
    int64_t index;
    while ((index = queue.claim(1)) >= 0) {
        claimed.insert(index);
        queue.complete(index, 1);
    }

    check(claimed.size() == static_cast<size_t>(itemCount), "single worker claims every item once");
    check(queue.allDone(), "queue reports all done after draining");
    check(queue.getDoneCount() == itemCount, "done count matches item count");
}

static void testQueueConcurrentClaims(const std::string &dir) {
    const std::string queuePath = dir + "/concurrent_queue";
    const int64_t itemCount = 200;
    const int childCount = 8;

    // Created up front so children do not race on initialisation semantics that
    // the resume test covers separately.
    {
        WorkQueue queue(queuePath, itemCount);
        check(queue.getDoneCount() == 0, "fresh queue starts empty");
    }

    for (int child = 0; child < childCount; child++) {
        pid_t pid = fork();
        if (pid == 0) {
            WorkQueue queue(queuePath, itemCount);
            std::vector<int64_t> mine;
            int64_t index;
            while ((index = queue.claim(child + 1)) >= 0) {
                mine.push_back(index);
                queue.complete(index, child + 1);
            }
            writeValues(dir + "/queue_child_" + std::to_string(child), mine);
            _exit(EXIT_SUCCESS);
        }
    }
    for (int child = 0; child < childCount; child++) {
        int status = 0;
        wait(&status);
    }

    std::vector<int64_t> all;
    for (int child = 0; child < childCount; child++) {
        std::vector<int64_t> values = readValues(dir + "/queue_child_" + std::to_string(child));
        all.insert(all.end(), values.begin(), values.end());
    }
    std::set<int64_t> unique(all.begin(), all.end());

    check(all.size() == static_cast<size_t>(itemCount), "no item was handed out twice");
    check(unique.size() == static_cast<size_t>(itemCount), "every item was handed out");

    WorkQueue queue(queuePath, itemCount);
    check(queue.allDone(), "queue is drained after concurrent workers finish");
}

// The property that makes a stage survive a 24 h walltime: a worker that dies
// holding a claim must not strand its item forever.
static void testLeaseExpiryRecovers(const std::string &dir) {
    const std::string queuePath = dir + "/lease_queue";
    WorkQueue queue(queuePath, 4);

    int64_t first = queue.claim(1, 1);
    check(first == 0, "first claim takes the lowest item");

    int64_t second = queue.claim(2, 60);
    check(second == 1, "a live claim is not stolen by another worker");

    sleep(2);

    int64_t reclaimed = queue.claim(3, 60);
    check(reclaimed == 0, "an expired claim is handed to another worker");

    // The original holder must not be able to mark it done any more -- but if it
    // does, the item stays done exactly once rather than double-counting.
    queue.complete(0, 3);
    queue.complete(0, 1);
    check(queue.getDoneCount() == 1, "completing twice counts once");
}

static void testReleaseRequeues(const std::string &dir) {
    const std::string queuePath = dir + "/release_queue";
    WorkQueue queue(queuePath, 4);

    int64_t index = queue.claim(1, 600);
    check(index == 0, "claimed the first item");
    queue.release(index, 1);

    int64_t again = queue.claim(2, 600);
    check(again == 0, "released item is immediately re-claimable");
}

// Reopening a queue must preserve progress; this is what lets a stage continue in
// a fresh Slurm job without redoing finished work.
static void testResumeKeepsProgress(const std::string &dir) {
    const std::string queuePath = dir + "/resume_queue";
    const int64_t itemCount = 16;

    {
        WorkQueue queue(queuePath, itemCount);
        for (int i = 0; i < 8; i++) {
            int64_t index = queue.claim(1);
            queue.complete(index, 1);
        }
        check(queue.getDoneCount() == 8, "half the queue is done before restart");
    }

    {
        WorkQueue queue(queuePath, itemCount);
        check(queue.getDoneCount() == 8, "reopened queue remembers finished work");

        std::set<int64_t> claimed;
        int64_t index;
        while ((index = queue.claim(2)) >= 0) {
            claimed.insert(index);
            queue.complete(index, 2);
        }
        check(claimed.size() == 8, "restart only redoes the unfinished items");
        check(queue.allDone(), "queue completes after restart");
    }
}

int main(int, const char**) {
    std::string dir = makeTempDir();

    testCounterAcrossProcesses(dir);
    testCounterAcrossThreads(dir);
    testQueueDrainsExactlyOnce(dir);
    testQueueConcurrentClaims(dir);
    testLeaseExpiryRecovers(dir);
    testReleaseRequeues(dir);
    testResumeKeepsProgress(dir);

    removeTempDir(dir);

    if (failures > 0) {
        fprintf(stderr, "\n%d check(s) failed\n", failures);
        return EXIT_FAILURE;
    }
    fprintf(stdout, "\nall checks passed\n");
    return EXIT_SUCCESS;
}
