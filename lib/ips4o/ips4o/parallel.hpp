/******************************************************************************
 * ips4o/parallel.hpp
 *
 * In-place Parallel Super Scalar Samplesort (IPS⁴o)
 *
 ******************************************************************************
 * BSD 2-Clause License
 *
 * Copyright © 2017, Michael Axtmann <michael.axtmann@kit.edu>
 * Copyright © 2017, Daniel Ferizovic <daniel.ferizovic@student.kit.edu>
 * Copyright © 2017, Sascha Witt <sascha.witt@kit.edu>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *****************************************************************************/

#pragma once
#if defined(_REENTRANT) || defined(_OPENMP)

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <functional>
#include <thread>
#include <utility>
#include <vector>

#include "ips4o_fwd.hpp"
#include "memory.hpp"
#include "sequential.hpp"
#include "partitioning.hpp"

namespace ips4o {
namespace detail {

/**
 * Processes sequential subtasks in the parallel algorithm.
 */
template <class Cfg>
void Sorter<Cfg>::processSmallTasks(const iterator begin, SharedData& shared) {
    std::size_t my_index = shared.small_task_index.fetch_add(1, std::memory_order_relaxed);
    while (my_index < shared.small_tasks.size()) {
        const auto my_task = shared.small_tasks[my_index];
        sequential(begin + my_task.begin, begin + my_task.end);
        my_index = shared.small_task_index.fetch_add(1, std::memory_order_relaxed);
    }
}

/**
 * Main loop for secondary threads in the parallel algorithm.
 */
template <class Cfg>
void Sorter<Cfg>::parallelSecondary(SharedData& shared, int id, int num_threads) {
    const auto begin = shared.begin_;
    shared.local[id] = &local_;
    do {
        const auto task = shared.big_tasks.back();
        partition<true>(begin + task.begin, begin + task.end, shared.bucket_start, &shared, id, num_threads);
        shared.sync.barrier();
    } while (!shared.big_tasks.empty());

    processSmallTasks(begin, shared);
}

/**
 * Main loop for the primary thread in the parallel algorithm.
 */
template <class Cfg>
template <class TaskSorter>
void Sorter<Cfg>::parallelPrimary(const iterator begin, const iterator end,
                                  SharedData& shared, const int num_threads,
                                  TaskSorter&& task_sorter) {
    const auto max_sequential_size = [&] {
        const auto div = (end - begin) / num_threads;
        const auto blocks = num_threads * Cfg::kMinParallelBlocksPerThread * Cfg::kBlockSize;
        return std::max(div, blocks);
    }();

    shared.small_tasks.clear();
    shared.small_task_index.store(0, std::memory_order_relaxed);
    // Queues a subtask either as a big task, a small task, or not at all, depending on the size
    const auto queueTask = [&shared, max_sequential_size](int i, std::ptrdiff_t offset, int level) {
        const auto start = offset + shared.bucket_start[i];
        const auto stop = offset + shared.bucket_start[i + 1];
        if (stop - start > max_sequential_size) {
            shared.big_tasks.push_back({start, stop, level});
        } else if (stop - start > 2 * Cfg::kBaseCaseSize) {
            shared.small_tasks.push_back({start, stop, level});
        }
    };

    do {
        // Do parallel partitioning
        const auto task = shared.big_tasks.back();
        const auto res = partition<true>(begin + task.begin, begin + task.end,
                                         shared.bucket_start, &shared, 0, num_threads);
        const int num_buckets = std::get<0>(res);
        const bool equal_buckets = std::get<1>(res);
        shared.big_tasks.pop_back();

        // Queue subtasks if we didn't reach the last level yet
        const bool is_last_level = end_ - begin_ <= Cfg::kSingleLevelThreshold;
        if (!is_last_level) {
            // Skip equality buckets
            for (int i = 0; i < num_buckets; i += 1 + equal_buckets)
                queueTask(i, task.begin, task.level + 1);
            if (equal_buckets)
                queueTask(num_buckets - 1, task.begin, task.level + 1);
        }
        if (shared.big_tasks.empty()) {
            // Sort small tasks by size, larger ones first
            task_sorter(shared.small_tasks.begin(), shared.small_tasks.end());
        }

        shared.reset();
        shared.sync.barrier();
    } while (!shared.big_tasks.empty());

    // Process remaining small tasks
    processSmallTasks(begin, shared);
}

}  // namespace detail

/**
 * Reusable parallel sorter.
 */
template <class Cfg>
class ParallelSorter {
    using Sorter = detail::Sorter<Cfg>;
    using iterator = typename Cfg::iterator;

 public:
    /**
     * Construct the sorter. Thread pool may be passed by reference.
     */
    ParallelSorter(typename Cfg::less comp, typename Cfg::ThreadPool thread_pool)
            : thread_pool_(std::forward<typename Cfg::ThreadPool>(thread_pool))
            , shared_ptr_(Cfg::kDataAlignment, std::move(comp), thread_pool_.sync(), thread_pool_.numThreads())
            , buffer_storage_(thread_pool_.numThreads())
            , local_ptrs_(new detail::AlignedPtr<typename Sorter::LocalData>[thread_pool_.numThreads()])
            , task_sorter_(false, {}, buffer_storage_.forThread(0))
    {
        // Allocate local data
        thread_pool_([this](int my_id, int) {
            auto& shared = this->shared_ptr_.get();
            this->local_ptrs_[my_id] = detail::AlignedPtr<typename Sorter::LocalData>(
                            Cfg::kDataAlignment, shared.classifier.getComparator(), buffer_storage_.forThread(my_id));
            shared.local[my_id] = &this->local_ptrs_[my_id].get();
        });
    }

    /**
     * Sort in parallel.
     */
    void operator()(iterator begin, iterator end) {
        // Sort small input sequentially
        const int num_threads = Cfg::numThreadsFor(begin, end, thread_pool_.numThreads());
        if (num_threads < 2) {
            Sorter(local_ptrs_[0].get()).sequential(std::move(begin), std::move(end));
            return;
        }

        // Set up base data before switching to parallel mode
        auto& shared = shared_ptr_.get();
        shared.begin_ = begin;
        shared.big_tasks.push_back({0, end - begin, 1});

        // Execute in parallel
        thread_pool_([this, begin, end](int my_id, int num_threads) {
            auto& shared = this->shared_ptr_.get();
            Sorter sorter(*shared.local[my_id]);
            if (my_id == 0)
                sorter.parallelPrimary(begin, end, shared, num_threads, task_sorter_);
            else
                sorter.parallelSecondary(shared, my_id, num_threads);
        }, num_threads);
    }

 private:
    typename Cfg::ThreadPool thread_pool_;
    detail::AlignedPtr<typename Sorter::SharedData> shared_ptr_;
    typename Sorter::BufferStorage buffer_storage_;
    std::unique_ptr<detail::AlignedPtr<typename Sorter::LocalData>[]> local_ptrs_;
    SequentialSorter<ExtendedConfig<
                std::vector<detail::ParallelTask>::iterator,
                std::greater<>, typename Cfg::BaseConfig
                       >> task_sorter_;
};

}  // namespace ips4o
#endif  // _REENTRANT || _OPENMP
