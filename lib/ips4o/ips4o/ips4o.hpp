/******************************************************************************
 * ips4o/ips4o.hpp
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

#include <functional>
#include <iterator>
#include <type_traits>

#include "ips4o_fwd.hpp"
#include "base_case.hpp"
#include "config.hpp"
#include "memory.hpp"
#include "parallel.hpp"
#include "sequential.hpp"

namespace ips4o {

/**
 * Helper function for creating a reusable sequential sorter.
 */
template <class It, class Cfg = Config<>, class Comp = std::less<>>
SequentialSorter<ExtendedConfig<It, Comp, Cfg>> make_sorter(Comp comp = Comp()) {
  return SequentialSorter<ExtendedConfig<It, Comp, Cfg>>{true, std::move(comp)};
}

/**
 * Configurable interface.
 */
template <class Cfg, class It, class Comp = std::less<>>
void sort(It begin, It end, Comp comp = Comp()) {
  if (detail::sortedCaseSort(begin, end, comp)) return;
  
  if ((end - begin) <= Cfg::kBaseCaseMultiplier * Cfg::kBaseCaseSize) {
    detail::baseCaseSort(std::move(begin), std::move(end), std::move(comp));
  } else {
    ips4o::SequentialSorter<ips4o::ExtendedConfig<It, Comp, Cfg>> sorter{false, std::move(comp)};
    sorter(std::move(begin), std::move(end));
  }
}

/**
 * Standard interface.
 */
template <class It, class Comp>
void sort(It begin, It end, Comp comp) {
    ips4o::sort<Config<>>(std::move(begin), std::move(end), std::move(comp));
}

template <class It>
void sort(It begin, It end) {
    ips4o::sort<Config<>>(std::move(begin), std::move(end), std::less<>());
}

#if defined(_REENTRANT) || defined(_OPENMP)
namespace parallel {

/**
 * Helper functions for creating a reusable parallel sorter.
 */
template <class It, class Cfg = Config<>, class ThreadPool, class Comp = std::less<>>
std::enable_if_t<std::is_class<std::remove_reference_t<ThreadPool>>::value,
                 ParallelSorter<ExtendedConfig<It, Comp, Cfg, ThreadPool>>>
make_sorter(ThreadPool&& thread_pool, Comp comp = Comp()) {
    return ParallelSorter<ExtendedConfig<It, Comp, Cfg, ThreadPool>>(
            std::move(comp), std::forward<ThreadPool>(thread_pool));
}

template <class It, class Cfg = Config<>, class Comp = std::less<>>
ParallelSorter<ExtendedConfig<It, Comp, Cfg>> make_sorter(
        int num_threads = DefaultThreadPool::maxNumThreads(), Comp comp = Comp()) {
    return ParallelSorter<ExtendedConfig<It, Comp, Cfg>>(
            std::move(comp), DefaultThreadPool(num_threads));
}

/**
 * Configurable interface.
 */
template <class Cfg = Config<>, class It, class Comp, class ThreadPool>
std::enable_if_t<std::is_class<std::remove_reference_t<ThreadPool>>::value>
sort(It begin, It end, Comp comp, ThreadPool&& thread_pool) {
    if (Cfg::numThreadsFor(begin, end, thread_pool.numThreads()) < 2)
        ips4o::sort<Cfg>(std::move(begin), std::move(end), std::move(comp));
    else
        ips4o::parallel::make_sorter<It, Cfg>(
                std::forward<ThreadPool>(thread_pool), std::move(comp))(std::move(begin),
                                                                        std::move(end));
}

template <class Cfg = Config<>, class It, class Comp>
void sort(It begin, It end, Comp comp, int num_threads) {
    num_threads = Cfg::numThreadsFor(begin, end, num_threads);
    if (num_threads < 2)
        ips4o::sort<Cfg>(std::move(begin), std::move(end), std::move(comp));
    else
        ips4o::parallel::make_sorter<It, Cfg>(num_threads, std::move(comp))
                (std::move(begin), std::move(end));
}

/**
 * Standard interface.
 */
template <class It, class Comp>
void sort(It begin, It end, Comp comp) {
    ips4o::parallel::sort<Config<>>(std::move(begin), std::move(end), std::move(comp),
                                    DefaultThreadPool::maxNumThreads());
}

template <class It>
void sort(It begin, It end) {
    ips4o::parallel::sort(std::move(begin), std::move(end), std::less<>());
}

}  // namespace parallel
#endif  // _REENTRANT || _OPENMP
}  // namespace ips4o
