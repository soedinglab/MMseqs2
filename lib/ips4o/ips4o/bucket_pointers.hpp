/******************************************************************************
 * ips4o/bucket_pointers.hpp
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

#include <atomic>
#include <climits>
#include <cstdint>
#include <utility>

#include "ips4o_fwd.hpp"

namespace ips4o {
namespace detail {

#if UINTPTR_MAX != UINT32_MAX && !defined(__SIZEOF_INT128__)
#error "Unsupported architecture"
#endif

template <class Cfg>
class Sorter<Cfg>::BucketPointers {
#if UINTPTR_MAX == UINT32_MAX
    using atomic_type = std::uint64_t;
#else
    using atomic_type = unsigned __int128;
#endif
    static constexpr const int kShift = sizeof(atomic_type) * CHAR_BIT / 2;
    static constexpr const atomic_type kMask = (static_cast<atomic_type>(1) << kShift) - 1;
    using diff_t = typename Cfg::difference_type;

 public:
    /**
     * Sets write/read pointers.
     */
    void set(diff_t w, diff_t r) {
        single_.w = w;
        single_.r = r;
        num_reading_.store(0, std::memory_order_relaxed);
    }

    /**
     * Gets the write pointer.
     */
    diff_t getWrite() const {
        return single_.w;
    }

    /**
     * Gets write/read pointers and increases the write pointer.
     */
    template <bool kAtomic>
    std::pair<diff_t, diff_t> incWrite() {
        if (kAtomic) {
            const auto p = __atomic_fetch_add(&all_, Cfg::kBlockSize, __ATOMIC_RELAXED);
            const diff_t w = p & kMask;
            const diff_t r = (p >> kShift);
            return {w, r};
        } else {
            const auto w = single_.w;
            single_.w += Cfg::kBlockSize;
            return {w, single_.r};
        }
    }

    /**
     * Gets write/read pointers, decreases the read pointer, and increases the read counter.
     */
    template <bool kAtomic>
    std::pair<diff_t, diff_t> decRead() {
        if (kAtomic) {
            // Must not be moved after the following fetch_sub, as that could lead to
            // another thread writing to our block, because isReading() returns false.
            num_reading_.fetch_add(1, std::memory_order_acquire);
            const auto p = __atomic_fetch_sub(&all_,
                    static_cast<atomic_type>(Cfg::kBlockSize) << kShift, __ATOMIC_RELAXED);
            const diff_t w = p & kMask;
            const diff_t r = (p >> kShift) & ~(Cfg::kBlockSize - 1);
            return {w, r};
        } else {
            const auto r = single_.r;
            single_.r -= Cfg::kBlockSize;
            return {single_.w, r};
        }
    }

    /**
     * Decreases the read counter.
     */
    void stopRead() {
        // synchronizes with threads wanting to write to this bucket
        num_reading_.fetch_sub(1, std::memory_order_release);
    }

    /**
     * Returns true if any thread is currently reading from here.
     */
    bool isReading() {
        // synchronize with threads currently reading from this bucket
        return num_reading_.load(std::memory_order_acquire) != 0;
    }

 private:
    struct Pointers { diff_t w, r; };
    union {
        atomic_type all_;
        Pointers single_;
    };
    std::atomic_int num_reading_;
};

}  // namespace detail
}  // namespace ips4o
