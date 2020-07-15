/******************************************************************************
 * ips4o/classifier.hpp
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

#include <type_traits>
#include <utility>

#include "ips4o_fwd.hpp"
#include "utils.hpp"

namespace ips4o {
namespace detail {

/**
 * Branch-free classifier.
 */
template <class Cfg>
class Sorter<Cfg>::Classifier {
    using iterator = typename Cfg::iterator;
    using value_type = typename Cfg::value_type;
    using bucket_type = typename Cfg::bucket_type;
    using less = typename Cfg::less;

 public:
    Classifier(less comp) : comp_(std::move(comp)) {}

    ~Classifier() {
        if (log_buckets_) cleanup();
    }

    /**
     * Calls destructors on splitter elements.
     */
    void reset() {
        if (log_buckets_) cleanup();
    }

    /**
     * The sorted array of splitters, to be filled externally.
     */
    value_type* getSortedSplitters() {
        return static_cast<value_type*>(static_cast<void*>(sorted_storage_));
    }

    /**
     * The comparison operator.
     */
    less getComparator() const { return comp_; }

    /**
     * Builds the tree from the sorted splitters.
     */
    void build(int log_buckets) {
        log_buckets_ = log_buckets;
        num_buckets_ = 1 << log_buckets;
        const auto num_splitters = (1 << log_buckets) - 1;
        IPS4O_ASSUME_NOT(getSortedSplitters() + num_splitters == nullptr);
        new (getSortedSplitters() + num_splitters) value_type(getSortedSplitters()[num_splitters - 1]);
        build(getSortedSplitters(), getSortedSplitters() + num_splitters, 1);
    }

    /**
     * Classifies a single element.
     */
    template <bool kEqualBuckets>
    bucket_type classify(const value_type& value) const {
        const int log_buckets = log_buckets_;
        const bucket_type num_buckets = num_buckets_;
        IPS4O_ASSUME_NOT(log_buckets < 1);
        IPS4O_ASSUME_NOT(log_buckets > 9);

        bucket_type b = 1;
        for (int l = 0; l < log_buckets; ++l)
            b = 2 * b + comp_(splitter(b), value);
        if (kEqualBuckets)
            b = 2 * b + !comp_(value, sortedSplitter(b - num_buckets));
        return b - (kEqualBuckets ? 2 * num_buckets : num_buckets);
    }

    /**
     * Classifies all elements using a callback.
     */
    template <bool kEqualBuckets, class Yield>
    void classify(iterator begin, iterator end, Yield&& yield) const {
        switch (log_buckets_) {
            case 1: classifyUnrolled<kEqualBuckets, 1>(begin, end, std::forward<Yield>(yield)); break;
            case 2: classifyUnrolled<kEqualBuckets, 2>(begin, end, std::forward<Yield>(yield)); break;
            case 3: classifyUnrolled<kEqualBuckets, 3>(begin, end, std::forward<Yield>(yield)); break;
            case 4: classifyUnrolled<kEqualBuckets, 4>(begin, end, std::forward<Yield>(yield)); break;
            case 5: classifyUnrolled<kEqualBuckets, 5>(begin, end, std::forward<Yield>(yield)); break;
            case 6: classifyUnrolled<kEqualBuckets, 6>(begin, end, std::forward<Yield>(yield)); break;
            case 7: classifyUnrolled<kEqualBuckets, 7>(begin, end, std::forward<Yield>(yield)); break;
            case 8: classifyUnrolled<kEqualBuckets, 8>(begin, end, std::forward<Yield>(yield)); break;
            default: IPS4O_ASSUME_NOT(true);
        }
    }

    /**
     * Classifies all elements using a callback.
     */
    template <bool kEqualBuckets, int kLogBuckets, class Yield>
    void classifyUnrolled(iterator begin, const iterator end, Yield&& yield) const {
        constexpr const bucket_type kNumBuckets = 1l << (kLogBuckets + kEqualBuckets);
        constexpr const int kUnroll = Cfg::kUnrollClassifier;
        IPS4O_ASSUME_NOT(begin >= end);
        IPS4O_ASSUME_NOT(begin > (end - kUnroll));

        bucket_type b[kUnroll];
        for (auto cutoff = end - kUnroll; begin <= cutoff; begin += kUnroll) {
            for (int i = 0; i < kUnroll; ++i)
                b[i] = 1;

            for (int l = 0; l < kLogBuckets; ++l)
                for (int i = 0; i < kUnroll; ++i)
                    b[i] = 2 * b[i] + comp_(splitter(b[i]), begin[i]);

            if (kEqualBuckets)
                for (int i = 0; i < kUnroll; ++i)
                    b[i] = 2 * b[i] + !comp_(begin[i], sortedSplitter(b[i] - kNumBuckets / 2));

            for (int i = 0; i < kUnroll; ++i)
                yield(b[i] - kNumBuckets, begin + i);
        }

        IPS4O_ASSUME_NOT(begin > end);
        for (; begin != end; ++begin) {
            bucket_type b = 1;
            for (int l = 0; l < kLogBuckets; ++l)
                b = 2 * b + comp_(splitter(b), *begin);
            if (kEqualBuckets)
                b = 2 * b + !comp_(*begin, sortedSplitter(b - kNumBuckets / 2));
            yield(b - kNumBuckets, begin);
        }
    }

 private:
    const value_type& splitter(bucket_type i) const {
        return static_cast<const value_type*>(static_cast<const void*>(storage_))[i];
    }

    const value_type& sortedSplitter(bucket_type i) const {
        return static_cast<const value_type*>(static_cast<const void*>(sorted_storage_))[i];
    }

    value_type* data() {
        return static_cast<value_type*>(static_cast<void*>(storage_));
    }

    /**
     * Recursively builds the tree.
     */
    void build(const value_type* const left, const value_type* const right, const bucket_type pos) {
        const auto mid = left + (right - left) / 2;
        IPS4O_ASSUME_NOT(data() + pos == nullptr);
        new (data() + pos) value_type(*mid);
        if (2 * pos < num_buckets_) {
            build(left, mid, 2 * pos);
            build(mid, right, 2 * pos + 1);
        }
    }

    /**
     * Destructs splitters.
     */
    void cleanup() {
        auto p = data() + 1;
        auto q = getSortedSplitters();
        for (int i = num_buckets_ - 1; i; --i) {
            p++->~value_type();
            q++->~value_type();
        }
        q->~value_type();
        log_buckets_ = 0;
    }

    // Filled from 1 to num_buckets_
    std::aligned_storage_t<sizeof(value_type), alignof(value_type)> storage_[Cfg::kMaxBuckets / 2];
    // Filled from 0 to num_buckets_, last one is duplicated
    std::aligned_storage_t<sizeof(value_type), alignof(value_type)> sorted_storage_[Cfg::kMaxBuckets / 2];
    int log_buckets_ = 0;
    bucket_type num_buckets_ = 0;
    less comp_;
};

}  // namespace detail
}  // namespace ips4o
