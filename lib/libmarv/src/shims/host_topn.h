#pragma once
#include <algorithm>
#include <utility>
#include <vector>

namespace libmarv {
// Heap selection of the best topN (score desc, id asc) pairs; scores/ids must be host-visible.
template<class Id>
void selectHostTopN(int* scores, Id* ids, int count, int topN) {
    const int keep = std::min(count, topN);
    auto better = [](const auto& a, const auto& b) {
        return a.first != b.first ? a.first > b.first : a.second < b.second;
    };
    std::vector<std::pair<int, Id>> heap;
    heap.reserve(keep);
    for (int i = 0; i < keep; ++i) heap.emplace_back(scores[i], ids[i]);
    std::make_heap(heap.begin(), heap.end(), better);
    for (int i = keep; i < count; ++i) {
        const std::pair<int, Id> candidate{scores[i], ids[i]};
        if (keep > 0 && better(candidate, heap.front())) {
            std::pop_heap(heap.begin(), heap.end(), better);
            heap.back() = candidate;
            std::push_heap(heap.begin(), heap.end(), better);
        }
    }
    std::sort_heap(heap.begin(), heap.end(), better);
    for (int i = 0; i < keep; ++i) {
        scores[i] = heap[i].first;
        ids[i] = heap[i].second;
    }
}
} // namespace libmarv
