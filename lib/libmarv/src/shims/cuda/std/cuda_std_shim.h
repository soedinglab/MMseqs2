// cuda::std replacement: the standard library plus a minimal mdspan with operator().
#ifndef CUDA_STD_SHIM_H
#define CUDA_STD_SHIM_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iterator>
#include <span>
#include <tuple>
#include <type_traits>
#include <utility>

namespace cuda {
namespace std {

using ::std::span;
using ::std::array;
using ::std::tuple;
using ::std::get;
using ::std::plus;
using ::std::distance;
using ::std::iterator_traits;
using ::std::random_access_iterator_tag;
using ::std::true_type;
using ::std::false_type;
using ::std::uint8_t;

template<class IndexType, ::std::size_t Rank>
struct dextents {
    using index_type = IndexType;
    ::std::array<IndexType, Rank> ext{};
    static constexpr ::std::size_t rank() { return Rank; }
    IndexType extent(::std::size_t r) const { return ext[r]; }
};

template<class T, class Extents>
struct mdspan {
    using element_type = T;
    using index_type = typename Extents::index_type;
    T* ptr = nullptr;
    Extents ext{};

    mdspan() = default;
    template<class... Ix>
    mdspan(T* p, Ix... e) : ptr(p) {
        static_assert(sizeof...(Ix) == Extents::rank(), "wrong number of extents");
        index_type tmp[] = {static_cast<index_type>(e)...};
        for (::std::size_t i = 0; i < Extents::rank(); i++) ext.ext[i] = tmp[i];
    }
    static constexpr ::std::size_t rank() { return Extents::rank(); }
    index_type extent(::std::size_t r) const { return ext.extent(r); }
    T* data_handle() const { return ptr; }
    ::std::size_t size() const { ::std::size_t n = 1; for (::std::size_t i = 0; i < Extents::rank(); i++) n *= (::std::size_t) ext.extent(i); return n; }
    template<class... Ix>
    T& operator()(Ix... idx) const {
        static_assert(sizeof...(Ix) == Extents::rank(), "wrong number of indices");
        index_type tmp[] = {static_cast<index_type>(idx)...};
        ::std::size_t off = 0;
        for (::std::size_t i = 0; i < Extents::rank(); i++) off = off * (::std::size_t) ext.extent(i) + (::std::size_t) tmp[i];
        return ptr[off];
    }
    T& operator[](index_type i) const { return ptr[i]; }
};

}

template<class T, class U>
constexpr auto ceil_div(T a, U b) { return (a + b - 1) / b; }
template<class T, class U>
constexpr auto round_up(T a, U b) { return ((a + b - 1) / b) * b; }

template<class R, class F>
constexpr F proclaim_return_type(F f) { return f; }

}

#endif
