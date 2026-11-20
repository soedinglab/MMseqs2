// CPU implementation of the thrust subset used by libmarv's host code; execution policies wait for pending device work.
#ifndef METAL_THRUST_IMPL_H
#define METAL_THRUST_IMPL_H

#include <cuda_runtime.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iterator>
#include <new>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace thrust {

struct execution_policy_tag {};
struct seq_t : execution_policy_tag {};
inline constexpr seq_t seq{};

struct device_t : execution_policy_tag {
    const device_t& on(cudaStream_t) const { marvBackendSync(); return *this; }
    template<class A> const device_t& operator()(const A&) const { return *this; }
};
inline constexpr device_t device{};

namespace cuda {
    struct par_t : execution_policy_tag {
        const par_t& on(cudaStream_t) const { marvBackendSync(); return *this; }
        template<class A> const par_t& operator()(const A&) const { return *this; }
    };
    inline constexpr par_t par{};
    inline constexpr par_t par_nosync{};
}

template<class P>
using is_exec = std::is_base_of<execution_policy_tag, std::decay_t<P>>;
#define THRUST_POL(P) class P, class = std::enable_if_t<::thrust::is_exec<P>::value>

template<class T = void> using greater = std::greater<T>;
using std::get;
using std::tuple;
using std::make_tuple;

template<class It>
typename std::iterator_traits<It>::difference_type distance(It a, It b) { return std::distance(a, b); }

template<class T>
class device_ptr {
    T* p_ = nullptr;
public:
    using value_type = std::remove_cv_t<T>;
    using difference_type = std::ptrdiff_t;
    using pointer = T*;
    using reference = T&;
    using iterator_category = std::random_access_iterator_tag;

    device_ptr() = default;
    explicit device_ptr(T* p) : p_(p) {}
    template<class U, class = std::enable_if_t<std::is_convertible<U*, T*>::value>>
    device_ptr(const device_ptr<U>& o) : p_(o.get()) {}

    T* get() const { return p_; }
    T& operator*() const { return *p_; }
    T* operator->() const { return p_; }
    T& operator[](difference_type i) const { return p_[i]; }
    device_ptr& operator++() { ++p_; return *this; }
    device_ptr operator++(int) { device_ptr t = *this; ++p_; return t; }
    device_ptr& operator--() { --p_; return *this; }
    device_ptr& operator+=(difference_type n) { p_ += n; return *this; }
    device_ptr& operator-=(difference_type n) { p_ -= n; return *this; }
    device_ptr operator+(difference_type n) const { return device_ptr(p_ + n); }
    device_ptr operator-(difference_type n) const { return device_ptr(p_ - n); }
    difference_type operator-(const device_ptr& o) const { return p_ - o.p_; }
    bool operator==(const device_ptr& o) const { return p_ == o.p_; }
    bool operator!=(const device_ptr& o) const { return p_ != o.p_; }
    bool operator<(const device_ptr& o) const { return p_ < o.p_; }
    explicit operator bool() const { return p_ != nullptr; }
};

template<class T> device_ptr<T> device_pointer_cast(T* p) { return device_ptr<T>(p); }
template<class T> device_ptr<T> device_pointer_cast(device_ptr<T> p) { return p; }
template<class T> T* raw_pointer_cast(device_ptr<T> p) { return p.get(); }
template<class T> T* raw_pointer_cast(T* p) { return p; }

template<class T>
struct device_malloc_allocator {
    using value_type = T;
    using pointer = device_ptr<T>;
    using const_pointer = device_ptr<const T>;
    using reference = T&;
    using const_reference = const T&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    template<class U> struct rebind { using other = device_malloc_allocator<U>; };

    device_malloc_allocator() = default;
    template<class U> device_malloc_allocator(const device_malloc_allocator<U>&) {}

    pointer allocate(size_type n) {
        T* p = nullptr;
        if (cudaMalloc((void**) &p, n * sizeof(T)) != cudaSuccess) throw std::bad_alloc();
        return pointer(p);
    }
    void deallocate(pointer p, size_type) { cudaFree(p.get()); }
};

struct no_init_t {};
inline constexpr no_init_t no_init{};

template<class T, class Alloc = device_malloc_allocator<T>>
class device_vector {
    Alloc alloc_;
    T* data_ = nullptr;
    std::size_t size_ = 0;
    std::size_t capacity_ = 0;

    void reallocate(std::size_t newCapacity) {
        T* fresh = newCapacity ? raw_pointer_cast(alloc_.allocate(newCapacity)) : nullptr;
        if (size_ && fresh) cudaMemcpy(fresh, data_, size_ * sizeof(T), cudaMemcpyDeviceToDevice);
        if (data_) alloc_.deallocate(typename Alloc::pointer(data_), capacity_);
        data_ = fresh; capacity_ = newCapacity;
    }
public:
    using value_type = T;
    using allocator_type = Alloc;
    using size_type = std::size_t;
    using pointer = device_ptr<T>;
    using iterator = T*;
    using const_iterator = const T*;
    using reference = T&;

    device_vector() : alloc_() {}
    explicit device_vector(const Alloc& a) : alloc_(a) {}
    explicit device_vector(size_type n) : alloc_() { resize(n); }
    device_vector(const device_vector&) = delete;
    device_vector& operator=(const device_vector&) = delete;
    device_vector(device_vector&& o) noexcept : alloc_(o.alloc_), data_(o.data_), size_(o.size_), capacity_(o.capacity_) { o.data_ = nullptr; o.size_ = o.capacity_ = 0; }
    device_vector& operator=(device_vector&& o) noexcept {
        if (this != &o) { release(); alloc_ = o.alloc_; data_ = o.data_; size_ = o.size_; capacity_ = o.capacity_; o.data_ = nullptr; o.size_ = o.capacity_ = 0; }
        return *this;
    }
    ~device_vector() { release(); }

    void release() { if (data_) { marvBackendSync(); alloc_.deallocate(typename Alloc::pointer(data_), capacity_); } data_ = nullptr; size_ = capacity_ = 0; }
    void reserve(size_type n) { if (n > capacity_) reallocate(n); }
    void resize(size_type n, no_init_t) { reserve(n); size_ = n; }
    void resize(size_type n) {
        const size_type old = size_;
        resize(n, no_init);
        if (n > old) cudaMemset(data_ + old, 0, (n - old) * sizeof(T));
    }
    size_type size() const { return size_; }
    pointer data() { return pointer(data_); }
    device_ptr<const T> data() const { return device_ptr<const T>(data_); }
    iterator begin() { return data_; }
    iterator end() { return data_ + size_; }
    const_iterator begin() const { return data_; }
    const_iterator end() const { return data_ + size_; }
    T& operator[](size_type i) { return data_[i]; }
    const T& operator[](size_type i) const { return data_[i]; }
};

template<class T>
class counting_iterator {
    T v_{};
public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = void;
    using reference = T;
    counting_iterator() = default;
    explicit counting_iterator(T v) : v_(v) {}
    T operator*() const { return v_; }
    T operator[](difference_type n) const { return (T)(v_ + n); }
    counting_iterator& operator++() { ++v_; return *this; }
    counting_iterator operator++(int) { auto t = *this; ++v_; return t; }
    counting_iterator& operator--() { --v_; return *this; }
    counting_iterator& operator+=(difference_type n) { v_ = (T)(v_ + n); return *this; }
    counting_iterator operator+(difference_type n) const { return counting_iterator((T)(v_ + n)); }
    counting_iterator operator-(difference_type n) const { return counting_iterator((T)(v_ - n)); }
    difference_type operator-(const counting_iterator& o) const { return (difference_type) v_ - (difference_type) o.v_; }
    bool operator==(const counting_iterator& o) const { return v_ == o.v_; }
    bool operator!=(const counting_iterator& o) const { return v_ != o.v_; }
    bool operator<(const counting_iterator& o) const { return v_ < o.v_; }
};
template<class T> counting_iterator<T> make_counting_iterator(T v) { return counting_iterator<T>(v); }

template<class T>
class constant_iterator {
    T v_{}; std::ptrdiff_t i_ = 0;
public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = void;
    using reference = T;
    constant_iterator() = default;
    explicit constant_iterator(T v) : v_(v) {}
    T operator*() const { return v_; }
    T operator[](difference_type) const { return v_; }
    constant_iterator& operator++() { ++i_; return *this; }
    constant_iterator operator+(difference_type n) const { auto t = *this; t.i_ += n; return t; }
    difference_type operator-(const constant_iterator& o) const { return i_ - o.i_; }
    bool operator==(const constant_iterator& o) const { return i_ == o.i_; }
    bool operator!=(const constant_iterator& o) const { return i_ != o.i_; }
};
template<class T> constant_iterator<T> make_constant_iterator(T v) { return constant_iterator<T>(v); }

struct discard_iterator {
    struct sink { template<class U> const sink& operator=(const U&) const { return *this; } };
    std::ptrdiff_t i_ = 0;
    using iterator_category = std::random_access_iterator_tag;
    using value_type = sink;
    using difference_type = std::ptrdiff_t;
    using pointer = void;
    using reference = sink;
    sink operator*() const { return sink{}; }
    sink operator[](difference_type) const { return sink{}; }
    discard_iterator& operator++() { ++i_; return *this; }
    discard_iterator operator++(int) { auto t = *this; ++i_; return t; }
    discard_iterator operator+(difference_type n) const { auto t = *this; t.i_ += n; return t; }
    difference_type operator-(const discard_iterator& o) const { return i_ - o.i_; }
    bool operator==(const discard_iterator& o) const { return i_ == o.i_; }
    bool operator!=(const discard_iterator& o) const { return i_ != o.i_; }
};
inline discard_iterator make_discard_iterator() { return discard_iterator{}; }

template<class It, class F>
class transform_iterator {
    It it_{}; mutable F f_{};
public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = std::decay_t<decltype(std::declval<F&>()(*std::declval<It&>()))>;
    using difference_type = std::ptrdiff_t;
    using pointer = void;
    using reference = value_type;
    transform_iterator() = default;
    transform_iterator(It it, F f) : it_(it), f_(f) {}
    reference operator*() const { return f_(*it_); }
    reference operator[](difference_type n) const { return f_(it_[n]); }
    transform_iterator& operator++() { ++it_; return *this; }
    transform_iterator operator++(int) { auto t = *this; ++it_; return t; }
    transform_iterator& operator--() { --it_; return *this; }
    transform_iterator& operator+=(difference_type n) { it_ += n; return *this; }
    transform_iterator operator+(difference_type n) const { return transform_iterator(it_ + n, f_); }
    transform_iterator operator-(difference_type n) const { return transform_iterator(it_ - n, f_); }
    difference_type operator-(const transform_iterator& o) const { return it_ - o.it_; }
    bool operator==(const transform_iterator& o) const { return it_ == o.it_; }
    bool operator!=(const transform_iterator& o) const { return it_ != o.it_; }
    bool operator<(const transform_iterator& o) const { return it_ < o.it_; }
    const It& base() const { return it_; }
    const F& functor() const { return f_; }
};
template<class It, class F>
transform_iterator<It, F> make_transform_iterator(It it, F f) { return transform_iterator<It, F>(it, f); }

template<class It, class F>
class transform_output_iterator {
    It it_{}; F f_{};
public:
    struct proxy {
        It it; mutable F f;
        template<class U> const proxy& operator=(const U& v) const { *it = f(v); return *this; }
    };
    using iterator_category = std::output_iterator_tag;
    using value_type = void;
    using difference_type = std::ptrdiff_t;
    using pointer = void;
    using reference = proxy;
    transform_output_iterator() = default;
    transform_output_iterator(It it, F f) : it_(it), f_(f) {}
    proxy operator*() const { return proxy{it_, f_}; }
    proxy operator[](difference_type n) const { return proxy{it_ + n, f_}; }
    transform_output_iterator& operator++() { ++it_; return *this; }
    transform_output_iterator operator++(int) { auto t = *this; ++it_; return t; }
    transform_output_iterator operator+(difference_type n) const { return transform_output_iterator(it_ + n, f_); }
};
template<class It, class F>
transform_output_iterator<It, F> make_transform_output_iterator(It it, F f) { return transform_output_iterator<It, F>(it, f); }

template<class... Its>
class zip_iterator {
    std::tuple<Its...> its_;
    using seq_t = std::index_sequence_for<Its...>;
    template<std::size_t... I> auto deref(std::index_sequence<I...>) const {
        return std::tuple<typename std::iterator_traits<Its>::reference...>(*std::get<I>(its_)...);
    }
    template<std::size_t... I> void adv(std::ptrdiff_t n, std::index_sequence<I...>) {
        int dummy[] = {(std::advance(std::get<I>(its_), n), 0)...}; (void) dummy;
    }
public:
    using iterator_category = std::random_access_iterator_tag;
    using reference = std::tuple<typename std::iterator_traits<Its>::reference...>;
    using value_type = std::tuple<typename std::iterator_traits<Its>::value_type...>;
    using difference_type = std::ptrdiff_t;
    using pointer = void;
    zip_iterator() = default;
    explicit zip_iterator(std::tuple<Its...> t) : its_(t) {}
    reference operator*() const { return deref(seq_t{}); }
    reference operator[](difference_type n) const { auto t = *this; t += n; return *t; }
    zip_iterator& operator++() { adv(1, seq_t{}); return *this; }
    zip_iterator operator++(int) { auto t = *this; adv(1, seq_t{}); return t; }
    zip_iterator& operator--() { adv(-1, seq_t{}); return *this; }
    zip_iterator& operator+=(difference_type n) { adv(n, seq_t{}); return *this; }
    zip_iterator& operator-=(difference_type n) { adv(-n, seq_t{}); return *this; }
    zip_iterator operator+(difference_type n) const { auto t = *this; t += n; return t; }
    zip_iterator operator-(difference_type n) const { auto t = *this; t -= n; return t; }
    difference_type operator-(const zip_iterator& o) const { return std::get<0>(its_) - std::get<0>(o.its_); }
    bool operator==(const zip_iterator& o) const { return std::get<0>(its_) == std::get<0>(o.its_); }
    bool operator!=(const zip_iterator& o) const { return !(*this == o); }
    bool operator<(const zip_iterator& o) const { return std::get<0>(its_) < std::get<0>(o.its_); }
    const std::tuple<Its...>& get_iterator_tuple() const { return its_; }
};
template<class... Its> zip_iterator<Its...> make_zip_iterator(std::tuple<Its...> t) { return zip_iterator<Its...>(t); }
template<class... Its> zip_iterator<Its...> make_zip_iterator(Its... its) { return zip_iterator<Its...>(std::make_tuple(its...)); }

template<class It, class T> void fill(It f, It l, const T& v) { std::fill(f, l, v); }
template<THRUST_POL(P), class It, class T> void fill(P, It f, It l, const T& v) { std::fill(f, l, v); }

template<class It, class T> void sequence(It f, It l, T init, T step) { for (T v = init; f != l; ++f, v = (T)(v + step)) *f = v; }
template<class It, class T> void sequence(It f, It l, T init) { thrust::sequence(f, l, init, T(1)); }
template<class It> void sequence(It f, It l) { using V = typename std::iterator_traits<It>::value_type; thrust::sequence(f, l, V(0), V(1)); }
template<THRUST_POL(P), class It> void sequence(P, It f, It l) { thrust::sequence(f, l); }
template<THRUST_POL(P), class It, class T> void sequence(P, It f, It l, T init) { thrust::sequence(f, l, init); }
template<THRUST_POL(P), class It, class T> void sequence(P, It f, It l, T init, T step) { thrust::sequence(f, l, init, step); }

template<class It, class F> F for_each(It f, It l, F fn) { return std::for_each(f, l, fn); }
template<THRUST_POL(P), class It, class F> F for_each(P, It f, It l, F fn) { return std::for_each(f, l, fn); }

template<class It, class T> It lower_bound(It f, It l, const T& v) { return std::lower_bound(f, l, v); }
template<class It, class T, class C> It lower_bound(It f, It l, const T& v, C c) { return std::lower_bound(f, l, v, c); }
template<THRUST_POL(P), class It, class T> It lower_bound(P, It f, It l, const T& v) { return std::lower_bound(f, l, v); }
template<THRUST_POL(P), class It, class T, class C> It lower_bound(P, It f, It l, const T& v, C c) { return std::lower_bound(f, l, v, c); }

template<class In, class Out, class Op> Out inclusive_scan(In f, In l, Out o, Op op) {
    using V = typename std::iterator_traits<Out>::value_type;
    bool first = true; V acc = V();
    for (; f != l; ++f, ++o) { acc = first ? (V)(*f) : (V) op(acc, (V)(*f)); first = false; *o = acc; }
    return o;
}
template<class In, class Out> Out inclusive_scan(In f, In l, Out o) { using V = typename std::iterator_traits<Out>::value_type; return thrust::inclusive_scan(f, l, o, std::plus<V>()); }
template<THRUST_POL(P), class In, class Out> Out inclusive_scan(P, In f, In l, Out o) { return thrust::inclusive_scan(f, l, o); }
template<THRUST_POL(P), class In, class Out, class Op> Out inclusive_scan(P, In f, In l, Out o, Op op) { return thrust::inclusive_scan(f, l, o, op); }

template<class It, class T, class Op> T reduce(It f, It l, T init, Op op) { for (; f != l; ++f) init = op(init, (T)(*f)); return init; }
template<class It, class T> T reduce(It f, It l, T init) { return thrust::reduce(f, l, init, std::plus<T>()); }
template<class It> typename std::iterator_traits<It>::value_type reduce(It f, It l) { using V = typename std::iterator_traits<It>::value_type; return thrust::reduce(f, l, V()); }
template<THRUST_POL(P), class It> auto reduce(P, It f, It l) { return thrust::reduce(f, l); }
template<THRUST_POL(P), class It, class T> T reduce(P, It f, It l, T init) { return thrust::reduce(f, l, init); }
template<THRUST_POL(P), class It, class T, class Op> T reduce(P, It f, It l, T init, Op op) { return thrust::reduce(f, l, init, op); }

template<class In, class Map, class Out> void scatter(In f, In l, Map map, Out out) { for (; f != l; ++f, ++map) out[*map] = *f; }
template<THRUST_POL(P), class In, class Map, class Out> void scatter(P, In f, In l, Map map, Out out) { thrust::scatter(f, l, map, out); }
template<class Map, class In, class Out> Out gather(Map mf, Map ml, In in, Out out) { for (; mf != ml; ++mf, ++out) *out = in[*mf]; return out; }
template<THRUST_POL(P), class Map, class In, class Out> Out gather(P, Map mf, Map ml, In in, Out out) { return thrust::gather(mf, ml, in, out); }

template<class In1, class In2, class Out> Out merge(In1 af, In1 al, In2 bf, In2 bl, Out o) { return std::merge(af, al, bf, bl, o); }
template<class In1, class In2, class Out, class C> Out merge(In1 af, In1 al, In2 bf, In2 bl, Out o, C c) { return std::merge(af, al, bf, bl, o, c); }
template<THRUST_POL(P), class In1, class In2, class Out> Out merge(P, In1 af, In1 al, In2 bf, In2 bl, Out o) { return std::merge(af, al, bf, bl, o); }
template<THRUST_POL(P), class In1, class In2, class Out, class C> Out merge(P, In1 af, In1 al, In2 bf, In2 bl, Out o, C c) { return std::merge(af, al, bf, bl, o, c); }

template<class It> void sort(It f, It l) { std::sort(f, l); }
template<class It, class C> void sort(It f, It l, C c) { std::sort(f, l, c); }
template<THRUST_POL(P), class It> void sort(P, It f, It l) { std::sort(f, l); }
template<THRUST_POL(P), class It, class C> void sort(P, It f, It l, C c) { std::sort(f, l, c); }

// zip iterators yield proxies: sort through an index permutation
template<class... Its, class C>
void sort(zip_iterator<Its...> f, zip_iterator<Its...> l, C comp) {
    const std::size_t n = (std::size_t)(l - f);
    if (n < 2) return;
    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), (std::size_t) 0);
    std::stable_sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return comp(f[a], f[b]); });
    using V = typename zip_iterator<Its...>::value_type;
    std::vector<V> tmp; tmp.reserve(n);
    for (std::size_t i = 0; i < n; ++i) tmp.push_back(V(f[idx[i]]));
    for (std::size_t i = 0; i < n; ++i) f[i] = tmp[i];
}
template<class... Its> void sort(zip_iterator<Its...> f, zip_iterator<Its...> l) {
    using V = typename zip_iterator<Its...>::value_type;
    sort(f, l, [](const V& a, const V& b) { return a < b; });
}
template<THRUST_POL(P), class... Its, class C> void sort(P, zip_iterator<Its...> f, zip_iterator<Its...> l, C c) { thrust::sort(f, l, c); }
template<THRUST_POL(P), class... Its> void sort(P, zip_iterator<Its...> f, zip_iterator<Its...> l) { thrust::sort(f, l); }

namespace detail {

// stable LSD radix sort of 4-byte keys; Desc flips the key order
template<class K, class V, bool Desc>
void radix_sort_by_key(K* keys, std::size_t n, V* vals) {
    if (n < 2) return;
    constexpr int Bits = 8, Passes = 4, Buckets = 1 << Bits;
    constexpr std::uint32_t Mask = Buckets - 1;
    auto enc = [](K k) -> std::uint32_t {
        std::uint32_t u; std::memcpy(&u, &k, 4);
        if constexpr (std::is_floating_point<K>::value) u ^= (u >> 31) ? 0xFFFFFFFFu : 0x80000000u;
        else if constexpr (std::is_signed<K>::value) u ^= 0x80000000u;
        if constexpr (Desc) u = ~u;
        return u;
    };
    std::vector<K> kbuf(n); std::vector<V> vbuf(n);
    K* ks = keys; K* kd = kbuf.data();
    V* vs = vals; V* vd = vbuf.data();
    std::vector<std::size_t> hist((std::size_t) Passes * Buckets, 0);
    for (std::size_t i = 0; i < n; ++i) {
        const std::uint32_t u = enc(ks[i]);
        for (int p = 0; p < Passes; ++p) hist[(std::size_t) p * Buckets + ((u >> (Bits * p)) & Mask)]++;
    }
    bool skip[Passes];
    for (int p = 0; p < Passes; ++p) {
        std::size_t sum = 0; bool sk = false;
        for (int b = 0; b < Buckets; ++b) { std::size_t c = hist[(std::size_t) p * Buckets + b]; if (c == n) sk = true; hist[(std::size_t) p * Buckets + b] = sum; sum += c; }
        skip[p] = sk;
    }
    for (int p = 0; p < Passes; ++p) {
        if (skip[p]) continue;
        const int shift = Bits * p;
        std::size_t* h = hist.data() + (std::size_t) p * Buckets;
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t pos = h[(enc(ks[i]) >> shift) & Mask]++;
            kd[pos] = ks[i]; vd[pos] = vs[i];
        }
        std::swap(ks, kd); std::swap(vs, vd);
    }
    if (ks != keys) { std::copy(ks, ks + n, keys); std::copy(vs, vs + n, vals); }
}

template<class KeyIt, class ValIt, class Comp, bool Stable>
void sort_by_key_impl(KeyIt kf, KeyIt kl, ValIt vf, Comp comp) {
    using K = typename std::iterator_traits<KeyIt>::value_type;
    using V = typename std::iterator_traits<ValIt>::value_type;
    const std::size_t n = (std::size_t) std::distance(kf, kl);
    constexpr bool isGreater = std::is_same<Comp, std::greater<K>>::value || std::is_same<Comp, std::greater<void>>::value;
    constexpr bool isLess = std::is_same<Comp, std::less<K>>::value || std::is_same<Comp, std::less<void>>::value;
    constexpr bool radixable = std::is_pointer<KeyIt>::value && std::is_pointer<ValIt>::value && std::is_arithmetic<K>::value && sizeof(K) == 4 && (isGreater || isLess);
    if constexpr (radixable) {
        if (n) radix_sort_by_key<K, V, isGreater>(&*kf, n, &*vf);
        return;
    }
    struct Rec { K key; V val; };
    std::vector<Rec> buf(n);
    for (std::size_t i = 0; i < n; ++i) buf[i] = Rec{kf[i], vf[i]};
    auto cmp = [&](const Rec& a, const Rec& b) { return comp(a.key, b.key); };
    if (Stable) std::stable_sort(buf.begin(), buf.end(), cmp); else std::sort(buf.begin(), buf.end(), cmp);
    for (std::size_t i = 0; i < n; ++i) { kf[i] = buf[i].key; vf[i] = buf[i].val; }
}

}

template<class KeyIt, class ValIt> void stable_sort_by_key(KeyIt kf, KeyIt kl, ValIt vf) { using K = typename std::iterator_traits<KeyIt>::value_type; detail::sort_by_key_impl<KeyIt, ValIt, std::less<K>, true>(kf, kl, vf, std::less<K>()); }
template<class KeyIt, class ValIt, class C> void stable_sort_by_key(KeyIt kf, KeyIt kl, ValIt vf, C c) { detail::sort_by_key_impl<KeyIt, ValIt, C, true>(kf, kl, vf, c); }
template<THRUST_POL(P), class KeyIt, class ValIt> void stable_sort_by_key(P, KeyIt kf, KeyIt kl, ValIt vf) { thrust::stable_sort_by_key(kf, kl, vf); }
template<THRUST_POL(P), class KeyIt, class ValIt, class C> void stable_sort_by_key(P, KeyIt kf, KeyIt kl, ValIt vf, C c) { thrust::stable_sort_by_key(kf, kl, vf, c); }

#undef THRUST_POL

}

#endif
