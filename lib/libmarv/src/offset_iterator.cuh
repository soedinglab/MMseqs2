#ifndef LIBMARV_OFFSET_ITERATOR_CUH
#define LIBMARV_OFFSET_ITERATOR_CUH

#include "namespace.hpp"
#include <cuda/std/iterator>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/counting_iterator.h>

LIBMARV_NAMESPACE_WITH_NESTING_BEGIN


template<class BaseOffsetIterator, class IndexType = int>
struct ReadOnlyCustomOffsetIteratorTransformOp{
    using value_type = typename cuda::std::iterator_traits<BaseOffsetIterator>::value_type;

    bool convertToZeroBased;
    BaseOffsetIterator parent;

    __host__ __device__
    value_type operator()(IndexType pos) const{
        value_type val = parent[pos];
        if(convertToZeroBased){
            val -= parent[0];
        }
        return val;
    }
};

template<class BaseOffsetIterator, class IndexType = int>
using ReadOnlyCustomOffsetIterator = decltype(
    thrust::make_transform_iterator(
        thrust::make_counting_iterator<IndexType>(0), 
        ReadOnlyCustomOffsetIteratorTransformOp<BaseOffsetIterator, IndexType>{}));

template<class BaseOffsetIterator, class IndexType = int>
__host__ __device__ 
auto makeCustomOffsetIterator(BaseOffsetIterator iter, bool convertToZeroBased){
    return thrust::make_transform_iterator(
        thrust::make_counting_iterator<IndexType>(0), 
        ReadOnlyCustomOffsetIteratorTransformOp<BaseOffsetIterator, IndexType>{convertToZeroBased, iter});
}


#if 0

template<class BaseOffsetIterator, class IndexType = int>
struct ReadOnlyCustomOffsetIterator{
    using value_type = typename cuda::std::iterator_traits<BaseOffsetIterator>::value_type;
    using difference_type = IndexType;
    using pointer = void;
    using reference = value_type&;
    using iterator_category = cuda::std::random_access_iterator_tag;


    bool convertToZeroBased = false;
    IndexType iteratorindex = 0;    
    BaseOffsetIterator parentIterator;

    ReadOnlyCustomOffsetIterator() = default;
    ReadOnlyCustomOffsetIterator(const ReadOnlyCustomOffsetIterator&) = default;
    ReadOnlyCustomOffsetIterator& operator=(const ReadOnlyCustomOffsetIterator&) = default;
    
    __host__ __device__
    value_type operator[](IndexType pos) const{
        value_type val = parentIterator[pos];
        if(convertToZeroBased){
            val -= parentIterator[0];
        }
        return val;
    }

    __host__ __device__
    value_type operator*() const{
        return operator[](iteratorindex);
    }

    __host__ __device__
    ReadOnlyCustomOffsetIterator operator+(IndexType steps) const{
        auto x = *this;
        x.iteratorindex += steps;
        return x;
    }

    __host__ __device__
    ReadOnlyCustomOffsetIterator operator-(IndexType steps) const{
        auto x = *this;
        x.iteratorindex -= steps;
        return x;
    }


    // prefix increment
    __host__ __device__
    ReadOnlyCustomOffsetIterator& operator++()
    {
        ++iteratorindex;
        return *this;
    }
    
    // postfix increment
    __host__ __device__
    ReadOnlyCustomOffsetIterator operator++(int)
    {
        auto old = *this;
        operator++();
        return old;
    }
    
    // prefix decrement
    __host__ __device__
    ReadOnlyCustomOffsetIterator& operator--()
    {
        --iteratorindex;
        return *this;
    }
    
    // postfix decrement
    __host__ __device__
    ReadOnlyCustomOffsetIterator operator--(int)
    {
        auto old = *this;
        operator--();
        return old;
    }

    __host__ __device__
    ReadOnlyCustomOffsetIterator& operator+=(IndexType steps){
        *this = operator+(steps);
        return *this;
    }

    __host__ __device__
    ReadOnlyCustomOffsetIterator& operator-=(IndexType steps){
        *this = operator-(steps);
        return *this;
    }

    __host__ __device__
    difference_type operator-(const ReadOnlyCustomOffsetIterator& rhs) const{
        return iteratorindex - rhs.iteratorindex;
    }


    __host__ __device__
    bool operator==(const ReadOnlyCustomOffsetIterator& rhs) const{
        if(convertToZeroBased == rhs.convertToZeroBased
            && iteratorindex == rhs.iteratorindex
            && parentIterator == rhs.parentIterator){

            return true;
        }
        return false;
    }

    __host__ __device__
    bool operator!=(const ReadOnlyCustomOffsetIterator& rhs) const{
        return !operator==(rhs);
    }
};

template<class BaseOffsetIterator, class IndexType = int>
__host__ __device__ 
auto makeCustomOffsetIterator(BaseOffsetIterator iter, bool convertToZeroBased){
    ReadOnlyCustomOffsetIterator<BaseOffsetIterator,IndexType> result;
    result.convertToZeroBased = convertToZeroBased;
    result.iteratorindex = 0;
    result.parentIterator = iter;

    return result;
}

#endif




LIBMARV_NAMESPACE_WITH_NESTING_END

#endif