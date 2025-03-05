#ifndef UTIL_CUH
#define UTIL_CUH

#include "config.hpp"
#include "hpc_helpers/all_helpers.cuh"

#include <thrust/device_malloc_allocator.h>
#include <thrust/sequence.h>
#include <thrust/fill.h>



namespace cudasw4{


template<class T, int numRows, int numColumns>
struct SharedPSSM_singletile{
    static_assert(16 % sizeof(T) == 0);
    //each row is padded to 16 bytes
    static constexpr int numPaddedColumns = SDIV(numColumns, 16/sizeof(T)) * 16/sizeof(T);
    alignas(16) T data[numRows][numPaddedColumns];
};

template<int groupsize, int factor_>
struct SmemIndexCalculator{
    static constexpr int factor = factor_;

    __device__
    int getIndex(int ithChunkOfFour){
        constexpr int groupsizeForSmem = factor*groupsize;
        const int groupLaneForSmem = threadIdx.x % groupsizeForSmem;
        return 4*(groupLaneForSmem+ithChunkOfFour*groupsizeForSmem);
    }
};


template<class T> struct Vectorized2;
template<> struct Vectorized2<int>{ using type = int2; };
template<> struct Vectorized2<float>{ using type = float2; };

template<class T> struct Vectorized4;
template<> struct Vectorized4<int>{ using type = int4; };
template<> struct Vectorized4<float>{ using type = float4; };

template<class Score, class Extra>
struct ScoreWithExtra{
    Score score;
    Extra extra;

    ScoreWithExtra() = default;

    __host__ __device__
    ScoreWithExtra(Score s, Extra e) : score(s), extra(e){}

    __host__ __device__
    Score getScore() const{
        return score;
    }

    __host__ __device__
    Extra getExtra() const{
        return extra;
    }
};






template <class T>
struct thrust_async_allocator : public thrust::device_malloc_allocator<T> {
public:
    using Base      = thrust::device_malloc_allocator<T>;
    using pointer   = typename Base::pointer;
    using size_type = typename Base::size_type;

    thrust_async_allocator(cudaStream_t stream_) : stream{stream_} {}

    pointer allocate(size_type num){
        //std::cout << "allocate " << num << "\n";
        T* result = nullptr;
        cudaError_t status = cudaMallocAsync(&result, sizeof(T) * num, stream);
        if(status != cudaSuccess){
            throw std::runtime_error("thrust_async_allocator error allocate");
        }
        return thrust::device_pointer_cast(result);
    }

    void deallocate(pointer ptr, size_type /*num*/){
        //std::cout << "deallocate \n";
        cudaError_t status = cudaFreeAsync(thrust::raw_pointer_cast(ptr), stream);
        if(status != cudaSuccess){
            throw std::runtime_error("thrust_async_allocator error deallocate");
        }
    }

private:
    cudaStream_t stream;
};

template <class T>
struct thrust_preallocated_single_allocator : public thrust::device_malloc_allocator<T> {
public:
    using Base      = thrust::device_malloc_allocator<T>;
    using pointer   = typename Base::pointer;
    using size_type = typename Base::size_type;

    thrust_preallocated_single_allocator(void* ptr, size_t size) : preallocated{ptr}, preallocatedSize{size} {}

    pointer allocate(size_type num){
        if(!free){
            throw std::runtime_error("thrust_async_allocator error allocate");
        }else{
            if(sizeof(T) * num <= preallocatedSize){
                T* result = (T*)preallocated;
                free = false;
                return thrust::device_pointer_cast(result);
            }else{
                throw std::runtime_error("thrust_async_allocator error allocate");
            }
        }
    }

    void deallocate(pointer ptr, size_type /*num*/){
        if(free){
            throw std::runtime_error("thrust_async_allocator error deallocate");
        }else{
            T* result = thrust::raw_pointer_cast(ptr);
            if((void*) result != preallocated){
                throw std::runtime_error("thrust_async_allocator error deallocate");
            }
            free = true;
        }
    }

private:
    bool free = true;
    void* preallocated;
    size_t preallocatedSize;
    cudaStream_t stream;
};

//Call cudaSetDevice on destruction
struct RevertDeviceId{
    RevertDeviceId(){
        cudaGetDevice(&id);
    }
    RevertDeviceId(int id_) : id(id_){}
    ~RevertDeviceId(){
        cudaSetDevice(id);
    }
    int id;
};


//template<size_t size>
struct TopNMaximaArray{
    struct Ref{
        size_t index;
        size_t indexOffset;
        float* d_scores;
        ReferenceIdT* d_indices;
        size_t size;

        __device__
        Ref& operator=(float newscore){     
            d_scores[index] = newscore;
            d_indices[index] = indexOffset + index;
            return *this;
        }
    };

    TopNMaximaArray(float* d_scores_, ReferenceIdT* d_indices_, size_t offset, size_t size_)
        : indexOffset(offset), d_scores(d_scores_), d_indices(d_indices_), size(size_){}

    template<class Index>
    __device__
    Ref operator[](Index index) const{
        Ref r;
        r.index = index;
        r.indexOffset = indexOffset;
        r.d_scores = d_scores;
        r.d_indices = d_indices;
        r.size = size;
        return r;
    }

    void setAllScoresToZero(cudaStream_t stream){
        thrust::fill(
            thrust::cuda::par_nosync.on(stream),
            d_scores,
            d_scores + size,
            0
        );
        thrust::sequence(
            thrust::cuda::par_nosync.on(stream),
            d_indices,
            d_indices + size,
            ReferenceIdT(0)
        );
    }

    size_t indexOffset = 0;
    float* d_scores;
    ReferenceIdT* d_indices;
    size_t size;
};


template<class ExtraData>
struct TopNMaximaArrayWithExtra{
    struct Ref{
        size_t index;
        size_t indexOffset;
        float* d_scores;
        ReferenceIdT* d_indices;
        ExtraData* d_extras;
        size_t size;

        template<class Payload>
        __device__
        Ref& operator=(const Payload& payload){     
            d_scores[index] = payload.getScore();
            d_indices[index] = indexOffset + index;
            d_extras[index] = payload.getExtra();
            return *this;
        }
    };

    TopNMaximaArrayWithExtra(float* d_scores_, ReferenceIdT* d_indices_, ExtraData* d_extras_, size_t offset, size_t size_)
        : indexOffset(offset), d_scores(d_scores_), d_indices(d_indices_), d_extras(d_extras_), size(size_){}

    template<class Index>
    __device__
    Ref operator[](Index index) const{
        Ref r;
        r.index = index;
        r.indexOffset = indexOffset;
        r.d_scores = d_scores;
        r.d_indices = d_indices;
        r.d_extras = d_extras;
        r.size = size;
        return r;
    }

    void setAllScoresToZero(cudaStream_t stream){
        thrust::fill(
            thrust::cuda::par_nosync.on(stream),
            d_scores,
            d_scores + size,
            0
        );
        thrust::sequence(
            thrust::cuda::par_nosync.on(stream),
            d_indices,
            d_indices + size,
            ReferenceIdT(0)
        );
        thrust::fill(
            thrust::cuda::par_nosync.on(stream),
            d_extras,
            d_extras + size,
            ExtraData{}
        );
    }

    size_t indexOffset = 0;
    float* d_scores;
    ReferenceIdT* d_indices;
    ExtraData* d_extras;
    size_t size;
};

}

#endif