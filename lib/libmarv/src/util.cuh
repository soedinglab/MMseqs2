#ifndef UTIL_CUH
#define UTIL_CUH

#if defined(__HIPCC__)
    #include "hip/hip_runtime.h"
#endif    

#include "cuda_hip_rename.h"
#include "stream_ordered_allocation_helper.cuh"

#include "config.hpp"
#include "hpc_helpers/all_helpers.cuh"
#include "custom_score_types.cuh"

#include <thrust/device_vector.h>
#include <thrust/device_malloc_allocator.h>
#include <thrust/sequence.h>
#include <thrust/fill.h>

// #if defined(__CUDACC__)
#include <cuda/std/array>
#include <cuda/cmath>
// #elif defined(__HIPCC__)
// #include <hip/std/array>
// #endif



#include "namespace.hpp"
LIBMARV_NAMESPACE_WITH_NESTING_BEGIN


template<class T, T factor>
struct RoundToNextMultiple{
    __host__ __device__ 
    T operator()(const T& value){
        return SDIV(value, factor) * factor;
    }        
};


template<class T, int numRows_, int numColumns_>
struct SharedPSSM_singletile{
    using value_type = T;

    static_assert(16 % sizeof(T) == 0);
    static constexpr int numRows = numRows_;
    static constexpr int numColumns = numColumns_;
    //each row is padded to 16 bytes
    static constexpr int numPaddedColumns = SDIV(numColumns, 16/sizeof(T)) * 16/sizeof(T);

    alignas(16) T data[numRows][numPaddedColumns];
};

template<int groupsize, class PssmScoreType, int elementsPerLoad>
struct PssmReplicationFactor{
    static_assert(32 % groupsize == 0);

    __host__ __device__
    static constexpr int value(){
        constexpr int bytesPerLoadPerThread = sizeof(PssmScoreType) * elementsPerLoad;
        constexpr int bytesPerLoadPerGroup = groupsize * bytesPerLoadPerThread;
        constexpr int smemTransactionSize = 128;
        static_assert(bytesPerLoadPerGroup > smemTransactionSize || smemTransactionSize % bytesPerLoadPerGroup == 0);

        return SDIV(smemTransactionSize, bytesPerLoadPerGroup);
    };
};

// template<int groupsize, int factor_>
// struct SmemIndexCalculator{
//     static constexpr int factor = factor_;

//     __device__
//     int getIndex(int ithChunkOfFour){
//         constexpr int groupsizeForSmem = factor*groupsize;
//         const int groupLaneForSmem = threadIdx.x % groupsizeForSmem;
//         return 4*(groupLaneForSmem+ithChunkOfFour*groupsizeForSmem);
//     }
// };

template<int groupsize, int replicationFactor, int elementsPerVectorLoad = 4>
struct PssmSmemIndexCalculator{
    static constexpr int factor = replicationFactor;

    __device__
    int getIndexOfFirstElementInVectorType(int ithGroupwideVectorLoad){
        constexpr int groupsizeForSmem = factor*groupsize;
        const int groupLaneForSmem = threadIdx.x % groupsizeForSmem;
        return elementsPerVectorLoad*(groupLaneForSmem+ithGroupwideVectorLoad*groupsizeForSmem);
    }

    __host__ __device__
    cuda::std::array<int, replicationFactor> computeReplicatedColumnIndices(int inputColumnIndex){
        const int vectorIndex = inputColumnIndex / elementsPerVectorLoad;
        const int elementIndexWithinVector = inputColumnIndex % elementsPerVectorLoad;
        const int ithGroupwideVectorLoad = vectorIndex / groupsize;
        const int vectorPositionInGroupwideVectorLoad = vectorIndex % groupsize;

        cuda::std::array<int, replicationFactor> result{};

        for(int i = 0; i < replicationFactor; i++){
            const int ith_outputVectorIndex = (ithGroupwideVectorLoad*replicationFactor*groupsize + i*groupsize) + vectorPositionInGroupwideVectorLoad;
            const int ith_outputElementIndex = elementsPerVectorLoad * ith_outputVectorIndex + elementIndexWithinVector;
            result[i] = ith_outputElementIndex;
        }

        return result;
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

    template <class Other>
    struct rebind {
        using other = thrust_async_allocator<Other>;
    };

    thrust_async_allocator(cudaStream_t stream_) : stream{stream_} {}

    template <class Other>
    thrust_async_allocator(thrust_async_allocator<Other> other) : stream{other.stream} {}


    pointer allocate(size_type num){
        //std::cout << "allocate " << num << "\n";
        T* result = nullptr;

        cudaError_t status = helper_gpuMallocAsync((void**)&result, sizeof(T) * num, stream);

        if(status != cudaSuccess){
            std::string msg = "thrust_async_allocator error allocate: ";
            msg += cudaGetErrorString(status);
            throw std::runtime_error(msg);
        }
        return thrust::device_pointer_cast(result);
    }

    void deallocate(pointer ptr, size_type /*num*/){
        //std::cout << "deallocate \n";
        
        cudaError_t status = helper_gpuFreeAsync((void*)thrust::raw_pointer_cast(ptr), stream);

        if(status != cudaSuccess){
            throw std::runtime_error("thrust_async_allocator error deallocate");
        }
    }

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

#if 0

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

#endif

struct PositionsIterator{
    ReferenceIdT* ptr;
    ReferenceIdT offset;

    static PositionsIterator fromPointer(ReferenceIdT* p){
        return {p, 0};
    }

    static PositionsIterator fromCountingIterator(ReferenceIdT startOffset){
        return {nullptr, startOffset};
    }

    __host__ __device__
    ReferenceIdT operator[](size_t i) const{
        return ptr ? ptr[i] : offset + static_cast<ReferenceIdT>(i);
    }

    __host__ __device__
    PositionsIterator operator+(size_t n) const{
        if(ptr) return {ptr + n, 0};
        return {nullptr, offset + static_cast<ReferenceIdT>(n)};
    }
};



template<class ScoreType> struct ScalarScoreType{};
template<> struct ScalarScoreType<half>{ using type = half; };
template<> struct ScalarScoreType<short>{ using type = short; };
template<> struct ScalarScoreType<int>{ using type = int; };
template<> struct ScalarScoreType<float>{ using type = float; };
template<> struct ScalarScoreType<char2>{ using type = char; };
template<> struct ScalarScoreType<half2>{ using type = half; };
template<> struct ScalarScoreType<short2>{ using type = short; };
template<> struct ScalarScoreType<int2>{ using type = int; };
template<> struct ScalarScoreType<float2>{ using type = float; };
template<> struct ScalarScoreType<char4>{ using type = char; };
template<> struct ScalarScoreType<short4>{ using type = short; };
template<> struct ScalarScoreType<int4>{ using type = int; };
template<> struct ScalarScoreType<float4>{ using type = float; };
template<> struct ScalarScoreType<ScoreType_u8x4>{ using type = cuda::std::uint8_t; };
template<> struct ScalarScoreType<__nv_fp8x2_e4m3>{ using type = __nv_fp8_e4m3; };

template<class T>
struct IsScalarType : std::false_type {};
template<> struct IsScalarType<float> : std::true_type {};
template<> struct IsScalarType<int> : std::true_type {};
template<> struct IsScalarType<half> : std::true_type {};
template<> struct IsScalarType<short> : std::true_type {};
template<> struct IsScalarType<__nv_fp8_e4m3> : std::true_type {};
template<> struct IsScalarType<__nv_fp8_e5m2> : std::true_type {};
template<> struct IsScalarType<char> : std::true_type {};
template<> struct IsScalarType<int8_t> : std::true_type {};
template<> struct IsScalarType<uint8_t> : std::true_type {};

template<class T>
struct IsVec2Type : cuda::std::false_type {};
template<> struct IsVec2Type<float2> : cuda::std::true_type {};
template<> struct IsVec2Type<int2> : cuda::std::true_type {};
template<> struct IsVec2Type<half2> : cuda::std::true_type {};
template<> struct IsVec2Type<short2> : cuda::std::true_type {};
template<> struct IsVec2Type<char2> : cuda::std::true_type {};

template<class T>
struct IsVec4Type : cuda::std::false_type {};
template<> struct IsVec4Type<float4> : cuda::std::true_type {};
template<> struct IsVec4Type<int4> : cuda::std::true_type {};
template<> struct IsVec4Type<char4> : cuda::std::true_type {};
template<> struct IsVec4Type<short4> : cuda::std::true_type {};
template<> struct IsVec4Type<ScoreType_u8x4> : cuda::std::true_type {};


enum class Approach : int{
    hardcodedzero = 0,
    kernelparamzero = 1,
    unused = 2,
    half2fmarelu = 3,
    floatloadhalf = 4,
};
enum class Datatype : int{
    Half2 = 0,
    Short2 = 1,
    UInt8x4 = 2,
    Float = 3,
    Int = 4,
    Invalid = 9999
};

template<class Type>
constexpr Datatype mapDatatypeToDatatypeEnum(){
    if constexpr(std::is_same_v<Type, half2>){
        return Datatype::Half2;
    }else if constexpr(std::is_same_v<Type, short2>){
        return Datatype::Short2;
    }else if constexpr(std::is_same_v<Type, ScoreType_u8x4>){
        return Datatype::UInt8x4;
    }else if constexpr(std::is_same_v<Type, float>){
        return Datatype::Float;
    }else if constexpr(std::is_same_v<Type, int>){
        return Datatype::Int;
    }else{
        return Datatype::Invalid;
    }
}

template<Datatype datatype>
struct ScoreTypeFromDatatypeE{
    using type = typename std::conditional<
        datatype == Datatype::Half2,
        half2,
        typename std::conditional<
            datatype == Datatype::Short2,
            short2,
            typename std::conditional<
                datatype == Datatype::UInt8x4,
                ScoreType_u8x4,
                typename std::conditional<
                    datatype == Datatype::Float,
                    float,
                    typename std::conditional<
                        datatype == Datatype::Int,
                        int,
                        void
                    >::type
                >::type
            >::type
        >::type
    >::type;
} ;

template<Approach approach_, Datatype datatype_>
struct ApproachAndType{
    static constexpr Approach approach = approach_;
    static constexpr Datatype datatype = datatype_;

    using ScoreType = typename ScoreTypeFromDatatypeE<datatype_>::type;
};

__forceinline__
std::string to_string(Approach approach){
    switch(approach){
        case Approach::hardcodedzero: return "hardcodedzero";
        case Approach::kernelparamzero: return "kernelparamzero";
        case Approach::unused: return "unused";
        case Approach::half2fmarelu: return "half2fmarelu";
        case Approach::floatloadhalf: return "floatloadhalf";
        
    }
    return "to_string: missing case for Approach";
}

__forceinline__
std::string to_string(Datatype datatype){
    switch(datatype){
        case Datatype::Half2: return "Half2";
        case Datatype::Short2: return "Short2";
        case Datatype::UInt8x4: return "UInt8x4";
        case Datatype::Float: return "Float";
        case Datatype::Int: return "Int";
        case Datatype::Invalid: return "Invalid";
    }
    return "to_string: missing case for Datatype";
}





template<class OutputT, class InputT>
struct ConvertType{
    __host__ __device__
    OutputT operator()(InputT input) const {
        return static_cast<OutputT>(input);
    }
};

template<>
struct ConvertType<__nv_fp8_e4m3, int>{
    __host__ __device__
    __nv_fp8_e4m3 operator()(int input) const {
        return __nv_fp8_e4m3(half(input));
    }
};

template<class OutputT, class InputT>
__host__ __device__
OutputT convert(InputT input){
    return ConvertType<OutputT,InputT>{}(input);
}


//stream-ordered uninitialized allocation
inline
thrust::device_vector<char, thrust_async_allocator<char>> allocateTempBuffer(size_t bytes, size_t minimumBytesForRetry, cudaStream_t stream){

    thrust::device_vector<char, thrust_async_allocator<char>> result{thrust_async_allocator<char>(stream)};

    
        // // cudaDeviceSynchronize();
        // int deviceId;
        // cudaGetDevice(&deviceId);
        // cudaMemPool_t mempool;
        // cudaDeviceGetDefaultMemPool(&mempool, deviceId);

        // size_t reservedpoolmem = 0;
        // cudaMemPoolGetAttribute(mempool, cudaMemPoolAttrReservedMemCurrent, &reservedpoolmem);
        // size_t usedpoolmem = 0;
        // cudaMemPoolGetAttribute(mempool, cudaMemPoolAttrUsedMemCurrent, &usedpoolmem);

        // // std::cerr << "reservedpoolmem " << reservedpoolmem << ", usedpoolmem " << usedpoolmem << "\n";

        // size_t free,total;
        // cudaMemGetInfo(&free, &total);
        
        // if(free > (64 << 20)){
        //     free -= (64 << 20);
        // }
        // std::cerr << "bytes " << bytes << ", free " << free << ", total " << total << ", updatesbytes " << std::min(free, bytes) << "\n";
        // bytes = std::min(free, bytes);
    

    do{
        // std::cerr << "try " << bytes << "\n";
        try{
            result.resize(bytes, thrust::no_init);
            // cudaError_t status = cudaPeekAtLastError();
            
            // cudaMemPoolGetAttribute(mempool, cudaMemPoolAttrReservedMemCurrent, &reservedpoolmem);
            // cudaMemPoolGetAttribute(mempool, cudaMemPoolAttrUsedMemCurrent, &usedpoolmem);
            // cudaMemGetInfo(&free, &total);

            // std::cerr << cudaGetErrorString(status) << ", reservedpoolmem " << reservedpoolmem << ", usedpoolmem " << usedpoolmem
            //     << ", free " << free << ", total " << total << "\n";
            // std::cerr << cudaGetErrorString(status) << "\n";
            return result;
        }catch(...){
            bytes /= 1.5;
            cudaGetLastError();
        }
    }
    while(bytes >= minimumBytesForRetry);

    throw std::bad_alloc();
}

//stream-ordered uninitialized allocation
inline
thrust::device_vector<char, thrust_async_allocator<char>> allocateTempBuffer(size_t bytes, cudaStream_t stream){
    thrust::device_vector<char, thrust_async_allocator<char>> result{thrust_async_allocator<char>(stream)};
    result.resize(bytes, thrust::no_init);
    return result;
}


//Buffers that contain sequence data will be padded to simplify vectorized memory accesses
inline
size_t getSequenceCharsAllocationSize(size_t bytes){
    return cuda::round_up(bytes, 32);
}


LIBMARV_NAMESPACE_WITH_NESTING_END

#endif