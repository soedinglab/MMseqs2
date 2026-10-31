#ifndef CUDA_HIP_COMPATIBILITY_
#define CUDA_HIP_COMPATIBILITY_

#if defined(__HIPCC__)
    #include <hip/hip_fp16.h>
    #include <hipcub/hipcub.hpp>
#endif

#if defined(__CUDACC__) || defined(__HIP_PLATFORM_NVIDIA__)

#include <cooperative_groups/reduce.h>
// #include <cooperative_groups/memcpy_async.h>

#endif

namespace compatibility{

//HALF2 MATH

#if defined(__CUDACC__) || defined(__HIP_PLATFORM_NVIDIA__)

__device__ __forceinline__
half2 hmax2(half2 a, half2 b){
    return __hmax2(a,b);
}

__device__ __forceinline__
half2 hmin2(half2 a, half2 b){
    return __hmin2(a,b);
}

#elif defined(__HIPCC__)

__device__
half2 hmax2(half2 a, half2 b){
    half2 res;

    asm volatile ("v_pk_max_f16 %0, %1, %2\n\t" : "=v" (res) : "v" (a) , "v" (b));

    return res;
}

__device__
half2 hmin2(half2 a, half2 b){
    half2 res;

    asm volatile ("v_pk_min_f16 %0, %1, %2\n\t" : "=v" (res) : "v" (a) , "v" (b));

    return res;
}

#endif









//GROUP SHFL DOWN

template<class Group, class T>
__device__
T group_shfl_down(const Group& group, T item, int offset){
    constexpr int numInts = (sizeof(T) + 3) / 4;
    unsigned int tmp[numInts];
    memcpy(&tmp, &item, sizeof(T));

    #pragma unroll
    for(int i = 0; i < numInts; i++){
        tmp[i] = group.shfl_down(tmp[i], offset);
    }

    memcpy(&item, &tmp, sizeof(T));
    return item;

}

//GROUP SHFL XOR

template<class Group, class T>
__device__
T group_shfl_xor(const Group& group, T item, int offset){
    constexpr int numInts = (sizeof(T) + 3) / 4;
    unsigned int tmp[numInts];
    memcpy(&tmp, &item, sizeof(T));

    #pragma unroll
    for(int i = 0; i < numInts; i++){
        tmp[i] = group.shfl_xor(tmp[i], offset);
    }

    memcpy(&item, &tmp, sizeof(T));
    return item;

}



// GROUP REDUCE

#if defined(__CUDACC__) || defined(__HIP_PLATFORM_NVIDIA__)

template<class Group, class T, class Func>
__device__
T group_reduce(const Group& group, T val, Func func){
    return cooperative_groups::reduce(group, val, func);
}

template<class Group>
__device__
half2 group_reduce_max(const Group& group, half2 val){
    return cooperative_groups::reduce(group, val, [](const auto& l, const auto& r){return __hmax2(l,r);});
}

template<class Group>
__device__
float group_reduce_max(const Group& group, float val){
    return cooperative_groups::reduce(group, val, cooperative_groups::greater<float>{});
}

template<class Group>
__device__
short2 group_reduce_max(const Group& group, short2 val){
    static_assert(sizeof(int) == sizeof(short2));
    unsigned int tmp; memcpy(&tmp, &val, sizeof(short2));
    tmp = cooperative_groups::reduce(group, tmp, [](const auto& l, const auto& r){return __vmaxs2(l,r);});
    memcpy(&val, &tmp, sizeof(short2));
    return val;
}

template<class Group>
__device__
int group_reduce_max(const Group& group, int val){
    return cooperative_groups::reduce(group, val, cooperative_groups::greater<int>{});
}

namespace cooperative_groups{

    using namespace ::cooperative_groups;

    template<int /*groupsize*/, int /*blocksize*/, class Group, class T>
    T exclusive_scan(const Group& group, T val){
        return cooperative_groups::exclusive_scan(group, val);
    }

}


#elif defined(__HIP_PLATFORM_AMD__)


template<class Group, class T, class Func>
__device__
T group_reduce(const Group& group, T val, Func func){
    for(int stride = 1; stride < group.size(); stride *= 2){
        val = func(val, group_shfl_xor(group, val, stride));
    }
    return val;
}

template<class Group>
__device__
half2 group_reduce_max(const Group& group, half2 val){
    return group_reduce(group, val, [](const auto& l, const auto& r){return hmax2(l,r);});
}

template<class Group>
__device__
float group_reduce_max(const Group& group, float val){
    return group_reduce(group, val, [](const auto& l, const auto& r){return max(l,r);});
}

template<class Group>
__device__
short2 group_reduce_max(const Group& group, short2 val){
    static_assert(sizeof(int) == sizeof(short2));
    unsigned int tmp; memcpy(&tmp, &val, sizeof(short2));
    tmp = group_reduce(group, tmp, [](const auto& l, const auto& r){return __vmaxs2(l,r);});
    memcpy(&val, &tmp, sizeof(short2));
    return val;
}

template<class Group>
__device__
int group_reduce_max(const Group& group, int val){
    return group_reduce(group, val, [](const auto& l, const auto& r){return max(l,r);});
}

namespace cooperative_groups{

template<class T>
struct plus{
    __device__
    T operator()(const T& l, const T& r) const{
        return l + r;
    }
};

template<class T>
struct bit_or{
    __device__
    T operator()(const T& l, const T& r) const{
        return l | r;
    }
};

template<class Group, class T, class Func>
__device__
T reduce(const Group& group, T val, Func func){
    for(int stride = 1; stride < group.size(); stride *= 2){
        val = func(val, group_shfl_xor(group, val, stride));
    }
    return val;
}

template<int groupsize, int blocksize, class Group, class T>
__device__
T exclusive_scan(const Group& group, T val){
    static_assert(blocksize % groupsize == 0);
    constexpr int numGroupsInBlock = blocksize / groupsize;

    using warp_scan = typename rocprim::warp_scan<T, 1>;

    __shared__ typename warp_scan::storage_type temp[numGroupsInBlock];

    const int logical_warp_id = threadIdx.x/groupsize;
    warp_scan().exclusive_scan(
        val, // input
        val, // output
        temp[logical_warp_id]
    );

    return val;
}

} //namespace cooperative_groups


#endif













} //namespace compatibility



#endif
