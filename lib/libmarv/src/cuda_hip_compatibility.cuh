#ifndef CUDA_HIP_COMPATIBILITY_
#define CUDA_HIP_COMPATIBILITY_

#if defined(__HIPCC__)
    #include <hip/hip_fp16.h>
#endif

#if defined(__CUDACC__) || defined(__HIP_PLATFORM_NVIDIA__)

#include <cooperative_groups/reduce.h>

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
T group_shfl_down(Group group, T item, int offset){
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
T group_shfl_xor(Group group, T item, int offset){
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
T group_reduce(Group group, T val, Func func){
    return cooperative_groups::reduce(group, val, func);
}

template<class Group>
__device__
half2 group_reduce_max(Group group, half2 val){
    return cooperative_groups::reduce(group, val, [](const auto& l, const auto& r){return __hmax2(l,r);});
}

template<class Group>
__device__
float group_reduce_max(Group group, float val){
    return cooperative_groups::reduce(group, val, cooperative_groups::greater<float>{});
}

template<class Group>
__device__
short2 group_reduce_max(Group group, short2 val){
    static_assert(sizeof(int) == sizeof(short2));
    unsigned int tmp; memcpy(&tmp, &val, sizeof(short2));
    tmp = cooperative_groups::reduce(group, tmp, [](const auto& l, const auto& r){return __vmaxs2(l,r);});
    memcpy(&val, &tmp, sizeof(short2));
    return val;
}

template<class Group>
__device__
int group_reduce_max(Group group, int val){
    return cooperative_groups::reduce(group, val, cooperative_groups::greater<int>{});
}

#elif defined(__HIP_PLATFORM_AMD__)


template<class Group, class T, class Func>
__device__
T group_reduce(Group group, T val, Func func){
    for(int stride = 1; stride < group.size(); stride *= 2){
        val = func(val, group_shfl_xor(group, val, stride));
    }
    return val;
}

template<class Group>
__device__
half2 group_reduce_max(Group group, half2 val){
    return group_reduce(group, val, [](const auto& l, const auto& r){return hmax2(l,r);});
}

template<class Group>
__device__
float group_reduce_max(Group group, float val){
    return group_reduce(group, val, [](const auto& l, const auto& r){return max(l,r);});
}

template<class Group>
__device__
short2 group_reduce_max(Group group, short2 val){
    static_assert(sizeof(int) == sizeof(short2));
    unsigned int tmp; memcpy(&tmp, &val, sizeof(short2));
    tmp = group_reduce(group, tmp, [](const auto& l, const auto& r){return __vmaxs2(l,r);});
    memcpy(&val, &tmp, sizeof(short2));
    return val;
}

template<class Group>
__device__
int group_reduce_max(Group group, int val){
    return group_reduce(group, val, [](const auto& l, const auto& r){return max(l,r);});
}
#endif













} //namespace compatibility



#endif
