#ifndef LIBMARV_NAMESPACE_HPP
#define LIBMARV_NAMESPACE_HPP

#if defined(__CUDACC__)
#include <cuda/std/__cccl/preprocessor.h>
#endif


#define LIBMARV_NAMESPACE_BEGIN namespace libmarv{ 

#define LIBMARV_NAMESPACE_END }



#if defined(__CUDACC__)

    #define LIBMARV_INLINE_NAMESPACE_BEGIN \
        inline namespace _CCCL_PP_SPLICE_WITH(_, SM, __CUDA_ARCH_LIST__){

    #define LIBMARV_INLINE_NAMESPACE_END }

#else


    #ifdef AMD_INLINE_ARCH_NAMESPACE_NAME
        #define LIBMARV_INLINE_NAMESPACE_BEGIN \
            inline namespace AMD_INLINE_ARCH_NAMESPACE_NAME{

    #else
    #error AMD_INLINE_ARCH_NAMESPACE_NAME not defined. This could lead to unexpected errors at runtime when compiling for multiple AMD architectures!
    #endif

    
    #define LIBMARV_INLINE_NAMESPACE_END }

#endif


#define LIBMARV_NAMESPACE_WITH_NESTING_BEGIN \
    LIBMARV_NAMESPACE_BEGIN \
    LIBMARV_INLINE_NAMESPACE_BEGIN

#define LIBMARV_NAMESPACE_WITH_NESTING_END \
    LIBMARV_INLINE_NAMESPACE_END \
    LIBMARV_NAMESPACE_END
    


#endif