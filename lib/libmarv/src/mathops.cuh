#ifndef MATH_OPS_CUH
#define MATH_OPS_CUH

#include <cuda_fp16.h>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>

// from https://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__INTRINSIC__HALF__CONSTANTS.html
#ifndef CUDART_ZERO_FP16
#define CUDART_ZERO_FP16 __ushort_as_half((unsigned short)0x0000U)
#endif

namespace cudasw4{

    template<class T>
    struct MathOps{};

    template<>
    struct MathOps<half2>{
        using Type = half;
        using VecType = half2;

        __host__ __device__
        static VecType zero_score(){
            return 	(make_half2(0,0));
        }

        __device__
        static VecType add(const VecType& a, const VecType& b){
            return __hadd2(a,b);
        }

        //max(a,b)
        __device__
        static Type max(const Type& a, const Type& b){
            return __hmax(a, b);
        }

        //max(a,b)
        __device__
        static VecType max(const VecType& a, const VecType& b){
            return __hmax2(a, b);
        }

        //max(a,b)
        __device__
        static VecType max(const VecType& a, const VecType& b, bool* a_xIsMax, bool* a_yIsMax){
            #if 0
                //uses prmt to extract x and y
                VecType result;
                if(a.x >= b.x){
                    result.x = a.x;
                    *a_xIsMax = true;
                }else{
                    result.x = b.x;
                    *a_xIsMax = false;
                }
                if(a.y >= b.y){
                    result.y = a.y;
                    *a_yIsMax = true;
                }else{
                    result.y = b.y;
                    *a_yIsMax = false;
                }
                return result;
            #else
                VecType result = max(a,b);
                if(a.x >= b.x){
                    *a_xIsMax = true;
                }else{
                    *a_xIsMax = false;
                }
                if(a.y >= b.y){
                    *a_yIsMax = true;
                }else{
                    *a_yIsMax = false;
                }
                return result;
            #endif
        }


        //max(a,b,c)
        __device__
        static VecType max3(const VecType& a, const VecType& b, const VecType& c){
            return __hmax2(a, __hmax2(b,c));
        }

        //max(a+b,c)
        __device__
        static VecType add_max(const VecType& a, const VecType& b, const VecType& c){
            return __hmax2(__hadd2(a,b), c);
        }

        //max(a+b,c)
        __device__
        static VecType add_max(const VecType& a, const VecType& b, const VecType& c, bool* sum_xIsMax, bool* sum_yIsMax){
            VecType sum = add(a,b);
            return max(sum, c, sum_xIsMax, sum_yIsMax);
        }

        //max(a+b,0)
        __device__
        static VecType add_relu(const VecType& a, const VecType& b){
            return add_max(a,b, make_half2(0,0));
        }

        //max(a+b,0)
        __device__
        static VecType add_relu(const VecType& a, const VecType& b, const VecType& zero){
            return add_max(a,b, zero);
        }

        //max(max(a + b, c), 0)
        __device__
        static VecType add_max_relu(const VecType& a, const VecType& b, const VecType& c){
            return max(add_max(a,b,c), make_half2(0,0));
        }

        template<class Group>
        __device__
        static VecType reduce_max(Group& group, VecType val){
            return cooperative_groups::reduce(group, val, [](const auto& l, const auto& r){return __hmax2(l,r);});
            //return cooperative_groups::reduce(group, val, cooperative_groups::greater<VecType>{});
        }
    };

    template<>
    struct MathOps<short2>{
        using Type = short;
        using VecType = short2;

        __host__ __device__
        static VecType zero_score(){
            return 	(make_short2(0,0));
        }

        __device__
        static unsigned int asUint(const short2& s){
            unsigned int u;
            memcpy(&u, &s, sizeof(unsigned int));
            return u;
        }

        __device__
        static VecType asVec(unsigned int u){
            VecType v;
            memcpy(&v, &u, sizeof(unsigned int));
            return v;
        }

        __device__
        static VecType add(const VecType& a, const VecType& b){
            return asVec(__vadd2(asUint(a), asUint(b)));
        }

        //max(a,b)
        __device__
        static Type max(const Type& a, const Type& b){
            return ::max(a, b);
        }

        //max(a,b)
        __device__
        static VecType max(const VecType& a, const VecType& b){
            return 	asVec(__vmaxs2(asUint(a), asUint(b)));
        }

        //max(a,b)
        __device__
        static VecType max(const VecType& a, const VecType& b, bool* a_xIsMax, bool* a_yIsMax){
            return  asVec(__vibmax_u16x2(asUint(a), asUint(b), a_xIsMax, a_yIsMax));
        }

        //max(a,b,c)
        __device__
        static VecType max3(const VecType& a, const VecType& b, const VecType& c){
            return 	asVec(__vimax3_s16x2(asUint(a), asUint(b), asUint(c)));
        }

        //max(a+b,c)
        __device__
        static VecType add_max(const VecType& a, const VecType& b, const VecType& c){
            return 	asVec(__viaddmax_s16x2(asUint(a), asUint(b), asUint(c)));
        }

        //max(a+b,c)
        __device__
        static VecType add_max(const VecType& a, const VecType& b, const VecType& c, bool* sum_xIsMax, bool* sum_yIsMax){
            VecType sum = add(a,b);
            return max(sum, c, sum_xIsMax, sum_yIsMax);
        }

        //max(a+b,0)
        __device__
        static VecType add_relu(const VecType& a, const VecType& b){
            return add_max(a,b, make_short2(0,0));
        }

        //max(a+b,0)
        __device__
        static VecType add_relu(const VecType& a, const VecType& b, const VecType& zero){
            return asVec(__viaddmax_s16x2(asUint(a), asUint(b), asUint(zero)));
        }

        //max(max(a + b, c), 0)
        __device__
        static VecType add_max_relu(const VecType& a, const VecType& b, const VecType& c){
            return asVec(__viaddmax_s16x2_relu(asUint(a), asUint(b), asUint(c)));
        }

        template<class Group>
        __device__
        static VecType reduce_max(Group& group, VecType val){
            return asVec(cooperative_groups::reduce(group, asUint(val), [](const auto& l, const auto& r){return __vmaxs2(l,r);}));
            //return cooperative_groups::reduce(group, val, cooperative_groups::greater<VecType>{});
        }
    };


    template<>
    struct MathOps<float>{
        using Type = float;

        __device__
        static Type add(const Type& a, const Type& b){
            return a + b;
        }

        //max(a,b)
        __device__
        static Type max(const Type& a, const Type& b){
            return ::max(a,b);
        }

        //max(a,b)
        __device__
        static Type max(const Type& a, const Type& b, bool* firstIsMax){
            if(a >= b){
                *firstIsMax = true;
                return a;
            }else{
                *firstIsMax = false;
                return b;
            }
        }

        //max(a,b,c)
        __device__
        static Type max3(const Type& a, const Type& b, const Type& c){
            return max(a,max(b,c));
        }

        //max(a+b,c)
        __device__
        static Type add_max(const Type& a, const Type& b, const Type& c){
            return max(a+b,c);
        }

        //max(a+b,c)
        __device__
        static Type add_max(const Type& a, const Type& b, const Type& c, bool* firstIsMax){
            if(a+b >= c){
                *firstIsMax = true;
                return a+b;
            }else{
                *firstIsMax = false;
                return c;
            }
        }

        //max(a+b,0)
        __device__
        static Type add_relu(const Type& a, const Type& b){
            return add_max(a,b,0.0f);
        }

        //max(max(a + b, c), 0)
        __device__
        static Type add_max_relu(const Type& a, const Type& b, const Type& c){
            return max(add_max(a,b,c), 0.0f);
        }

        template<class Group>
        __device__
        static Type reduce_max(Group& group, Type val){
            return cooperative_groups::reduce(group, val, cooperative_groups::greater<Type>{});
        }
    };


    template<>
    struct MathOps<int>{
        using Type = int;

        __device__
        static Type add(const Type& a, const Type& b){
            return a + b;
        }

        //max(a,b)
        __device__
        static Type max(const Type& a, const Type& b){
            return ::max(a,b);
        }

        //max(a,b)
        __device__
        static Type max(const Type& a, const Type& b, bool* firstIsMax){
            return __vibmax_s32(a, b, firstIsMax);
        } 

        //max(a,b,c)
        __device__
        static Type max3(const Type& a, const Type& b, const Type& c){
            return __vimax3_s32(a,b,c);
        }

        //max(a+b,c)
        __device__
        static Type add_max(const Type& a, const Type& b, const Type& c){
            return __viaddmax_s32(a,b,c);
        }

        //max(a+b,c)
        __device__
        static Type add_max(const Type& a, const Type& b, const Type& c, bool* firstIsMax){
            return __vibmax_s32(a+b, c, firstIsMax);
        }

        //max(a+b,0)
        __device__
        static Type add_relu(const Type& a, const Type& b){
            return add_max(a,b, 0);
        }

        //max(max(a + b, c), 0)
        __device__
        static Type add_max_relu(const Type& a, const Type& b, const Type& c){
            return __viaddmax_s32_relu(a,b,c);
        }    

        template<class Group>
        __device__
        static Type reduce_max(Group& group, Type val){
            return cooperative_groups::reduce(group, val, cooperative_groups::greater<Type>{});
        }
    };

    template<>
    struct MathOps<half>{
        using Type = half;

        //max(a,b)
        __device__
        static Type max(const Type& a, const Type& b){
            return __hmax(a,b);
        }

        //max(a,b,c)
        __device__
        static Type max3(const Type& a, const Type& b, const Type& c){
            return max(a,max(b,c));
        }

        //max(a+b,c)
        __device__
        static Type add_max(const Type& a, const Type& b, const Type& c){
            return max(a+b,c);
        }

        //max(a+b,0)
        __device__
        static Type add_relu(const Type& a, const Type& b){
            return add_max(a,b, CUDART_ZERO_FP16);
        }

        //max(max(a + b, c), 0)
        __device__
        static Type add_max_relu(const Type& a, const Type& b, const Type& c){
            return max(add_max(a,b,c), CUDART_ZERO_FP16);
        }

        template<class Group>
        __device__
        static Type reduce_max(Group& group, Type val){
            return cooperative_groups::reduce(group, val, cooperative_groups::greater<Type>{});
        }
    };

    template<>
    struct MathOps<short>{
        using Type = short;

        //max(a,b)
        __device__
        static Type max(const Type& a, const Type& b){
            return ::max(a,b);
        }

        //max(a,b,c)
        __device__
        static Type max3(const Type& a, const Type& b, const Type& c){
            return max(a,max(b,c));
        }

        //max(a+b,c)
        __device__
        static Type add_max(const Type& a, const Type& b, const Type& c){
            return max(a+b,c);
        }

        //max(a+b,0)
        __device__
        static Type add_relu(const Type& a, const Type& b){
            return add_max(a,b,0);
        }

        //max(max(a + b, c), 0)
        __device__
        static Type add_max_relu(const Type& a, const Type& b, const Type& c){
            return max(add_max(a,b,c), 0);
        }

        template<class Group>
        __device__
        static Type reduce_max(Group& group, Type val){
            return cooperative_groups::reduce(group, val, cooperative_groups::greater<Type>{});
        }
    };


} //namespace cudasw4

#endif