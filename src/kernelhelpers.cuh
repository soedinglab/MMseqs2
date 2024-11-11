#include <type_traits>

namespace cudasw4{

    template<class SequenceLengthT>
    __host__ __device__
    SequenceLengthT getPaddedQueryLength(SequenceLengthT queryLength){
        //pad query length to char4, add warpsize char4 border.
        return ((queryLength + 4 - 1) / 4) * 4 + sizeof(char4) * 32;
    }

    typedef union
    {
        int32_t   i;
        uint32_t  u;
        short2   s2;
    } data_pack;

    typedef short2	  score_type;
    typedef data_pack score_pack;

    template<class T, std::enable_if_t<(sizeof(T) <= 4), bool> = true>
    __device__
    T warp_max_reduce_broadcast(unsigned int mask, T val){
        #if __CUDA_ARCH__ >= 800
            return __reduce_max_sync(mask, val);
        #else
            for (int offset = 16; offset > 0; offset /= 2){
                T tmp = __shfl_down_sync(mask, val, offset);
                val = tmp > val ? tmp : val;
            }
            return __shfl_sync(mask, val, 0);
        #endif
    }

    template<class T, std::enable_if_t<(sizeof(T) > 4), bool> = true>
    __device__
    T warp_max_reduce_broadcast(unsigned int mask, T val){
        for (int offset = 16; offset > 0; offset /= 2){
            T tmp = __shfl_down_sync(mask, val, offset);
            val = tmp > val ? tmp : val;
        }
        return __shfl_sync(mask, val, 0);
    }

    inline __device__ short2 viaddmax(const short2 a_in, const short2 b_in, const short2 c_in) {
        score_pack a, b, c, d;
        a.s2 = a_in;
        b.s2 = b_in;
        c.s2 = c_in;
        d.u = __viaddmax_s16x2(a.u, b.u, c.u);
        return(d.s2);
    }

    inline __device__ short2 viadd(const short2 a_in, const short2 b_in) {
        score_pack a, b, d;
        a.s2 = a_in;
        b.s2 = b_in;
        d.u = __vadd2(a.u, b.u);
        return(d.s2);
    }

    inline __device__ short2 vimax(const short2 a_in, const short2 b_in) {
        score_pack a, b, d;
        a.s2 = a_in;
        b.s2 = b_in;
        d.u = __vmaxs2(a.u, b.u);
        return(d.s2);
    }

    inline __device__ short2 vimax3(const short2 a_in, const short2 b_in, const short2 c_in) {
        score_pack a, b, c, d;
        a.s2 = a_in;
        b.s2 = b_in;
        c.s2 = c_in;
        d.u = __vimax3_s16x2_relu(a.u, b.u, c.u);
        return(d.s2);
    }

    inline __device__ short2 shfl_up_2xint16(const uint32_t bitmap, const short2 value, const int lane, const int group_size) {
        score_pack v, res;
        v.s2 = value;
        res.u =__shfl_up_sync(bitmap, v.u, lane, group_size);
        return(res.s2);
    }

    inline __device__ short2 shfl_down_2xint16(const uint32_t bitmap, const short2 value, const int lane, const int group_size) {
        score_pack v, res;
        v.s2 = value;
        res.u =__shfl_down_sync(bitmap, v.u, lane, group_size);
        return(res.s2);
    }

} //namespace cudasw4