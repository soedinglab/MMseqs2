#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <cstdint>
#include <type_traits>

#include "namespace.hpp"
LIBMARV_NAMESPACE_BEGIN

//MODIFY AT OWN RISK

//data type to enumerate all sequences in the database
using ReferenceIdT = std::int32_t;

//data type for length of of both query sequences and databases sequences
using SequenceLengthT = std::int32_t;

static_assert(std::is_same_v<ReferenceIdT, std::int32_t>, "unexpected reference type");
static_assert(std::is_same_v<SequenceLengthT, std::int32_t>, "unexpected sequence length type");

struct MaxSequencesInDB{
    static constexpr ReferenceIdT value(){
        return std::numeric_limits<ReferenceIdT>::max() - 1;
    }
};

struct MaxSequenceLength{
    static constexpr SequenceLengthT value(){
        return std::numeric_limits<SequenceLengthT>::max() - 128 - 4;
    }
};

struct MaxNumberOfResults{
    static constexpr int value(){
        return 512*1024;
    }
};

struct alignas(8) AlignmentEndPosition{
    int queryEndExclusive;
    int subjectEndExclusive;

    AlignmentEndPosition() = default;
    AlignmentEndPosition(int subjectEndExcl, int queryEndExcl)
        : queryEndExclusive(queryEndExcl), subjectEndExclusive(subjectEndExcl){}

    #ifdef __CUDACC__
    __host__ __device__
    #endif
    int getQueryEndInclusive() const{
        return queryEndExclusive-1;
    }

    #ifdef __CUDACC__
    __host__ __device__
    #endif
    int getSubjectEndInclusive() const{
        return subjectEndExclusive-1;
    }

    #ifdef __CUDACC__
    __host__ __device__
    #endif
    int getQueryEndExclusive() const{
        return queryEndExclusive;
    }

    #ifdef __CUDACC__
    __host__ __device__
    #endif
    int getSubjectEndExclusive() const{
        return subjectEndExclusive;
    }
};


LIBMARV_NAMESPACE_END


#endif