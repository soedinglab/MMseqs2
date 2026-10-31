#ifndef LIBMARV_ALIGNMENT_INTERFACE_DATA_CUH
#define LIBMARV_ALIGNMENT_INTERFACE_DATA_CUH

#include "../offset_iterator.cuh"
#include "../util.cuh"
#include <cstdint>
#include <array>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN


struct OneToAllInputDataPSSM{
    public:
    
        OneToAllInputDataPSSM() = default;
        
        __host__ __device__
        int getNumAlignments() const{
            return numAlignments;
        }
    
        __host__ __device__
        const std::int8_t* getSubject(ReferenceIdT i) const{
            const auto index = indexIndirection[i];
            return d_subjects + d_subjectOffsets[index];
        }
    
        __host__ __device__
        SequenceLengthT getSubjectLength(ReferenceIdT i) const{
            const auto index = indexIndirection[i];
            return d_subjectLengths[index];
        }
    
        __host__ __device__
        SequenceLengthT getQueryLength() const{
            return queryLength;
        }
    
        __host__ __device__
        int getMaximumSubjectLength() const{
            return maximumSubjectLength;
        }   
    
    public: // API does not care about those
        const std::int8_t* d_subjects{};
        ReadOnlyCustomOffsetIterator<size_t*> d_subjectOffsets{};
        const SequenceLengthT* d_subjectLengths{};
        SequenceLengthT queryLength{};
        int numAlignments{};
        int maximumSubjectLength{};
        PositionsIterator indexIndirection = PositionsIterator::fromCountingIterator(0);
        const unsigned char* letterTranslation{}; //optional, maps stored bytes to alphabet codes (Metal and SIMD backends)
    };

    enum class AlignmentAlgorithmE{
        Gapless,
        SW
    };

    struct GapScoreArgs{
        int gapopenscore;
        int gapextendscore;
    };

    template<
        int alphabetSize_,
        Approach approach_,
        Datatype datatype_,
        int blocksize_,
        int groupsize_,
        int numRegs_
    >
    struct GaplessConfig{
        static constexpr AlignmentAlgorithmE algorithm = AlignmentAlgorithmE::Gapless;

        static constexpr int alphabetSize = alphabetSize_;
        static constexpr Approach approach = approach_;
        static constexpr Datatype datatype = datatype_;
        static constexpr int blocksize = blocksize_;
        static constexpr int groupsize = groupsize_;
        static constexpr int numRegs = numRegs_;

        using ScoreType = typename ScoreTypeFromDatatypeE<datatype>::type;
        static constexpr int tileSize = groupsize * numRegs * (
            IsScalarType<ScoreType>::value ? 1 : 
                (IsVec2Type<ScoreType>::value ? 2 : 
                    (IsVec4Type<ScoreType>::value ? 4 : 0)));
        
    };

    template<
        int alphabetSize_,
        Approach approach_,
        Datatype datatype_,
        int blocksize_,
        int groupsize_,
        int numRegs_
    >
    struct SWConfig{
        static constexpr AlignmentAlgorithmE algorithm = AlignmentAlgorithmE::SW;

        static constexpr int alphabetSize = alphabetSize_;
        static constexpr Approach approach = approach_;
        static constexpr Datatype datatype = datatype_;
        static constexpr int blocksize = blocksize_;
        static constexpr int groupsize = groupsize_;
        static constexpr int numRegs = numRegs_;

        static constexpr int tileSize = groupsize * numRegs;
    };

    template<class ConfigTuple>
    struct ArrayOfTileSizes;

    template<class... Configs>
    struct ArrayOfTileSizes<std::tuple<Configs...>>{
        static constexpr std::array<int, sizeof...(Configs)> array{Configs::tileSize ...};
    };

    template <typename T, std::size_t N, std::size_t... I>
    constexpr bool array_is_sorted_impl(const std::array<T, N>& a,
                                  std::index_sequence<I...>)
    {
        // pairwise: a[I] < a[I+1]
        return ((I + 1 < N ? a[I] < a[I + 1] : true) && ...);
    }
    
    template <typename T, std::size_t N>
    constexpr bool array_is_sorted(const std::array<T, N>& a)
    {
        return array_is_sorted_impl(a, std::make_index_sequence<N>{});
    }

    template<class ConfigTuple>
    constexpr bool configsAreSortedByTileSize(){
        constexpr auto tileSizeArray = ArrayOfTileSizes<ConfigTuple>::array;

        return array_is_sorted(tileSizeArray);
    }

LIBMARV_NAMESPACE_END


#endif //LIBMARV_ALIGNMENT_INTERFACE_DATA_CUH