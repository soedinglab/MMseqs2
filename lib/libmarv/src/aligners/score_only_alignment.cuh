#ifndef LIBMARV_SCORE_ONLY_ALIGNMENT_CUH
#define LIBMARV_SCORE_ONLY_ALIGNMENT_CUH

#include "../cuda_hip_rename.h"

#if defined(__HIPCC__)
    #include "hip/hip_runtime.h"
#endif

#include <cuda/cmath>


#include "../pssm.cuh"
#include "../util.cuh"
#include "../cuda_errorcheck.cuh"
#include "score_only_alignment_interfaces.cuh"
#include "../alignment_algorithms/alignment_interface_data.cuh"
#include "../offset_iterator.cuh"

#include "../alignment_algorithms/gapless/pssmkernels_gapless.cuh"
#include "../alignment_algorithms/gapless/pssmkernels_gapless_float.cuh"
#include "../alignment_algorithms/gapless/pssmkernels_gapless_int8.cuh"
#include "../alignment_algorithms/smithwaterman/pssmkernels_smithwaterman.cuh"

#include <map>
#include <cuda/std/mdspan>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_WITH_NESTING_BEGIN



    template<
        int alphabetSize,
        class ApproachAndType,
        int blocksize,
        int groupsize,
        int numRegs,
        bool enableSingleTile,
        bool enableMultiTile
    >
    struct GaplessAlignment_scalarscore : public ScoreOnlyAlignmentInterface{
        public:
        int deviceId;
        int maxSharedMemoryPerBlockOptin;
        int numSMs;

        using ScoreType = typename ApproachAndType::ScoreType;
        static_assert(IsScalarType<ScoreType>::value);

        using PssmScoreType = typename std::conditional<
            std::is_same_v<ScoreType, float>,
            typename std::conditional<
                ApproachAndType::approach == Approach::floatloadhalf,
                half,
                float
            >::type,
            ScoreType
        >::type;


        static_assert(enableSingleTile || enableMultiTile);

        static constexpr int numVectorLanes = 1;

        static constexpr int tileSize = groupsize * numRegs * numVectorLanes;

        static constexpr bool subjectIsCaseSensitive = true;
        static constexpr bool withEndPosition = false;

        GaplessAlignment_scalarscore(){
            CUDACHECK(cudaGetDevice(&deviceId));
            CUDACHECK(cudaDeviceGetAttribute(&maxSharedMemoryPerBlockOptin, cudaDevAttrMaxSharedMemoryPerBlockOptin, deviceId));
            CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));
        }

        int getTileSize() const override{
            return tileSize;
        }

        AlignmentAlgorithmE getAlignmentAlgorithm() const override{
            return AlignmentAlgorithmE::Gapless;
        }

        Datatype getScoreDatatype() const override{
            return mapDatatypeToDatatypeEnum<ScoreType>();
        }

        bool isSingleTile(const OneToAllInputDataPSSM& inputData) const override{
            if constexpr(enableSingleTile){
                return inputData.getQueryLength() <= tileSize;
            }else{
                return false;
            }
        }

        void makeGpuPssm(
            GpuConvertedPSSM& result,
            const cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>> hostPssmView,
            cudaStream_t stream
        ) override{
            makeConvertedGpuPssm_gapless<ScoreType, PssmScoreType, std::int8_t, groupsize, numRegs>(
                result,
                hostPssmView,
                stream
            );
        }


        size_t getMinimumSuggestedTempBytes_multiTile(const OneToAllInputDataPSSM& inputData) const override{
            constexpr int numGroupsPerBlock = blocksize / groupsize;
            constexpr int alignmentsPerBlock = numGroupsPerBlock;
            const size_t tileTempBytesPerGroup = getTileTempBytesPerGroup(inputData);
            const int numBlocks = cuda::ceil_div(inputData.getNumAlignments(), alignmentsPerBlock);
            const size_t tempBytes1BlockPerSM = tileTempBytesPerGroup * numGroupsPerBlock * std::min(numSMs, numBlocks);
            return tempBytes1BlockPerSM;
        }

      

        void scoreOnly_singleTile(
            int* const devAlignmentScores,
            const OneToAllInputDataPSSM& inputData,
            const GpuConvertedPSSM& strided_PSSM,
            GapScoreArgs /*gapArgs*/,
            cudaStream_t stream
        ) override{
            if(strided_PSSM.kindOfPssm != GpuConvertedPSSM::Kind::Gapless){
                throw std::runtime_error("invalid pssm kind");
            }
            if constexpr(!enableSingleTile){
                throw std::runtime_error("attempted to call scoreOnly_singleTile, but singleTile was not enabled");
            }else{

                if(strided_PSSM.groupsize != groupsize) throw std::runtime_error("pssm groupsize does not match aligner groupsize");
                if(strided_PSSM.numRegs != numRegs) throw std::runtime_error("pssm numRegs does not match aligner numRegs");

                const int numAlignments = inputData.getNumAlignments();

                constexpr int groupsPerBlock = blocksize / groupsize;
                constexpr int alignmentsPerGroup = 1;
                constexpr int alignmentsPerBlock = groupsPerBlock * alignmentsPerGroup;
                // std::cout << "blocksize " << blocksize << ", groupsize " << groupsize 
                //     << ", alignmentsPerBlock " << alignmentsPerBlock << ", numAlignments " << numAlignments << "\n";

                constexpr int relaxChunkSize = 16 / sizeof(PssmScoreType);
                constexpr int pssmReplicationFactor = PssmReplicationFactor<groupsize, PssmScoreType, relaxChunkSize>::value();

                constexpr int numRowsPSSM = 21;
                constexpr int pssmTileSize = groupsize * relaxChunkSize * SDIV(numRegs, relaxChunkSize);
                constexpr int numColumnsPSSM = pssmTileSize * pssmReplicationFactor;

                
                using SharedPSSM = SharedPSSM_singletile<PssmScoreType, numRowsPSSM, numColumnsPSSM>;
                // std::cout << "host: pssmTileSize = " << pssmTileSize << ", numColumnsPSSM = " << numColumnsPSSM << "\n";
                // std::cout << "host: smem pssm takes " << sizeof(SharedPSSM) << " bytes\n";

                
                int smem = sizeof(SharedPSSM);
                auto kernel = [&](){
                    return  GaplessAlignment_PSSM_scalar_singletile_kernel<
                        ScoreType,
                        blocksize, 
                        groupsize, 
                        numRegs, 
                        withEndPosition,
                        subjectIsCaseSensitive,
                        int*,
                        PssmScoreType,
                        OneToAllInputDataPSSM>;                    
                }();

                auto setSmemKernelAttribute = [&](){
                    static std::map<int, bool> isSet;
                    if(smem > 48*1024){
                        int deviceId;
                        CUDACHECK(cudaGetDevice(&deviceId));
                        if(!isSet[deviceId]){
                            cudaError_t status = cudaFuncSetAttribute(reinterpret_cast<const void*>(kernel), cudaFuncAttributeMaxDynamicSharedMemorySize, smem);
                            if(status != cudaSuccess){
                                cudaGetLastError();
                                throw std::runtime_error("Unable to set cudaFuncAttributeMaxDynamicSharedMemorySize to " + std::to_string(smem));
                            }
                            isSet[deviceId] = true;
                        }
                    }
                };
                setSmemKernelAttribute();

                dim3 grid = (numAlignments + alignmentsPerBlock - 1) / alignmentsPerBlock;

                kernel<<<grid, blocksize, smem, stream>>>(
                    devAlignmentScores,
                    nullptr,
                    nullptr,
                    inputData,
                    getPssmView(strided_PSSM)
                ); CUERR;
                
            }
        }

        void scoreOnly_multiTile(
            char* d_temp,
            size_t tempBytes,
            int* const devAlignmentScores,
            const OneToAllInputDataPSSM& inputData,
            const GpuConvertedPSSM& strided_PSSM,
            GapScoreArgs /*gapArgs*/,
            cudaStream_t stream
        ) override{
            if(strided_PSSM.kindOfPssm != GpuConvertedPSSM::Kind::Gapless){
                throw std::runtime_error("invalid pssm kind");
            }
            if constexpr(!enableMultiTile){
                throw std::runtime_error("attempted to call scoreOnly_multiTile, but multiTile was not enabled");
            }else{
                if(strided_PSSM.groupsize != groupsize) throw std::runtime_error("pssm groupsize does not match aligner groupsize");
                if(strided_PSSM.numRegs != numRegs) throw std::runtime_error("pssm numRegs does not match aligner numRegs");

                const int numAlignments = inputData.getNumAlignments();

                constexpr int relaxChunkSize = 16 / sizeof(PssmScoreType);
                constexpr int pssmReplicationFactor = PssmReplicationFactor<groupsize, PssmScoreType, relaxChunkSize>::value();

                constexpr int numRowsPSSM = 21;
                constexpr int pssmTileSize = groupsize * relaxChunkSize * SDIV(numRegs, relaxChunkSize);
                constexpr int numColumnsPSSM = pssmTileSize * pssmReplicationFactor;
                using SharedPSSM = SharedPSSM_singletile<PssmScoreType, numRowsPSSM, numColumnsPSSM>;

                int smem = sizeof(SharedPSSM);
                auto kernel = [&](){
                    return GaplessAlignment_PSSM_scalar_multitile_kernel<ScoreType, blocksize, groupsize, numRegs, 
                        withEndPosition,
                        subjectIsCaseSensitive,
                        int*, PssmScoreType,
                        OneToAllInputDataPSSM>;
                }();

                auto setSmemKernelAttribute = [&](){
                    static std::map<int, bool> isSet;
                    if(smem > 48*1024){
                        int deviceId;
                        CUDACHECK(cudaGetDevice(&deviceId));
                        if(!isSet[deviceId]){
                            cudaError_t status = cudaFuncSetAttribute(reinterpret_cast<const void*>(kernel), cudaFuncAttributeMaxDynamicSharedMemorySize, smem);
                            if(status != cudaSuccess){
                                cudaGetLastError();
                                throw std::runtime_error("Unable to set cudaFuncAttributeMaxDynamicSharedMemorySize to " + std::to_string(smem));
                            }
                            isSet[deviceId] = true;
                        }
                    }
                };
                setSmemKernelAttribute();

                int maxBlocksPerSM = 0;
                CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
                    &maxBlocksPerSM,
                    kernel,
                    blocksize, 
                    smem
                ));

                const size_t tileTempBytesPerGroup = getTileTempBytesPerGroup(inputData);
                constexpr int groupsPerBlock = (blocksize / groupsize);
                constexpr int alignmentsPerBlock = groupsPerBlock;
                const int maxNumBlocksByInputSize = (numAlignments + alignmentsPerBlock - 1) / alignmentsPerBlock;
                const int maxNumBlocksByOccupancy = maxBlocksPerSM * numSMs;
                const int maxNumBlocksByTempBytes = tempBytes / (tileTempBytesPerGroup * groupsPerBlock);

                const int numBlocks = std::min(maxNumBlocksByTempBytes, std::min(maxNumBlocksByInputSize, maxNumBlocksByOccupancy));
                if(numBlocks <= 0){
                    throw std::runtime_error("attempted to launch multitile kernel with 0 blocks");
                }

                dim3 grid = numBlocks;


                kernel<<<grid, blocksize, smem, stream>>>(
                    devAlignmentScores,
                    nullptr,
                    nullptr,
                    inputData,
                    getPssmView(strided_PSSM),
                    (float4*)d_temp,
                    tileTempBytesPerGroup / sizeof(float4)
                ); CUERR;                
            }
        }

        private:
        size_t getTileTempBytesPerGroup(const OneToAllInputDataPSSM& inputData) const{
            static_assert(sizeof(ScoreType) == 4); //e.g. 1 float

            const int maximumSequenceLengthPadded = inputData.getMaximumSubjectLength() + groupsize;

            //1 float per position, 4 floats are stored vectorized to temp storage, round up appropriatly to ensure pointer alignment
            const size_t tileTempBytesPerGroup = cuda::round_up(4ull * maximumSequenceLengthPadded, 16ull);
            return tileTempBytesPerGroup;
        }

        PSSM_2D_View<PssmScoreType> getPssmView(const GpuPermutedPSSMforGapless& pssm) const{
            return pssm.makeScalarView<PssmScoreType>();
        }

        PSSM_2D_View<PssmScoreType> getPssmView(const GpuConvertedPSSM& pssm) const{
            return pssm.makeScalarView<PssmScoreType>();
        }
    };


    
    template<
        int alphabetSize,
        class ApproachAndType,
        int blocksize,
        int groupsize,
        int numRegs,
        bool enableSingleTile,
        bool enableMultiTile
    >
    struct GaplessAlignment_vec2score : public ScoreOnlyAlignmentInterface{
    public:
        int deviceId;
        int maxSharedMemoryPerBlockOptin;
        int numSMs;

        using ScoreType = typename ApproachAndType::ScoreType;
        static_assert(IsVec2Type<ScoreType>::value);

        using PssmScoreType = typename std::conditional<
            ApproachAndType::approach == Approach::half2fmarelu,
            __nv_fp8x2_e4m3,
            ScoreType
        >::type;
        

        static_assert(enableSingleTile || enableMultiTile);

        static constexpr int numVectorLanes = 2;

        static constexpr int tileSize = groupsize * numRegs * numVectorLanes;

        static constexpr bool subjectIsCaseSensitive = true;

        GaplessAlignment_vec2score(){
            CUDACHECK(cudaGetDevice(&deviceId));
            CUDACHECK(cudaDeviceGetAttribute(&maxSharedMemoryPerBlockOptin, cudaDevAttrMaxSharedMemoryPerBlockOptin, deviceId));
            CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));
        }

        int getTileSize() const override{
            return tileSize;
        }

        AlignmentAlgorithmE getAlignmentAlgorithm() const override{
            return AlignmentAlgorithmE::Gapless;
        }

        Datatype getScoreDatatype() const override{
            return mapDatatypeToDatatypeEnum<ScoreType>();
        }

        bool isSingleTile(const OneToAllInputDataPSSM& inputData) const override{
            if constexpr(enableSingleTile){
                return inputData.getQueryLength() <= tileSize;
            }else{
                return false;
            }
        }

        void makeGpuPssm(
            GpuConvertedPSSM& result,
            const cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>> hostPssmView,
            cudaStream_t stream
        ) override{
            makeConvertedGpuPssm_gapless<ScoreType, PssmScoreType, std::int8_t, groupsize, numRegs>(
                result,
                hostPssmView,
                stream
            );
        }


        size_t getMinimumSuggestedTempBytes_multiTile(const OneToAllInputDataPSSM& inputData) const override{
            constexpr int numGroupsPerBlock = blocksize / groupsize;
            constexpr int alignmentsPerBlock = numGroupsPerBlock;
            const size_t tileTempBytesPerGroup = getTileTempBytesPerGroup(inputData);
            const int numBlocks = cuda::ceil_div(inputData.getNumAlignments(), alignmentsPerBlock);
            const size_t tempBytes1BlockPerSM = tileTempBytesPerGroup * numGroupsPerBlock * std::min(numSMs, numBlocks);
            return tempBytes1BlockPerSM;
        }


        void scoreOnly_singleTile(
            int* const devAlignmentScores,
            const OneToAllInputDataPSSM& inputData,
            const GpuConvertedPSSM& strided_PSSM,
            GapScoreArgs /*gapArgs*/,
            cudaStream_t stream
        ) override{
            if(strided_PSSM.kindOfPssm != GpuConvertedPSSM::Kind::Gapless){
                throw std::runtime_error("invalid pssm kind");
            }
            if constexpr(!enableSingleTile){
                throw std::runtime_error("attempted to call scoreOnly_singleTile, but singleTile was not enabled");
            }else{

                if(strided_PSSM.groupsize != groupsize) throw std::runtime_error("pssm groupsize does not match aligner groupsize");
                if(strided_PSSM.numRegs != numRegs) throw std::runtime_error("pssm numRegs does not match aligner numRegs");

                const int numAlignments = inputData.getNumAlignments();

                constexpr int groupsPerBlock = blocksize / groupsize;
                constexpr int alignmentsPerGroup = 1;
                constexpr int alignmentsPerBlock = groupsPerBlock * alignmentsPerGroup;
                // std::cout << "blocksize " << blocksize << ", groupsize " << groupsize 
                //     << ", alignmentsPerBlock " << alignmentsPerBlock << ", numAlignments " << numAlignments << "\n";
       
                constexpr int relaxChunkSize = 16 / sizeof(PssmScoreType);
                constexpr int pssmReplicationFactor = PssmReplicationFactor<groupsize, PssmScoreType, relaxChunkSize>::value();
        
                constexpr int numRowsPSSM = 21;
                constexpr int pssmTileSize = groupsize * relaxChunkSize * SDIV(numRegs, relaxChunkSize);
                constexpr int numColumnsPSSM = pssmTileSize * pssmReplicationFactor;

                
                using SharedPSSM = SharedPSSM_singletile<PssmScoreType, numRowsPSSM, numColumnsPSSM>;
                // std::cout << "host: pssmTileSize = " << pssmTileSize << ", numColumnsPSSM = " << numColumnsPSSM << "\n";
                // std::cout << "host: smem pssm takes " << sizeof(SharedPSSM) << " bytes\n";

                
                int smem = sizeof(SharedPSSM);
                auto kernel = [&](){
                    if constexpr(ApproachAndType::approach == Approach::hardcodedzero){
                        return hardcodedzero::GaplessFilter_strided_PSSM_singletile_kernel<
                            ScoreType,
                            blocksize, 
                            groupsize, 
                            numRegs, 
                            subjectIsCaseSensitive,
                            int*,
                            PssmScoreType, 
                            OneToAllInputDataPSSM>;
                    }else if constexpr(ApproachAndType::approach == Approach::kernelparamzero){
                        return kernelparamzero::GaplessFilter_strided_PSSM_singletile_kernel<
                            ScoreType,
                            blocksize, 
                            groupsize, 
                            numRegs, 
                            subjectIsCaseSensitive,
                            int*,
                            PssmScoreType,
                            OneToAllInputDataPSSM>;
                    }else if constexpr(ApproachAndType::approach == Approach::half2fmarelu){
                        return half2fmarelu::GaplessFilter_strided_PSSM_singletile_kernel<
                            ScoreType,
                            blocksize, 
                            groupsize, 
                            numRegs, 
                            subjectIsCaseSensitive,
                            int*,
                            PssmScoreType,
                            OneToAllInputDataPSSM>;
                    }
                }();
        
                auto setSmemKernelAttribute = [&](){
                    static std::map<int, bool> isSet;
                    if(smem > 48*1024){
                        int deviceId;
                        CUDACHECK(cudaGetDevice(&deviceId));
                        if(!isSet[deviceId]){
                            cudaError_t status = cudaFuncSetAttribute(reinterpret_cast<const void*>(kernel), cudaFuncAttributeMaxDynamicSharedMemorySize, smem);
                            if(status != cudaSuccess){
                                cudaGetLastError();
                                throw std::runtime_error("Unable to set cudaFuncAttributeMaxDynamicSharedMemorySize to " + std::to_string(smem));
                            }
                            isSet[deviceId] = true;
                        }
                    }
                };
                setSmemKernelAttribute();
        
                dim3 grid = (numAlignments + alignmentsPerBlock - 1) / alignmentsPerBlock;
        
                if constexpr(ApproachAndType::approach == Approach::hardcodedzero){
                    kernel<<<grid, blocksize, smem, stream>>>(
                        devAlignmentScores,
                        inputData,
                        getPssmView(strided_PSSM)
                    ); CUERR;
                }else if constexpr(ApproachAndType::approach == Approach::kernelparamzero){
                    ScoreType zero = MathOps<ScoreType>::zero_score();
                    
                    kernel<<<grid, blocksize, smem, stream>>>(
                        devAlignmentScores,
                        inputData,
                        getPssmView(strided_PSSM),
                        zero
                    ); CUERR;
                }else if constexpr(ApproachAndType::approach == Approach::half2fmarelu){
                    
                    kernel<<<grid, blocksize, smem, stream>>>(
                        devAlignmentScores,
                        inputData,
                        getPssmView(strided_PSSM)
                    ); CUERR;
                }
            }
        }

        void scoreOnly_multiTile(
            char* d_temp,
            size_t tempBytes,
            int* const devAlignmentScores,
            const OneToAllInputDataPSSM& inputData,
            const GpuConvertedPSSM& strided_PSSM,
            GapScoreArgs /*gapArgs*/,
            cudaStream_t stream
        ) override{
            if(strided_PSSM.kindOfPssm != GpuConvertedPSSM::Kind::Gapless){
                throw std::runtime_error("invalid pssm kind");
            }
            if constexpr(!enableMultiTile){
                throw std::runtime_error("attempted to call scoreOnly_multiTile, but multiTile was not enabled");
            }else{
                if(strided_PSSM.groupsize != groupsize) throw std::runtime_error("pssm groupsize does not match aligner groupsize");
                if(strided_PSSM.numRegs != numRegs) throw std::runtime_error("pssm numRegs does not match aligner numRegs");

                const int numAlignments = inputData.getNumAlignments();

                constexpr int relaxChunkSize = 16 / sizeof(PssmScoreType);
                constexpr int pssmReplicationFactor = PssmReplicationFactor<groupsize, PssmScoreType, relaxChunkSize>::value();
        
                constexpr int numRowsPSSM = 21;
                constexpr int pssmTileSize = groupsize * relaxChunkSize * SDIV(numRegs, relaxChunkSize);
                constexpr int numColumnsPSSM = pssmTileSize * pssmReplicationFactor;
                using SharedPSSM = SharedPSSM_singletile<PssmScoreType, numRowsPSSM, numColumnsPSSM>;

                int smem = sizeof(SharedPSSM);
                auto kernel = [&](){
                    if constexpr(ApproachAndType::approach == Approach::hardcodedzero){
                        return hardcodedzero::GaplessFilter_strided_PSSM_multitile_kernel<ScoreType, blocksize, groupsize, numRegs, subjectIsCaseSensitive,
                            int*, PssmScoreType, OneToAllInputDataPSSM>;
                    }else if constexpr(ApproachAndType::approach == Approach::kernelparamzero){
                        return kernelparamzero::GaplessFilter_strided_PSSM_multitile_kernel<ScoreType, blocksize, groupsize, numRegs, subjectIsCaseSensitive,
                            int*, PssmScoreType, OneToAllInputDataPSSM>;
                    }else if constexpr(ApproachAndType::approach == Approach::half2fmarelu){
                        return half2fmarelu::GaplessFilter_strided_PSSM_multitile_kernel<ScoreType, blocksize, groupsize, numRegs, subjectIsCaseSensitive,
                            int*, PssmScoreType, OneToAllInputDataPSSM>;
                    }
                }();

                auto setSmemKernelAttribute = [&](){
                    static std::map<int, bool> isSet;
                    if(smem > 48*1024){
                        int deviceId;
                        CUDACHECK(cudaGetDevice(&deviceId));
                        if(!isSet[deviceId]){
                            cudaError_t status = cudaFuncSetAttribute(reinterpret_cast<const void*>(kernel), cudaFuncAttributeMaxDynamicSharedMemorySize, smem);
                            if(status != cudaSuccess){
                                cudaGetLastError();
                                throw std::runtime_error("Unable to set cudaFuncAttributeMaxDynamicSharedMemorySize to " + std::to_string(smem));
                            }
                            isSet[deviceId] = true;
                        }
                    }
                };
                setSmemKernelAttribute();

                int maxBlocksPerSM = 0;
                CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
                    &maxBlocksPerSM,
                    kernel,
                    blocksize, 
                    smem
                ));

                const size_t tileTempBytesPerGroup = getTileTempBytesPerGroup(inputData);
                constexpr int groupsPerBlock = (blocksize / groupsize);
                constexpr int alignmentsPerBlock = groupsPerBlock;
                const int maxNumBlocksByInputSize = (numAlignments + alignmentsPerBlock - 1) / alignmentsPerBlock;
                const int maxNumBlocksByOccupancy = maxBlocksPerSM * numSMs;
                const int maxNumBlocksByTempBytes = tempBytes / (tileTempBytesPerGroup * groupsPerBlock);
        
                const int numBlocks = std::min(maxNumBlocksByTempBytes, std::min(maxNumBlocksByInputSize, maxNumBlocksByOccupancy));
                if(numBlocks <= 0){
                    throw std::runtime_error("attempted to launch multitile kernel with 0 blocks");
                }

                dim3 grid = numBlocks;

                if constexpr(ApproachAndType::approach == Approach::hardcodedzero){

                    kernel<<<grid, blocksize, smem, stream>>>(
                        devAlignmentScores,
                        inputData,
                        getPssmView(strided_PSSM),
                        (float2*)d_temp,
                        tileTempBytesPerGroup / sizeof(float2)
                    ); CUERR;
                }else if constexpr(ApproachAndType::approach == Approach::kernelparamzero){
                    ScoreType zero = MathOps<ScoreType>::zero_score();

                    kernel<<<grid, blocksize, smem, stream>>>(
                        devAlignmentScores,
                        inputData,
                        getPssmView(strided_PSSM),
                        (float2*)d_temp,
                        tileTempBytesPerGroup / sizeof(float2),
                        zero
                    ); CUERR;
                }else if constexpr(ApproachAndType::approach == Approach::half2fmarelu){

                    kernel<<<grid, blocksize, smem, stream>>>(
                        devAlignmentScores,
                        inputData,
                        getPssmView(strided_PSSM),
                        (float2*)d_temp,
                        tileTempBytesPerGroup / sizeof(float2)
                    ); CUERR;
                }
            }
        }

    private:
        size_t getTileTempBytesPerGroup(const OneToAllInputDataPSSM& inputData) const{
            static_assert(sizeof(ScoreType) == 4); //e.g. 2 halfs

            const int maximumSequenceLengthPadded = inputData.getMaximumSubjectLength() + groupsize;

            //1 half per position, 4 halfs are stored vectorized to temp storage, round up appropriatly to ensure pointer alignment
            const size_t tileTempBytesPerGroup = cuda::round_up(2 * maximumSequenceLengthPadded, 8);
            return tileTempBytesPerGroup;
        }

        PSSM_2D_View<PssmScoreType> getPssmView(const GpuConvertedPSSM& pssm) const{
            return pssm.makeVec2View<PssmScoreType>();
        }
    };

    template<
        int alphabetSize,
        class ApproachAndType,
        int blocksize,
        int groupsize,
        int numRegs,
        bool enableSingleTile,
        bool enableMultiTile
    >
    struct GaplessAlignment_vec4score : public ScoreOnlyAlignmentInterface{
    public:
        int deviceId;
        int maxSharedMemoryPerBlockOptin;
        int numSMs;

        using ScoreType = typename ApproachAndType::ScoreType;
        static_assert(IsVec4Type<ScoreType>::value);

        using PssmScoreType = ScoreType_u8x4;

        static_assert(enableSingleTile || enableMultiTile);

        static constexpr int numVectorLanes = 4;

        static constexpr int tileSize = groupsize * numRegs * numVectorLanes;

        static constexpr bool subjectIsCaseSensitive = true;

        GaplessAlignment_vec4score(){
            CUDACHECK(cudaGetDevice(&deviceId));
            CUDACHECK(cudaDeviceGetAttribute(&maxSharedMemoryPerBlockOptin, cudaDevAttrMaxSharedMemoryPerBlockOptin, deviceId));
            CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));
        }

        int getTileSize() const override{
            return tileSize;
        }

        AlignmentAlgorithmE getAlignmentAlgorithm() const override{
            return AlignmentAlgorithmE::Gapless;
        }

        Datatype getScoreDatatype() const override{
            return mapDatatypeToDatatypeEnum<ScoreType>();
        }

        bool isSingleTile(const OneToAllInputDataPSSM& inputData) const override{
            if constexpr(enableSingleTile){
                return inputData.getQueryLength() <= tileSize;
            }else{
                return false;
            }
        }

        void makeGpuPssm(
            GpuConvertedPSSM& result,
            const cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>> hostPssmView,
            cudaStream_t stream
        ) override{
            makeConvertedGpuPssm_gapless<ScoreType, PssmScoreType, std::int8_t, groupsize, numRegs>(
                result,
                hostPssmView,
                stream
            );
        }


        size_t getMinimumSuggestedTempBytes_multiTile(const OneToAllInputDataPSSM& inputData) const override{
            constexpr int numGroupsPerBlock = blocksize / groupsize;
            constexpr int alignmentsPerBlock = numGroupsPerBlock;
            const size_t tileTempBytesPerGroup = getTileTempBytesPerGroup(inputData);
            const int numBlocks = cuda::ceil_div(inputData.getNumAlignments(), alignmentsPerBlock);
            const size_t tempBytes1BlockPerSM = tileTempBytesPerGroup * numGroupsPerBlock * std::min(numSMs, numBlocks);
            return tempBytes1BlockPerSM;
        }

        void scoreOnly_singleTile(
            int* const devAlignmentScores,
            const OneToAllInputDataPSSM& inputData,
            const GpuConvertedPSSM& strided_PSSM,
            GapScoreArgs /*gapArgs*/,
            cudaStream_t stream
        ) override{
            if(strided_PSSM.kindOfPssm != GpuConvertedPSSM::Kind::Gapless){
                throw std::runtime_error("invalid pssm kind");
            }
            if constexpr(!enableSingleTile){
                throw std::runtime_error("attempted to call scoreOnly_singleTile, but singleTile was not enabled");
            }else{
                if(strided_PSSM.groupsize != groupsize) throw std::runtime_error("pssm groupsize does not match aligner groupsize");
                if(strided_PSSM.numRegs != numRegs) throw std::runtime_error("pssm numRegs does not match aligner numRegs");
                if(!strided_PSSM.substitutionScoreBias_optional.has_value()) throw std::runtime_error("pssm substitution score bias is unknown");

                const int numAlignments = inputData.getNumAlignments();

                const int substitutionScoreBias = strided_PSSM.substitutionScoreBias_optional.value();
                if(substitutionScoreBias > 255){
                    throw std::runtime_error("substitutionScoreBias cannot be represented with 8 bits.");
                }
                const ScoreType substitutionScoreBias_vec4(
                    substitutionScoreBias,
                    substitutionScoreBias,
                    substitutionScoreBias,
                    substitutionScoreBias
                );

                constexpr int groupsPerBlock = blocksize / groupsize;
                constexpr int alignmentsPerGroup = 1;
                constexpr int alignmentsPerBlock = groupsPerBlock * alignmentsPerGroup;
                // std::cout << "blocksize " << blocksize << ", groupsize " << groupsize 
                //     << ", alignmentsPerBlock " << alignmentsPerBlock << ", numAlignments " << numAlignments << "\n";

                constexpr int relaxChunkSize = 16 / sizeof(PssmScoreType);
                constexpr int pssmReplicationFactor = PssmReplicationFactor<groupsize, PssmScoreType, relaxChunkSize>::value();
        
                constexpr int numRowsPSSM = 21;
                constexpr int pssmTileSize = groupsize * relaxChunkSize * SDIV(numRegs, relaxChunkSize);
                constexpr int numColumnsPSSM = pssmTileSize * pssmReplicationFactor;
                using SharedPSSM = SharedPSSM_singletile<PssmScoreType, numRowsPSSM, numColumnsPSSM>;
                
                int smem = sizeof(SharedPSSM);
                auto kernel = [&](){
                    return uint8x4::GaplessFilter_strided_PSSM_singletile_uint8x4_kernel<
                        ScoreType,
                        blocksize, 
                        groupsize, 
                        numRegs, 
                        false, //false = score only
                        subjectIsCaseSensitive,
                        int*,
                        PssmScoreType, 
                        OneToAllInputDataPSSM>;
                }();

                auto setSmemKernelAttribute = [&](){
                    static std::map<int, bool> isSet;
                    if(smem > 48*1024){
                        int deviceId;
                        CUDACHECK(cudaGetDevice(&deviceId));
                        if(!isSet[deviceId]){
                            cudaError_t status = cudaFuncSetAttribute(reinterpret_cast<const void*>(kernel), cudaFuncAttributeMaxDynamicSharedMemorySize, smem);
                            if(status != cudaSuccess){
                                cudaGetLastError();
                                throw std::runtime_error("Unable to set cudaFuncAttributeMaxDynamicSharedMemorySize to " + std::to_string(smem));
                            }
                            isSet[deviceId] = true;
                        }
                    }
                };
                setSmemKernelAttribute();

                dim3 grid = (numAlignments + alignmentsPerBlock - 1) / alignmentsPerBlock;
                // std::cout << "GaplessFilter_strided_PSSM_singletile_uint8x4_kernel " << grid.x << " " << blocksize << " " << smem << ", numAlignments " << numAlignments << "\n";
                kernel<<<grid, blocksize, smem, stream>>>(
                    devAlignmentScores,
                    nullptr,
                    nullptr,
                    inputData,
                    getPssmView(strided_PSSM),
                    substitutionScoreBias_vec4
                ); CUERR;
            }
        }

        void scoreOnly_multiTile(
            char* d_temp,
            size_t tempBytes,
            int* const devAlignmentScores,
            const OneToAllInputDataPSSM& inputData,
            const GpuConvertedPSSM& strided_PSSM,
            GapScoreArgs /*gapArgs*/,
            cudaStream_t stream
        ) override{
            if(strided_PSSM.kindOfPssm != GpuConvertedPSSM::Kind::Gapless){
                throw std::runtime_error("invalid pssm kind");
            }
            if constexpr(!enableMultiTile){
                throw std::runtime_error("attempted to call scoreOnly_multiTile, but multiTile was not enabled");
            }else{
                if(strided_PSSM.groupsize != groupsize) throw std::runtime_error("pssm groupsize does not match aligner groupsize");
                if(strided_PSSM.numRegs != numRegs) throw std::runtime_error("pssm numRegs does not match aligner numRegs");
                if(!strided_PSSM.substitutionScoreBias_optional.has_value()) throw std::runtime_error("pssm substitution score bias is unknown");

                const int numAlignments = inputData.getNumAlignments();

                const int substitutionScoreBias = strided_PSSM.substitutionScoreBias_optional.value();
                if(substitutionScoreBias > 255){
                    throw std::runtime_error("substitutionScoreBias cannot be represented with 8 bits.");
                }
                const ScoreType substitutionScoreBias_vec4(
                    substitutionScoreBias,
                    substitutionScoreBias,
                    substitutionScoreBias,
                    substitutionScoreBias
                );

                constexpr int relaxChunkSize = 16 / sizeof(PssmScoreType);
                constexpr int pssmReplicationFactor = PssmReplicationFactor<groupsize, PssmScoreType, relaxChunkSize>::value();
        
                constexpr int numRowsPSSM = 21;
                constexpr int pssmTileSize = groupsize * relaxChunkSize * SDIV(numRegs, relaxChunkSize);
                constexpr int numColumnsPSSM = pssmTileSize * pssmReplicationFactor;
                using SharedPSSM = SharedPSSM_singletile<PssmScoreType, numRowsPSSM, numColumnsPSSM>;

                int smem = sizeof(SharedPSSM);
                auto kernel = [&](){
                    return uint8x4::GaplessFilter_strided_PSSM_multitile_uint8x4_kernel<ScoreType, blocksize, groupsize, 
                        numRegs, 
                        false, //false = score only
                        subjectIsCaseSensitive,
                        int*, PssmScoreType, OneToAllInputDataPSSM>;
                }();

                auto setSmemKernelAttribute = [&](){
                    static std::map<int, bool> isSet;
                    if(smem > 48*1024){
                        int deviceId;
                        CUDACHECK(cudaGetDevice(&deviceId));
                        if(!isSet[deviceId]){
                            cudaError_t status = cudaFuncSetAttribute(reinterpret_cast<const void*>(kernel), cudaFuncAttributeMaxDynamicSharedMemorySize, smem);
                            if(status != cudaSuccess){
                                cudaGetLastError();
                                throw std::runtime_error("Unable to set cudaFuncAttributeMaxDynamicSharedMemorySize to " + std::to_string(smem));
                            }
                            isSet[deviceId] = true;
                        }
                    }
                };
                setSmemKernelAttribute();

                int maxBlocksPerSM = 0;
                CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
                    &maxBlocksPerSM,
                    kernel,
                    blocksize, 
                    smem
                ));

                const size_t tileTempBytesPerGroup = getTileTempBytesPerGroup(inputData);
                constexpr int groupsPerBlock = (blocksize / groupsize);
                constexpr int alignmentsPerBlock = groupsPerBlock;
                const int maxNumBlocksByInputSize = (numAlignments + alignmentsPerBlock - 1) / alignmentsPerBlock;
                const int maxNumBlocksByOccupancy = maxBlocksPerSM * numSMs;
                const int maxNumBlocksByTempBytes = tempBytes / (tileTempBytesPerGroup * groupsPerBlock);
                // std::cout << "tileTempBytesPerGroup " << tileTempBytesPerGroup << ", "
                //     << "groupsPerBlock " << groupsPerBlock << ", "
                //     << "alignmentsPerBlock " << alignmentsPerBlock << ", "
                //     << "maxNumBlocksByInputSize " << maxNumBlocksByInputSize << ", "
                //     << "maxNumBlocksByOccupancy " << maxNumBlocksByOccupancy << ", "
                //     << "maxNumBlocksByTempBytes " << maxNumBlocksByTempBytes << "\n";

                const int numBlocks = std::min(maxNumBlocksByTempBytes, std::min(maxNumBlocksByInputSize, maxNumBlocksByOccupancy));
                if(numBlocks <= 0){
                    throw std::runtime_error("attempted to launch multitile kernel with 0 blocks");
                }

                dim3 grid = numBlocks;

                kernel<<<grid, blocksize, smem, stream>>>(
                    devAlignmentScores,
                    nullptr,
                    nullptr,
                    inputData,
                    getPssmView(strided_PSSM),
                    substitutionScoreBias_vec4,
                    (std::uint32_t*)d_temp,
                    tileTempBytesPerGroup / sizeof(std::uint32_t)
                ); CUERR;
            }
        }

    private:
        size_t getTileTempBytesPerGroup(const OneToAllInputDataPSSM& inputData) const{
            static_assert(sizeof(ScoreType) == 4); //e.g. 1 uint8x4

            const int maximumSequenceLengthPadded = inputData.getMaximumSubjectLength() + groupsize;

            //1 int8 per position, 4 int8 are stored vectorized to temp storage, round up appropriatly to ensure pointer alignment
            const size_t tileTempBytesPerGroup = cuda::round_up(1ull * maximumSequenceLengthPadded, 4ull);
            return tileTempBytesPerGroup;
        }


        PSSM_2D_View<PssmScoreType> getPssmView(const GpuConvertedPSSM& pssm) const{
            return pssm.makeVec4View<PssmScoreType>();
        }
    };




    template<
        int alphabetSize,
        class ApproachAndType,
        int blocksize,
        int groupsize,
        int numRegs,
        bool enableSingleTile,
        bool enableMultiTile
    >
    struct SWAlignment_scalarscore : public ScoreOnlyAlignmentInterface{
    public:
        int deviceId;
        int maxSharedMemoryPerBlockOptin;
        int numSMs;

        using ScoreType = typename ApproachAndType::ScoreType;

        using PssmScoreType = ScoreType;
        

        static_assert(enableSingleTile || enableMultiTile);

        static constexpr int tileSize = groupsize * numRegs;

        static constexpr bool subjectIsCaseSensitive = true;
        static constexpr bool withEndPosition = false;

        SWAlignment_scalarscore(){
            CUDACHECK(cudaGetDevice(&deviceId));
            CUDACHECK(cudaDeviceGetAttribute(&maxSharedMemoryPerBlockOptin, cudaDevAttrMaxSharedMemoryPerBlockOptin, deviceId));
            CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));
        }

        int getTileSize() const override{
            return tileSize;
        }

        AlignmentAlgorithmE getAlignmentAlgorithm() const override{
            return AlignmentAlgorithmE::SW;
        }

        Datatype getScoreDatatype() const override{
            return mapDatatypeToDatatypeEnum<ScoreType>();
        }

        bool isSingleTile(const OneToAllInputDataPSSM& inputData) const override{
            if constexpr(enableSingleTile){
                return inputData.getQueryLength() <= tileSize;
            }else{
                return false;
            }
        }


        void makeGpuPssm(
            GpuConvertedPSSM& result,
            const cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>> hostPssmView,
            cudaStream_t stream
        ) override{
            makeConvertedGpuPssm_sw<ScoreType, PssmScoreType, std::int8_t, groupsize, numRegs>(
                result,
                hostPssmView,
                stream
            );
        }

        size_t getMinimumSuggestedTempBytes_multiTile(const OneToAllInputDataPSSM& inputData) const override{
            constexpr int numGroupsPerBlock = blocksize / groupsize;
            constexpr int alignmentsPerBlock = numGroupsPerBlock;
            const size_t tileTempBytesPerGroup = getTileTempBytesPerGroup(inputData);
            const int numBlocks = cuda::ceil_div(inputData.getNumAlignments(), alignmentsPerBlock);
            const size_t tempBytes1BlockPerSM = tileTempBytesPerGroup * numGroupsPerBlock * std::min(numSMs, numBlocks);
            return tempBytes1BlockPerSM;
        }

        void scoreOnly_singleTile(
            int* const devAlignmentScores,
            const OneToAllInputDataPSSM& inputData,
            const GpuConvertedPSSM& strided_PSSM,
            GapScoreArgs gapArgs,
            cudaStream_t stream
        ) override{
            if(strided_PSSM.kindOfPssm != GpuConvertedPSSM::Kind::SW){
                throw std::runtime_error("invalid pssm kind");
            }
            if constexpr(!enableSingleTile){
                throw std::runtime_error("attempted to call scoreOnly_singleTile, but singleTile was not enabled");
            }else{

                if(strided_PSSM.groupsize != groupsize) throw std::runtime_error("pssm groupsize does not match aligner groupsize");
                if(strided_PSSM.numRegs != numRegs) throw std::runtime_error("pssm numRegs does not match aligner numRegs");

                const int numAlignments = inputData.getNumAlignments();

                constexpr int groupsPerBlock = blocksize / groupsize;
                constexpr int alignmentsPerGroup = 1;
                constexpr int alignmentsPerBlock = groupsPerBlock * alignmentsPerGroup;
                // std::cout << "blocksize " << blocksize << ", groupsize " << groupsize 
                //     << ", alignmentsPerBlock " << alignmentsPerBlock << ", numAlignments " << numAlignments << "\n";
               
                constexpr int numRowsPSSM = 21;
                constexpr int pssmTileSize = groupsize * numRegs;
                constexpr int numColumnsPSSM = pssmTileSize;

                
                using SharedPSSM = SharedPSSM_singletile<PssmScoreType, numRowsPSSM, numColumnsPSSM>;
                // std::cout << "host: pssmTileSize = " << pssmTileSize << ", numColumnsPSSM = " << numColumnsPSSM << "\n";
                // std::cout << "host: smem pssm takes " << sizeof(SharedPSSM) << " bytes\n";

                
                int smem = sizeof(SharedPSSM);
                auto kernel = [&](){
                    return amino_gpu_localAlignmentKernel_affinegap_floatOrInt_pssm_singletile<
                        ScoreType,
                        blocksize, 
                        groupsize, 
                        numRegs, 
                        withEndPosition,
                        subjectIsCaseSensitive,
                        OneToAllInputDataPSSM
                    >;
                }();
        
                auto setSmemKernelAttribute = [&](){
                    static std::map<int, bool> isSet;
                    if(smem > 48*1024){
                        int deviceId;
                        CUDACHECK(cudaGetDevice(&deviceId));
                        if(!isSet[deviceId]){
                            cudaError_t status = cudaFuncSetAttribute(reinterpret_cast<const void*>(kernel), cudaFuncAttributeMaxDynamicSharedMemorySize, smem);
                            if(status != cudaSuccess){
                                cudaGetLastError();
                                throw std::runtime_error("Unable to set cudaFuncAttributeMaxDynamicSharedMemorySize to " + std::to_string(smem));
                            }
                            isSet[deviceId] = true;
                        }
                    }
                };
                setSmemKernelAttribute();
        
                dim3 grid = (numAlignments + alignmentsPerBlock - 1) / alignmentsPerBlock;        

                int* subjectEndPositions_exclusive = nullptr;
                int* queryEndPositions_exclusive = nullptr;

                kernel<<<grid, blocksize, smem, stream>>>(
                    devAlignmentScores,
                    subjectEndPositions_exclusive,
                    queryEndPositions_exclusive,
                    inputData,
                    getPssmView(strided_PSSM),
                    gapArgs.gapopenscore, 
                    gapArgs.gapextendscore
                ); CUERR;

            }
        }

        void scoreOnly_multiTile(
            char* d_temp,
            size_t tempBytes,
            int* const devAlignmentScores,
            const OneToAllInputDataPSSM& inputData,
            const GpuConvertedPSSM& strided_PSSM,
            GapScoreArgs gapArgs,
            cudaStream_t stream
        ) override{
            if(strided_PSSM.kindOfPssm != GpuConvertedPSSM::Kind::SW){
                throw std::runtime_error("invalid pssm kind");
            }
            if constexpr(!enableMultiTile){
                throw std::runtime_error("attempted to call scoreOnly_multiTile, but multiTile was not enabled");
            }else{
                if(strided_PSSM.groupsize != groupsize) throw std::runtime_error("pssm groupsize does not match aligner groupsize");
                if(strided_PSSM.numRegs != numRegs) throw std::runtime_error("pssm numRegs does not match aligner numRegs");

                const int numAlignments = inputData.getNumAlignments();
        
                constexpr int numRowsPSSM = 21;
                constexpr int pssmTileSize = groupsize * numRegs;
                constexpr int numColumnsPSSM = pssmTileSize;
                using SharedPSSM = SharedPSSM_singletile<PssmScoreType, numRowsPSSM, numColumnsPSSM>;

                int smem = sizeof(SharedPSSM);
                auto kernel = [&](){
                    return amino_gpu_localAlignmentKernel_affinegap_floatOrInt_pssm_multitile<
                        ScoreType,
                        blocksize, 
                        groupsize, 
                        numRegs, 
                        withEndPosition,
                        subjectIsCaseSensitive,
                        OneToAllInputDataPSSM
                    >;
                }();

                auto setSmemKernelAttribute = [&](){
                    static std::map<int, bool> isSet;
                    if(smem > 48*1024){
                        int deviceId;
                        CUDACHECK(cudaGetDevice(&deviceId));
                        if(!isSet[deviceId]){
                            cudaError_t status = cudaFuncSetAttribute(reinterpret_cast<const void*>(kernel), cudaFuncAttributeMaxDynamicSharedMemorySize, smem);
                            if(status != cudaSuccess){
                                cudaGetLastError();
                                throw std::runtime_error("Unable to set cudaFuncAttributeMaxDynamicSharedMemorySize to " + std::to_string(smem));
                            }
                            isSet[deviceId] = true;
                        }
                    }
                };
                setSmemKernelAttribute();

                int maxBlocksPerSM = 0;
                CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
                    &maxBlocksPerSM,
                    kernel,
                    blocksize, 
                    smem
                ));

                const size_t tileTempBytesPerGroup = getTileTempBytesPerGroup(inputData);
                constexpr int groupsPerBlock = (blocksize / groupsize);
                constexpr int alignmentsPerBlock = groupsPerBlock;
                const int maxNumBlocksByInputSize = (numAlignments + alignmentsPerBlock - 1) / alignmentsPerBlock;
                const int maxNumBlocksByOccupancy = maxBlocksPerSM * numSMs;
                const int maxNumBlocksByTempBytes = tempBytes / (tileTempBytesPerGroup * groupsPerBlock);
        
                const int numBlocks = std::min(maxNumBlocksByTempBytes, std::min(maxNumBlocksByInputSize, maxNumBlocksByOccupancy));
                if(numBlocks <= 0){
                    throw std::runtime_error("attempted to launch multitile kernel with 0 blocks");
                }

                dim3 grid = numBlocks;

                int* subjectEndPositions_exclusive = nullptr;
                int* queryEndPositions_exclusive = nullptr;
                kernel<<<grid, blocksize, smem, stream>>>(
                    devAlignmentScores,
                    subjectEndPositions_exclusive,
                    queryEndPositions_exclusive,
                    inputData,
                    getPssmView(strided_PSSM),
                    gapArgs.gapopenscore, 
                    gapArgs.gapextendscore,
                    d_temp,
                    tileTempBytesPerGroup
                ); CUERR;
            }
        }

    private:
        size_t getTileTempBytesPerGroup(const OneToAllInputDataPSSM& inputData) const{
            using TempStorageDataType = float2; //two floats per position (affine scores)

            const int maximumSequenceLengthPadded = inputData.getMaximumSubjectLength() + groupsize;
            const size_t tileTempBytesPerGroup = sizeof(TempStorageDataType) * maximumSequenceLengthPadded;
            return tileTempBytesPerGroup;
        }

        PSSM_2D_View<PssmScoreType> getPssmView(const GpuConvertedPSSM& pssm) const{
            return pssm.makeScalarView<PssmScoreType>();
        }
    };






    template<class SingleTileConfigs, class MultiTileConfigs>
    struct MultiConfigScoreOnlyAligner : public MultiConfigScoreOnlyAlignerInterface{

        static_assert(configsAreSortedByTileSize<typename SingleTileConfigs::type>());
        static_assert(configsAreSortedByTileSize<typename MultiTileConfigs::type>());

        std::vector<std::unique_ptr<ScoreOnlyAlignmentInterface>> singleTileAlignerVector;
        std::vector<std::unique_ptr<ScoreOnlyAlignmentInterface>> multiTileAlignerVector;

        std::vector<int> performanceHintsSingleTile;
        std::vector<int> performanceHintsMultiTile;

        AlignmentAlgorithmE algorithm;

        MultiConfigScoreOnlyAligner(AlignmentAlgorithmE algorithm_) : algorithm(algorithm_){

            auto addSingleTileAligner = [&](auto config, auto& vector){
                using Config = decltype(config);
                using ApproachAndType = ApproachAndType<Config::approach, Config::datatype>;

                if(Config::algorithm != algorithm){
                    throw std::runtime_error("Aligner algorithm misconfigured");
                }

                if constexpr(Config::algorithm == AlignmentAlgorithmE::Gapless){

                    
                    if constexpr(IsScalarType<typename ApproachAndType::ScoreType>::value){
                        using AlignerClass = GaplessAlignment_scalarscore<
                            Config::alphabetSize, 
                            ApproachAndType,
                            Config::blocksize, Config::groupsize, Config::numRegs,
                            true,
                            false
                        >;
                        vector.push_back(std::make_unique<AlignerClass>());    
                    }else if constexpr(IsVec2Type<typename ApproachAndType::ScoreType>::value){
                        using AlignerClass = GaplessAlignment_vec2score<
                            Config::alphabetSize, 
                            ApproachAndType,
                            Config::blocksize, Config::groupsize, Config::numRegs,
                            true,
                            false
                        >;
                        vector.push_back(std::make_unique<AlignerClass>());    
                    }else if constexpr(IsVec4Type<typename ApproachAndType::ScoreType>::value){

                        using AlignerClass = GaplessAlignment_vec4score<
                            Config::alphabetSize, 
                            ApproachAndType,
                            Config::blocksize, Config::groupsize, Config::numRegs,
                            true,
                            false
                        >;
                        vector.push_back(std::make_unique<AlignerClass>());
                    }else{
                        throw std::runtime_error("unexpected score type");
                    }
                }else if constexpr(Config::algorithm == AlignmentAlgorithmE::SW){
                    using AlignerClass = SWAlignment_scalarscore<
                        Config::alphabetSize, 
                        ApproachAndType,
                        Config::blocksize, Config::groupsize, Config::numRegs,
                        true,
                        false
                    >;
                    vector.push_back(std::make_unique<AlignerClass>());  
                }else{
                    throw std::runtime_error("unexpected config kind");
                }
            };

            auto addMultiTileAligner = [&](auto config, auto& vector){
                using Config = decltype(config);
                using ApproachAndType = ApproachAndType<Config::approach, Config::datatype>;

                if(Config::algorithm != algorithm){
                    throw std::runtime_error("Aligner algorithm misconfigured");
                }

                if constexpr(Config::algorithm == AlignmentAlgorithmE::Gapless){
                
                    if constexpr(IsScalarType<typename ApproachAndType::ScoreType>::value){
                        using AlignerClass = GaplessAlignment_scalarscore<
                            Config::alphabetSize, 
                            ApproachAndType,
                            Config::blocksize, Config::groupsize, Config::numRegs,
                            false,
                            true
                        >;
                        vector.push_back(std::make_unique<AlignerClass>());    
                    }else if constexpr(IsVec2Type<typename ApproachAndType::ScoreType>::value){
                        using AlignerClass = GaplessAlignment_vec2score<
                            Config::alphabetSize, 
                            ApproachAndType,
                            Config::blocksize, Config::groupsize, Config::numRegs,
                            false,
                            true
                        >;
                        vector.push_back(std::make_unique<AlignerClass>());    
                    }else if constexpr(IsVec4Type<typename ApproachAndType::ScoreType>::value){

                        using AlignerClass = GaplessAlignment_vec4score<
                            Config::alphabetSize, 
                            ApproachAndType,
                            Config::blocksize, Config::groupsize, Config::numRegs,
                            false,
                            true
                        >;
                        vector.push_back(std::make_unique<AlignerClass>());
                    }else{
                        throw std::runtime_error("unexpected score type");
                    }
                }else if constexpr(Config::algorithm == AlignmentAlgorithmE::SW){
                    using AlignerClass = SWAlignment_scalarscore<
                        Config::alphabetSize, 
                        ApproachAndType,
                        Config::blocksize, Config::groupsize, Config::numRegs,
                        false,
                        true
                    >;
                    vector.push_back(std::make_unique<AlignerClass>());  
                }else{
                    throw std::runtime_error("unexpected config kind");
                }
            };

            typename SingleTileConfigs::type singleTileConfigs;
            std::apply(
                [&](auto& ... x){
                    (addSingleTileAligner(x, singleTileAlignerVector), ...);
                },
                singleTileConfigs
            );

            typename MultiTileConfigs::type multiTileConfigs;
            std::apply(
                [&](auto& ... x){
                    (addMultiTileAligner(x, multiTileAlignerVector), ...);
                },
                multiTileConfigs
            );

            if(singleTileAlignerVector.size() == 0 && multiTileAlignerVector.size() == 0){
                throw std::runtime_error("MultiConfigScoreOnlyAligner requires at least one config");
            }

            if constexpr (SingleTileConfigs::hasPerformanceHints()){
                auto hints = SingleTileConfigs::getPerformanceHints();
                performanceHintsSingleTile.insert(performanceHintsSingleTile.end(), 
                    hints.begin(),
                    hints.end());
            }

            if constexpr (MultiTileConfigs::hasPerformanceHints()){
                auto hints = MultiTileConfigs::getPerformanceHints();
                performanceHintsMultiTile.insert(performanceHintsMultiTile.end(), 
                    hints.begin(),
                    hints.end());
            }

        }

        bool isSingleTile(const OneToAllInputDataPSSM& inputData) const override{
            return isSingleTile(inputData.getQueryLength());
        }

        int getTileSize(const OneToAllInputDataPSSM& inputData) const override{
            const auto* aligner = getBestAligner(inputData.getQueryLength());
            return aligner->getTileSize();
        }

        AlignmentAlgorithmE getAlignmentAlgorithm() const override{
            return algorithm;
        }

        Datatype getScoreDatatype(int queryLength) const override{
            const auto* aligner = getBestAligner(queryLength);
            return aligner->getScoreDatatype();
        }

        void makeGpuPssm(
            GpuConvertedPSSM& result,
            const cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>> hostPssmView,
            cudaStream_t stream
        ) override{
            const int queryLength = hostPssmView.extent(1);
            auto* aligner = getBestAligner(queryLength);
            aligner->makeGpuPssm(result, hostPssmView, stream);
        }

        size_t getMinimumSuggestedTempBytes_multiTile(const OneToAllInputDataPSSM& inputData) const override{
            auto* aligner = getBestMultiTileAligner(inputData.getQueryLength());
            return aligner->getMinimumSuggestedTempBytes_multiTile(inputData);
        }

        void scoreOnly_singleTile(
            int* const devAlignmentScores,
            const OneToAllInputDataPSSM& inputData,
            const GpuConvertedPSSM& strided_PSSM,
            GapScoreArgs gapArgs,
            cudaStream_t stream
        ) override{
            assert(inputData.getQueryLength() == getQueryLengthFromPssm(strided_PSSM));

            auto* aligner = getBestSingleTileAligner(getQueryLengthFromPssm(strided_PSSM));

            aligner->scoreOnly_singleTile(
                devAlignmentScores,
                inputData,
                strided_PSSM,
                gapArgs,
                stream
            );
        }

        void scoreOnly_multiTile(
            char* d_temp,
            size_t tempBytes,
            int* const devAlignmentScores,
            const OneToAllInputDataPSSM& inputData,
            const GpuConvertedPSSM& strided_PSSM,
            GapScoreArgs gapArgs,
            cudaStream_t stream
        ) override{
            assert(inputData.getQueryLength() == getQueryLengthFromPssm(strided_PSSM));

            auto* aligner = getBestMultiTileAligner(getQueryLengthFromPssm(strided_PSSM));

            aligner->scoreOnly_multiTile(
                d_temp,
                tempBytes,
                devAlignmentScores,
                inputData,
                strided_PSSM,
                gapArgs,
                stream
            );
        }

    private:
        ScoreOnlyAlignmentInterface* getBestSingleTileAligner(SequenceLengthT queryLength) const {
            if(singleTileAlignerVector.size() == 0){
                throw std::runtime_error("no single tile aligner available");
            }
            ScoreOnlyAlignmentInterface* bestAligner = std::lower_bound(
                singleTileAlignerVector.begin(),
                singleTileAlignerVector.end(),
                queryLength,
                [](const auto& aligner, SequenceLengthT length){
                    return aligner->getTileSize() < length;
                }
            )->get();
            return bestAligner;
        }

        ScoreOnlyAlignmentInterface* getBestMultiTileAligner(SequenceLengthT queryLength) const{
            if(multiTileAlignerVector.size() == 0){
                throw std::runtime_error("no multi tile aligner available");
            }

            ScoreOnlyAlignmentInterface* bestAligner = multiTileAlignerVector[0].get();

            if(performanceHintsMultiTile.size() == 0){
                const int remainderInLastTile0 = queryLength % bestAligner->getTileSize();
                double utilization = remainderInLastTile0 == 0 ? 1.0 : double(remainderInLastTile0) / bestAligner->getTileSize();
                for(size_t i = 1; i < multiTileAlignerVector.size(); i++){
                    ScoreOnlyAlignmentInterface* newAligner = multiTileAlignerVector[i].get();

                    const int remainderInLastTile = queryLength % newAligner->getTileSize();
                    const double newUtilization = remainderInLastTile == 0 ? 1.0 : double(remainderInLastTile) / newAligner->getTileSize();
                    if(newUtilization >= utilization){
                        utilization = newUtilization;
                        bestAligner = newAligner;
                    }

                    // std::cout << "selected based on utilization: tilesize " << bestAligner->getTileSize() << "\n";
                }
            }else{
                const int remainderInLastTile0 = queryLength % bestAligner->getTileSize();
                int numTiles = (queryLength + bestAligner->getTileSize() - 1) / bestAligner->getTileSize();
                double utilization = remainderInLastTile0 == 0 ? 1.0 : double(remainderInLastTile0) / bestAligner->getTileSize();
                double estimatedPerf = utilization * performanceHintsMultiTile[0];
                for(size_t i = 1; i < multiTileAlignerVector.size(); i++){
                    ScoreOnlyAlignmentInterface* newAligner = multiTileAlignerVector[i].get();

                    const int remainderInLastTile = queryLength % newAligner->getTileSize();
                    int newNumTiles = (queryLength + newAligner->getTileSize() - 1) / newAligner->getTileSize();
                    const double newUtilization = remainderInLastTile == 0 ? 1.0 : double(remainderInLastTile) / newAligner->getTileSize();
                    const double newEstimatedPerf = newUtilization * performanceHintsMultiTile[i];
                    if(newEstimatedPerf >= estimatedPerf){
                        estimatedPerf = newEstimatedPerf;
                        bestAligner = newAligner;
                    }
                }

                //std::cout << "selected based on hints: tilesize " << bestAligner->getTileSize() << "\n";
            }

            return bestAligner;
        }

        ScoreOnlyAlignmentInterface* getBestAligner(SequenceLengthT queryLength){
            if(isSingleTile(queryLength)){
                return getBestSingleTileAligner(queryLength);
            }else{                
                return getBestMultiTileAligner(queryLength);
            }
        }

        const ScoreOnlyAlignmentInterface* getBestAligner(SequenceLengthT queryLength) const{
            if(isSingleTile(queryLength)){
                return getBestSingleTileAligner(queryLength);
            }else{                
                return getBestMultiTileAligner(queryLength);
            }
        }

        int getQueryLengthFromPssm(const GpuConvertedPSSM& pssm){
            return pssm.queryLength;
        }

        bool isSingleTile(int queryLength) const{
            if(singleTileAlignerVector.size() == 0){
                return false;
            }else{
                return queryLength <= singleTileAlignerVector.back()->getTileSize();
            }
        }
    };

LIBMARV_NAMESPACE_WITH_NESTING_END


#endif //LIBMARV_SCORE_ONLY_ALIGNMENT_CUH