#include "score_only_foldseek_search.cuh"

#include <cuda/functional>
#include <cuda/cmath>
#include <cuda/std/span>
#include <cuda/std/tuple>

#include <vector>
#include <memory>
#include <algorithm>
#include <numeric>

#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/iterator/discard_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/transform_output_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/logical.h>
#include <thrust/merge.h>

#if defined(__CUDACC__)
#include <cub/cub.cuh>
#else
#include <hipcub/hipcub.hpp>
namespace cub = hipcub;
#endif

#include "../aligners/score_only_alignment_interfaces.cuh"
#include "../aligners/score_with_endpos_alignment_interfaces.cuh"
#include "../cuda_errorcheck.cuh"
#include "../config.hpp"
#include "../pssm.cuh"
#include "../gpudatabase.cuh"
#include "search_common.cuh"
#include "../sequence_masking.cuh"
#include "../diagonal_rescore.cuh"
#include "../stream_ordered_allocation_helper.cuh"

#include "../hpc_helpers/nvtx_markers.cuh"
#include "../hpc_helpers/simple_allocation.cuh"
#include "../hpc_helpers/cuda_raiiwrappers.cuh"


#include "../namespace.hpp"
LIBMARV_NAMESPACE_WITH_NESTING_BEGIN




    template<int maxMaskingThreshold, class InputOffsetIterator, class InputTransformOp>
    void mask3diSequences_maxthreshold(
        char* __restrict__ d_outputData_3di,
        const int* __restrict__ d_inputLengths,
        InputOffsetIterator d_inputOffsets,
        const char* __restrict__ d_inputData_combined,
        int numSequences,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp convertInputTo3di,
        int sequenceLengthUpperBound,
        cudaStream_t stream
    ){
        #define CALLSINGLETHREADSINGLETILE(blocksize, maxtilesize){ \
            maskingsinglethread::callSingleTileWithPlainFlagsKernel<blocksize, maxMaskingThreshold, maxtilesize>( \
                d_outputData_3di, \
                d_inputLengths, \
                d_inputOffsets, \
                d_inputData_combined, \
                numSequences, \
                maskingThreshold, \
                maskingLetter, \
                applyExtraPaddingLetter, \
                paddingLetter, \
                convertInputTo3di, \
                stream \
            ); \
        }
        #define CALLSINGLEGROUPSINGLETILE(blocksize, groupsize, maxtilesize){ \
            maskingsinglegroup::callSingleTileWithPlainFlagsKernel<blocksize, groupsize, maxMaskingThreshold, maxtilesize>( \
                d_outputData_3di, \
                d_inputLengths, \
                d_inputOffsets, \
                d_inputData_combined, \
                numSequences, \
                maskingThreshold, \
                maskingLetter, \
                applyExtraPaddingLetter, \
                paddingLetter, \
                convertInputTo3di, \
                stream \
            ); \
        }
        #define CALLSINGLEGROUPMULTITILE(blocksize, groupsize, maxtilesize){ \
            maskingsinglegroup::callMultiTileWithPlainFlagsKernel<blocksize, groupsize, maxMaskingThreshold, maxtilesize>( \
                d_outputData_3di, \
                d_inputLengths, \
                d_inputOffsets, \
                d_inputData_combined, \
                numSequences, \
                maskingThreshold, \
                maskingLetter, \
                applyExtraPaddingLetter, \
                paddingLetter, \
                convertInputTo3di, \
                stream \
            ); \
        }
        #define CALLSINGLEBLOCKMULTITILE(blocksize, maxtilesize){ \
            maskingsingleblock::callMultiTileWithPlainFlagsKernel<blocksize, maxMaskingThreshold, maxtilesize>( \
                d_outputData_3di, \
                d_inputLengths, \
                d_inputOffsets, \
                d_inputData_combined, \
                numSequences, \
                maskingThreshold, \
                maskingLetter, \
                applyExtraPaddingLetter, \
                paddingLetter, \
                convertInputTo3di, \
                stream \
            ); \
        }
        if(sequenceLengthUpperBound <= 16){
            CALLSINGLETHREADSINGLETILE(128, 16)
        }else if(sequenceLengthUpperBound <= 32){
            CALLSINGLETHREADSINGLETILE(128, 16)
        }else if(sequenceLengthUpperBound <= 32){
            CALLSINGLETHREADSINGLETILE(512, 48)
        }else if(sequenceLengthUpperBound <= 64){
            CALLSINGLETHREADSINGLETILE(512, 64)
        }else if(sequenceLengthUpperBound <= 80){
            CALLSINGLETHREADSINGLETILE(512, 80)
        }else if(sequenceLengthUpperBound <= 96){
            CALLSINGLETHREADSINGLETILE(512, 96)
        }else if(sequenceLengthUpperBound <= 112){
            CALLSINGLETHREADSINGLETILE(256, 112)
        }else if(sequenceLengthUpperBound <= 128){
            CALLSINGLETHREADSINGLETILE(256, 128)
        }else if(sequenceLengthUpperBound <= 160){
            CALLSINGLETHREADSINGLETILE(256, 160)
        }else if(sequenceLengthUpperBound <= 192){
            CALLSINGLETHREADSINGLETILE(256, 192)
        }else if(sequenceLengthUpperBound <= 256){
            CALLSINGLEGROUPSINGLETILE(512, 16, 256)
        }else if(sequenceLengthUpperBound <= 320){
            CALLSINGLEGROUPSINGLETILE(512, 16, 320)
        }else if(sequenceLengthUpperBound <= 384){
            CALLSINGLEGROUPSINGLETILE(512, 16, 384)
        }else if(sequenceLengthUpperBound <= 448){
            CALLSINGLEGROUPSINGLETILE(512, 16, 448)
        }else if(sequenceLengthUpperBound <= 512){
            CALLSINGLEGROUPSINGLETILE(256, 32, 512)
        }else if(sequenceLengthUpperBound <= 640){
            CALLSINGLEGROUPSINGLETILE(256, 32, 640)
        }else if(sequenceLengthUpperBound <= 768){
            CALLSINGLEGROUPSINGLETILE(256, 32, 768)
        }else if(sequenceLengthUpperBound <= 896){
            CALLSINGLEGROUPSINGLETILE(256, 32, 896)
        }else if(sequenceLengthUpperBound <= 1024){
            CALLSINGLEGROUPSINGLETILE(256, 32, 1024)
        }else{
            if(sequenceLengthUpperBound < 8192 || numSequences >= 10000){
                CALLSINGLEGROUPMULTITILE(256, 32, 512)
            }else{
                CALLSINGLEBLOCKMULTITILE(128, 1024)
            }
        }

        #undef CALLSINGLETHREADSINGLETILE
        #undef CALLSINGLEGROUPSINGLETILE
        #undef CALLSINGLEGROUPMULTITILE
        #undef CALLSINGLEBLOCKMULTITILE
    }

    template<class InputOffsetIterator, class InputTransformOp>
    void mask3diSequences(
        char* __restrict__ d_outputData_3di,
        const int* __restrict__ d_inputLengths,
        InputOffsetIterator d_inputOffsets,
        const char* __restrict__ d_inputData_combined,
        int numSequences,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp convertInputTo3di,
        int sequenceLengthUpperBound,
        cudaStream_t stream
    ){
        #define CALLMASKINGMAXTHRESHOLD(maxthreshold){ \
            mask3diSequences_maxthreshold<maxthreshold>( \
                d_outputData_3di, \
                d_inputLengths, \
                d_inputOffsets, \
                d_inputData_combined, \
                numSequences, \
                maskingThreshold, \
                maskingLetter, \
                applyExtraPaddingLetter, \
                paddingLetter, \
                convertInputTo3di, \
                sequenceLengthUpperBound, \
                stream \
            ); \
        }
        if(maskingThreshold <= 6){
            CALLMASKINGMAXTHRESHOLD(6)
        }else if(maskingThreshold <= 16){
            CALLMASKINGMAXTHRESHOLD(16)
        }else{
            throw std::runtime_error("masking threshold > 16 not supported");
        }

        #undef CALLMASKINGMAXTHRESHOLD
    }



    struct ScoreOnlyFoldseekSearch::GpuWorkingSet{
        int maxNumQueries;
        int kernelOutputArraySize;
        std::vector<int> kernelOutputArrayUsedSize;
        int topNSize;

        MyDeviceBuffer<int> d_kernelOutputArrayScores;
        MyDeviceBuffer<ReferenceIdT> d_kernelOutputArrayIndices;
        MyDeviceBuffer<int> d_kernelOutputArraySubjectEndPositionsExcl;
        MyDeviceBuffer<int> d_kernelOutputArrayQueryEndPositionsExcl;

        std::vector<MyDeviceBuffer<int>> d_topN_scores_vec;
        std::vector<MyDeviceBuffer<ReferenceIdT>> d_topN_refIds_vec;
        std::vector<MyDeviceBuffer<int>> d_topN_scores_tmp_vec;
        std::vector<MyDeviceBuffer<ReferenceIdT>> d_topN_refIds_tmp_vec;
    
        // MyDeviceBuffer<char> d_query;
        // MyDeviceBuffer<char> d_tempStorage;

        CudaEvent forkStreamEvent;

        int* getKernelOutputArrayScores(int queryIndex){
            return d_kernelOutputArrayScores.data() + size_t(kernelOutputArraySize) * queryIndex;
        }
        ReferenceIdT* getKernelOutputArrayIndices(int queryIndex){
            return d_kernelOutputArrayIndices.data() + size_t(kernelOutputArraySize) * queryIndex;
        }
        int* getKernelOutputArraySubjectEndPositionsExcl(int queryIndex){
            return d_kernelOutputArraySubjectEndPositionsExcl.data() + size_t(kernelOutputArraySize) * queryIndex;
        }
        int* getKernelOutputArrayQueryEndPositionsExcl(int queryIndex){
            return d_kernelOutputArrayQueryEndPositionsExcl.data() + size_t(kernelOutputArraySize) * queryIndex;
        }
        MyDeviceBuffer<int>& getTopN_scores(int queryIndex){
            return d_topN_scores_vec[queryIndex];
        }
        MyDeviceBuffer<ReferenceIdT>& getTopN_refIds(int queryIndex){
            return d_topN_refIds_vec[queryIndex];
        }
        MyDeviceBuffer<int>& getTopN_scores_tmp(int queryIndex){
            return d_topN_scores_tmp_vec[queryIndex];
        }
        MyDeviceBuffer<ReferenceIdT>& getTopN_refIds_tmp(int queryIndex){
            return d_topN_refIds_tmp_vec[queryIndex];
        }

        GpuWorkingSet(int kernelOutputArraySize_, int topNSize_, int numQueries_)
        : maxNumQueries(numQueries_), kernelOutputArraySize(kernelOutputArraySize_), 
            kernelOutputArrayUsedSize(maxNumQueries, 0), topNSize(topNSize_){
            nvtx::ScopedRange sr("GpuWorkingSet constructor", 0);
            d_kernelOutputArrayScores.resize(kernelOutputArraySize * maxNumQueries);
            d_kernelOutputArrayIndices.resize(kernelOutputArraySize * maxNumQueries);
            d_kernelOutputArraySubjectEndPositionsExcl.resize(kernelOutputArraySize * maxNumQueries);
            d_kernelOutputArrayQueryEndPositionsExcl.resize(kernelOutputArraySize * maxNumQueries);

            for(int i = 0; i < maxNumQueries; i++){
                //the topN results are stored at the begin of buffers.
                //2*topNSize is necessary because thrust::merge will produce 2*topNSize elements
                d_topN_scores_vec.emplace_back(2*topNSize);
                d_topN_refIds_vec.emplace_back(2*topNSize);
                d_topN_scores_tmp_vec.emplace_back(2*topNSize);
                d_topN_refIds_tmp_vec.emplace_back(2*topNSize);
            }


            // d_tempStorage.resize(1ull << 32);

            forkStreamEvent = CudaEvent{cudaEventDisableTiming};
            CUDACHECKASYNC;
        }

    };


    ScoreOnlyFoldseekSearch::ScoreOnlyFoldseekSearch(bool isAA20db, //for testing with a standard protein db
        GpuDatabase& gpuDatabase_, int kernelOutputArraySize_, int topNSize_, int maxNumQueries_) 
        : kernelOutputArraySize(kernelOutputArraySize_), 
        topNSize(int(std::min(size_t(topNSize_), gpuDatabase_.getNumSequences()))), 
        maxNumQueries(maxNumQueries_),
        deviceIds(gpuDatabase_.getDeviceIds()), gpuDatabase(&gpuDatabase_),
        isOrdinaryAA20db(isAA20db)
    {
        initWorkingSets();

        const int numGpus = deviceIds.size();

        processedSequencesPerGpu.resize(numGpus, 0);
        h_finalAlignmentScoresPerQuery.resize(maxNumQueries);
        h_finalReferenceIdsPerQuery.resize(maxNumQueries);
        for(int q = 0; q < maxNumQueries; q++){
            h_finalAlignmentScoresPerQuery[q].resize(topNSize);
            h_finalReferenceIdsPerQuery[q].resize(topNSize);
        }

        resultSizePerQueryPerGpu.resize(maxNumQueries, std::vector<int>(numGpus));
        h_numOverflowsPerGpuPerQuery.resize(maxNumQueries * numGpus);

        RevertDeviceId rd{};
        CUDACHECK(cudaSetDevice(deviceIds[0]));
        scanTimer = std::make_unique<helpers::GpuTimer>("Scan");

        streamsPerGpuPerQuery.reserve(numGpus * maxNumQueries);
        for(int gpu = 0; gpu < numGpus; gpu++){
            CUDACHECK(cudaSetDevice(deviceIds[gpu]));
            for(int queryIndex = 0; queryIndex < maxNumQueries; queryIndex++){
                streamsPerGpuPerQuery.emplace_back(cudaStreamNonBlocking);
            }
        }
        
    }

    ScoreOnlyFoldseekSearch::~ScoreOnlyFoldseekSearch(){
        destroyWorkingSets();
        destroyPssms_3di();
        destroyPssms_aa12();
    }


    void ScoreOnlyFoldseekSearch::reset(){
        nvtx::ScopedRange sr("ScoreOnlyFoldseekSearch::reset", 0);

        clearProcessedSequencesAndCachedResults();

        // const int numGpus = deviceIds.size();
        // for(int gpu = 0; gpu < numGpus; gpu++){
        //     CUDACHECK(cudaSetDevice(deviceIds[gpu]));
        //     auto& ws = *workingSets[gpu];
        //     cudaStream_t stream = cudaStreamPerThread;
        //     thrust::fill(thrust::cuda::par_nosync.on(stream),
        //         ws.d_topN_scores.data(),
        //         ws.d_topN_scores.data() + ws.d_topN_scores.size(),
        //         -1
        //     );
        //     CUDACHECK(cudaMemsetAsync(ws.d_topN_refIds.data(), 0, sizeof(ReferenceIdT) * ws.d_topN_refIds.size(), stream));
        //     thrust::fill(thrust::cuda::par_nosync.on(stream),
        //         ws.d_kernelOutputArrayScores.data(),
        //         ws.d_kernelOutputArrayScores.data() + kernelOutputArraySize,
        //         -1
        //     );
        //     CUDACHECK(cudaMemsetAsync(ws.d_kernelOutputArrayIndices.data(), 0, sizeof(ReferenceIdT) * kernelOutputArraySize, stream));
        // }

        // for(int gpu = 0; gpu < numGpus; gpu++){
        //     CUDACHECK(cudaSetDevice(deviceIds[gpu]));
        //     CUDACHECK(cudaStreamSynchronize(cudaStreamPerThread));
        // }
    }



    void ScoreOnlyFoldseekSearch::destroyWorkingSets(){
        const int numGpus = deviceIds.size();
        for(int i = 0; i < numGpus; i++){
            CUDACHECK(cudaSetDevice(deviceIds[i]));
            workingSets[i].reset();
        }
    }

    void ScoreOnlyFoldseekSearch::destroyPssms_3di(){
        const int numGpus = deviceIds.size();
        const int numPssms = queryPssmPerGpuPerQuery_aa12.size();

        for(int i = 0; i < numGpus; i++){
            CUDACHECK(cudaSetDevice(deviceIds[i]));

            for(int queryIndex = 0; queryIndex < numPssms; queryIndex++){
                queryPssmPerGpuPerQuery_forStandardScoreOnly[queryIndex][i].reset();
                if(scoreOnlyOverflowSettings.checkForOverflow){
                    queryPssmPerGpuPerQuery_forOverflowScoreOnly[queryIndex][i].reset();
                }
                queryPssmPerGpuPerQuery_forStandardEndPos[queryIndex][i].reset();
                if(endposOverflowSettings.checkForOverflow){
                    queryPssmPerGpuPerQuery_forOverflowEndPos[queryIndex][i].reset();
                }
            }
        }

        queryPssmPerGpuPerQuery_forStandardScoreOnly.clear();
        queryPssmPerGpuPerQuery_forOverflowScoreOnly.clear();
        queryPssmPerGpuPerQuery_forStandardEndPos.clear();
        queryPssmPerGpuPerQuery_forOverflowEndPos.clear();
    }

    void ScoreOnlyFoldseekSearch::destroyPssms_aa12(){
        const int numGpus = deviceIds.size();
        const int numPssms = queryPssmPerGpuPerQuery_aa12.size();

        for(int i = 0; i < numGpus; i++){
            CUDACHECK(cudaSetDevice(deviceIds[i]));

            for(int queryIndex = 0; queryIndex < numPssms; queryIndex++){
                queryPssmPerGpuPerQuery_aa12[queryIndex][i].reset();
            }
        }

        queryPssmPerGpuPerQuery_aa12.clear();
    }


    std::vector<SearchResult> ScoreOnlyFoldseekSearch::execute(
        cuda::std::span<MultiConfigScoreOnlyAlignerInterface*> alignersScoreOnly,
        const OverflowSettings& scoreOnlyOverflowSettings_,
        cuda::std::span<MultiConfigScoreWithEndPosAlignerInterface*> alignersEndPos,
        const OverflowSettings& endposOverflowSettings_,
        cuda::std::span<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>> queryPssmViews_3di,
        cuda::std::span<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>> queryPssmViews_aa12,
        MaskingOptions maskingOptions_,
        int gapopenscore,
        int gapextendscore
    ){
        nvtx::ScopedRange sr("ScoreOnlyFoldseekSearch::execute", 0);
        helpers::GpuTimer timer("Total");
        clearProcessedSequencesAndCachedResults();

        setAligner(alignersScoreOnly, scoreOnlyOverflowSettings_, alignersEndPos, endposOverflowSettings_);
        setPssms_3di(queryPssmViews_3di);
        setPssms_aa12(queryPssmViews_aa12);

        gop = gapopenscore;
        gex = gapextendscore;
        maskingOptions = maskingOptions_;

        if(workingSets.size() == 0){
            initWorkingSets();
        }

        const int numQueries = queryPssmViews_3di.size();

        std::vector<cudaStream_t> streams(deviceIds.size(), cudaStreamPerThread);
        cuda::std::span<cudaStream_t> streamsSpan(streams.data(), streams.size());

        scanTimer->reset();
        scanTimer->start();

        constexpr bool processCachedDBInSeparateBatches = true;
        // score only search in complete database
        gpuDatabase->forEachBatchOfSequencesInDatabase(
            std::bind(&ScoreOnlyFoldseekSearch::forEachCallback_fullscan, this, std::placeholders::_1), 
            streamsSpan,
            processCachedDBInSeparateBatches
        );

        finalizeFullScanResults(streamsSpan);

        scanTimer->stop();

        double firstPassSeconds = scanTimer->elapsed() / 1000.0;

        const int sumOfQueryLengths = std::reduce(currentQueryLengths.begin(), currentQueryLengths.end());
        const auto& sequenceLengthStatistics = gpuDatabase->getSequenceLengthStatistics();
        size_t computedCellsFirstPass = size_t(sequenceLengthStatistics.sumOfLengths) * sumOfQueryLengths;

        std::vector<SearchResult> allResults(numQueries);

        for(int q = 0; q < numQueries; q++){
            clearProcessedSequencesAndCachedResults(q);

            std::vector<ReferenceIdT> referencesToConsider(h_finalReferenceIdsPerQuery[q].begin(), h_finalReferenceIdsPerQuery[q].begin() + topNSize);
            cuda::std::span<ReferenceIdT> referenceIdsSpan(referencesToConsider.begin(), referencesToConsider.size());
            // std::vector<int> firstPassScores(h_finalAlignmentScoresPerQuery[q].begin(), h_finalAlignmentScoresPerQuery[q].begin() + topNSize);
            // std::cout << "first pass scores\n";
            // for(auto x : firstPassScores){
            //     std::cout << x << " ";
            // }
            // std::cout << "\n";
            // std::cout << "first pass refids\n";
            // for(auto x : referenceIdsSpan){
            //     std::cout << x << " ";
            // }
            // std::cout << "\n";

            //for debugging
            std::fill(h_finalReferenceIdsPerQuery[q].begin(), h_finalReferenceIdsPerQuery[q].begin() + topNSize, -1);

            scanTimer->reset();
            scanTimer->start();

            //for each score-only result, compute score along best diagonal and add it to gapless score
            gpuDatabase->forEachBatchOfSequencesInDatabaseSubset(
                std::bind(&ScoreOnlyFoldseekSearch::forEachCallback_subset, this, std::placeholders::_1, q), 
                referenceIdsSpan, 
                streamsSpan
            );

            finalizeSubScanResults(q, streamsSpan);

            scanTimer->stop();

            double secondPassSeconds = scanTimer->elapsed() / 1000.0;

            int resultTopNSize = std::min(int(referencesToConsider.size()), topNSize);

            SearchResult results;
            results.scores.insert(results.scores.end(), h_finalAlignmentScoresPerQuery[q].begin(), h_finalAlignmentScoresPerQuery[q].begin() + resultTopNSize);
            results.referenceIds.insert(results.referenceIds.end(), h_finalReferenceIdsPerQuery[q].begin(), h_finalReferenceIdsPerQuery[q].begin() + resultTopNSize);

            //results' referenceIds are in the range [0, referencesToConsider.size()-1]
            //convert them into the correct global ids
            for(auto& refId : results.referenceIds){
                assert(refId != -1);
                refId = referencesToConsider[refId];
            }

            results.stats = std::make_unique<BenchmarkStats>();
            results.stats->seconds = (firstPassSeconds + secondPassSeconds);
            results.stats->gcups = computedCellsFirstPass / 1000. / 1000. / 1000. / results.stats->seconds;

            allResults[q] = std::move(results);
        }

        return allResults;
    }


    void ScoreOnlyFoldseekSearch::forEachCallback_fullscan(cuda::std::span<GpuDatabase::CallbackInputPerGpu> callbackInput){
        assert(workingSets.size() == deviceIds.size());

        const int numGpus = deviceIds.size();
        const int numInputGpus = callbackInput.size();
        assert(numGpus == numInputGpus);
        for(int gpu = 0; gpu < numGpus; gpu++){
            assert(deviceIds[gpu] == callbackInput[gpu].deviceId);
            if(cudaStreamPerThread != callbackInput[gpu].stream){
                auto& ws = *workingSets[gpu];
                CUDACHECK(cudaEventRecord(ws.forkStreamEvent, callbackInput[gpu].stream));
                CUDACHECK(cudaStreamWaitEvent(cudaStreamPerThread, ws.forkStreamEvent, 0));
            }
        }

        std::vector<char*> d_masked3di_sequences_vec(numGpus);

        {
            nvtx::ScopedRange sr("convert_to_3di_and_masking_full", 0);
            //convert to 3di and apply masking
            for(int gpu = 0; gpu < numGpus; gpu++){
                if(callbackInput[gpu].numInputSequences > 0){
                    CUDACHECK(cudaSetDevice(callbackInput[gpu].deviceId));
                    cudaStream_t stream = cudaStreamPerThread;
                    auto& ws = *workingSets[gpu];

                    char* d_masked3di_ptr;
                    CUDACHECK(helper_gpuMallocAsync((void**)&d_masked3di_ptr, sizeof(char) * getSequenceCharsAllocationSize(callbackInput[gpu].numInputChars), stream));

                    thrust::device_vector<size_t, thrust_async_allocator<size_t>> d_plainOffsets{thrust_async_allocator<size_t>(stream)};
                    d_plainOffsets.resize(callbackInput[gpu].numInputSequences, thrust::no_init);

                    thrust::copy(
                        thrust::cuda::par_nosync.on(stream),
                        callbackInput[gpu].d_inputOffsets,
                        callbackInput[gpu].d_inputOffsets + callbackInput[gpu].numInputSequences,
                        d_plainOffsets.begin()
                    );

                    if(isOrdinaryAA20db){
                        mask3diSequences(
                            d_masked3di_ptr,
                            callbackInput[gpu].d_inputLengths,
                            //callbackInput[gpu].d_inputOffsets,
                            d_plainOffsets.data().get(),
                            callbackInput[gpu].d_inputChars, 
                            callbackInput[gpu].numInputSequences,
                            maskingOptions.maskingThreshold, //masking threshold
                            maskingOptions.maskingLetter, //maskingLetter,
                            true, // applyExtraPaddingLetter,
                            20, //padding letter
                            ConvertCombinedTo3diDummy{},
                            callbackInput[gpu].sequenceLengthUpperBound,
                            stream
                        );
                    }else{
                        mask3diSequences(
                            d_masked3di_ptr,
                            callbackInput[gpu].d_inputLengths,
                            //callbackInput[gpu].d_inputOffsets,
                            d_plainOffsets.data().get(),
                            callbackInput[gpu].d_inputChars, 
                            callbackInput[gpu].numInputSequences,
                            maskingOptions.maskingThreshold, //masking threshold
                            maskingOptions.maskingLetter, //maskingLetter,
                            true, // applyExtraPaddingLetter,
                            20, //padding letter
                            ConvertCombinedTo3di{},
                            callbackInput[gpu].sequenceLengthUpperBound,
                            stream
                        );
                    }
                    
                    d_masked3di_sequences_vec[gpu] = d_masked3di_ptr;                  
                }
            }

        }


        //determine maximum number of sequences to process over all gpus
        size_t maxNumSequencesInBatchForGpus = 0;
        for(int gpu = 0; gpu < numGpus; gpu++){
            maxNumSequencesInBatchForGpus = std::max(maxNumSequencesInBatchForGpus, callbackInput[gpu].numInputSequences);
        }
        // std::cout << maxNumSequencesInBatchForGpus << "\n";
        const size_t seqsPerPass = kernelOutputArraySize;
        const int numQueries = currentQueryLengths.size();
        for(size_t sequencePassOffset = 0; sequencePassOffset < maxNumSequencesInBatchForGpus; sequencePassOffset += seqsPerPass){


            // check if there is enough room in kernel output arrays. if not, combine them with topN output
            for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
                for(int gpu = 0; gpu < numGpus; gpu++){
                    if(sequencePassOffset < callbackInput[gpu].numInputSequences){
                        CUDACHECK(cudaSetDevice(callbackInput[gpu].deviceId));
                        cudaStream_t stream = cudaStreamPerThread;
                        auto& ws = *workingSets[gpu];
                        const int numInPass = std::min(callbackInput[gpu].numInputSequences - sequencePassOffset, seqsPerPass);
                        // std::cout << "numInPass " << numInPass << ",  ws.kernelOutputArrayUsedSize[queryIndex] " <<  ws.kernelOutputArrayUsedSize[queryIndex] << "\n";
                        if(numInPass + ws.kernelOutputArrayUsedSize[queryIndex] > kernelOutputArraySize){
                            consumeKernelOutputArraysFullScan(queryIndex, gpu, stream);
                        }
                    }
                }
            }


            //score only alignment of queries against DB
            for(int gpu = 0; gpu < numGpus; gpu++){                
                if(sequencePassOffset < callbackInput[gpu].numInputSequences){
                    CUDACHECK(cudaSetDevice(callbackInput[gpu].deviceId));
                    cudaStream_t stream = cudaStreamPerThread;

                    auto& ws = *workingSets[gpu];

                    PositionsIterator d_selectedPositions = PositionsIterator::fromCountingIterator(sequencePassOffset);
                    const size_t numInPass = std::min(callbackInput[gpu].numInputSequences - sequencePassOffset, seqsPerPass);


                    for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
                        nvtx::ScopedRange sr("align_query_score_only", 0);
                        OneToAllInputDataPSSM inputData;
                        inputData.d_subjects = reinterpret_cast<const int8_t*>(d_masked3di_sequences_vec[gpu]);
                        inputData.d_subjectOffsets = callbackInput[gpu].d_inputOffsets;
                        inputData.d_subjectLengths = callbackInput[gpu].d_inputLengths;
                        inputData.queryLength = currentQueryLengths[queryIndex];
                        inputData.numAlignments = numInPass;
                        inputData.maximumSubjectLength = callbackInput[gpu].sequenceLengthUpperBound;
                        inputData.indexIndirection = d_selectedPositions;

                        if(standardAlignerPerGpu[gpu]->isSingleTile(inputData)){
                            standardAlignerPerGpu[gpu]->scoreOnly_singleTile(
                                ws.getKernelOutputArrayScores(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex],
                                inputData,
                                *queryPssmPerGpuPerQuery_forStandardScoreOnly[queryIndex][gpu],
                                GapScoreArgs{.gapopenscore = gop, .gapextendscore = gex},
                                stream
                            );
                        }else{
                            thrust::device_vector<char, thrust_async_allocator<char>> d_tempStorage
                                = allocateTempBuffer(8ull << 30, 256ull << 20, stream);

                            standardAlignerPerGpu[gpu]->scoreOnly_multiTile(
                                d_tempStorage.data().get(),
                                d_tempStorage.size(),
                                ws.getKernelOutputArrayScores(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex],
                                inputData,
                                *queryPssmPerGpuPerQuery_forStandardScoreOnly[queryIndex][gpu],
                                GapScoreArgs{.gapopenscore = gop, .gapextendscore = gex},
                                stream
                            );
                        }

                        thrust::sequence(
                            thrust::cuda::par_nosync.on(stream),
                            ws.getKernelOutputArrayIndices(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex], 
                            ws.getKernelOutputArrayIndices(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex] + numInPass, 
                            ReferenceIdT(processedSequencesPerGpu[gpu])
                        );
                    }
                }
            }

            //handle alignment score overflows
            if(scoreOnlyOverflowSettings.checkForOverflow){
                
                // std::vector<std::unique_ptr<thrust::device_vector<char, thrust_async_allocator<char>>>> d_cub_temp_vec(numGpus);
                // for(int gpu = 0; gpu < numGpus; gpu++){
                //     cudaStream_t stream = cudaStreamPerThread;
                //     CUDACHECK(cudaSetDevice(callbackInput[gpu].deviceId));
                //     d_cub_temp_vec[gpu] = std::make_unique<thrust::device_vector<char, thrust_async_allocator<char>>>(thrust_async_allocator<char>(stream));
                // }
                    
                { 
                    nvtx::ScopedRange srov("overflow_correction_count_overflows", 2);
                    std::fill_n(h_numOverflowsPerGpuPerQuery.data(), numGpus * numQueries, 0);

                    for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
                        for(int gpu = 0; gpu < numGpus; gpu++){
                            if(sequencePassOffset < callbackInput[gpu].numInputSequences){
                                auto scoreDatatype = standardAlignerPerGpu[gpu]->getScoreDatatype(currentQueryLengths[queryIndex]);
                                const int overflowThreshold = getOverflowThreshold(queryIndex, scoreDatatype);
                                auto overflowDatatype = scoreOnlyOverflowSettings.overflowAlignersScoreOnly[gpu]->getScoreDatatype(currentQueryLengths[queryIndex]);
                                const int newThreshold = getOverflowThreshold(queryIndex, overflowDatatype);
                                
                                if(newThreshold > overflowThreshold){
                                    CUDACHECK(cudaSetDevice(callbackInput[gpu].deviceId));
                                    cudaStream_t stream = cudaStreamPerThread;

                                    auto& ws = *workingSets[gpu];

                                    const size_t numInPass = std::min(callbackInput[gpu].numInputSequences - sequencePassOffset, seqsPerPass);

                                    size_t cubBytes = 0;
                                    CUDACHECK(cub::DeviceSelect::If(
                                        nullptr, 
                                        cubBytes, 
                                        ws.getKernelOutputArrayScores(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex], 
                                        thrust::make_discard_iterator(),
                                        h_numOverflowsPerGpuPerQuery.data() + gpu * numQueries + queryIndex, 
                                        numInPass,
                                        IsOverflowScore{overflowThreshold},
                                        stream
                                    ));
                                    thrust::device_vector<char, thrust_async_allocator<char>> d_tempStorage = allocateTempBuffer(cubBytes, stream);
                                    CUDACHECK(cub::DeviceSelect::If(
                                        d_tempStorage.data().get(), 
                                        cubBytes, 
                                        ws.getKernelOutputArrayScores(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex], 
                                        thrust::make_discard_iterator(),
                                        h_numOverflowsPerGpuPerQuery.data() + gpu * numQueries + queryIndex, 
                                        numInPass,
                                        IsOverflowScore{overflowThreshold},
                                        stream
                                    ));
                                }
                            }
                        }
                    }
                    for(int gpu = 0; gpu < numGpus; gpu++){
                        cudaStream_t stream = cudaStreamPerThread;
                        CUDACHECK(cudaStreamSynchronize(stream))
                    }
                }

                {
                    nvtx::ScopedRange srov("overflow_correction_correct_overflows", 3);
                    for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
                        for(int gpu = 0; gpu < numGpus; gpu++){
                            if(sequencePassOffset < callbackInput[gpu].numInputSequences){
                                auto scoreDatatype = standardAlignerPerGpu[gpu]->getScoreDatatype(currentQueryLengths[queryIndex]);
                                const int overflowThreshold = getOverflowThreshold(queryIndex, scoreDatatype);
                                const int numOverflows = h_numOverflowsPerGpuPerQuery[gpu * numQueries + queryIndex];

                                if(numOverflows > 0){
                                    CUDACHECK(cudaSetDevice(callbackInput[gpu].deviceId));
                                    cudaStream_t stream = streamsPerGpuPerQuery[gpu * numQueries + queryIndex];

                                    auto& ws = *workingSets[gpu];

                                    const size_t numInPass = std::min(callbackInput[gpu].numInputSequences - sequencePassOffset, seqsPerPass);

                                    auto selectop = cuda::proclaim_return_type<bool>(IsIndexToOverflowScore{overflowThreshold, ws.getKernelOutputArrayScores(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex]});
                                    
                                    thrust::device_vector<int, thrust_async_allocator<int>> d_arrays{thrust_async_allocator<int>(stream)};
                                    const int numOverflowsRounded = cuda::ceil_div(numOverflows, 128) * 128;
                                    d_arrays.resize(3*numOverflowsRounded, thrust::no_init);
                                    int* const d_overflowIndices = d_arrays.data().get();
                                    int* const d_overflowIndices_global = d_overflowIndices + numOverflowsRounded;
                                    int* const d_correctedScores = d_overflowIndices_global + numOverflowsRounded;                                    

                                    auto selectIfOutputBuffer = thrust::make_transform_output_iterator(
                                        d_overflowIndices_global,
                                        cuda::proclaim_return_type<int>(
                                            [sequencePassOffset] __device__ (int i){
                                                return int(sequencePassOffset + i);
                                            }
                                        )
                                    );

                                    size_t cubBytes = 0;
                                    CUDACHECK(cub::DeviceSelect::If(
                                        nullptr, 
                                        cubBytes, 
                                        thrust::make_counting_iterator<int>(0), 
                                        // d_overflowIndices, 
                                        selectIfOutputBuffer,
                                        thrust::make_discard_iterator(),
                                        numInPass,
                                        selectop,
                                        stream
                                    ));
                                    thrust::device_vector<char, thrust_async_allocator<char>> d_cubTemp{thrust_async_allocator<char>(stream)};
                                    d_cubTemp.resize(cubBytes, thrust::no_init);
                                    CUDACHECK(cub::DeviceSelect::If(
                                        d_cubTemp.data().get(), 
                                        cubBytes, 
                                        thrust::make_counting_iterator<int>(0), 
                                        // d_overflowIndices, 
                                        selectIfOutputBuffer,
                                        thrust::make_discard_iterator(),
                                        numInPass,
                                        selectop,
                                        stream
                                    ));

                                    //the indices computed are local to the current pass, but for the alignment call we require the global indices within the batch
                                    // thrust::transform(
                                    //     thrust::cuda::par_nosync.on(stream),
                                    //     d_overflowIndices,
                                    //     d_overflowIndices + numOverflows,
                                    //     thrust::make_constant_iterator(sequencePassOffset),
                                    //     d_overflowIndices_global,
                                    //     cuda::std::plus<int>{}
                                    // );

                                    // cudaStream_t alignmentStream = streamsPerGpuPerQuery[gpu * numQueries + queryIndex];
                                    // CUDACHECK(cudaEventRecord(ws.forkStreamEvent, stream));
                                    // CUDACHECK(cudaStreamWaitEvent(alignmentStream, ws.forkStreamEvent, 0));

                                    PositionsIterator d_overflowPositionsIterator = PositionsIterator::fromPointer(d_overflowIndices_global);

                                    //sequence data is same as before, but we have to pass the indices to overflow alignments
                                    OneToAllInputDataPSSM inputData;
                                    inputData.d_subjects = reinterpret_cast<const int8_t*>(d_masked3di_sequences_vec[gpu]);
                                    inputData.d_subjectOffsets = callbackInput[gpu].d_inputOffsets;
                                    inputData.d_subjectLengths = callbackInput[gpu].d_inputLengths;
                                    inputData.queryLength = currentQueryLengths[queryIndex];
                                    inputData.numAlignments = numOverflows;
                                    inputData.maximumSubjectLength = callbackInput[gpu].sequenceLengthUpperBound;
                                    inputData.indexIndirection = d_overflowPositionsIterator;

                                    if(scoreOnlyOverflowSettings.overflowAlignersScoreOnly[gpu]->isSingleTile(inputData)){
                                        scoreOnlyOverflowSettings.overflowAlignersScoreOnly[gpu]->scoreOnly_singleTile(
                                            d_correctedScores,
                                            inputData,
                                            *queryPssmPerGpuPerQuery_forOverflowScoreOnly[queryIndex][gpu],
                                            GapScoreArgs{.gapopenscore = gop, .gapextendscore = gex},
                                            stream
                                        );
                                    }else{
                                        size_t suggestedBytes = scoreOnlyOverflowSettings.overflowAlignersScoreOnly[gpu]->getMinimumSuggestedTempBytes_multiTile(
                                            inputData
                                        );
                                        thrust::device_vector<char, thrust_async_allocator<char>> d_tempStorage
                                            = allocateTempBuffer(suggestedBytes, 256ull << 20, stream);
                                        scoreOnlyOverflowSettings.overflowAlignersScoreOnly[gpu]->scoreOnly_multiTile(
                                            d_tempStorage.data().get(),
                                            d_tempStorage.size(),
                                            d_correctedScores,
                                            inputData,
                                            *queryPssmPerGpuPerQuery_forOverflowScoreOnly[queryIndex][gpu],
                                            GapScoreArgs{.gapopenscore = gop, .gapextendscore = gex},
                                            stream
                                        );
                                    }

                                    thrust::scatter(
                                        thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                                        d_correctedScores,
                                        d_correctedScores + numOverflows,
                                        //d_overflowIndices,
                                        thrust::make_transform_iterator(d_overflowIndices_global,
                                            cuda::proclaim_return_type<int>(
                                                [sequencePassOffset] __device__ (int i){
                                                    return int(i - sequencePassOffset);
                                                }
                                            )
                                        ),
                                        ws.getKernelOutputArrayScores(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex]
                                    );
                                }
                            }
                        }
                    }

                    for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
                        for(int gpu = 0; gpu < numGpus; gpu++){
                            CUDACHECK(cudaSetDevice(deviceIds[gpu]));
                            auto& ws = *workingSets[gpu];
                            cudaStream_t stream = cudaStreamPerThread;
                            cudaStream_t alignmentStream = streamsPerGpuPerQuery[gpu * numQueries + queryIndex];
                            CUDACHECK(cudaEventRecord(ws.forkStreamEvent, alignmentStream));
                            CUDACHECK(cudaStreamWaitEvent(stream, ws.forkStreamEvent, 0));
                        }
                    }
                }

            }



            for(int gpu = 0; gpu < numGpus; gpu++){             
                auto& ws = *workingSets[gpu];   
                if(sequencePassOffset < callbackInput[gpu].numInputSequences){
                    const size_t numInPass = std::min(callbackInput[gpu].numInputSequences - sequencePassOffset, seqsPerPass);
                    for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
                        ws.kernelOutputArrayUsedSize[queryIndex] += numInPass;
                    }
                    processedSequencesPerGpu[gpu] += numInPass;
                }
            }


        } // while sequencePassOffset < maxNumSequencesInBatchForGpus

        //free the masked 3di sequences
        for(int gpu = 0; gpu < numGpus; gpu++){
            if(callbackInput[gpu].numInputSequences > 0){
                cudaStream_t stream = cudaStreamPerThread;
                CUDACHECK(cudaSetDevice(callbackInput[gpu].deviceId));
                CUDACHECK(helper_gpuFreeAsync((void*)d_masked3di_sequences_vec[gpu], stream));
            }
        }
    }

    void ScoreOnlyFoldseekSearch::consumeKernelOutputArraysFullScan(int queryIndex, int gpu, cudaStream_t stream){
        const int numGpus = deviceIds.size();
            
        auto& ws = *workingSets[gpu];
        if(ws.kernelOutputArrayUsedSize[queryIndex] > 0){
            CUDACHECK(cudaSetDevice(deviceIds[gpu]));

            //db sequences are processed in ascending order. stable sort ensures that sequences with same score are sorted by ascending id without a custom comparator
            thrust::stable_sort_by_key(
                thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                ws.getKernelOutputArrayScores(queryIndex),
                ws.getKernelOutputArrayScores(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex],
                ws.getKernelOutputArrayIndices(queryIndex),
                thrust::greater<int>()
            );

            // {
            //     std::cout << "consumeKernelOutputArraysFullScan queryIndex " << queryIndex << "\n";
            //     std::cout << "input\n";
            //     std::vector<int> hscores(ws.kernelOutputArrayUsedSize[queryIndex]);
            //     CUDACHECK(cudaMemcpy(hscores.data(), ws.getKernelOutputArrayScores(queryIndex), sizeof(int) * hscores.size(), cudaMemcpyDeviceToHost));
            //     std::vector<ReferenceIdT> hids(ws.kernelOutputArrayUsedSize[queryIndex]);
            //     CUDACHECK(cudaMemcpy(hids.data(), ws.getKernelOutputArrayIndices(queryIndex), sizeof(int) * hids.size(), cudaMemcpyDeviceToHost));
            //     for(int i = 0; i < std::min(int(hscores.size()), 10); i++){
            //         std::cout << hscores[i] << " " << hids[i] << "\n";
            //     }
            // }

            if(resultSizePerQueryPerGpu[queryIndex][gpu] > 0){
                //there are already topN results, combine with kernel outputs
                auto mergeInput1 = thrust::make_zip_iterator(
                    ws.getKernelOutputArrayScores(queryIndex),
                    ws.getKernelOutputArrayIndices(queryIndex)
                );
                const int input1Size = std::min(ws.kernelOutputArrayUsedSize[queryIndex], topNSize);

                auto mergeInput2 = thrust::make_zip_iterator(
                    ws.getTopN_scores(queryIndex).data(), 
                    ws.getTopN_refIds(queryIndex).data()
                );
                const int input2Size = std::min(resultSizePerQueryPerGpu[queryIndex][gpu], topNSize);

                auto mergeOutput = thrust::make_zip_iterator(
                    ws.getTopN_scores_tmp(queryIndex).data(), 
                    ws.getTopN_refIds_tmp(queryIndex).data()
                );
                thrust::merge(
                    thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                    mergeInput1,
                    mergeInput1 + input1Size,
                    mergeInput2,
                    mergeInput2 + input2Size,
                    mergeOutput,
                    CompareScoresDescendingRefIdsAscending{}
                );
                resultSizePerQueryPerGpu[queryIndex][gpu] = std::min(input1Size + input2Size, topNSize);

                std::swap(ws.getTopN_scores(queryIndex), ws.getTopN_scores_tmp(queryIndex));
                std::swap(ws.getTopN_refIds(queryIndex), ws.getTopN_refIds_tmp(queryIndex));
            }else{
                //it is the first set of topN results, simple copy
                const int numResults = std::min(ws.kernelOutputArrayUsedSize[queryIndex], topNSize);
                CUDACHECK(cudaMemcpyAsync(
                    ws.getTopN_scores(queryIndex).data(),
                    ws.getKernelOutputArrayScores(queryIndex), 
                    sizeof(int) * numResults,
                    cudaMemcpyDeviceToDevice,
                    stream
                ));
                CUDACHECK(cudaMemcpyAsync(
                    ws.getTopN_refIds(queryIndex).data(),
                    ws.getKernelOutputArrayIndices(queryIndex), 
                    sizeof(ReferenceIdT) * numResults,
                    cudaMemcpyDeviceToDevice,
                    stream
                ));
                resultSizePerQueryPerGpu[queryIndex][gpu] = numResults;
            }

            ws.kernelOutputArrayUsedSize[queryIndex] = 0;
        }
    }

    void ScoreOnlyFoldseekSearch::finalizeFullScanResults(cuda::std::span<cudaStream_t> streams){
        const int numGpus = deviceIds.size();
        const int numQueries = currentQueryLengths.size();


        // combine any leftover kernel outputs with topN output
        for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
            for(int gpu = 0; gpu < numGpus; gpu++){
                consumeKernelOutputArraysFullScan(queryIndex, gpu, streams[gpu]);
            }
        }

        //copy final topN results to host

        //Note: after full scan, a subset search is executed which requires reference ids to be sorted in ascending order
        //Thus, host results will not be sorted by descending score, but by ascending reference id

        if(numGpus == 1){
            //straight-forward copy
            CUDACHECK(cudaSetDevice(deviceIds[0]));
            auto& ws = *workingSets[0];
            for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
                //sort tops by reference id
                thrust::sort_by_key(
                    thrust::cuda::par_nosync(thrust_async_allocator<char>(streams[0])).on(streams[0]),
                    ws.getTopN_refIds(queryIndex).data(),
                    ws.getTopN_refIds(queryIndex).data() + topNSize,
                    ws.getTopN_scores(queryIndex).data()
                );
                CUDACHECK(cudaMemcpyAsync(
                    h_finalAlignmentScoresPerQuery[queryIndex].data(), 
                    ws.getTopN_scores(queryIndex).data(),
                    sizeof(int) * topNSize, 
                    cudaMemcpyDeviceToHost, 
                    streams[0]
                ));
                CUDACHECK(cudaMemcpyAsync(
                    h_finalReferenceIdsPerQuery[queryIndex].data(), 
                    ws.getTopN_refIds(queryIndex).data(),
                    sizeof(ReferenceIdT) * topNSize, 
                    cudaMemcpyDeviceToHost, 
                    streams[0]
                ));
            }
            CUDACHECK(cudaStreamSynchronize(streams[0]));
        }else{
            CUDACHECK(cudaSetDevice(deviceIds[0]));
            
            std::array<size_t, 2> bytes;
            bytes[0] = SDIV(sizeof(int) * numGpus * topNSize, 256) * 256;
            bytes[1] = SDIV(sizeof(ReferenceIdT) * numGpus * topNSize, 256) * 256;
            std::array<size_t, 2> bytesInclPS;
            std::inclusive_scan(bytes.begin(), bytes.end(), bytesInclPS.begin());
            size_t maxBytesForAllResultsPerQuery = bytesInclPS.back();

            //standard async mempool allocation is only accessible by current device so it cannot be used for peer copies, use plain allocation instead
            //(we could create a different pool for this purpose)
            thrust::device_vector<char> d_tempStorage;
            d_tempStorage.resize(maxBytesForAllResultsPerQuery, thrust::no_init);

            int* const d_finalAlignmentScores_allGpus = reinterpret_cast<int*>(d_tempStorage.data().get());
            ReferenceIdT* const d_finalReferenceIds_allGpus = reinterpret_cast<ReferenceIdT*>(d_tempStorage.data().get() + bytesInclPS[0]);


            //combine per-gpu results, then copy final topN results to host
            for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
                CUDACHECK(cudaSetDevice(deviceIds[0]));
                CUDACHECK(cudaEventRecord(workingSets[0]->forkStreamEvent, streams[0]));
                for(int gpu = 1; gpu < numGpus; gpu++){
                    CUDACHECK(cudaSetDevice(deviceIds[gpu]));
                    CUDACHECK(cudaStreamWaitEvent(streams[gpu], workingSets[0]->forkStreamEvent));
                }

                std::vector<int> gatherDstOffsets(numGpus);
                std::exclusive_scan(
                    resultSizePerQueryPerGpu[queryIndex].begin(), resultSizePerQueryPerGpu[queryIndex].end(), 
                    gatherDstOffsets.begin(),
                    0
                );
                const int totalNumResultsAllGpus = gatherDstOffsets.back() + resultSizePerQueryPerGpu[queryIndex].back();



                for(int gpu = 0; gpu < numGpus; gpu++){
                    const int numResults = resultSizePerQueryPerGpu[queryIndex][gpu];
                    if(numResults > topNSize){
                        throw std::runtime_error("error numResults A");
                    }
                    if(numResults > 0){
                        CUDACHECK(cudaSetDevice(deviceIds[gpu]));
                        auto& ws = *workingSets[gpu];            

                        gpuDatabase->convertLocalIndicesToGlobalIndices(gpu, ws.getTopN_refIds(queryIndex).data(), numResults, streams[gpu]);

                        CUDACHECK(cudaMemcpyAsync(
                            d_finalAlignmentScores_allGpus + gatherDstOffsets[gpu],
                            ws.getTopN_scores(queryIndex).data(),
                            sizeof(int) * numResults,
                            cudaMemcpyDeviceToDevice,
                            streams[gpu]
                        ));
                        CUDACHECK(cudaMemcpyAsync(
                            d_finalReferenceIds_allGpus + gatherDstOffsets[gpu],
                            ws.getTopN_refIds(queryIndex).data(),
                            sizeof(ReferenceIdT) * numResults,
                            cudaMemcpyDeviceToDevice,
                            streams[gpu]
                        ));            

                        CUDACHECK(cudaEventRecord(ws.forkStreamEvent, streams[gpu]));

                        CUDACHECK(cudaSetDevice(deviceIds[0]));
                        CUDACHECK(cudaStreamWaitEvent(streams[0], ws.forkStreamEvent, 0));
                    }
                }

                CUDACHECK(cudaSetDevice(deviceIds[0]));

                //sort per-gpu top results to find overall top results
                auto sortInput = thrust::make_zip_iterator(
                    d_finalAlignmentScores_allGpus,
                    d_finalReferenceIds_allGpus
                );
                thrust::sort(
                    thrust::cuda::par_nosync(thrust_async_allocator<char>(streams[0])).on(streams[0]),
                    sortInput,
                    sortInput + totalNumResultsAllGpus,
                    CompareScoresDescendingRefIdsAscending{}
                );

                //sort tops by reference id for subsequent subset search pass
                thrust::sort_by_key(
                    thrust::cuda::par_nosync(thrust_async_allocator<char>(streams[0])).on(streams[0]),
                    d_finalReferenceIds_allGpus,
                    d_finalReferenceIds_allGpus + topNSize,
                    d_finalAlignmentScores_allGpus
                );

                CUDACHECK(cudaMemcpyAsync(
                    h_finalAlignmentScoresPerQuery[queryIndex].data(), 
                    d_finalAlignmentScores_allGpus, 
                    sizeof(int) * topNSize, 
                    cudaMemcpyDeviceToHost, 
                    streams[0]
                ));
                CUDACHECK(cudaMemcpyAsync(
                    h_finalReferenceIdsPerQuery[queryIndex].data(), 
                    d_finalReferenceIds_allGpus, 
                    sizeof(ReferenceIdT) * topNSize, 
                    cudaMemcpyDeviceToHost, 
                    streams[0]
                ));
            }
            CUDACHECK(cudaSetDevice(deviceIds[0]));
            CUDACHECK(cudaStreamSynchronize(streams[0]));
        }


    }

    void ScoreOnlyFoldseekSearch::consumeKernelOutputArraysSubScan(int queryIndex, int gpu, cudaStream_t stream){
        const int numGpus = deviceIds.size();
            
        auto& ws = *workingSets[gpu];
        // std::cout << "consumeKernelOutputArraysSubScan, kernel output size " << ws.kernelOutputArrayUsedSize[queryIndex] << ", processedSequencesPerGpu[gpu] " << processedSequencesPerGpu[gpu] << "\n";
        if(ws.kernelOutputArrayUsedSize[queryIndex] > 0){
            CUDACHECK(cudaSetDevice(deviceIds[gpu]));

            //db sequences are processed in ascending order. stable sort ensures that sequences with same score are sorted by ascending id without a custom comparator
            thrust::stable_sort_by_key(
                thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                ws.getKernelOutputArrayScores(queryIndex),
                ws.getKernelOutputArrayScores(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex],
                ws.getKernelOutputArrayIndices(queryIndex),
                thrust::greater<int>()
            );

            // {
            //     std::cout << "consumeKernelOutputArraysSubScan\n";
            //     std::cout << "input\n";
            //     std::vector<int> hscores(ws.kernelOutputArrayUsedSize[queryIndex]);
            //     CUDACHECK(cudaMemcpy(hscores.data(), ws.getKernelOutputArrayScores(queryIndex), sizeof(int) * hscores.size(), cudaMemcpyDeviceToHost));
            //     std::vector<ReferenceIdT> hids(ws.kernelOutputArrayUsedSize[queryIndex]);
            //     CUDACHECK(cudaMemcpy(hids.data(), ws.getKernelOutputArrayIndices(queryIndex), sizeof(int) * hids.size(), cudaMemcpyDeviceToHost));
            //     for(int i = 0; i < std::min(int(hscores.size()), 10); i++){
            //         std::cout << hscores[i] << " " << hids[i] << "\n";
            //     }
            // }

            if(resultSizePerQueryPerGpu[queryIndex][gpu] > 0){
                //there are already topN results, combine with kernel outputs
                auto mergeInput1 = thrust::make_zip_iterator(
                    ws.getKernelOutputArrayScores(queryIndex),
                    ws.getKernelOutputArrayIndices(queryIndex)
                );
                const int input1Size = std::min(ws.kernelOutputArrayUsedSize[queryIndex], topNSize);

                auto mergeInput2 = thrust::make_zip_iterator(
                    ws.getTopN_scores(queryIndex).data(), 
                    ws.getTopN_refIds(queryIndex).data()
                );
                const int input2Size = std::min(resultSizePerQueryPerGpu[queryIndex][gpu], topNSize);

                auto mergeOutput = thrust::make_zip_iterator(
                    ws.getTopN_scores_tmp(queryIndex).data(), 
                    ws.getTopN_refIds_tmp(queryIndex).data()
                );
                thrust::merge(
                    thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                    mergeInput1,
                    mergeInput1 + input1Size,
                    mergeInput2,
                    mergeInput2 + input2Size,
                    mergeOutput,
                    CompareScoresDescendingRefIdsAscending{}
                );
                resultSizePerQueryPerGpu[queryIndex][gpu] = std::min(input1Size + input2Size, topNSize);

                std::swap(ws.getTopN_scores(queryIndex), ws.getTopN_scores_tmp(queryIndex));
                std::swap(ws.getTopN_refIds(queryIndex), ws.getTopN_refIds_tmp(queryIndex));
            }else{
                //it is the first set of topN results, simple copy
                const int numResults = std::min(ws.kernelOutputArrayUsedSize[queryIndex], topNSize);
                CUDACHECK(cudaMemcpyAsync(
                    ws.getTopN_scores(queryIndex).data(),
                    ws.getKernelOutputArrayScores(queryIndex), 
                    sizeof(int) * numResults,
                    cudaMemcpyDeviceToDevice,
                    stream
                ));
                CUDACHECK(cudaMemcpyAsync(
                    ws.getTopN_refIds(queryIndex).data(),
                    ws.getKernelOutputArrayIndices(queryIndex), 
                    sizeof(ReferenceIdT) * numResults,
                    cudaMemcpyDeviceToDevice,
                    stream
                ));
                resultSizePerQueryPerGpu[queryIndex][gpu] = numResults;
                // {
                //     std::vector<ReferenceIdT> vec(topNSize);
                //     CUDACHECK(cudaMemcpy(vec.data(), ws.getTopN_refIds(queryIndex).data(), sizeof(ReferenceIdT) * topNSize, cudaMemcpyDeviceToHost));
                //     std::cout << "topn refids new:\n";
                //     for(int i = 0; i < topNSize; i++){
                //         std::cout << vec[i] << "\n";
                //     }
                //     std::cout << "\n";
                // }
            }



            ws.kernelOutputArrayUsedSize[queryIndex] = 0;
        }
    }

    void ScoreOnlyFoldseekSearch::finalizeSubScanResults(int queryIndex, cuda::std::span<cudaStream_t> streams){
        const int numGpus = deviceIds.size();
        const int numQueries = currentQueryLengths.size();


        // combine any leftover kernel outputs with topN output
        for(int gpu = 0; gpu < numGpus; gpu++){
            consumeKernelOutputArraysSubScan(queryIndex, gpu, streams[gpu]);
        }

        //copy final topN results to host
        if(numGpus == 1){
            //straight-forward copy
            CUDACHECK(cudaSetDevice(deviceIds[0]));
            auto& ws = *workingSets[0];
            CUDACHECK(cudaMemcpyAsync(
                h_finalAlignmentScoresPerQuery[queryIndex].data(), 
                ws.getTopN_scores(queryIndex).data(),
                sizeof(int) * topNSize, 
                cudaMemcpyDeviceToHost, 
                streams[0]
            ));
            CUDACHECK(cudaMemcpyAsync(
                h_finalReferenceIdsPerQuery[queryIndex].data(), 
                ws.getTopN_refIds(queryIndex).data(),
                sizeof(ReferenceIdT) * topNSize, 
                cudaMemcpyDeviceToHost, 
                streams[0]
            ));
            CUDACHECK(cudaStreamSynchronize(streams[0]));
        }else{
            CUDACHECK(cudaSetDevice(deviceIds[0]));
            
            std::array<size_t, 2> bytes;
            bytes[0] = SDIV(sizeof(int) * numGpus * topNSize, 256) * 256;
            bytes[1] = SDIV(sizeof(ReferenceIdT) * numGpus * topNSize, 256) * 256;
            std::array<size_t, 2> bytesInclPS;
            std::inclusive_scan(bytes.begin(), bytes.end(), bytesInclPS.begin());
            size_t maxBytesForAllResultsPerQuery = bytesInclPS.back();

            //standard async mempool allocation is only accessible by current device so it cannot be used for peer copies, use plain allocation instead
            //(we could create a different pool for this purpose)
            thrust::device_vector<char> d_tempStorage;
            d_tempStorage.resize(maxBytesForAllResultsPerQuery, thrust::no_init);

            int* const d_finalAlignmentScores_allGpus = reinterpret_cast<int*>(d_tempStorage.data().get());
            ReferenceIdT* const d_finalReferenceIds_allGpus = reinterpret_cast<ReferenceIdT*>(d_tempStorage.data().get() + bytesInclPS[0]);


            //combine per-gpu results, then copy final topN results to host
            std::vector<int> gatherDstOffsets(numGpus);
            std::exclusive_scan(
                resultSizePerQueryPerGpu[queryIndex].begin(), resultSizePerQueryPerGpu[queryIndex].end(), 
                gatherDstOffsets.begin(),
                0
            );
            const int totalNumResultsAllGpus = gatherDstOffsets.back() + resultSizePerQueryPerGpu[queryIndex].back();

            for(int gpu = 0; gpu < numGpus; gpu++){
                const int numResults = resultSizePerQueryPerGpu[queryIndex][gpu];
                if(numResults > 0){
                    CUDACHECK(cudaSetDevice(deviceIds[gpu]));
                    auto& ws = *workingSets[gpu];            

                    gpuDatabase->convertLocalIndicesToGlobalIndices(gpu, ws.getTopN_refIds(queryIndex).data(), numResults, streams[gpu]);

                    CUDACHECK(cudaMemcpyAsync(
                        d_finalAlignmentScores_allGpus + gatherDstOffsets[gpu],
                        ws.getTopN_scores(queryIndex).data(),
                        sizeof(int) * numResults,
                        cudaMemcpyDeviceToDevice,
                        streams[gpu]
                    ));
                    CUDACHECK(cudaMemcpyAsync(
                        d_finalReferenceIds_allGpus + gatherDstOffsets[gpu],
                        ws.getTopN_refIds(queryIndex).data(),
                        sizeof(ReferenceIdT) * numResults,
                        cudaMemcpyDeviceToDevice,
                        streams[gpu]
                    ));            

                    CUDACHECK(cudaEventRecord(ws.forkStreamEvent, streams[gpu]));

                    CUDACHECK(cudaSetDevice(deviceIds[0]));
                    CUDACHECK(cudaStreamWaitEvent(streams[0], ws.forkStreamEvent, 0));
                }
            }

            CUDACHECK(cudaSetDevice(deviceIds[0]));

            //sort per-gpu top results to find overall top results
            auto sortInput = thrust::make_zip_iterator(
                d_finalAlignmentScores_allGpus,
                d_finalReferenceIds_allGpus
            );
            thrust::sort(
                thrust::cuda::par_nosync(thrust_async_allocator<char>(streams[0])).on(streams[0]),
                sortInput,
                sortInput + totalNumResultsAllGpus,
                CompareScoresDescendingRefIdsAscending{}
            );

            CUDACHECK(cudaMemcpyAsync(
                h_finalAlignmentScoresPerQuery[queryIndex].data(), 
                d_finalAlignmentScores_allGpus, 
                sizeof(int) * topNSize, 
                cudaMemcpyDeviceToHost, 
                streams[0]
            ));
            CUDACHECK(cudaMemcpyAsync(
                h_finalReferenceIdsPerQuery[queryIndex].data(), 
                d_finalReferenceIds_allGpus, 
                sizeof(ReferenceIdT) * topNSize, 
                cudaMemcpyDeviceToHost, 
                streams[0]
            ));

            CUDACHECK(cudaSetDevice(deviceIds[0]));
            CUDACHECK(cudaStreamSynchronize(streams[0]));
        }


    }


    void ScoreOnlyFoldseekSearch::forEachCallback_subset(cuda::std::span<GpuDatabase::CallbackInputPerGpu> callbackInput, int queryIndex){
        assert(workingSets.size() == deviceIds.size());

        const int numGpus = deviceIds.size();
        const int numInputGpus = callbackInput.size();
        assert(numGpus == numInputGpus);
        for(int gpu = 0; gpu < numGpus; gpu++){
            assert(deviceIds[gpu] == callbackInput[gpu].deviceId);
            if(cudaStreamPerThread != callbackInput[gpu].stream){
                auto& ws = *workingSets[gpu];
                CUDACHECK(cudaEventRecord(ws.forkStreamEvent, callbackInput[gpu].stream));
                CUDACHECK(cudaStreamWaitEvent(cudaStreamPerThread, ws.forkStreamEvent, 0));
            }
        }


        for(int gpu = 0; gpu < numGpus; gpu++){
            //subset search searches for topN sequences
            if(int(processedSequencesPerGpu[gpu] + callbackInput[gpu].numInputSequences) > topNSize){
                throw std::runtime_error("invalid subset search numInputSequences");
            }
        }

        std::vector<char*> d_masked3di_sequences_vec(numGpus);
        std::vector<char*> d_aa12_sequences_vec(numGpus);
        std::vector<size_t*> d_compactedSequenceOffsets_vec(numGpus);

        {
            nvtx::ScopedRange sr("convert_to_3di_and_masking_subset", 0);

            //convert to 3di and apply masking
            for(int gpu = 0; gpu < numGpus; gpu++){
                if(callbackInput[gpu].numInputSequences > 0){
                    CUDACHECK(cudaSetDevice(callbackInput[gpu].deviceId));
                    cudaStream_t stream = cudaStreamPerThread;
                    auto& ws = *workingSets[gpu];

                    char* d_masked3di_ptr;
                    CUDACHECK(helper_gpuMallocAsync((void**)&d_masked3di_ptr, sizeof(char) * getSequenceCharsAllocationSize(callbackInput[gpu].numInputChars), stream));

                    char* d_compactedSequences;
                    CUDACHECK(helper_gpuMallocAsync((void**)&d_compactedSequences, sizeof(char) * getSequenceCharsAllocationSize(callbackInput[gpu].numInputChars), stream));
                    
                    size_t* d_compactedSequenceOffsets;
                    CUDACHECK(helper_gpuMallocAsync((void**)&d_compactedSequenceOffsets, sizeof(size_t) * callbackInput[gpu].numInputSequences, stream));

                    // {
                    //     std::cout << "subset callbackInput[gpu].numInputChars " << callbackInput[gpu].numInputChars << "\n";
                    //     std::cout << "subset callbackInput[gpu].numInputSequences " << callbackInput[gpu].numInputSequences << "\n";

                    //     thrust::device_vector<size_t> dfoo(callbackInput[gpu].numInputSequences);
                    //     thrust::copy(thrust::device, callbackInput[gpu].d_inputLengths, callbackInput[gpu].d_inputLengths + callbackInput[gpu].numInputSequences, dfoo.begin());
                    //     std::cout << "subset input lengths\n";
                    //     for(int i = 0; i < int(dfoo.size()); i++){
                    //         std::cout << dfoo[i] << " ";
                    //     }
                    //     std::cout << "\n";
                    //     thrust::device_vector<size_t> dfoo2(callbackInput[gpu].numInputSequences);
                    //     thrust::copy(thrust::device, callbackInput[gpu].d_inputOffsets, callbackInput[gpu].d_inputOffsets + callbackInput[gpu].numInputSequences, dfoo2.begin());
                    //     std::cout << "subset input offsets\n";
                    //     for(int i = 0; i < int(dfoo2.size()); i++){
                    //         std::cout << dfoo2[i] << " ";
                    //     }
                    //     std::cout << "\n";
                    // }

                    //for subset search, input sequences could be stored non-contiguous in the buffer,
                    //compute offsets for contiguous storage
                    auto d_paddedLengths = thrust::make_transform_iterator(
                        callbackInput[gpu].d_inputLengths,
                        RoundToNextMultiple<size_t, 4>{}
                    );
                    thrust::exclusive_scan(
                        thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                        d_paddedLengths,
                        d_paddedLengths + callbackInput[gpu].numInputSequences,
                        d_compactedSequenceOffsets
                    );
                    d_compactedSequenceOffsets_vec[gpu] = d_compactedSequenceOffsets;

                    //copy sequences to contiguous offsets
                    {
                        auto inputSeqPointers = thrust::make_transform_iterator(
                            thrust::make_counting_iterator(0),
                            cuda::proclaim_return_type<char*>(
                                [baseptr = callbackInput[gpu].d_inputChars, offsets = callbackInput[gpu].d_inputOffsets] __device__ (int index){
                                    return const_cast<char*>(baseptr + offsets[index]);
                                }
                            )
                        );
                        auto outputSeqPointers = thrust::make_transform_iterator(
                            thrust::make_counting_iterator(0),
                            cuda::proclaim_return_type<char*>(
                                [d_compactedSequences, d_compactedSequenceOffsets] __device__ (int index){
                                    return d_compactedSequences + d_compactedSequenceOffsets[index];
                                }
                            )
                        );
                        
                        size_t cubBytes = 0;
                        CUDACHECK(cub::DeviceMemcpy::Batched(
                            nullptr,
                            cubBytes,
                            inputSeqPointers,
                            outputSeqPointers,
                            d_paddedLengths,
                            callbackInput[gpu].numInputSequences,
                            stream
                        ));
                        thrust::device_vector<char, thrust_async_allocator<char>> d_tempStorage = allocateTempBuffer(cubBytes, stream);
                        CUDACHECK(cub::DeviceMemcpy::Batched(
                            d_tempStorage.data().get(),
                            cubBytes,
                            inputSeqPointers,
                            outputSeqPointers,
                            d_paddedLengths,
                            callbackInput[gpu].numInputSequences,
                            stream
                        ));

                    }

                    if(isOrdinaryAA20db){
                        mask3diSequences(
                            d_masked3di_ptr,
                            callbackInput[gpu].d_inputLengths,
                            d_compactedSequenceOffsets,
                            d_compactedSequences, 
                            callbackInput[gpu].numInputSequences,
                            maskingOptions.maskingThreshold, //masking threshold
                            maskingOptions.maskingLetter, //maskingLetter,
                            true, // applyExtraPaddingLetter,
                            20, //padding letter
                            ConvertCombinedTo3diDummy{},
                            callbackInput[gpu].sequenceLengthUpperBound,
                            stream
                        );
                    }else{
                        mask3diSequences(
                            d_masked3di_ptr,
                            callbackInput[gpu].d_inputLengths,
                            d_compactedSequenceOffsets,
                            d_compactedSequences, 
                            callbackInput[gpu].numInputSequences,
                            maskingOptions.maskingThreshold, //masking threshold
                            maskingOptions.maskingLetter, //maskingLetter,
                            true, // applyExtraPaddingLetter,
                            20, //padding letter
                            ConvertCombinedTo3di{},
                            callbackInput[gpu].sequenceLengthUpperBound,
                            stream
                        );
                    }

                    d_masked3di_sequences_vec[gpu] = d_masked3di_ptr;

                    //transform combined to aa12 (inplace)
                    char* const d_aa12_ptr = d_compactedSequences;
                    thrust::transform(
                        thrust::cuda::par_nosync.on(stream),
                        d_compactedSequences,
                        d_compactedSequences + callbackInput[gpu].numInputChars,
                        d_aa12_ptr,
                        ConvertCombinedToAA12{}
                    );

                    d_aa12_sequences_vec[gpu] = d_aa12_ptr;
                }
            }
        }

        //determine maximum number of sequences to process over all gpus
        size_t maxNumSequencesInBatchForGpus = 0;
        for(int gpu = 0; gpu < numGpus; gpu++){
            maxNumSequencesInBatchForGpus = std::max(maxNumSequencesInBatchForGpus, callbackInput[gpu].numInputSequences);
        }
        // std::cout << maxNumSequencesInBatchForGpus << "\n";
        const size_t seqsPerPass = kernelOutputArraySize;
        for(size_t sequencePassOffset = 0; sequencePassOffset < maxNumSequencesInBatchForGpus; sequencePassOffset += seqsPerPass){

            // check if there is enough room in kernel output arrays. if not, combine them with topN output
            for(int gpu = 0; gpu < numGpus; gpu++){                
                if(sequencePassOffset < callbackInput[gpu].numInputSequences){
                    CUDACHECK(cudaSetDevice(callbackInput[gpu].deviceId));
                    cudaStream_t stream = cudaStreamPerThread;
                    auto& ws = *workingSets[gpu];
                    const int numInPass = std::min(callbackInput[gpu].numInputSequences - sequencePassOffset, seqsPerPass);
                    // std::cout << "numInPass " << numInPass << ",  ws.kernelOutputArrayUsedSize[queryIndex] " <<  ws.kernelOutputArrayUsedSize[queryIndex] << "\n";
                    if(numInPass + ws.kernelOutputArrayUsedSize[queryIndex] > kernelOutputArraySize){
                        consumeKernelOutputArraysSubScan(queryIndex, gpu, stream);
                    }
                }
            }


            for(int gpu = 0; gpu < numGpus; gpu++){                
                if(sequencePassOffset < callbackInput[gpu].numInputSequences){
                    CUDACHECK(cudaSetDevice(callbackInput[gpu].deviceId));
                    cudaStream_t stream = cudaStreamPerThread;

                    auto& ws = *workingSets[gpu];

                    PositionsIterator d_selectedPositions = PositionsIterator::fromCountingIterator(sequencePassOffset);
                    const size_t numInPass = std::min(callbackInput[gpu].numInputSequences - sequencePassOffset, seqsPerPass);


                    {
                        nvtx::ScopedRange sr("align_endpos", 0);

                        OneToAllInputDataPSSM inputData;
                        inputData.d_subjects = reinterpret_cast<const int8_t*>(d_masked3di_sequences_vec[gpu]);
                        inputData.d_subjectOffsets = makeCustomOffsetIterator(d_compactedSequenceOffsets_vec[gpu], false);
                        inputData.d_subjectLengths = callbackInput[gpu].d_inputLengths;
                        inputData.queryLength = currentQueryLengths[queryIndex];
                        inputData.numAlignments = numInPass;
                        inputData.maximumSubjectLength = callbackInput[gpu].sequenceLengthUpperBound;
                        inputData.indexIndirection = d_selectedPositions;

                        if(standardEndposAlignerPerGpu[gpu]->isSingleTile(inputData)){
                            standardEndposAlignerPerGpu[gpu]->scoreEndPos_singleTile(
                                ws.getKernelOutputArrayScores(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex],
                                ws.getKernelOutputArraySubjectEndPositionsExcl(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex],
                                ws.getKernelOutputArrayQueryEndPositionsExcl(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex],
                                inputData,
                                *queryPssmPerGpuPerQuery_forStandardEndPos[queryIndex][gpu],
                                GapScoreArgs{.gapopenscore = gop, .gapextendscore = gex},
                                stream
                            );
                        }else{
                            thrust::device_vector<char, thrust_async_allocator<char>> d_tempStorage
                                = allocateTempBuffer(8ull << 30, 256ull << 20, stream);
                            standardEndposAlignerPerGpu[gpu]->scoreEndPos_multiTile(
                                d_tempStorage.data().get(),
                                d_tempStorage.size(),
                                ws.getKernelOutputArrayScores(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex],
                                ws.getKernelOutputArraySubjectEndPositionsExcl(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex],
                                ws.getKernelOutputArrayQueryEndPositionsExcl(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex],
                                inputData,
                                *queryPssmPerGpuPerQuery_forStandardEndPos[queryIndex][gpu],
                                GapScoreArgs{.gapopenscore = gop, .gapextendscore = gex},
                                stream
                            );
                        }


                    }

                    if(endposOverflowSettings.checkForOverflow){
                        nvtx::ScopedRange sr("endpos_overflow_correction", 2);

                        auto scoreDatatype = standardEndposAlignerPerGpu[gpu]->getScoreDatatype(currentQueryLengths[queryIndex]);
                        const int overflowThreshold = getOverflowThreshold(queryIndex, scoreDatatype);
                        auto overflowDatatype = endposOverflowSettings.overflowAlignersScoreWithEndpos[gpu]->getScoreDatatype(currentQueryLengths[queryIndex]);
                        const int newThreshold = getOverflowThreshold(queryIndex, overflowDatatype);

                        if(newThreshold > overflowThreshold){
                            int* const d_scores = ws.getKernelOutputArrayScores(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex];
                            int* const d_subjectEndPos = ws.getKernelOutputArraySubjectEndPositionsExcl(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex];
                            int* const d_queryEndPos = ws.getKernelOutputArrayQueryEndPositionsExcl(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex];

                            thrust::device_vector<int, thrust_async_allocator<int>> d_numOverflows{thrust_async_allocator<int>(stream)};
                            d_numOverflows.resize(1, thrust::no_init);
                            {
                                size_t cubBytes = 0;
                                CUDACHECK(cub::DeviceSelect::If(nullptr, cubBytes, d_scores, thrust::make_discard_iterator(), d_numOverflows.data().get(), numInPass, IsOverflowScore{overflowThreshold}, stream));
                                thrust::device_vector<char, thrust_async_allocator<char>> d_cubTemp{thrust_async_allocator<char>(stream)};
                                d_cubTemp.resize(cubBytes, thrust::no_init);
                                CUDACHECK(cub::DeviceSelect::If(d_cubTemp.data().get(), cubBytes, d_scores, thrust::make_discard_iterator(), d_numOverflows.data().get(), numInPass, IsOverflowScore{overflowThreshold}, stream));
                            }
                            int numOverflows = 0;
                            CUDACHECK(cudaMemcpyAsync(&numOverflows, d_numOverflows.data().get(), sizeof(int), cudaMemcpyDeviceToHost, stream));
                            CUDACHECK(cudaStreamSynchronize(stream));

                            if(numOverflows > 0){
                                auto selectop = cuda::proclaim_return_type<bool>(IsIndexToOverflowScore{overflowThreshold, d_scores});

                                thrust::device_vector<int, thrust_async_allocator<int>> d_arrays{thrust_async_allocator<int>(stream)};
                                const int numOverflowsRounded = cuda::ceil_div(numOverflows, 128) * 128;
                                d_arrays.resize(4*numOverflowsRounded, thrust::no_init);
                                int* const d_overflowIndices_global = d_arrays.data().get();
                                int* const d_correctedScores = d_overflowIndices_global + numOverflowsRounded;
                                int* const d_correctedSubjectPos = d_correctedScores + numOverflowsRounded;
                                int* const d_correctedQueryPos = d_correctedSubjectPos + numOverflowsRounded;

                                auto selectIfOutputBuffer = thrust::make_transform_output_iterator(
                                    d_overflowIndices_global,
                                    cuda::proclaim_return_type<int>(
                                        [sequencePassOffset] __device__ (int i){
                                            return int(sequencePassOffset + i);
                                        }
                                    )
                                );

                                {
                                    size_t cubBytes = 0;
                                    CUDACHECK(cub::DeviceSelect::If(nullptr, cubBytes, thrust::make_counting_iterator<int>(0), selectIfOutputBuffer, thrust::make_discard_iterator(), numInPass, selectop, stream));
                                    thrust::device_vector<char, thrust_async_allocator<char>> d_cubTemp{thrust_async_allocator<char>(stream)};
                                    d_cubTemp.resize(cubBytes, thrust::no_init);
                                    CUDACHECK(cub::DeviceSelect::If(d_cubTemp.data().get(), cubBytes, thrust::make_counting_iterator<int>(0), selectIfOutputBuffer, thrust::make_discard_iterator(), numInPass, selectop, stream));
                                }

                                PositionsIterator d_overflowPositionsIterator = PositionsIterator::fromPointer(d_overflowIndices_global);

                                OneToAllInputDataPSSM overflowInputData;
                                overflowInputData.d_subjects = reinterpret_cast<const int8_t*>(d_masked3di_sequences_vec[gpu]);
                                overflowInputData.d_subjectOffsets = makeCustomOffsetIterator(d_compactedSequenceOffsets_vec[gpu], false);
                                overflowInputData.d_subjectLengths = callbackInput[gpu].d_inputLengths;
                                overflowInputData.queryLength = currentQueryLengths[queryIndex];
                                overflowInputData.numAlignments = numOverflows;
                                overflowInputData.maximumSubjectLength = callbackInput[gpu].sequenceLengthUpperBound;
                                overflowInputData.indexIndirection = d_overflowPositionsIterator;

                                if(endposOverflowSettings.overflowAlignersScoreWithEndpos[gpu]->isSingleTile(overflowInputData)){
                                    endposOverflowSettings.overflowAlignersScoreWithEndpos[gpu]->scoreEndPos_singleTile(
                                        d_correctedScores,
                                        d_correctedSubjectPos,
                                        d_correctedQueryPos,
                                        overflowInputData,
                                        *queryPssmPerGpuPerQuery_forOverflowEndPos[queryIndex][gpu],
                                        GapScoreArgs{.gapopenscore = gop, .gapextendscore = gex},
                                        stream
                                    );
                                }else{
                                    thrust::device_vector<char, thrust_async_allocator<char>> d_tempStorage
                                        = allocateTempBuffer(8ull << 30, 256ull << 20, stream);
                                    endposOverflowSettings.overflowAlignersScoreWithEndpos[gpu]->scoreEndPos_multiTile(
                                        d_tempStorage.data().get(),
                                        d_tempStorage.size(),
                                        d_correctedScores,
                                        d_correctedSubjectPos,
                                        d_correctedQueryPos,
                                        overflowInputData,
                                        *queryPssmPerGpuPerQuery_forOverflowEndPos[queryIndex][gpu],
                                        GapScoreArgs{.gapopenscore = gop, .gapextendscore = gex},
                                        stream
                                    );
                                }

                                auto scatterInput = thrust::make_zip_iterator(
                                    d_correctedScores,
                                    d_correctedSubjectPos,
                                    d_correctedQueryPos
                                );
                                auto scatterOutput = thrust::make_zip_iterator(
                                    d_scores,
                                    d_subjectEndPos,
                                    d_queryEndPos
                                );
                                thrust::scatter(
                                    thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                                    scatterInput,
                                    scatterInput + numOverflows,
                                    thrust::make_transform_iterator(d_overflowIndices_global,
                                        cuda::proclaim_return_type<int>(
                                            [sequencePassOffset] __device__ (int i){
                                                return int(i - sequencePassOffset);
                                            }
                                        )
                                    ),
                                    scatterOutput
                                );
                            }
                        }
                    }

                    {
                        nvtx::ScopedRange sr("diagonal_rescore", 0);

                        assert(queryPssmPerGpuPerQuery_aa12[queryIndex][gpu]->numRows == 13);
                        call_diagonal_rescore_simple_kernel(
                            ws.getKernelOutputArrayScores(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex], //result[i] will be added to output[i]
                            reinterpret_cast<const std::int8_t*>(d_aa12_sequences_vec[gpu]), // must be in [0,11]
                            callbackInput[gpu].d_inputLengths,
                            d_compactedSequenceOffsets_vec[gpu],
                            numInPass,
                            ws.getKernelOutputArraySubjectEndPositionsExcl(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex],
                            ws.getKernelOutputArrayQueryEndPositionsExcl(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex],
                            queryPssmPerGpuPerQuery_aa12[queryIndex][gpu]->data<const half*>(), //12+1 rows, queryLength columns. last row is filled with 0
                            currentQueryLengths[queryIndex],
                            stream
                        );
                    }

                    thrust::sequence(
                        thrust::cuda::par_nosync.on(stream),
                        ws.getKernelOutputArrayIndices(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex], 
                        ws.getKernelOutputArrayIndices(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex] + numInPass, 
                        ReferenceIdT(processedSequencesPerGpu[gpu])
                    );

                    ws.kernelOutputArrayUsedSize[queryIndex] += numInPass;
                    processedSequencesPerGpu[gpu] += numInPass;
                }
            }


        } // while sequencePassOffset < maxNumSequencesInBatchForGpus

        //free the masked 3di sequences
        for(int gpu = 0; gpu < numGpus; gpu++){
            if(callbackInput[gpu].numInputSequences > 0){
                cudaStream_t stream = cudaStreamPerThread;
                CUDACHECK(cudaSetDevice(callbackInput[gpu].deviceId));
                CUDACHECK(helper_gpuFreeAsync((void*)d_masked3di_sequences_vec[gpu], stream));
                CUDACHECK(helper_gpuFreeAsync((void*)d_aa12_sequences_vec[gpu], stream));
                CUDACHECK(helper_gpuFreeAsync((void*)d_compactedSequenceOffsets_vec[gpu], stream));
                
            }
        }
    }

    void ScoreOnlyFoldseekSearch::initWorkingSets(){
        const int numGpus = deviceIds.size();
        for(int i = 0; i < numGpus; i++){
            CUDACHECK(cudaSetDevice(deviceIds[i]));
            workingSets.emplace_back(
                std::make_unique<GpuWorkingSet>(
                    kernelOutputArraySize,
                    topNSize,
                    maxNumQueries
                )
            );
        }
    }

    void ScoreOnlyFoldseekSearch::setAligner(
        cuda::std::span<MultiConfigScoreOnlyAlignerInterface*> alignersScoreOnly,
        const OverflowSettings& scoreOnlyOverflowSettings_,
        cuda::std::span<MultiConfigScoreWithEndPosAlignerInterface*> alignersEndPos,
        const OverflowSettings& endposOverflowSettings_
    ){
        assert(alignersScoreOnly.size() == deviceIds.size());
        assert(alignersEndPos.size() == deviceIds.size());
        standardAlignerPerGpu = alignersScoreOnly;
        scoreOnlyOverflowSettings = scoreOnlyOverflowSettings_;
        standardEndposAlignerPerGpu = alignersEndPos;
        endposOverflowSettings = endposOverflowSettings_;
    }


    int ScoreOnlyFoldseekSearch::getOverflowThreshold(int queryIndex, Datatype datatype){
        switch(datatype){
            case Datatype::Float: return 16777216; // 2^24
            case Datatype::Int: return std::numeric_limits<int>::max() - std::abs(std::min(0, currentMinMaxValuesInPssms[queryIndex].first));
            case Datatype::Half2: return 2048; // 2^11
            case Datatype::Short2: return std::numeric_limits<short>::max() - std::abs(std::min(0, currentMinMaxValuesInPssms[queryIndex].first));
            case Datatype::UInt8x4: return std::numeric_limits<uint8_t>::max() - std::abs(std::min(0, currentMinMaxValuesInPssms[queryIndex].first));
            default: throw std::runtime_error("unhandled case");
        }   
    }

    void ScoreOnlyFoldseekSearch::setPssms_3di(cuda::std::span<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>> hostQueryPssmViews_3di_){
        destroyPssms_3di();

        hostQueryPssmViews_3di = hostQueryPssmViews_3di_;
        const int numQueries = hostQueryPssmViews_3di.size();
        if(numQueries > maxNumQueries){
            throw std::runtime_error("too many queries for setPssm\n");
        }
        currentQueryLengths.resize(numQueries);
        currentMinMaxValuesInPssms.resize(numQueries);

        // std::cout << "setPssms_3di: " << "numQueries = " << numQueries << "\n"; 


        const int numGpus = deviceIds.size();

        for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
            currentQueryLengths[queryIndex] = hostQueryPssmViews_3di[queryIndex].extent(1);

            int minValInPssm = std::numeric_limits<int>::max();
            int maxValInPssm = std::numeric_limits<int>::min();
            for(int r = 0; r < hostQueryPssmViews_3di[queryIndex].extent(0); r++){
                for(int c = 0; c < hostQueryPssmViews_3di[queryIndex].extent(1); c++){
                    int val = hostQueryPssmViews_3di[queryIndex](r,c);
                    minValInPssm = std::min(minValInPssm, val);
                    maxValInPssm = std::max(maxValInPssm, val);
                }
            }
            currentMinMaxValuesInPssms[queryIndex] = std::make_pair(minValInPssm, maxValInPssm);

            std::vector<std::unique_ptr<libmarv::GpuConvertedPSSM>> gpuStandardPssmScoreOnly_vec;
            std::vector<std::unique_ptr<libmarv::GpuConvertedPSSM>> gpuOverflowPssm_vec;
            std::vector<std::unique_ptr<libmarv::GpuConvertedPSSM>> gpuStandardPssmEndPos_vec;
            std::vector<std::unique_ptr<libmarv::GpuConvertedPSSM>> gpuOverflowPssmEndPos_vec;

            for(int gpu = 0; gpu < numGpus; gpu++){
                CUDACHECK(cudaSetDevice(deviceIds[gpu]));

                
    
                auto gpuStandardPssm = std::make_unique<libmarv::GpuConvertedPSSM>();    
                standardAlignerPerGpu[gpu]->makeGpuPssm(
                    *gpuStandardPssm,
                    hostQueryPssmViews_3di[queryIndex],
                    cudaStreamPerThread
                );    
                gpuStandardPssmScoreOnly_vec.push_back(std::move(gpuStandardPssm));

                if(scoreOnlyOverflowSettings.checkForOverflow){
                    assert(int(scoreOnlyOverflowSettings.overflowAlignersScoreOnly.size()) > gpu);

                    auto gpuOverflowPssm = std::make_unique<libmarv::GpuConvertedPSSM>();    
                    scoreOnlyOverflowSettings.overflowAlignersScoreOnly[gpu]->makeGpuPssm(
                        *gpuOverflowPssm,
                        hostQueryPssmViews_3di[queryIndex],
                        cudaStreamPerThread
                    );    
                    gpuOverflowPssm_vec.push_back(std::move(gpuOverflowPssm));
                }

                
                auto gpuStandardPssmEndpos = std::make_unique<libmarv::GpuConvertedPSSM>();
                standardEndposAlignerPerGpu[gpu]->makeGpuPssm(
                    *gpuStandardPssmEndpos,
                    hostQueryPssmViews_3di[queryIndex],
                    cudaStreamPerThread
                );
                gpuStandardPssmEndPos_vec.push_back(std::move(gpuStandardPssmEndpos));

                if(endposOverflowSettings.checkForOverflow){
                    assert(int(endposOverflowSettings.overflowAlignersScoreWithEndpos.size()) > gpu);

                    auto gpuOverflowPssmEndpos = std::make_unique<libmarv::GpuConvertedPSSM>();
                    endposOverflowSettings.overflowAlignersScoreWithEndpos[gpu]->makeGpuPssm(
                        *gpuOverflowPssmEndpos,
                        hostQueryPssmViews_3di[queryIndex],
                        cudaStreamPerThread
                    );
                    gpuOverflowPssmEndPos_vec.push_back(std::move(gpuOverflowPssmEndpos));
                }
            }
            queryPssmPerGpuPerQuery_forStandardScoreOnly.push_back(std::move(gpuStandardPssmScoreOnly_vec));
            if(scoreOnlyOverflowSettings.checkForOverflow){
                queryPssmPerGpuPerQuery_forOverflowScoreOnly.push_back(std::move(gpuOverflowPssm_vec));
            }
            queryPssmPerGpuPerQuery_forStandardEndPos.push_back(std::move(gpuStandardPssmEndPos_vec));
            if(endposOverflowSettings.checkForOverflow){
                queryPssmPerGpuPerQuery_forOverflowEndPos.push_back(std::move(gpuOverflowPssmEndPos_vec));
            }
        }
    }

    
    void ScoreOnlyFoldseekSearch::setPssms_aa12(cuda::std::span<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>> hostQueryPssmViews_aa12_){
        destroyPssms_aa12();

        hostQueryPssmViews_aa12 = hostQueryPssmViews_aa12_;
        const int numQueries = hostQueryPssmViews_aa12.size();
        if(numQueries > maxNumQueries){
            throw std::runtime_error("too many queries for setPssm\n");
        }
        assert(numQueries == currentQueryLengths.size());

        const int numGpus = deviceIds.size();

        for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
            auto hostPssmView = hostQueryPssmViews_aa12_[queryIndex];
            const int numRows = hostPssmView.extent(0);
            const int queryLength = hostPssmView.extent(1);
            assert(numRows == 12);
            assert(queryLength == currentQueryLengths[queryIndex]);
            

            std::vector<half> h_pssm(numRows * queryLength);
            for(int s = 0; s < numRows; s++){
                for(int p = 0; p < queryLength; p++){
                    h_pssm[s * queryLength + p] = hostPssmView(s,p);
                }
            }


            std::vector<std::unique_ptr<libmarv::GpuPlainPSSM>> gpuPssms_vec;
    
            for(int gpu = 0; gpu < numGpus; gpu++){
                CUDACHECK(cudaSetDevice(deviceIds[gpu]));
                auto gpuPssm = std::make_unique<libmarv::GpuPlainPSSM>();    
                gpuPssm->resize<half>(numRows+1, queryLength, cudaStreamPerThread);
                //set last row to zero
                CUDACHECK(cudaMemsetAsync(
                    gpuPssm->data<half*>() + (numRows) * queryLength, 
                    0,
                    sizeof(half) * queryLength, 
                    cudaStreamPerThread
                ));
                gpuPssms_vec.push_back(std::move(gpuPssm));
            }

            //broadcast h_pssm to all gpus
            for(int gpu = 0; gpu < numGpus; gpu++){
                CUDACHECK(cudaSetDevice(deviceIds[gpu]));            
                CUDACHECK(cudaMemcpyAsync(
                    gpuPssms_vec[gpu]->data<half*>(), 
                    h_pssm.data(), 
                    sizeof(half) * (numRows) * queryLength, 
                    cudaMemcpyHostToDevice, 
                    cudaStreamPerThread
                ));
                //last row was already set to zero
            }
            queryPssmPerGpuPerQuery_aa12.push_back(std::move(gpuPssms_vec));
        }
    }

    void ScoreOnlyFoldseekSearch::clearProcessedSequencesAndCachedResults(){
        const int numGpus = deviceIds.size();
        for(int gpu = 0; gpu < numGpus; gpu++){
            auto& ws = *workingSets[gpu];
            std::fill(ws.kernelOutputArrayUsedSize.begin(), ws.kernelOutputArrayUsedSize.end(), 0);
        }

        for(auto& vec : resultSizePerQueryPerGpu){
            std::fill(vec.begin(), vec.end(), 0);
        }
        std::fill(processedSequencesPerGpu.begin(), processedSequencesPerGpu.end(), 0);
    }

    void ScoreOnlyFoldseekSearch::clearProcessedSequencesAndCachedResults(int queryIndex){
        const int numGpus = deviceIds.size();
        for(int gpu = 0; gpu < numGpus; gpu++){
            auto& ws = *workingSets[gpu];
            ws.kernelOutputArrayUsedSize[queryIndex] = 0;
        }

        std::fill(resultSizePerQueryPerGpu[queryIndex].begin(), resultSizePerQueryPerGpu[queryIndex].end(), 0);

        std::fill(processedSequencesPerGpu.begin(), processedSequencesPerGpu.end(), 0);
    }


LIBMARV_NAMESPACE_WITH_NESTING_END

