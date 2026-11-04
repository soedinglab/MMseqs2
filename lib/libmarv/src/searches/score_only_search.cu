#include "score_only_search.cuh"
#include "topn_select.cuh"

#include <cuda/functional>
#include <numeric>
#include <cuda/cmath>
#include <cuda/std/span>
#include <cuda/std/mdspan>
#include <cuda/std/tuple>


#include <vector>
#include <memory>
#include <algorithm>

#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/iterator/transform_output_iterator.h>
#include <thrust/iterator/discard_iterator.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
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
#include "../cuda_errorcheck.cuh"
#include "../config.hpp"
#include "../pssm.cuh"
#include "../gpudatabase.cuh"
#include "search_common.cuh"

#include "../hpc_helpers/nvtx_markers.cuh"
#include "../hpc_helpers/simple_allocation.cuh"
#include "../hpc_helpers/cuda_raiiwrappers.cuh"


#include "../namespace.hpp"
LIBMARV_NAMESPACE_WITH_NESTING_BEGIN



    struct ScoreOnlySearch::GpuWorkingSet{
        int maxNumQueries;
        int kernelOutputArraySize;
        std::vector<int> kernelOutputArrayUsedSize;
        int topNSize;

        MyDeviceBuffer<int> d_kernelOutputArrayScores;
        MyDeviceBuffer<ReferenceIdT> d_kernelOutputArrayIndices;

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
            d_kernelOutputArrayScores.resize(size_t(kernelOutputArraySize) * maxNumQueries);
            d_kernelOutputArrayIndices.resize(size_t(kernelOutputArraySize) * maxNumQueries);

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


    ScoreOnlySearch::ScoreOnlySearch(GpuDatabase& gpuDatabase_, int kernelOutputArraySize_, int topNSize_, int maxNumQueries_) 
        // a small database cannot fill the pass buffer
        : kernelOutputArraySize(int(std::min(size_t(kernelOutputArraySize_),
                                            std::max(size_t(1), gpuDatabase_.getNumSequences())))),
        topNSize(int(std::min(size_t(topNSize_), gpuDatabase_.getNumSequences()))), 
        maxTopNSize(topNSize), maxNumQueries(maxNumQueries_),
        deviceIds(gpuDatabase_.getDeviceIds()), gpuDatabase(&gpuDatabase_)
    {
        initWorkingSets();

        const int numGpus = deviceIds.size();

        processedSequencesPerGpu.resize(deviceIds.size(), 0);
        h_finalAlignmentScoresPerQuery.resize(maxNumQueries);
        h_finalReferenceIdsPerQuery.resize(maxNumQueries);
        for(int q = 0; q < maxNumQueries; q++){
            h_finalAlignmentScoresPerQuery[q].resize(topNSize);
            h_finalReferenceIdsPerQuery[q].resize(topNSize);
        }

        resultSizePerQueryPerGpu.resize(maxNumQueries, std::vector<int>(deviceIds.size()));
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

    ScoreOnlySearch::~ScoreOnlySearch(){
        releaseQueryResources();
        destroyWorkingSets();
    }

    size_t ScoreOnlySearch::workspaceBytes() const {
        auto bytes = [](const auto& buffer) {
            return marvBackendAllocationSize(buffer.capacity() * sizeof(*buffer.data()));
        };
        size_t total = bytes(h_numOverflowsPerGpuPerQuery);
        for (const auto& b : h_finalAlignmentScoresPerQuery) total += bytes(b);
        for (const auto& b : h_finalReferenceIdsPerQuery) total += bytes(b);
        for (const auto& ws : workingSets) {
            total += bytes(ws->d_kernelOutputArrayScores) + bytes(ws->d_kernelOutputArrayIndices);
            for (const auto& b : ws->d_topN_scores_vec) total += bytes(b);
            for (const auto& b : ws->d_topN_refIds_vec) total += bytes(b);
            for (const auto& b : ws->d_topN_scores_tmp_vec) total += bytes(b);
            for (const auto& b : ws->d_topN_refIds_tmp_vec) total += bytes(b);
        }
        return total;
    }

    void ScoreOnlySearch::releaseQueryResources() {
        RevertDeviceId rd{};
        // drain the retry streams before profiles or buffers are released, also on exceptions
        if (!queryPssmPerGpuPerQuery_forStandard.empty() || !queryPssmPerGpuPerQuery_forOverflow.empty()) {
            for (size_t gpu = 0; gpu < deviceIds.size(); ++gpu) {
                CUDACHECK(cudaSetDevice(deviceIds[gpu]));
                for (int q = 0; q < maxNumQueries; ++q)
                    CUDACHECK(cudaStreamSynchronize(streamsPerGpuPerQuery[gpu * maxNumQueries + q]));
                CUDACHECK(cudaStreamSynchronize(cudaStreamPerThread));
            }
        }
        destroyPssms();
        hostQueryPssmViews = {};
        standardAlignerPerGpu = {};
        overflowSettings = OverflowSettings::withoutOverflowCorrection();
        currentQueryLengths.clear();
        currentMinMaxValuesInPssms.clear();
    }


    void ScoreOnlySearch::reset(){
        nvtx::ScopedRange sr("ScoreOnlySearch::reset", 0);

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



    void ScoreOnlySearch::destroyWorkingSets(){
        const int numGpus = deviceIds.size();
        for(int i = 0; i < numGpus; i++){
            CUDACHECK(cudaSetDevice(deviceIds[i]));
            workingSets[i].reset();
        }
    }

    void ScoreOnlySearch::destroyPssms(){
        for (size_t gpu = 0; gpu < deviceIds.size(); ++gpu) {
            CUDACHECK(cudaSetDevice(deviceIds[gpu]));
            for (auto& profiles : queryPssmPerGpuPerQuery_forStandard)
                if (gpu < profiles.size()) profiles[gpu].reset();
            for (auto& profiles : queryPssmPerGpuPerQuery_forOverflow)
                if (gpu < profiles.size()) profiles[gpu].reset();
        }
        queryPssmPerGpuPerQuery_forStandard.clear();
        queryPssmPerGpuPerQuery_forOverflow.clear();
    }


    /*
        Search in complete database
    */
    std::vector<SearchResult> ScoreOnlySearch::execute(
        cuda::std::span<MultiConfigScoreOnlyAlignerInterface*> aligners,
        const OverflowSettings& overflowSettings_,
        cuda::std::span<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>> queryPssmViews,
        int gapopenscore,
        int gapextendscore
    ){
        helpers::GpuTimer timer("Total");
        clearProcessedSequencesAndCachedResults();

        setAligner(aligners, overflowSettings_);
        setPssms(queryPssmViews);
        gop = gapopenscore;
        gex = gapextendscore;

        if(workingSets.size() == 0){
            initWorkingSets();
        }

        const int numQueries = queryPssmViews.size();

        std::vector<cudaStream_t> streams(deviceIds.size(), cudaStreamPerThread);
        cuda::std::span<cudaStream_t> streamsSpan(streams.data(), streams.size());

        scanTimer->reset();
        scanTimer->start();

        // TODO: for multi-tile queries, it could improve performance to process in separate batches 
        // because it would enable more efficient temp storage usage
        constexpr bool processCachedDBInSeparateBatches = false;
        gpuDatabase->forEachBatchOfSequencesInDatabase(
            std::bind(&ScoreOnlySearch::forEachCallback, this, std::placeholders::_1), 
            streamsSpan, 
            processCachedDBInSeparateBatches
        );

        finalizeResults(streamsSpan);

        scanTimer->stop();

        const int sumOfQueryLengths = std::accumulate(currentQueryLengths.begin(), currentQueryLengths.end(), 0);
        
        std::vector<SearchResult> allResults(numQueries);
        for(int q = 0; q < numQueries; q++){
            SearchResult results;
            results.scores.insert(results.scores.end(), h_finalAlignmentScoresPerQuery[q].begin(), h_finalAlignmentScoresPerQuery[q].begin() + topNSize);
            results.referenceIds.insert(results.referenceIds.end(), h_finalReferenceIdsPerQuery[q].begin(), h_finalReferenceIdsPerQuery[q].begin() + topNSize);
            
            const auto& sequenceLengthStatistics = gpuDatabase->getSequenceLengthStatistics();
            size_t computedCells = sequenceLengthStatistics.sumOfLengths * sumOfQueryLengths;
            results.stats = std::make_unique<BenchmarkStats>();
            results.stats->seconds = scanTimer->elapsed() / 1000;
            results.stats->gcups = computedCells / 1000. / 1000. / 1000. / results.stats->seconds;
            
            allResults[q] = std::move(results);
        }

        // for(int i = 0; i < 3; i++){
        //     std::cout << h_finalAlignmentScoresPerQuery[i] << " " << "referenceId " << h_finalReferenceIdsPerQuery[i] << "\n";
        // }
        timer.print();
        return allResults;
    }


    /*
        Search only in subset of database given by reference ids
    */
    std::vector<SearchResult> ScoreOnlySearch::execute(
        cuda::std::span<MultiConfigScoreOnlyAlignerInterface*> aligners,
        const OverflowSettings& overflowSettings_,
        cuda::std::span<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>> queryPssmViews,
        std::vector<ReferenceIdT> referencesToConsider,
        int gapopenscore,
        int gapextendscore
    ){
        helpers::GpuTimer timer("Total");
        clearProcessedSequencesAndCachedResults();

        setAligner(aligners, overflowSettings_);
        setPssms(queryPssmViews);
        gop = gapopenscore;
        gex = gapextendscore;

        if(workingSets.size() == 0){
            initWorkingSets();
        }

        std::sort(referencesToConsider.begin(), referencesToConsider.end());

        const int numQueries = queryPssmViews.size();

        std::vector<cudaStream_t> streams(deviceIds.size(), cudaStreamPerThread);
        cuda::std::span<cudaStream_t> streamsSpan(streams.data(), streams.size());
        cuda::std::span<ReferenceIdT> referenceIdsSpan(referencesToConsider.data(), referencesToConsider.size());

        scanTimer->reset();
        scanTimer->start();

        gpuDatabase->forEachBatchOfSequencesInDatabaseSubset(
            std::bind(&ScoreOnlySearch::forEachCallback, this, std::placeholders::_1), 
            referenceIdsSpan,
            streamsSpan
        );

        finalizeResults(streamsSpan);

        scanTimer->stop();

        const int sumOfQueryLengths = std::accumulate(currentQueryLengths.begin(), currentQueryLengths.end(), 0);

        int resultTopNSize = std::min(int(referencesToConsider.size()), topNSize);
        
        std::vector<SearchResult> allResults(numQueries);
        for(int q = 0; q < numQueries; q++){
            SearchResult results;
            results.scores.insert(results.scores.end(), h_finalAlignmentScoresPerQuery[q].begin(), h_finalAlignmentScoresPerQuery[q].begin() + resultTopNSize);
            results.referenceIds.insert(results.referenceIds.end(), h_finalReferenceIdsPerQuery[q].begin(), h_finalReferenceIdsPerQuery[q].begin() + resultTopNSize);

            //results' referenceIds are in the range [0, referencesToConsider.size()-1]
            //convert them into the correct global ids
            for(auto& refId : results.referenceIds){
                refId = referencesToConsider[refId];
            }
            
            const auto& sequenceLengthStatistics = gpuDatabase->getSequenceLengthStatistics();
            size_t computedCells = sequenceLengthStatistics.sumOfLengths * sumOfQueryLengths;
            results.stats = std::make_unique<BenchmarkStats>();
            results.stats->seconds = scanTimer->elapsed() / 1000;
            results.stats->gcups = computedCells / 1000. / 1000. / 1000. / results.stats->seconds;
            
            allResults[q] = std::move(results);
        }

        // for(int i = 0; i < 3; i++){
        //     std::cout << h_finalAlignmentScoresPerQuery[i] << " " << "referenceId " << h_finalReferenceIdsPerQuery[i] << "\n";
        // }
        timer.print();
        return allResults;
    }


    void ScoreOnlySearch::forEachCallback(cuda::std::span<GpuDatabase::CallbackInputPerGpu> callbackInput){
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

        #if 0
        {
            static int callbackCount = 0;

            if(callbackCount == 0){
                std::vector<char> h_inputChars(callbackInput[0].numInputChars);
                CUDACHECK(cudaMemcpy(h_inputChars.data(), callbackInput[0].d_inputChars, sizeof(char) * callbackInput[0].numInputChars, cudaMemcpyDeviceToHost));

                std::vector<int> h_inputLengths(callbackInput[0].numInputSequences);
                CUDACHECK(cudaMemcpy(h_inputLengths.data(), callbackInput[0].d_inputLengths, sizeof(int) * callbackInput[0].numInputSequences, cudaMemcpyDeviceToHost));

                std::vector<size_t> h_inputOffsets(callbackInput[0].numInputSequences);
                thrust::device_vector<size_t> materializedOffsets(callbackInput[0].d_inputOffsets, callbackInput[0].d_inputOffsets + callbackInput[0].numInputSequences);

                // thrust::device_vector<size_t> materializedOffsets(callbackInput[0].numInputSequences, -13);
                // thrust::copy(thrust::device, callbackInput[0].d_inputOffsets, callbackInput[0].d_inputOffsets + callbackInput[0].numInputSequences, materializedOffsets.begin());
                CUDACHECK(cudaMemcpy(h_inputOffsets.data(), materializedOffsets.data().get(), sizeof(size_t) * callbackInput[0].numInputSequences, cudaMemcpyDeviceToHost));

                std::ofstream batchfile("batchfile" + std::to_string(callbackCount) + ".bin", std::ios::binary);
                size_t a = callbackInput[0].numInputChars;
                batchfile.write((const char*)&a, sizeof(size_t));
                a = callbackInput[0].numInputSequences;
                batchfile.write((const char*)&a, sizeof(size_t));
                batchfile.write((const char*)h_inputChars.data(), sizeof(char) * h_inputChars.size());
                batchfile.write((const char*)h_inputLengths.data(), sizeof(int) * h_inputLengths.size());
                batchfile.write((const char*)h_inputOffsets.data(), sizeof(size_t) * h_inputOffsets.size());
            }
            callbackCount++;
        }
        #endif


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
                            consumeKernelOutputArrays(queryIndex, gpu, stream);
                        }
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


                    for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
                        thrust::sequence(
                            thrust::cuda::par_nosync.on(stream),
                            ws.getKernelOutputArrayIndices(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex],
                            ws.getKernelOutputArrayIndices(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex] + numInPass,
                            ReferenceIdT(processedSequencesPerGpu[gpu])
                        );
                    }

                    std::vector<ScoreOnlyBatchQuery> batchQueries(numQueries);
                    for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
                        OneToAllInputDataPSSM& inputData = batchQueries[queryIndex].input;
                        inputData.d_subjects = reinterpret_cast<const int8_t*>(callbackInput[gpu].d_inputChars);
                        inputData.d_subjectOffsets = callbackInput[gpu].d_inputOffsets;
                        inputData.d_subjectLengths = callbackInput[gpu].d_inputLengths;
                        inputData.queryLength = currentQueryLengths[queryIndex];
                        inputData.numAlignments = numInPass;
                        inputData.maximumSubjectLength = callbackInput[gpu].sequenceLengthUpperBound;
                        inputData.indexIndirection = d_selectedPositions;
                        inputData.letterTranslation = gpuDatabase->getLetterTranslation();
                        batchQueries[queryIndex].scores = ws.getKernelOutputArrayScores(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex];
                        batchQueries[queryIndex].pssm = queryPssmPerGpuPerQuery_forStandard[queryIndex][gpu].get();
                    }

                    const bool scoredBatch = standardAlignerPerGpu[gpu]->scoreOnly_batch(
                        cuda::std::span<const ScoreOnlyBatchQuery>(batchQueries.data(), batchQueries.size()),
                        GapScoreArgs{.gapopenscore = gop, .gapextendscore = gex},
                        stream
                    );
                    for(int queryIndex = 0; !scoredBatch && queryIndex < numQueries; queryIndex++){
                        const OneToAllInputDataPSSM& inputData = batchQueries[queryIndex].input;

                        if(standardAlignerPerGpu[gpu]->isSingleTile(inputData)){
                            standardAlignerPerGpu[gpu]->scoreOnly_singleTile(
                                ws.getKernelOutputArrayScores(queryIndex) + ws.kernelOutputArrayUsedSize[queryIndex],
                                inputData,
                                *queryPssmPerGpuPerQuery_forStandard[queryIndex][gpu],
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
                                *queryPssmPerGpuPerQuery_forStandard[queryIndex][gpu],
                                GapScoreArgs{.gapopenscore = gop, .gapextendscore = gex},
                                stream
                            );
                        }

                    }
                }
            }

            if(overflowSettings.checkForOverflow){

                { 
                    nvtx::ScopedRange srov("overflow_correction_count_overflows", 2);
                    std::fill_n(h_numOverflowsPerGpuPerQuery.data(), numGpus * numQueries, 0);

                    for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
                        for(int gpu = 0; gpu < numGpus; gpu++){
                            if(sequencePassOffset < callbackInput[gpu].numInputSequences){
                                auto scoreDatatype = standardAlignerPerGpu[gpu]->getScoreDatatype(currentQueryLengths[queryIndex]);
                                const int overflowThreshold = getOverflowThreshold(queryIndex, scoreDatatype);
                                auto overflowDatatype = overflowSettings.overflowAlignersScoreOnly[gpu]->getScoreDatatype(currentQueryLengths[queryIndex]);
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
                                    thrust::device_vector<char, thrust_async_allocator<char>> d_tempStorage
                                        = allocateTempBuffer(cubBytes, stream);
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
                        CUDACHECK(cudaSetDevice(deviceIds[gpu]));
                        cudaStream_t stream = cudaStreamPerThread;
                        CUDACHECK(cudaStreamSynchronize(stream))
                    }

                    // std::cout << "num overflows\n";
                    // for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
                    //     std::cout << h_numOverflowsPerGpuPerQuery[0 * numQueries + queryIndex] << " ";
                    // }
                    // std::cout << "\n";
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
                                    cudaStream_t stream = streamsPerGpuPerQuery[gpu * maxNumQueries + queryIndex];

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
                                    inputData.d_subjects = reinterpret_cast<const int8_t*>(callbackInput[gpu].d_inputChars);
                                    inputData.d_subjectOffsets = callbackInput[gpu].d_inputOffsets;
                                    inputData.d_subjectLengths = callbackInput[gpu].d_inputLengths;
                                    inputData.queryLength = currentQueryLengths[queryIndex];
                                    inputData.numAlignments = numOverflows;
                                    inputData.maximumSubjectLength = callbackInput[gpu].sequenceLengthUpperBound;
                                    inputData.indexIndirection = d_overflowPositionsIterator;
                                    inputData.letterTranslation = gpuDatabase->getLetterTranslation();

                                    if(overflowSettings.overflowAlignersScoreOnly[gpu]->isSingleTile(inputData)){
                                        overflowSettings.overflowAlignersScoreOnly[gpu]->scoreOnly_singleTile(
                                            d_correctedScores,
                                            inputData,
                                            *queryPssmPerGpuPerQuery_forOverflow[queryIndex][gpu],
                                            GapScoreArgs{.gapopenscore = gop, .gapextendscore = gex},
                                            stream
                                        );
                                    }else{
                                        size_t suggestedBytes = overflowSettings.overflowAlignersScoreOnly[gpu]->getMinimumSuggestedTempBytes_multiTile(
                                            inputData
                                        );
                                        thrust::device_vector<char, thrust_async_allocator<char>> d_tempStorage
                                            = allocateTempBuffer(suggestedBytes, 256ull << 20, stream);
                                        overflowSettings.overflowAlignersScoreOnly[gpu]->scoreOnly_multiTile(
                                            d_tempStorage.data().get(),
                                            d_tempStorage.size(),
                                            d_correctedScores,
                                            inputData,
                                            *queryPssmPerGpuPerQuery_forOverflow[queryIndex][gpu],
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
                            const int numOverflows = h_numOverflowsPerGpuPerQuery[gpu * numQueries + queryIndex];
                            if(numOverflows > 0){
                                CUDACHECK(cudaSetDevice(deviceIds[gpu]));
                                auto& ws = *workingSets[gpu];
                                cudaStream_t stream = cudaStreamPerThread;
                                cudaStream_t alignmentStream = streamsPerGpuPerQuery[gpu * maxNumQueries + queryIndex];
                                CUDACHECK(cudaEventRecord(ws.forkStreamEvent, alignmentStream));
                                CUDACHECK(cudaStreamWaitEvent(stream, ws.forkStreamEvent, 0));
                            }
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

        for(int gpu = 0; gpu < numGpus; gpu++){
            totalNumberOfProcessedSequences += callbackInput[gpu].numInputSequences;
        }
    }

    void ScoreOnlySearch::consumeKernelOutputArrays(int queryIndex, int gpu, cudaStream_t stream){
        const int numGpus = deviceIds.size();
            
        auto& ws = *workingSets[gpu];
        if(ws.kernelOutputArrayUsedSize[queryIndex] > 0){
            CUDACHECK(cudaSetDevice(deviceIds[gpu]));

            //db sequences are processed in ascending order. stable sort ensures that sequences with same score are sorted by ascending id without a custom comparator
            selectTopN(
                ws.getKernelOutputArrayScores(queryIndex),
                ws.getKernelOutputArrayIndices(queryIndex),
                ws.kernelOutputArrayUsedSize[queryIndex],
                topNSize,
                stream
            );

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

    void ScoreOnlySearch::finalizeResults(cuda::std::span<cudaStream_t> streams){
        const int numGpus = deviceIds.size();
        const int numQueries = currentQueryLengths.size();


        // combine any leftover kernel outputs with topN output
        for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
            for(int gpu = 0; gpu < numGpus; gpu++){
                consumeKernelOutputArrays(queryIndex, gpu, streams[gpu]);
            }
        }

        //copy final topN results to host
        if(numGpus == 1){
            //straight-forward copy
            CUDACHECK(cudaSetDevice(deviceIds[0]));
            auto& ws = *workingSets[0];
            for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
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
            std::partial_sum(bytes.begin(), bytes.end(), bytesInclPS.begin());
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
                std::partial_sum( // exclusive scan, gatherDstOffsets[0] stays 0
                    resultSizePerQueryPerGpu[queryIndex].begin(), resultSizePerQueryPerGpu[queryIndex].end() - 1, 
                    gatherDstOffsets.begin() + 1
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

    void ScoreOnlySearch::initWorkingSets(){
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

    void ScoreOnlySearch::setAligner(
        cuda::std::span<MultiConfigScoreOnlyAlignerInterface*> aligners,
        const OverflowSettings& overflowSettings_
    ){
        assert(aligners.size() == deviceIds.size());
        standardAlignerPerGpu = aligners;
        overflowSettings = overflowSettings_;
    }


    int ScoreOnlySearch::getOverflowThreshold(int queryIndex, Datatype datatype){
        switch(datatype){
            case Datatype::Float: return 16777216; // 2^24
            case Datatype::Int: return std::numeric_limits<int>::max() - std::abs(std::min(0, currentMinMaxValuesInPssms[queryIndex].first));
            case Datatype::Half2: return 2048; // 2^11
            // the packed signed add wraps: stop one positive substitution below 32767
            case Datatype::Short2: return std::numeric_limits<short>::max() - std::max(0, currentMinMaxValuesInPssms[queryIndex].second);
            case Datatype::UInt8x4: return std::numeric_limits<uint8_t>::max() - std::abs(std::min(0, currentMinMaxValuesInPssms[queryIndex].first));
            default: throw std::runtime_error("unhandled case");
        }   
    }

    void ScoreOnlySearch::setPssms(cuda::std::span<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>> hostQueryPssmViews_){
        destroyPssms();

        hostQueryPssmViews = hostQueryPssmViews_;
        const int numQueries = hostQueryPssmViews.size();
        if(numQueries > maxNumQueries){
            throw std::runtime_error("too many queries for setPssm\n");
        }
        currentQueryLengths.resize(numQueries);
        currentMinMaxValuesInPssms.resize(numQueries);


        const int numGpus = deviceIds.size();

        for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
            currentQueryLengths[queryIndex] = hostQueryPssmViews[queryIndex].extent(1);

            int minValInPssm = std::numeric_limits<int>::max();
            int maxValInPssm = std::numeric_limits<int>::min();
            for(int r = 0; r < hostQueryPssmViews[queryIndex].extent(0); r++){
                for(int c = 0; c < hostQueryPssmViews[queryIndex].extent(1); c++){
                    int val = hostQueryPssmViews[queryIndex](r,c);
                    minValInPssm = std::min(minValInPssm, val);
                    maxValInPssm = std::max(maxValInPssm, val);
                }
            }
            currentMinMaxValuesInPssms[queryIndex] = std::make_pair(minValInPssm, maxValInPssm);

            std::vector<std::unique_ptr<libmarv::GpuConvertedPSSM>> gpuStandardPssm_vec;
            std::vector<std::unique_ptr<libmarv::GpuConvertedPSSM>> gpuOverflowPssm_vec;
    
            for(int gpu = 0; gpu < numGpus; gpu++){
                CUDACHECK(cudaSetDevice(deviceIds[gpu]));
    
                auto gpuStandardPssm = std::make_unique<libmarv::GpuConvertedPSSM>();    
                standardAlignerPerGpu[gpu]->makeGpuPssm(
                    *gpuStandardPssm,
                    hostQueryPssmViews[queryIndex],
                    cudaStreamPerThread
                );    
                gpuStandardPssm_vec.push_back(std::move(gpuStandardPssm));

                if(overflowSettings.checkForOverflow){
                    assert(int(overflowSettings.overflowAlignersScoreOnly.size()) > gpu);

                    auto gpuOverflowPssm = std::make_unique<libmarv::GpuConvertedPSSM>();    
                    overflowSettings.overflowAlignersScoreOnly[gpu]->makeGpuPssm(
                        *gpuOverflowPssm,
                        hostQueryPssmViews[queryIndex],
                        cudaStreamPerThread
                    );    
                    gpuOverflowPssm_vec.push_back(std::move(gpuOverflowPssm));
                }
            }
            queryPssmPerGpuPerQuery_forStandard.push_back(std::move(gpuStandardPssm_vec));
            if(overflowSettings.checkForOverflow){
                queryPssmPerGpuPerQuery_forOverflow.push_back(std::move(gpuOverflowPssm_vec));
            }
        }
    }

    void ScoreOnlySearch::clearProcessedSequencesAndCachedResults(){
        totalNumberOfProcessedSequences = 0;
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

    void ScoreOnlySearch::clearProcessedSequencesAndCachedResults(int queryIndex){
        const int numGpus = deviceIds.size();
        for(int gpu = 0; gpu < numGpus; gpu++){
            auto& ws = *workingSets[gpu];
            ws.kernelOutputArrayUsedSize[queryIndex] = 0;
        }

        std::fill(resultSizePerQueryPerGpu[queryIndex].begin(), resultSizePerQueryPerGpu[queryIndex].end(), 0);

        std::fill(processedSequencesPerGpu.begin(), processedSequencesPerGpu.end(), 0);
    }



LIBMARV_NAMESPACE_WITH_NESTING_END
