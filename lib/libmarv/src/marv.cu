#include <algorithm>
#include <exception>
#include <memory>
#include <string>
#include <vector>
#include <string_view>

#include "hpc_helpers/all_helpers.cuh"
#include "hpc_helpers/peer_access.cuh"

#include "dbdata.hpp"
#include "config.hpp"
#include "gpudatabaseallocation.cuh"

#include "gpudatabase.cuh"
#include "searches/score_only_search.cuh"
#include "searches/score_endpos_search.cuh"
#include "searches/score_only_foldseek_search.cuh"
#include "aligners/make_aligner_headers.cuh"

namespace b64 {
#include "base64.h"
}
#include "marv.h"

// score-only search workspace retained between calls, owned by the Marv instance
struct SearchWorkspaceCache {
    std::unique_ptr<libmarv::ScoreOnlySearch> search;
    size_t maxBytes = 0;

    std::unique_ptr<libmarv::ScoreOnlySearch> acquire(libmarv::GpuDatabase& db, int outputSize, int topN, int queries) {
        const int effectiveOutput = int(std::min(size_t(outputSize), std::max(size_t(1), db.getNumSequences())));
        const int effectiveTopN = int(std::min(size_t(topN), db.getNumSequences()));
        if (search && search->kernelOutputArraySize == effectiveOutput
            && search->maxTopNSize >= effectiveTopN && search->maxNumQueries >= queries) {
            search->gpuDatabase = &db;
            search->topNSize = effectiveTopN;
            return std::move(search);
        }
        search.reset();
        return std::make_unique<libmarv::ScoreOnlySearch>(db, outputSize, topN, queries);
    }
    void recycle(std::unique_ptr<libmarv::ScoreOnlySearch> s) {
        s->releaseQueryResources();
        if (maxBytes && s->workspaceBytes() <= maxBytes) search = std::move(s);
    }
};

std::vector<Marv::Stats> convertSearchResultToMarvResults(
    Marv::Result* results,
    const std::vector<libmarv::SearchResult>& searchResultPerQuery,
    int topNSize,
    libmarv::ScanType scanType
){
    const int numQueries = searchResultPerQuery.size();
    std::vector<Marv::Stats> marvStatsVector(numQueries);

    const int hasEndPositions = (scanType == libmarv::ScanType::Gapless_Endpos 
                || scanType == libmarv::ScanType::SW_Endpos 
                || scanType == libmarv::ScanType::GaplessPlusSW_Endpos);

    for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
        Marv::Result* marvresults = results + size_t(queryIndex) * topNSize;
        const auto& searchResult = searchResultPerQuery[queryIndex];

        for (size_t i = 0; i < searchResult.scores.size(); i++) {
            marvresults[i] = Marv::Result(
                searchResult.referenceIds[i],
                searchResult.scores[i],
                hasEndPositions ? searchResult.endPositions[i].getQueryEndInclusive() : -1,
                hasEndPositions ? searchResult.endPositions[i].getSubjectEndInclusive() : -1
            );
        }

        Marv::Stats stats;
        stats.results = searchResult.scores.size();
        if(searchResult.stats != nullptr){
            stats.seconds = searchResult.stats->seconds;
            stats.gcups = searchResult.stats->gcups;
        }else{
            stats.seconds = 0;
            stats.gcups = 0;
        }
        marvStatsVector[queryIndex] = stats;
    }

    return marvStatsVector;
}



std::vector<libmarv::SearchResult> do_search_scoreOnly_impl(
    libmarv::GpuDatabase& gpuDatabase, 
    SearchWorkspaceCache& workspaceCache,
    const std::vector<std::string_view>& queries,
    const std::vector<int8_t*>& queryPssms,
    int topNSize, 
    libmarv::ScanType scanType, 
    int gop, 
    int gex,
    bool int8IsAllowed,
    bool handleOverflows,
    std::vector<libmarv::ReferenceIdT>* subsetIdsPtr, //nullptr if unused
    int pssmRows = 21 //alphabet size of the pssms
){
    const std::vector<int>& deviceIds = gpuDatabase.getDeviceIds();
    
    std::vector<cudaStream_t> streams{cudaStreamPerThread};
    cuda::std::span<cudaStream_t> streamsSpan(streams.data(), streams.size());

    const int numQueries = queries.size();
    const int numGpus = deviceIds.size();

    std::vector<std::unique_ptr<libmarv::MultiConfigScoreOnlyAlignerInterface>> scoreOnlyAlignerVector;
    std::vector<libmarv::MultiConfigScoreOnlyAlignerInterface*> scoreOnlyAlignerPtrVector;    
    for(int i = 0; i < numGpus; i++){
        CUDACHECK(cudaSetDevice(deviceIds[i]));
        if(scanType == libmarv::ScanType::Gapless){
            if(int8IsAllowed){
                scoreOnlyAlignerVector.emplace_back(libmarv::make_multiconfig_gapless_aligner_8bit());
            }else{
                scoreOnlyAlignerVector.emplace_back(libmarv::make_multiconfig_gapless_aligner_16bit());
            }
        }else if(scanType == libmarv::ScanType::SW){
            scoreOnlyAlignerVector.emplace_back(libmarv::make_multiconfig_sw_aligner_32bit());
        }else{
            throw std::runtime_error("invalid scanType for score search");
        }
        scoreOnlyAlignerPtrVector.push_back(scoreOnlyAlignerVector.back().get());
    }
    cuda::std::span<libmarv::MultiConfigScoreOnlyAlignerInterface*> alignerSpan(scoreOnlyAlignerPtrVector.data(), scoreOnlyAlignerPtrVector.size());


    std::optional<cuda::std::span<libmarv::MultiConfigScoreOnlyAlignerInterface*>> overflowAlignerSpan_optional = std::nullopt;

    std::vector<std::unique_ptr<libmarv::MultiConfigScoreOnlyAlignerInterface>> overflowAlignerVector;
    std::vector<libmarv::MultiConfigScoreOnlyAlignerInterface*> overflowAlignerPtrVector;    
    // one correction pass: the 16-bit retry must be exact, FP16 aligners are not above 2048
    std::vector<int64_t> scoreBounds(numQueries, 0);
    if(int8IsAllowed && handleOverflows && scanType == libmarv::ScanType::Gapless){
        for(int q = 0; q < numQueries; ++q){
            const size_t length = queries[q].size();
            for(size_t c = 0; c < length; ++c){
                int columnMax = 0;
                for(int r = 0; r < pssmRows; ++r)
                    columnMax = std::max(columnMax, int(queryPssms[q][size_t(r) * length + c]));
                scoreBounds[q] += columnMax;
            }
        }
    }
    for(int i = 0; i < numGpus; i++){
        CUDACHECK(cudaSetDevice(deviceIds[i]));
        if(scanType == libmarv::ScanType::Gapless){
            if(int8IsAllowed){
                auto retry = libmarv::make_multiconfig_gapless_aligner_16bit();
                for(int q = 0; handleOverflows && q < numQueries; ++q){
                    const auto type = retry->getScoreDatatype(queries[q].size());
                    const int64_t limit = type == libmarv::Datatype::Half2 ? 2048 :
                        type == libmarv::Datatype::Short2 ? 32767 :
                        type == libmarv::Datatype::Float ? 16777216 : std::numeric_limits<int>::max();
                    if(scoreBounds[q] > limit){
                        retry = libmarv::make_multiconfig_gapless_aligner_32bit();
                        break;
                    }
                }
                overflowAlignerVector.emplace_back(std::move(retry));
            }else{
                overflowAlignerVector.emplace_back(libmarv::make_multiconfig_gapless_aligner_32bit());
            }
        }else if(scanType == libmarv::ScanType::SW){
            overflowAlignerVector.emplace_back(libmarv::make_multiconfig_sw_aligner_32bit()); //does not make sense to use the same aligner. only for testing
        }else{
            throw std::runtime_error("invalid scanType for score search");
        }
        overflowAlignerPtrVector.push_back(overflowAlignerVector.back().get());
    }
    overflowAlignerSpan_optional = cuda::std::span<libmarv::MultiConfigScoreOnlyAlignerInterface*>(overflowAlignerPtrVector.data(), overflowAlignerPtrVector.size());

    std::vector<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>> queryPssmViews;

    for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
        const auto& query = queries[queryIndex];
        const int queryLength = query.size();

        queryPssmViews.emplace_back(queryPssms[queryIndex], pssmRows, queryLength);
    }

    const int kernelOutputArraySize = scoreOnlyAlignerPtrVector[0]->preferredKernelOutputArraySize();
    auto search = workspaceCache.acquire(gpuDatabase, kernelOutputArraySize, topNSize, numQueries);
    struct RecycleSearch {
        SearchWorkspaceCache& cache;
        std::unique_ptr<libmarv::ScoreOnlySearch>& search;
        int exceptions = std::uncaught_exceptions();
        ~RecycleSearch() { // a failed search is discarded, a successful one retains its workspace
            if (exceptions == std::uncaught_exceptions()) cache.recycle(std::move(search));
        }
    } recycle{workspaceCache, search};
    auto& scoreOnlySearch = *search;

    auto overflowSettings = libmarv::OverflowSettings::withoutOverflowCorrection();
    if(handleOverflows){
        overflowSettings = libmarv::OverflowSettings::withOverflowCorrection(overflowAlignerSpan_optional.value());
    }

    if(subsetIdsPtr == nullptr){
        auto searchResultPerQuery = scoreOnlySearch.execute(
            alignerSpan, 
            overflowSettings,
            cuda::std::span<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>>(queryPssmViews.data(), queryPssmViews.size()), 
            gop, 
            gex
        );

        return searchResultPerQuery;
    }else{
        auto searchResultPerQuery = scoreOnlySearch.execute(
            alignerSpan, 
            overflowSettings,
            cuda::std::span<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>>(queryPssmViews.data(), queryPssmViews.size()), 
            *subsetIdsPtr,
            gop, 
            gex
        );

        return searchResultPerQuery;
    }
}


libmarv::SearchResult do_search_scoreOnly_impl(
    libmarv::GpuDatabase& gpuDatabase, 
    SearchWorkspaceCache& workspaceCache,
    std::string_view query, 
    int8_t* queryPssm,
    int topNSize, 
    libmarv::ScanType scanType, 
    int gop, 
    int gex,
    bool int8IsAllowed,
    bool handleOverflows,
    std::vector<libmarv::ReferenceIdT>* subsetIdsPtr //nullptr if unused
){
    auto resultvector = do_search_scoreOnly_impl(
        gpuDatabase, 
        workspaceCache,
        std::vector<std::string_view>{query}, std::vector<int8_t*>{queryPssm},
        topNSize, scanType, gop, gex, int8IsAllowed, handleOverflows, subsetIdsPtr
    );
    return std::move(resultvector[0]);
}





std::vector<libmarv::SearchResult> do_search_scoreEndpos_impl(
    libmarv::GpuDatabase& gpuDatabase, 
    const std::vector<std::string_view>& queries, 
    const std::vector<int8_t*> queryPssms,
    int topNSize, 
    libmarv::ScanType scanType, 
    int gop, 
    int gex,
    bool /*int8IsAllowed*/,
    bool handleOverflows,
    std::vector<libmarv::ReferenceIdT>* subsetIdsPtr //nullptr if unused
){
    const std::vector<int>& deviceIds = gpuDatabase.getDeviceIds();
    
    std::vector<cudaStream_t> streams{cudaStreamPerThread};
    cuda::std::span<cudaStream_t> streamsSpan(streams.data(), streams.size());

    const int numQueries = queries.size();
    const int numGpus = deviceIds.size();

    std::vector<std::unique_ptr<libmarv::MultiConfigScoreWithEndPosAlignerInterface>> scoreEndposAlignerVector;
    std::vector<libmarv::MultiConfigScoreWithEndPosAlignerInterface*> scoreEndposAlignerPtrVector;    
    for(int i = 0; i < numGpus; i++){
        CUDACHECK(cudaSetDevice(deviceIds[i]));
        if(scanType == libmarv::ScanType::Gapless_Endpos){
            scoreEndposAlignerVector.emplace_back(libmarv::make_multiconfig_gapless_endpos_aligner_16bit());
        }else if(scanType == libmarv::ScanType::SW_Endpos){
            scoreEndposAlignerVector.emplace_back(libmarv::make_multiconfig_sw_endpos_aligner_32bit());
        }else{
            throw std::runtime_error("invalid scanType for score endpos search");
        }
        scoreEndposAlignerPtrVector.emplace_back(scoreEndposAlignerVector.back().get());
    }
    cuda::std::span<libmarv::MultiConfigScoreWithEndPosAlignerInterface*> alignerSpan(scoreEndposAlignerPtrVector.data(), scoreEndposAlignerPtrVector.size());


    std::optional<cuda::std::span<libmarv::MultiConfigScoreWithEndPosAlignerInterface*>> overflowAlignerSpan_optional = std::nullopt;

    std::vector<std::unique_ptr<libmarv::MultiConfigScoreWithEndPosAlignerInterface>> overflowAlignerVector;
    std::vector<libmarv::MultiConfigScoreWithEndPosAlignerInterface*> overflowAlignerPtrVector;    
    for(int i = 0; i < numGpus; i++){
        CUDACHECK(cudaSetDevice(deviceIds[i]));
        if(scanType == libmarv::ScanType::Gapless_Endpos){
            overflowAlignerVector.emplace_back(libmarv::make_multiconfig_gapless_endpos_aligner_16bit());
        }else if(scanType == libmarv::ScanType::SW_Endpos){
            overflowAlignerVector.emplace_back(libmarv::make_multiconfig_sw_endpos_aligner_32bit()); //does not make sense to use the same aligner. only for testing
        }else{
            throw std::runtime_error("invalid scanType for score search");
        }
        overflowAlignerPtrVector.push_back(overflowAlignerVector.back().get());
    }
    overflowAlignerSpan_optional = cuda::std::span<libmarv::MultiConfigScoreWithEndPosAlignerInterface*>(overflowAlignerPtrVector.data(), overflowAlignerPtrVector.size());


    std::vector<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>> queryPssmViews;

    for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
        const auto& query = queries[queryIndex];
        const int queryLength = query.size();

        queryPssmViews.emplace_back(queryPssms[queryIndex], 21, queryLength);
    }

    constexpr int kernelOutputArraySize = 10'000'000;
    libmarv::ScoreEndposSearch scoreEndposSearch(gpuDatabase, kernelOutputArraySize, topNSize, numQueries);


    auto overflowSettings = libmarv::OverflowSettings::withoutOverflowCorrection();
    if(handleOverflows){
        overflowSettings = libmarv::OverflowSettings::withOverflowCorrection(overflowAlignerSpan_optional.value());
    }

    if(subsetIdsPtr == nullptr){    
        auto searchResultPerQuery = scoreEndposSearch.execute(
            alignerSpan, 
            overflowSettings,
            cuda::std::span<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>>(queryPssmViews.data(), queryPssmViews.size()), 
            gop, 
            gex
        );

        return searchResultPerQuery;
    }else{
        auto searchResultPerQuery = scoreEndposSearch.execute(
            alignerSpan, 
            overflowSettings,
            cuda::std::span<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>>(queryPssmViews.data(), queryPssmViews.size()), 
            *subsetIdsPtr,
            gop, 
            gex
        );

        return searchResultPerQuery;
    }
}

libmarv::SearchResult do_search_scoreEndpos_impl(
    libmarv::GpuDatabase& gpuDatabase, 
    std::string_view query, 
    int8_t* queryPssm,
    int topNSize, 
    libmarv::ScanType scanType, 
    int gop, 
    int gex,
    bool int8IsAllowed,
    bool handleOverflows,
    std::vector<libmarv::ReferenceIdT>* subsetIdsPtr //nullptr if unused
){
    auto resultvector = do_search_scoreEndpos_impl(
        gpuDatabase, 
        std::vector<std::string_view>{query}, std::vector<int8_t*>{queryPssm},
        topNSize, scanType, gop, gex, int8IsAllowed, handleOverflows, subsetIdsPtr
    );
    return std::move(resultvector[0]);
}


std::vector<libmarv::SearchResult> do_search_scoreOnly_foldseek_impl(
    libmarv::GpuDatabase& gpuDatabase, 
    const std::vector<int>& queryLengths,
    const std::vector<int8_t*>& queryPssms_3di,
    const std::vector<int8_t*>& queryPssms_aa12,
    int topNSize, 
    const libmarv::MaskingOptions& maskingOptions,
    bool int8IsAllowed,
    bool handleOverflows
){
#if defined(__METAL_BACKEND__) || defined(__SIMD_BACKEND__)
    throw std::runtime_error("libmarv/Metal: the Foldseek search is not implemented on the Metal backend");
#else
    const std::vector<int>& deviceIds = gpuDatabase.getDeviceIds();

    std::vector<cudaStream_t> streams{cudaStreamPerThread};
    cuda::std::span<cudaStream_t> streamsSpan(streams.data(), streams.size());

    assert(queryPssms_3di.size() == queryPssms_aa12.size());
    const int numQueries = queryPssms_3di.size();
    const int numGpus = deviceIds.size();

    std::vector<std::unique_ptr<libmarv::MultiConfigScoreOnlyAlignerInterface>> scoreOnlyAlignerVector;
    std::vector<libmarv::MultiConfigScoreOnlyAlignerInterface*> scoreOnlyAlignerPtrVector;    
    for(int i = 0; i < numGpus; i++){
        CUDACHECK(cudaSetDevice(deviceIds[i]));
        if(int8IsAllowed){
            scoreOnlyAlignerVector.emplace_back(libmarv::make_multiconfig_gapless_aligner_8bit());
        }else{
            scoreOnlyAlignerVector.emplace_back(libmarv::make_multiconfig_gapless_aligner_16bit());
        }
        scoreOnlyAlignerPtrVector.push_back(scoreOnlyAlignerVector.back().get());
    }
    cuda::std::span<libmarv::MultiConfigScoreOnlyAlignerInterface*> alignerSpan(scoreOnlyAlignerPtrVector.data(), scoreOnlyAlignerPtrVector.size());


    std::optional<cuda::std::span<libmarv::MultiConfigScoreOnlyAlignerInterface*>> overflowAlignerSpan_optional = std::nullopt;

    std::vector<std::unique_ptr<libmarv::MultiConfigScoreOnlyAlignerInterface>> overflowAlignerVector;
    std::vector<libmarv::MultiConfigScoreOnlyAlignerInterface*> overflowAlignerPtrVector;    
    for(int i = 0; i < numGpus; i++){
        CUDACHECK(cudaSetDevice(deviceIds[i]));
        if(int8IsAllowed){
            overflowAlignerVector.emplace_back(libmarv::make_multiconfig_gapless_aligner_16bit());
        }else{
            overflowAlignerVector.emplace_back(libmarv::make_multiconfig_gapless_aligner_32bit());
        }
        overflowAlignerPtrVector.push_back(overflowAlignerVector.back().get());
    }
    overflowAlignerSpan_optional = cuda::std::span<libmarv::MultiConfigScoreOnlyAlignerInterface*>(overflowAlignerPtrVector.data(), overflowAlignerPtrVector.size());


    std::vector<std::unique_ptr<libmarv::MultiConfigScoreWithEndPosAlignerInterface>> scoreEndposAlignerVector;
    std::vector<libmarv::MultiConfigScoreWithEndPosAlignerInterface*> scoreEndposAlignerPtrVector;    
    for(int i = 0; i < numGpus; i++){
        CUDACHECK(cudaSetDevice(deviceIds[i]));
        scoreEndposAlignerVector.emplace_back(libmarv::make_multiconfig_gapless_endpos_aligner_16bit());
        scoreEndposAlignerPtrVector.emplace_back(scoreEndposAlignerVector.back().get());
    }
    cuda::std::span<libmarv::MultiConfigScoreWithEndPosAlignerInterface*> alignerEndposSpan(scoreEndposAlignerPtrVector.data(), scoreEndposAlignerPtrVector.size());

    std::vector<std::unique_ptr<libmarv::MultiConfigScoreWithEndPosAlignerInterface>> overflowEndposAlignerVector;
    std::vector<libmarv::MultiConfigScoreWithEndPosAlignerInterface*> overflowEndposAlignerPtrVector;
    for(int i = 0; i < numGpus; i++){
        CUDACHECK(cudaSetDevice(deviceIds[i]));
        overflowEndposAlignerVector.emplace_back(libmarv::make_multiconfig_gapless_endpos_aligner_32bit());
        overflowEndposAlignerPtrVector.emplace_back(overflowEndposAlignerVector.back().get());
    }
    cuda::std::span<libmarv::MultiConfigScoreWithEndPosAlignerInterface*> overflowEndposSpan(overflowEndposAlignerPtrVector.data(), overflowEndposAlignerPtrVector.size());

    std::vector<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>> queryPssmViews_3di;
    std::vector<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>> queryPssmViews_aa12;

    for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
        const int queryLength = queryLengths[queryIndex];

        queryPssmViews_3di.emplace_back(queryPssms_3di[queryIndex], 21, queryLength);
        queryPssmViews_aa12.emplace_back(queryPssms_aa12[queryIndex], 12, queryLength);
    }

    constexpr int kernelOutputArraySize = 10'000'000;
    constexpr bool isAA20dbTesting = false;
    libmarv::ScoreOnlyFoldseekSearch scoreOnlyFoldseekSearch(isAA20dbTesting, gpuDatabase, kernelOutputArraySize, topNSize, numQueries);

    auto overflowSettings = libmarv::OverflowSettings::withoutOverflowCorrection();
    if(handleOverflows){
        overflowSettings = libmarv::OverflowSettings::withOverflowCorrection(overflowAlignerSpan_optional.value());
    }

    auto endposOverflowSettings = libmarv::OverflowSettings::withOverflowCorrection(overflowEndposSpan);

    auto searchResultPerQuery = scoreOnlyFoldseekSearch.execute(
        alignerSpan, 
        overflowSettings,
        alignerEndposSpan,
        endposOverflowSettings,
        cuda::std::span<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>>(queryPssmViews_3di.data(), queryPssmViews_3di.size()), 
        cuda::std::span<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>>(queryPssmViews_aa12.data(), queryPssmViews_aa12.size()), 
        maskingOptions
    );

    return searchResultPerQuery;

#endif
}

libmarv::SearchResult do_search_scoreOnly_foldseek_impl(
    libmarv::GpuDatabase& gpuDatabase, 
    int queryLength,
    int8_t* queryPssm_3di,
    int8_t* queryPssm_aa12,
    int topNSize, 
    const libmarv::MaskingOptions& maskingOptions,
    bool int8IsAllowed,
    bool handleOverflows
){
    auto resultvector = do_search_scoreOnly_foldseek_impl(
        gpuDatabase, 
        std::vector<int>{queryLength}, 
        std::vector<int8_t*>{queryPssm_3di}, 
        std::vector<int8_t*>{queryPssm_aa12}, 
        topNSize, maskingOptions, int8IsAllowed, handleOverflows
    );
    return std::move(resultvector[0]);
}



libmarv::ScanType mapCswScanType(Marv::AlignmentType type) {
    switch (type) {
        case Marv::AlignmentType::GAPLESS:
            return libmarv::ScanType::Gapless;
        case Marv::AlignmentType::SMITH_WATERMAN:
            return libmarv::ScanType::SW;
        case Marv::AlignmentType::GAPLESS_ENDPOS:
            return libmarv::ScanType::Gapless_Endpos;
        case Marv::AlignmentType::SMITH_WATERMAN_ENDPOS:
            return libmarv::ScanType::SW_Endpos;
        case Marv::AlignmentType::GAPLESS_THEN_SMITH_WATERMAN:
            return libmarv::ScanType::GaplessPlusSW_Endpos;
        case Marv::AlignmentType::FOLDSEEK_3DI_AA12:
            return libmarv::ScanType::Foldseek;
        default:
            return libmarv::ScanType::Gapless;
    }
}


Marv::Marv(size_t dbEntries, int alphabetSize) : Marv(dbEntries, alphabetSize, MemoryConfig{}) {}

Marv::Marv(size_t dbEntries, int alphabetSize, const MemoryConfig& config)
    : dbEntries(dbEntries), alphabetSize(alphabetSize), dbmanager(NULL) {
    std::vector<int> deviceIds = getDeviceIds();
    helpers::PeerAccess peerAccess(deviceIds, false);
    peerAccess.enableAllPeerAccesses();

    libmarv::GpuDatabase::MemoryConfig memoryConfig;
    memoryConfig.maxBatchBytes = config.maxBatchBytes;
    memoryConfig.maxBatchSequences = config.maxBatchSequences;
    memoryConfig.maxGpuMem = config.maxGpuMem;
    auto* workspaceCache = new SearchWorkspaceCache();
    workspaceCache->maxBytes = config.maxCachedQueryBytes;
    searchWorkspaceVoidPtr = static_cast<void*>(workspaceCache);

    const bool verbose = false;

    libmarv::GpuDatabase* gpuDatabase = new libmarv::GpuDatabase(
        deviceIds, 
        memoryConfig, 
        verbose
    );
    gpuDatabaseVoidPtr = static_cast<void*>(gpuDatabase);
}
std::vector<std::shared_ptr<libmarv::GpuDatabaseAllocationBase>> allocations_all;

Marv::~Marv() {
    allocations_all.clear();
    delete static_cast<SearchWorkspaceCache*>(searchWorkspaceVoidPtr);
    libmarv::GpuDatabase* gpuDatabase = static_cast<libmarv::GpuDatabase*>(gpuDatabaseVoidPtr);
    delete gpuDatabase;
}

void Marv::releaseSearchWorkspace() {
    static_cast<SearchWorkspaceCache*>(searchWorkspaceVoidPtr)->search.reset();
}

std::vector<int> Marv::getDeviceIds() {
    std::vector<int> deviceIds;
    int num = 0;
    cudaGetDeviceCount(&num); CUERR
    for (int i = 0; i < num; i++) {
        deviceIds.push_back(i);
    }
    return deviceIds;
}

void* Marv::loadDb(char* data, size_t* offset, int32_t* length, size_t dbByteSize) {
    return static_cast<void*>(new libmarv::MMseqsDB(libmarv::loadMMseqsDB(
        dbEntries, data, offset, length, dbByteSize
    )));
}

void* Marv::loadDb(char* /*data*/, size_t dbByteSize, void* otherdb) {
    libmarv::MMseqsDB* db = static_cast<libmarv::MMseqsDB*>(otherdb);
    const libmarv::DBdataMetaData& meta = db->getData().getMetaData();
    return static_cast<void*>(new libmarv::ExternalDB(libmarv::loadExternalDB(
        dbEntries, dbByteSize, meta
    )));
}

std::vector<std::string> split(const std::string &str, const std::string &sep) {
    std::vector<std::string> arr;

    char *cstr = strdup(str.c_str());
    const char* csep = sep.c_str();
    char *rest;
    char *current = strtok_r(cstr, csep, &rest);
    while (current != NULL) {
        arr.emplace_back(current);
        current = strtok_r(NULL, csep, &rest);
    }
    free(cstr);

    return arr;
}

void Marv::setDbWithAllocation(void* dbhandle, const std::string& allocationinfo) {
    auto parts = split(allocationinfo, ":");

    cudaIpcMemHandle_t h1, h2, h3;
    char* charData;
    libmarv::SequenceLengthT* lengths;
    size_t* offsets;
    size_t numChars, numSubjects;

    std::string decode;

    decode = b64::base64_decode(parts[0].data(), parts[0].length());
    memcpy((unsigned char *)(&h1), (unsigned char *)decode.data(), decode.length());
    cudaIpcOpenMemHandle((void **)&charData, h1, cudaIpcMemLazyEnablePeerAccess);
    CUERR
    decode = b64::base64_decode(parts[1].data(), parts[1].length());
    memcpy((unsigned char *)(&h2), (unsigned char *)decode.data(), decode.length());
    cudaIpcOpenMemHandle((void **)&lengths, h2, cudaIpcMemLazyEnablePeerAccess);
    CUERR
    decode = b64::base64_decode(parts[2].data(), parts[2].length());
    memcpy((unsigned char *)(&h3), (unsigned char *)decode.data(), decode.length());
    cudaIpcOpenMemHandle((void **)&offsets, h3, cudaIpcMemLazyEnablePeerAccess);
    CUERR
    numChars = strtoull(parts[3].c_str(), NULL, 10);
    numSubjects = strtoull(parts[4].c_str(), NULL, 10);
    std::vector<std::shared_ptr<libmarv::GpuDatabaseAllocationBase>> allocations_remote;
    allocations_remote.emplace_back(std::make_shared<libmarv::GpuDatabaseAllocationView>(libmarv::GpuDatabaseAllocationView(charData, lengths, offsets, numChars, numSubjects)));

    libmarv::MMseqsDB* db = static_cast<libmarv::MMseqsDB*>(dbhandle);
    auto doNothingDeleter = [](libmarv::MMseqsDB*){ /* do nothing */ };
    std::shared_ptr<libmarv::MMseqsDB> dbPtr(static_cast<libmarv::MMseqsDB*>(db), doNothingDeleter);

    libmarv::GpuDatabase* gpuDatabase = static_cast<libmarv::GpuDatabase*>(gpuDatabaseVoidPtr);
    // OpaqueAllocationManager* manager = static_cast<OpaqueAllocationManager*>(allocationhandle);
    static_cast<SearchWorkspaceCache*>(searchWorkspaceVoidPtr)->search.reset();
    gpuDatabase->setDatabase(dbPtr, allocations_remote);
}

void Marv::setDb(void* dbhandle) {
    libmarv::MMseqsDB* db = static_cast<libmarv::MMseqsDB*>(dbhandle);
    auto doNothingDeleter = [](libmarv::MMseqsDB*){ /* do nothing */ };
    std::shared_ptr<libmarv::MMseqsDB> dbPtr(static_cast<libmarv::MMseqsDB*>(db), doNothingDeleter);
    libmarv::GpuDatabase* gpuDatabase = static_cast<libmarv::GpuDatabase*>(gpuDatabaseVoidPtr);
    static_cast<SearchWorkspaceCache*>(searchWorkspaceVoidPtr)->search.reset();
    gpuDatabase->setDatabase(dbPtr);
}

std::string Marv::getDbMemoryHandle() {
    libmarv::GpuDatabase* gpuDatabase = static_cast<libmarv::GpuDatabase*>(gpuDatabaseVoidPtr);
    // sw->printDBInfo();
    // sw->printDBLengthPartitions();
    gpuDatabase->prefetchDBToGpus();
    allocations_all = gpuDatabase->getFullGpuDBAllocations();
    cudaIpcMemHandle_t h1, h2, h3;
    // char* charData;
    // libmarv::SequenceLengthT* lengths;
    // size_t* offsets;
    size_t numChars, numSubjects;
    // for(const auto& alloc : allocations_all){
        const auto& alloc = allocations_all[0];
        cudaIpcGetMemHandle(&h1, alloc->getCharData());
        CUERR
        cudaIpcGetMemHandle(&h2, alloc->getLengthData());
        CUERR
        cudaIpcGetMemHandle(&h3, alloc->getOffsetData());
        CUERR
        // charData = alloc->getCharData();
        // lengths = alloc->getLengthData();
        // offsets = alloc->getOffsetData();
        numChars = alloc->getNumChars();
        numSubjects = alloc->getNumSubjects();
    // }

    std::vector<std::string> handles;
    std::string enc1 = b64::base64_encode(&h1, sizeof(cudaIpcMemHandle_t));
    std::string enc2 = b64::base64_encode(&h2, sizeof(cudaIpcMemHandle_t));
    std::string enc3 = b64::base64_encode(&h3, sizeof(cudaIpcMemHandle_t));

    std::string res;
    res.append(enc1);
    res.append(1, ':');
    res.append(enc2);
    res.append(1, ':');
    res.append(enc3);
    res.append(1, ':');
    res.append(std::to_string(numChars));
    res.append(1, ':');
    res.append(std::to_string(numSubjects));
    res.append(1, '\n');

    return res;
}

void Marv::printInfo() {
    libmarv::GpuDatabase* gpuDatabase = static_cast<libmarv::GpuDatabase*>(gpuDatabaseVoidPtr);
    gpuDatabase->printDBInfo();
    gpuDatabase->printDBLengthPartitions();
}

void Marv::setLetterTranslation(const unsigned char* table) {
#if defined(__METAL_BACKEND__) || defined(__SIMD_BACKEND__)
    static_cast<libmarv::GpuDatabase*>(gpuDatabaseVoidPtr)->setLetterTranslation(table);
#else
    (void)table;
    throw std::runtime_error("letter translation is only available on the Metal and SIMD backends; use a padded GPU database");
#endif
}

void Marv::prefetch() {
    libmarv::GpuDatabase* gpuDatabase = static_cast<libmarv::GpuDatabase*>(gpuDatabaseVoidPtr);
    gpuDatabase->prefetchDBToGpus();
}


//sequence must be encoded
// Marv::Stats Marv::scan(Marv::AlignmentType alignmentType, const char* sequence, size_t sequenceLength, int8_t* pssm, Result* results) {
//     libmarv::ScanType scanType = mapCswScanType(alignmentType);
    
//     libmarv::CudaSW4* sw = static_cast<libmarv::CudaSW4*>(cudasw);
//     libmarv::EncodedQueryView queryView(sequence, sequenceLength);
//     libmarv::ScanResult scanResult = sw->scan(queryView, pssm);
//     for (size_t i = 0; i < scanResult.scores.size(); i++) {
//         results[i] = Result(
//             scanResult.referenceIds[i],
//             scanResult.scores[i],
//             alignmentType != GAPLESS ? scanResult.endPositions[i].getQueryEndInclusive() : -1,
//             alignmentType != GAPLESS ? scanResult.endPositions[i].getSubjectEndInclusive() : -1
//         );
//     }
//     Stats stats;
//     stats.results = scanResult.scores.size();
//     stats.numOverflows = scanResult.stats.numOverflows;
//     stats.seconds = scanResult.stats.seconds;
//     stats.gcups = scanResult.stats.gcups;
//     return stats;
// }




std::vector<Marv::Stats> Marv::search_scoreOnly(
    Marv::Result* results,
    const std::vector<const char*>& queries,
    const std::vector<int>& queryLengths,
    const std::vector<int8_t*>& queryPssms,
    int topNSize, 
    Marv::AlignmentType alignmentType, 
    int gop, 
    int gex,
    bool int8IsAllowed,
    bool handleOverflows,
    std::vector<std::int32_t>* subsetIdsPtr_ //nullptr if unused
){
    libmarv::GpuDatabase* const gpuDatabase = static_cast<libmarv::GpuDatabase*>(gpuDatabaseVoidPtr);
    libmarv::ScanType scanType = mapCswScanType(alignmentType);
    assert(queries.size() == queryPssms.size());
    static_assert(std::is_same_v<libmarv::ReferenceIdT, std::int32_t>);
    std::vector<libmarv::ReferenceIdT>* subsetIdsPtr = reinterpret_cast<std::vector<libmarv::ReferenceIdT>*>(subsetIdsPtr_);
    const int numQueries = queries.size();

    std::vector<std::string_view> queryViews(numQueries);
    for(int i = 0; i < numQueries; i++){
        queryViews[i] = std::string_view(queries[i], queryLengths[i]);
    }

    auto searchResultPerQuery = do_search_scoreOnly_impl(
        *gpuDatabase,
        *static_cast<SearchWorkspaceCache*>(searchWorkspaceVoidPtr),
        queryViews,
        queryPssms,
        topNSize,
        scanType,
        gop,
        gex,
        int8IsAllowed,
        handleOverflows,
        subsetIdsPtr
#if defined(__METAL_BACKEND__) || defined(__SIMD_BACKEND__)
        , alphabetSize
#endif
    );

    std::vector<Marv::Stats> marvStatsVector = convertSearchResultToMarvResults(
        results,
        searchResultPerQuery,
        topNSize,
        scanType        
    );

    return marvStatsVector;
}


Marv::Stats Marv::search_scoreOnly(
    Marv::Result* results,
    const char* query,
    int queryLength,
    int8_t* queryPssm,
    int topNSize, 
    Marv::AlignmentType alignmentType, 
    int gop, 
    int gex,
    bool int8IsAllowed,
    bool handleOverflows,
    std::vector<std::int32_t>* subsetIdsPtr_ //nullptr if unused
){
    auto statsvector = search_scoreOnly(
        results,
        std::vector<const char*>{query},
        std::vector<int>{queryLength},
        std::vector<int8_t*>{queryPssm}, 
        topNSize, alignmentType, gop, gex, int8IsAllowed, handleOverflows, subsetIdsPtr_
    );
    return statsvector[0];
}



std::vector<Marv::Stats> Marv::search_scoreEndpos(
    Marv::Result* results,
    const std::vector<const char*>& queries,
    const std::vector<int>& queryLengths,
    const std::vector<int8_t*>& queryPssms,
    int topNSize, 
    Marv::AlignmentType alignmentType, 
    int gop, 
    int gex,
    bool int8IsAllowed,
    bool handleOverflows,
    std::vector<std::int32_t>* subsetIdsPtr_ //nullptr if unused
){
    libmarv::GpuDatabase* const gpuDatabase = static_cast<libmarv::GpuDatabase*>(gpuDatabaseVoidPtr);
    libmarv::ScanType scanType = mapCswScanType(alignmentType);
    assert(queries.size() == queryPssms.size());
    static_assert(std::is_same_v<libmarv::ReferenceIdT, std::int32_t>);
    std::vector<libmarv::ReferenceIdT>* subsetIdsPtr = reinterpret_cast<std::vector<libmarv::ReferenceIdT>*>(subsetIdsPtr_);
    const int numQueries = queries.size();

    std::vector<std::string_view> queryViews(numQueries);
    for(int i = 0; i < numQueries; i++){
        queryViews[i] = std::string_view(queries[i], queryLengths[i]);
    }

    auto searchResultPerQuery = do_search_scoreEndpos_impl(
        *gpuDatabase,
        queryViews,
        queryPssms,
        topNSize,
        scanType,
        gop,
        gex,
        int8IsAllowed,
        handleOverflows,
        subsetIdsPtr
    );

    std::vector<Marv::Stats> marvStatsVector = convertSearchResultToMarvResults(
        results,
        searchResultPerQuery,
        topNSize,
        scanType
    );

    return marvStatsVector;
}


Marv::Stats Marv::search_scoreEndpos(
    Marv::Result* results,
    const char* query,
    int queryLength,
    int8_t* queryPssm,
    int topNSize, 
    Marv::AlignmentType alignmentType, 
    int gop, 
    int gex,
    bool int8IsAllowed,
    bool handleOverflows,
    std::vector<std::int32_t>* subsetIdsPtr_ //nullptr if unused
){
    auto statsvector = search_scoreEndpos(
        results,
        std::vector<const char*>{query},
        std::vector<int>{queryLength},
        std::vector<int8_t*>{queryPssm}, 
        topNSize, alignmentType, gop, gex, int8IsAllowed, handleOverflows, subsetIdsPtr_
    );
    return statsvector[0];
}



std::vector<Marv::Stats> Marv::search_scoreOnly_then_scoreEndpos_for_tops(
    Marv::Result* results,
    const std::vector<const char*>& queries,
    const std::vector<int>& queryLengths,
    const std::vector<int8_t*>& queryPssms,
    int topNSize, 
    int gop, 
    int gex,
    Marv::AlignmentType alignmentType_scoreOnly, 
    bool int8IsAllowed_scoreOnly,
    bool handleOverflows_scoreOnly,
    Marv::AlignmentType alignmentType_endpos, 
    bool int8IsAllowed_endpos, 
    bool handleOverflows_endpos
){
    libmarv::GpuDatabase* const gpuDatabase = static_cast<libmarv::GpuDatabase*>(gpuDatabaseVoidPtr);
    libmarv::ScanType scanType_scoreOnly = mapCswScanType(alignmentType_scoreOnly);
    libmarv::ScanType scanType_endpos = mapCswScanType(alignmentType_endpos);
    assert(queries.size() == queryPssms.size());

    const int numQueries = queries.size();

    std::vector<std::string_view> queryViews(numQueries);
    for(int i = 0; i < numQueries; i++){
        queryViews[i] = std::string_view(queries[i], queryLengths[i]);
    }

    auto scoreOnlyResults = do_search_scoreOnly_impl(
        *gpuDatabase,
        *static_cast<SearchWorkspaceCache*>(searchWorkspaceVoidPtr),
        queryViews,
        queryPssms,
        topNSize,
        scanType_scoreOnly,
        gop,
        gex,
        int8IsAllowed_scoreOnly,
        handleOverflows_scoreOnly,
        nullptr
    );

    std::vector<libmarv::SearchResult> swEndposResults;

    for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
        swEndposResults.push_back(
            do_search_scoreEndpos_impl(
                *gpuDatabase,
                queryViews[queryIndex],
                queryPssms[queryIndex],
                scoreOnlyResults[queryIndex].referenceIds.size(), //topNSize
                scanType_endpos,
                gop,
                gex,
                int8IsAllowed_endpos,
                handleOverflows_endpos,
                &(scoreOnlyResults[queryIndex].referenceIds)
            )
        );
    }

    std::vector<Marv::Stats> marvStatsVector = convertSearchResultToMarvResults(
        results,
        swEndposResults,
        topNSize,
        libmarv::ScanType::SW_Endpos
    );

    return marvStatsVector;
}


Marv::Stats Marv::search_scoreOnly_then_scoreEndpos_for_tops(
    Marv::Result* results,
    const char* query,
    int queryLength,
    int8_t* queryPssm,
    int topNSize, 
    int gop, 
    int gex,
    Marv::AlignmentType alignmentType_scoreOnly, 
    bool int8IsAllowed_scoreOnly,
    bool handleOverflows_scoreOnly,
    Marv::AlignmentType alignmentType_endpos, 
    bool int8IsAllowed_endpos, 
    bool handleOverflows_endpos
){
    auto statsvector = search_scoreOnly_then_scoreEndpos_for_tops(
        results,
        std::vector<const char*>{query},
        std::vector<int>{queryLength},
        std::vector<int8_t*>{queryPssm}, 
        topNSize, gop, gex, 
        alignmentType_scoreOnly, int8IsAllowed_scoreOnly, handleOverflows_scoreOnly,
        alignmentType_endpos, int8IsAllowed_endpos, handleOverflows_endpos
    );
    return statsvector[0];
}


std::vector<Marv::Stats> Marv::search_scoreOnly_foldseek(
    Marv::Result* results,
    const std::vector<int>& queryLengths,
    const std::vector<int8_t*>& queryPssms_3di,
    const std::vector<int8_t*>& queryPssms_aa12,
    int topNSize, 
    bool int8IsAllowed,
    bool handleOverflows,
    char maskingLetter3di,
    int maskingThreshold
){
    libmarv::GpuDatabase* const gpuDatabase = static_cast<libmarv::GpuDatabase*>(gpuDatabaseVoidPtr);
    libmarv::MaskingOptions maskingOptions(maskingThreshold, maskingLetter3di);

    assert(queryLengths.size() == queryPssms_3di.size());
    assert(queryLengths.size() == queryPssms_aa12.size());
    
    const int numQueries = queryLengths.size();

    auto searchResultPerQuery = do_search_scoreOnly_foldseek_impl(
        *gpuDatabase,
        queryLengths,
        queryPssms_3di,
        queryPssms_aa12,
        topNSize,
        maskingOptions,
        int8IsAllowed,
        handleOverflows
    );

    std::vector<Marv::Stats> marvStatsVector = convertSearchResultToMarvResults(
        results,
        searchResultPerQuery,
        topNSize,
        libmarv::ScanType::Foldseek        
    );

    return marvStatsVector;
}


Marv::Stats Marv::search_scoreOnly_foldseek(
    Marv::Result* results,
    int queryLength,
    int8_t* queryPssm_3di,
    int8_t* queryPssm_aa12,
    int topNSize, 
    bool int8IsAllowed,
    bool handleOverflows,
    char maskingLetter3di,
    int maskingThreshold
){
    auto statsvector = search_scoreOnly_foldseek(
        results,
        std::vector<int>{queryLength}, 
        std::vector<int8_t*>{queryPssm_3di}, 
        std::vector<int8_t*>{queryPssm_aa12}, 
        topNSize, int8IsAllowed, handleOverflows,
        maskingLetter3di, maskingThreshold
    );
    return statsvector[0];
}




