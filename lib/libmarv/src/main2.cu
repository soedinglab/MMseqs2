

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "hpc_helpers/all_helpers.cuh"
#include "hpc_helpers/peer_access.cuh"

#include "kseqpp/kseqpp.hpp"
#include "sequence_io.h"
#include "options.hpp"
#include "dbdata.hpp"

#include "config.hpp"
#include "pssm.cuh"
#include "convert.cuh"
#include "target_subject_ids.cuh"



#include "searches/score_only_search.cuh"
#include "searches/score_endpos_search.cuh"
#include "searches/score_only_foldseek_search.cuh"
#include "gpudatabase.cuh"
#include "aligners/make_aligner_headers.cuh"
#include "alignment_algorithms/cpu_referencealignments.hpp"



std::vector<std::string> split(const std::string& str, char c){
	std::vector<std::string> result;

	std::stringstream ss(str);
	std::string s;

	while (std::getline(ss, s, c)) {
		result.emplace_back(s);
	}

	return result;
}

void printScanResultPlain(std::ostream& os, const libmarv::SearchResult& scanResult, const libmarv::GpuDatabase& gpuDatabase){
    const int n = scanResult.scores.size();
    for(int i = 0; i < n; i++){
        const auto referenceId = scanResult.referenceIds[i];
        os << "Result " << i << ".";
        os << " Score: " << scanResult.scores[i] << ".";
        os << " Length: " << gpuDatabase.getReferenceLength(referenceId) << ".";
        os << " Header " << gpuDatabase.getReferenceHeader(referenceId) << ".";
        os << " referenceId " << referenceId << ".";
        if(scanResult.endPositions.size() > 0){
            os << " Alignment_end_query " << scanResult.endPositions[i].getQueryEndInclusive() << ".";
            os << " Alignment_end_ref " << scanResult.endPositions[i].getSubjectEndInclusive() << ".";
        }
        os << "\n";
        //std::cout << " Sequence " << cudaSW4.getReferenceSequence(referenceId) << "\n";

    }
}

void printTSVHeader(std::ostream& os){
    constexpr char sep = '\t';

    os << "Query number" << sep 
        << "Query length" << sep 
        << "Query header" << sep
        << "Result number" << sep
        << "Result score" << sep
        << "Reference length" << sep
        << "Reference header" << sep
        << "Reference ID in DB" << sep
        << "Alignment_end_query" << sep
        << "Alignment_end_ref" << sep
        << "\n";
}

void printScanResultTSV(
    std::ostream& os, 
    const libmarv::SearchResult& scanResult, 
    const libmarv::GpuDatabase& gpuDatabase, 
    int64_t queryId,
    libmarv::SequenceLengthT queryLength,
    std::string_view queryHeader
){
    constexpr char sep = '\t';

    const int n = scanResult.scores.size();
    for(int i = 0; i < n; i++){
        const auto referenceId = scanResult.referenceIds[i];
        
        os << queryId << sep 
            << queryLength << sep
            << queryHeader << sep
            << i << sep
            << scanResult.scores[i] << sep
            << gpuDatabase.getReferenceLength(referenceId) << sep
            << gpuDatabase.getReferenceHeader(referenceId) << sep
            << referenceId << sep;
            if(scanResult.endPositions.size() > 0){
                os << scanResult.endPositions[i].getQueryEndInclusive() << sep
                    << scanResult.endPositions[i].getSubjectEndInclusive();
            }else{
                os << -1 << sep
                    << -1;
            }
            os << "\n";

        //std::cout << " Sequence " << cudaSW4.getReferenceSequence(referenceId) << "\n";
    }
}


std::vector<libmarv::SearchResult> do_search_scoreOnly(
    libmarv::GpuDatabase& gpuDatabase, 
    const std::vector<std::string_view>& queries, 
    int topNSize, 
    libmarv::ScanType scanType, 
    int gop, 
    int gex,
    bool int8IsAllowed,
    bool handleOverflows,
    std::vector<libmarv::ReferenceIdT>* subsetIdsPtr //nullptr if unused
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
    for(int i = 0; i < numGpus; i++){
        CUDACHECK(cudaSetDevice(deviceIds[i]));
        if(scanType == libmarv::ScanType::Gapless){
            if(int8IsAllowed){
                overflowAlignerVector.emplace_back(libmarv::make_multiconfig_gapless_aligner_16bit());
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

    

    auto blosum = libmarv::BLOSUM62_20::get2D();

    std::vector<std::vector<std::int8_t>> queryPssmVector;
    std::vector<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>> queryPssmViews;

    for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
        const auto& query = queries[queryIndex];
        const int queryLength = query.size();

        std::vector<char> currentQueryEncodedHost(queryLength);
        std::transform(
            query.data(),
            query.data() + queryLength,
            currentQueryEncodedHost.begin(),
            libmarv::ConvertAA_20{}
        );
        std::vector<std::int8_t> queryPssm(21 * queryLength);
        for(int s = 0; s < 21; s++){
            for(int p = 0; p < queryLength; p++){
                auto convert = libmarv::ConvertAA_20_mmseqs_to_ncbi{};
                queryPssm[s * queryLength + p] = blosum[convert(s)][convert(currentQueryEncodedHost[p])];
            }
        }
        queryPssmVector.push_back(std::move(queryPssm));
        queryPssmViews.emplace_back(queryPssmVector.back().data(), 21, queryLength); 
    }

    constexpr int kernelOutputArraySize = 10'000'000;
    libmarv::ScoreOnlySearch scoreOnlySearch(gpuDatabase, kernelOutputArraySize, topNSize, numQueries);

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


libmarv::SearchResult do_search_scoreOnly(
    libmarv::GpuDatabase& gpuDatabase, 
    std::string_view query, 
    int topNSize, 
    libmarv::ScanType scanType, 
    int gop, 
    int gex,
    bool int8IsAllowed,
    bool handleOverflows,
    std::vector<libmarv::ReferenceIdT>* subsetIdsPtr //nullptr if unused
){
    auto resultvector = do_search_scoreOnly(
        gpuDatabase, std::vector<std::string_view>{query}, topNSize, scanType, gop, gex, int8IsAllowed, handleOverflows, subsetIdsPtr
    );
    return std::move(resultvector[0]);
}





std::vector<libmarv::SearchResult> do_search_scoreEndpos(
    libmarv::GpuDatabase& gpuDatabase, 
    const std::vector<std::string_view>& queries, 
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
            overflowAlignerVector.emplace_back(libmarv::make_multiconfig_gapless_endpos_aligner_32bit());
        }else if(scanType == libmarv::ScanType::SW_Endpos){
            overflowAlignerVector.emplace_back(libmarv::make_multiconfig_sw_endpos_aligner_32bit()); //does not make sense to use the same aligner. only for testing
        }else{
            throw std::runtime_error("invalid scanType for score search");
        }
        overflowAlignerPtrVector.push_back(overflowAlignerVector.back().get());
    }
    overflowAlignerSpan_optional = cuda::std::span<libmarv::MultiConfigScoreWithEndPosAlignerInterface*>(overflowAlignerPtrVector.data(), overflowAlignerPtrVector.size());


    auto blosum = libmarv::BLOSUM62_20::get2D();

    std::vector<std::vector<std::int8_t>> queryPssmVector;
    std::vector<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>> queryPssmViews;

    for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
        const auto& query = queries[queryIndex];
        const int queryLength = query.size();

        std::vector<char> currentQueryEncodedHost(queryLength);
        std::transform(
            query.data(),
            query.data() + queryLength,
            currentQueryEncodedHost.begin(),
            libmarv::ConvertAA_20{}
        );
        std::vector<std::int8_t> queryPssm(21 * queryLength);
        for(int s = 0; s < 21; s++){
            for(int p = 0; p < queryLength; p++){
                auto convert = libmarv::ConvertAA_20_mmseqs_to_ncbi{};
                queryPssm[s * queryLength + p] = blosum[convert(s)][convert(currentQueryEncodedHost[p])];
            }
        }
        queryPssmVector.push_back(std::move(queryPssm));
        queryPssmViews.emplace_back(queryPssmVector.back().data(), 21, queryLength); 
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

libmarv::SearchResult do_search_scoreEndpos(
    libmarv::GpuDatabase& gpuDatabase, 
    std::string_view query, 
    int topNSize, 
    libmarv::ScanType scanType, 
    int gop, 
    int gex,
    bool int8IsAllowed,
    bool handleOverflows,
    std::vector<libmarv::ReferenceIdT>* subsetIdsPtr //nullptr if unused
){
    auto resultvector = do_search_scoreEndpos(
        gpuDatabase, std::vector<std::string_view>{query}, topNSize, scanType, gop, gex, int8IsAllowed, handleOverflows, subsetIdsPtr
    );
    return std::move(resultvector[0]);
}



std::vector<libmarv::SearchResult> do_search_GaplessScoreOnly_and_SW_Endpos_for_tops(
    libmarv::GpuDatabase& gpuDatabase, 
    const std::vector<std::string_view>& queries, 
    int topNSize, 
    int gop, 
    int gex,
    bool int8IsAllowed,
    bool handleOverflows
){


    std::vector<libmarv::SearchResult> scoreOnlyResults = do_search_scoreOnly(
        gpuDatabase,
        queries,
        topNSize,
        libmarv::ScanType::Gapless,
        gop,
        gex,
        int8IsAllowed,
        handleOverflows,
        nullptr
    );

    std::vector<libmarv::SearchResult> swEndposResults;
    const int numQueries = queries.size();
    for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
        swEndposResults.push_back(
            do_search_scoreEndpos(
                gpuDatabase,
                queries[queryIndex],
                scoreOnlyResults[queryIndex].referenceIds.size(), //topNSize
                libmarv::ScanType::SW_Endpos,
                gop,
                gex,
                int8IsAllowed,
                handleOverflows,
                &(scoreOnlyResults[queryIndex].referenceIds)
            )
        );
    }


    return swEndposResults;
}


libmarv::SearchResult do_search_GaplessScoreOnly_and_SW_Endpos_for_tops(
    libmarv::GpuDatabase& gpuDatabase, 
    std::string_view query, 
    int topNSize, 
    int gop, 
    int gex,
    bool int8IsAllowed,
    bool handleOverflows
){
    auto resultvector = do_search_GaplessScoreOnly_and_SW_Endpos_for_tops(
        gpuDatabase, std::vector<std::string_view>{query}, topNSize, gop, gex, int8IsAllowed, handleOverflows
    );
    return std::move(resultvector[0]);
}





std::vector<libmarv::SearchResult> do_search_scoreOnly_foldseek(
    libmarv::GpuDatabase& gpuDatabase, 
    const std::vector<std::string_view>& queries, 
    int topNSize, 
    const libmarv::MaskingOptions& maskingOptions,
    int gop, 
    int gex,
    bool int8IsAllowed,
    bool handleOverflows
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



    auto blosum = libmarv::BLOSUM62_20::get2D();

    std::vector<std::vector<std::int8_t>> queryPssmVector_3di;
    std::vector<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>> queryPssmViews_3di;

    std::vector<std::vector<std::int8_t>> queryPssmVector_aa12;
    std::vector<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>> queryPssmViews_aa12;

    for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
        const auto& query = queries[queryIndex];
        const int queryLength = query.size();

        std::vector<char> currentQueryEncodedHost(queryLength);
        std::transform(
            query.data(),
            query.data() + queryLength,
            currentQueryEncodedHost.begin(),
            libmarv::ConvertAA_20{}
        );
        std::vector<std::int8_t> queryPssm_3di(21 * queryLength);
        for(int s = 0; s < 21; s++){
            for(int p = 0; p < queryLength; p++){
                auto convert = libmarv::ConvertAA_20_mmseqs_to_ncbi{};
                queryPssm_3di[s * queryLength + p] = blosum[convert(s)][convert(currentQueryEncodedHost[p])];
            }
        }
        queryPssmVector_3di.push_back(std::move(queryPssm_3di));
        queryPssmViews_3di.emplace_back(queryPssmVector_3di.back().data(), 21, queryLength); 

        std::vector<std::int8_t> queryPssm_aa12(12 * queryLength);
        for(int s = 0; s < 12; s++){
            for(int p = 0; p < queryLength; p++){
                auto convert = libmarv::ConvertAA_20_mmseqs_to_ncbi{};
                queryPssm_aa12[s * queryLength + p] = blosum[convert(s)][convert(currentQueryEncodedHost[p])];
            }
        }
        queryPssmVector_aa12.push_back(std::move(queryPssm_aa12));
        queryPssmViews_aa12.emplace_back(queryPssmVector_aa12.back().data(), 12, queryLength); 
    }

    constexpr int kernelOutputArraySize = 10'000'000;
    constexpr bool isAA20dbTesting = true;
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
        maskingOptions,
        gop,
        gex
    );

    return searchResultPerQuery;

}

libmarv::SearchResult do_search_scoreOnly_foldseek(
    libmarv::GpuDatabase& gpuDatabase, 
    std::string_view query, 
    int topNSize, 
    const libmarv::MaskingOptions& maskingOptions,
    int gop, 
    int gex,
    bool int8IsAllowed,
    bool handleOverflows
){
    auto resultvector = do_search_scoreOnly_foldseek(
        gpuDatabase, std::vector<std::string_view>{query}, topNSize, maskingOptions, gop, gex, int8IsAllowed, handleOverflows
    );
    return std::move(resultvector[0]);
}



//return true if error
bool verifyResult(libmarv::GpuDatabase& gpuDatabase, const libmarv::SearchResult& searchResult, std::string_view query, libmarv::ScanType scanType, int gop, int gex){
    if(scanType == libmarv::ScanType::Foldseek) return false;

    size_t numResults = searchResult.scores.size();

    //simple structural checks
    if(!std::all_of(searchResult.scores.begin(), searchResult.scores.end(), [](int score){return score >= 0; })){
        std::cout << "there are scores < 0\n";
        return true;
    }
    std::map<int,int> frequencyOfIds;
    for(auto refId : searchResult.referenceIds){
        frequencyOfIds[refId]++;
    }
    for(const auto& p : frequencyOfIds){
        if(p.second != 1){
            std::cout << "refid " << p.first << ", frequency " << p.second << "\n";
            return false;
        }
    }

    std::vector<char> convertedQuery(query.size());
    auto convertQuery = [&](auto c){
        libmarv::ConvertAA_20 charToMMseqs;
        libmarv::ConvertAA_20_mmseqs_to_ncbi mmseqsToNcbi;
        return mmseqsToNcbi(charToMMseqs(c));
    };
    std::transform(
        query.data(),
        query.data() + query.size(),
        convertedQuery.data(),
        convertQuery
    );

    // for(size_t i = 0; i < std::min(size_t(10), numResults); i++){        
    //     int gpuScore = searchResult.scores[i];
    //     auto referenceId = searchResult.referenceIds[i];
    //     std::cout << gpuScore << " " << referenceId << "\n";
    // }

    std::atomic_bool hasError{false};

    #pragma omp parallel for
    for(size_t i = 0; i < numResults; i++){
        if(!hasError){
            int gpuScore = searchResult.scores[i];
            auto referenceId = searchResult.referenceIds[i];
            int subjectEndExcl = 0;
            int queryEndExcl = 0;
            if(searchResult.endPositions.size() > 0){
                subjectEndExcl = searchResult.endPositions[i].getSubjectEndExclusive();
                queryEndExcl = searchResult.endPositions[i].getQueryEndExclusive();
            }

            const auto referenceLength = gpuDatabase.getReferenceLength(referenceId);
            const char* referencePtr = gpuDatabase.getReferenceDatabaseSequenceHostPtr(referenceId);

            std::vector<char> convertedSubject(referenceLength);
            std::transform(
                referencePtr,
                referencePtr + referenceLength,
                convertedSubject.data(),
                libmarv::ConvertAA_20_mmseqs_to_ncbi{}
            );

            int cpuScore = 0;
            std::vector<std::pair<int,int>> cpuEndpositions;

            if(scanType == libmarv::ScanType::Gapless){
                cpuScore = libmarv::GaplessFilter_host_protein_converted_blosum62(
                    convertedQuery.data(),
                    convertedSubject.data(),
                    convertedQuery.size(),
                    referenceLength
                );
            }else if(scanType == libmarv::ScanType::SW){
                cpuScore = libmarv::affine_local_DP_host_protein_blosum62_converted(
                    convertedQuery.data(),
                    convertedSubject.data(),
                    convertedQuery.size(),
                    referenceLength,
                    gop,
                    gex
                );
            }else if(scanType == libmarv::ScanType::Gapless_Endpos){
                std::tie(cpuScore, cpuEndpositions) = libmarv::GaplessFilter_host_protein_converted_blosum62_endPos(
                    convertedQuery.data(),
                    convertedSubject.data(),
                    convertedQuery.size(),
                    referenceLength
                );
            }else if(scanType == libmarv::ScanType::SW_Endpos || libmarv::ScanType::GaplessPlusSW_Endpos){
                std::tie(cpuScore, cpuEndpositions) = libmarv::affine_local_DP_host_protein_blosum62_converted_endPos(
                    convertedQuery.data(),
                    convertedSubject.data(),
                    convertedQuery.size(),
                    referenceLength,
                    gop,
                    gex
                );
            }else{
                throw std::runtime_error("not implemented\n");
            }



            #pragma omp critical
            {
                if(!hasError){
                    if(cpuScore < 2048){
                        if(cpuScore != gpuScore){
                            hasError = true;
                            std::cout << "error. i " << i << ", sequence id " << referenceId 
                                << ", cpu score " << cpuScore << ", gpu score " << gpuScore 
                                //<< ", exact cpu score " << cpuScoresExact[refId] 
                                << ".";
                            std::cout << "Query:\n";
                            std::copy(query.begin(), query.end(), std::ostream_iterator<char>(std::cout, ""));
                            std::cout << "\n";
                            std::cout << "db sequence:\n";
                            std::cout << gpuDatabase.getReferenceSequence(referenceId) << "\n";
                        }

                        if(cpuEndpositions.size() > 0){
                            auto iterator = std::find(cpuEndpositions.begin(), cpuEndpositions.end(), std::make_pair(queryEndExcl, subjectEndExcl));
                            if(iterator == cpuEndpositions.end()){
                                hasError = true;

                                std::cout << "error. i " << i << ", sequence id " << referenceId 
                                    << ", cpu score " << cpuScore << ", gpu score " << gpuScore 
                                    << ", kernel query/subject end pos (" << queryEndExcl << "," << subjectEndExcl << ")\n";
                                std::cout << "cpu end positions:\n";
                                for(const auto& pair : cpuEndpositions){
                                    std::cout << "(" << pair.first << "," << pair.second << "), ";
                                }
                                std::cout << "\n";
                                std::cout << "Query:\n";
                                std::copy(query.begin(), query.end(), std::ostream_iterator<char>(std::cout, ""));
                                std::cout << "\n";
                                std::cout << "db sequence:\n";
                                std::cout << gpuDatabase.getReferenceSequence(referenceId) << "\n";
                            }
                        }
                    }
                }
            }
        }
    }

    return hasError;
}



int main(int argc, char* argv[])
{
    ProgramOptions options;
    bool parseSuccess = parseArgs(argc, argv, options);

    if(!parseSuccess || options.help){
        printHelp(argc, argv);
        return 0;
    }

    printOptions(options);

    std::vector<int> deviceIds;
    {
        int num = 0;
        cudaGetDeviceCount(&num); CUERR
        for(int i = 0; i < num; i++){
            deviceIds.push_back(i);
        }
        if(deviceIds.size() > 0){
            if(options.verbose){
                std::cout << "Will use GPU";
                for(auto x : deviceIds){
                    std::cout << " " << x;
                }
                std::cout << "\n";
            }
        }else{
            throw std::runtime_error("No GPU found");
        }
    }

    helpers::PeerAccess peerAccess(deviceIds, false);
    peerAccess.enableAllPeerAccesses();
 
    using MemoryConfig = libmarv::GpuDatabase::MemoryConfig;


    MemoryConfig memoryConfig;
    memoryConfig.maxBatchBytes = options.maxBatchBytes;
    memoryConfig.maxBatchSequences = options.maxBatchSequences;
    memoryConfig.maxGpuMem = options.maxGpuMem;


    std::shared_ptr<libmarv::TargetSubjectIds> targetSubjectIds;
    if(options.subjectIdsFilename.has_value()){
        targetSubjectIds = std::make_shared<libmarv::TargetSubjectIds>(options.subjectIdsFilename.value());
        // for(auto x : targetSubjectIds->subjectIds){
        //     std::cout << x << ", ";
        // }
        // std::cout << "\n";
    }

    std::ofstream outputfile(options.outputfile);
    if(!bool(outputfile)){
        throw std::runtime_error("Cannot open file " + options.outputfile);
    }
    // if(options.outputMode == ProgramOptions::OutputMode::TSV){
    //     printTSVHeader(outputfile);
    // }

    int numTopOutputs = options.numTopOutputs;
    // numTopOutputs = std::min(numTopOutputs, libmarv::MaxNumberOfResults::value());
    // std::cout << "Will output up to " << numTopOutputs << " results\n";

    libmarv::GpuDatabase gpuDatabase(
        deviceIds, 
        memoryConfig, 
        options.verbose
    );

    
    // gpuDatabase.allowInt8(options.allowInt8);

    if(!options.usePseudoDB){
        if(options.verbose){
            std::cout << "Reading Database: \n";
        }
        try{
            helpers::CpuTimer timer_read_db("Read DB");
            constexpr bool writeAccess = false;
            const bool prefetchSeq = options.prefetchDBFile;
            auto fullDB_tmp = std::make_shared<libmarv::DB>(libmarv::loadDB(options.dbPrefix, writeAccess, prefetchSeq));
            if(options.verbose){
                timer_read_db.print();
            }

            gpuDatabase.setDatabase(fullDB_tmp);
        }catch(libmarv::LoadDBException& ex){
            if(options.verbose){
                std::cout << "Failed to map db files. Using fallback db. Error message: " << ex.what() << "\n";
            }
            helpers::CpuTimer timer_read_db("Read DB");
            auto fullDB_tmp = std::make_shared<libmarv::DBWithVectors>(libmarv::loadDBWithVectors(options.dbPrefix));
            if(options.verbose){
                timer_read_db.print();
            }

            gpuDatabase.setDatabase(fullDB_tmp);
        }
    }else{
        if(options.verbose){
            std::cout << "Generating pseudo db\n";
        }
        helpers::CpuTimer timer_read_db("Generate DB");
        auto fullDB_tmp = std::make_shared<libmarv::PseudoDB>(libmarv::loadPseudoDB(
            options.pseudoDBSize, 
            options.pseudoDBLength,
            options.pseudoDBSameSequence
        ));
        if(options.verbose){
            timer_read_db.print();
        }
        
        gpuDatabase.setDatabase(fullDB_tmp);
    }

    if(options.verbose){
        gpuDatabase.printDBInfo();
        if(options.printLengthPartitions){
            gpuDatabase.printDBLengthPartitions();
        }
    }

    if(options.loadFullDBToGpu){
        gpuDatabase.prefetchDBToGpus();
    }

    {

        for(const auto& queryFile : options.queryFiles){
            std::cout << "Processing query file " << queryFile << "\n";

            kseqpp::KseqPP reader(queryFile);
            int64_t query_num = 0;

            // cudaSW4.totalTimerStart();

            std::vector<std::string> queries;
            std::vector<std::string> headers;

            std::vector<libmarv::ReferenceIdT>* subsetIdsPtr = nullptr;
            if(targetSubjectIds){
                subsetIdsPtr = &(targetSubjectIds->subjectIds);
            }

            auto processQueries = [&](){
                std::cout << "Processing queries " << query_num << " - " << (query_num + queries.size() - 1) << " ...\n";
                std::vector<std::string_view> queryViews(queries.begin(), queries.end());

                auto searchResultPerQuery = [&](){
                    if(options.scanType == libmarv::ScanType::Gapless || options.scanType == libmarv::ScanType::SW){
                        return do_search_scoreOnly(
                            gpuDatabase, queryViews, numTopOutputs, options.scanType, 
                            options.gop, options.gex, options.allowInt8, options.handleOverflows,
                            subsetIdsPtr
                        );
                    }else if(options.scanType == libmarv::ScanType::Gapless_Endpos || options.scanType == libmarv::ScanType::SW_Endpos){
                        return do_search_scoreEndpos(
                            gpuDatabase, queryViews, numTopOutputs, options.scanType, 
                            options.gop, options.gex, options.allowInt8, options.handleOverflows,
                            subsetIdsPtr
                        );
                    }else if(options.scanType == libmarv::ScanType::GaplessPlusSW_Endpos){
                        return do_search_GaplessScoreOnly_and_SW_Endpos_for_tops(
                            gpuDatabase, queryViews, numTopOutputs, options.gop, 
                            options.gex, options.allowInt8, options.handleOverflows
                        );
                    }else if(options.scanType == libmarv::ScanType::Foldseek){
                        return do_search_scoreOnly_foldseek(
                            gpuDatabase, queryViews, numTopOutputs, 
                            libmarv::MaskingOptions(options.maskingThreshold, options.maskingLetter),
                            options.gop, 
                            options.gex, 
                            options.allowInt8, 
                            options.handleOverflows
                        );
                    }else{
                        throw std::runtime_error("cannot process with scan type");
                    }
                }();
                
                
                if(options.verbose && searchResultPerQuery[0].stats != nullptr){
                    std::cout << "Done. Scan time: " << searchResultPerQuery[0].stats->seconds << " s, " << searchResultPerQuery[0].stats->gcups << " GCUPS\n";
                }else{
                    std::cout << "Done.\n";
                }

                const int numQueries = queries.size();
                for(int queryIndex = 0; queryIndex < numQueries; queryIndex++){
                    const auto& searchResult = searchResultPerQuery[queryIndex];
                    if(options.checkResults){
                        std::cout << "checking results of query " << query_num << "\n";
                        bool hasError = verifyResult(gpuDatabase, searchResult, queries[queryIndex], options.scanType, options.gop, options.gex);
                        if(!hasError){
                            std::cout << "results look ok\n";
                        }
                    }else{

                        if(numTopOutputs > 0){
                            if(options.outputMode == ProgramOptions::OutputMode::Plain){
                                outputfile << "Query " << query_num << ", header " <<  headers[queryIndex]
                                    << ", length " << queries[queryIndex].size()
                                    << "\n";

                                printScanResultPlain(outputfile, searchResult, gpuDatabase);
                            }else if(options.outputMode == ProgramOptions::OutputMode::TSV){
                                printScanResultTSV(outputfile, searchResult, gpuDatabase, query_num, queries[queryIndex].size(), headers[queryIndex]);
                            }
                            outputfile.flush();
                        }
                    }
                    query_num++;
                }
            };

            while(reader.next() >= 0){
                const std::string& header = reader.getCurrentHeader();
                const std::string& sequence = reader.getCurrentSequence();

                headers.push_back(header);
                queries.push_back(sequence);

                const int numQueries = queries.size();

                if(numQueries >= options.queryBatchsize){
                    processQueries();
                    queries.clear();
                    headers.clear();
                }
            }

            if(queries.size() > 0){
                processQueries();
                queries.clear();
                headers.clear();
            }

        }
    }

}
