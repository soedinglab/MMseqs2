
#include "benchmarkdata.hpp"

#include "../cuda_errorcheck.cuh"

#include "common.hpp"

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>


std::string to_string(AlignmentResultKind e){
    switch(e){
        case AlignmentResultKind::ScoreOnly: return "ScoreOnly";
        case AlignmentResultKind::ScoreEndpos: return "ScoreEndpos";
        default: return "Forgot to name AlignmentResultKind!";
    }
}

std::ostream& operator<<(std::ostream& os, const BenchmarkData& data){

    os << data.blocksize << " " << data.tilesize << " " << data.groupsize << " " << data.numItems 
        << " " << data.mingcups << " " << data.avggcups << " " << data.maxgcups
        << " " << data.datatypeString << " " << data.approachString << " " << data.resultkindString << " " << data.numTiles;
    return os;
}


void printBenchmarkData_full(std::vector<BenchmarkData> vec, char sep){
    //sort results for same tile size by average gcups
    std::sort(vec.begin(), vec.end(), [](const auto& l, const auto& r){
        if(l.tilesize < r.tilesize) return true;
        if(l.tilesize > r.tilesize) return false;
        if(l.numTiles < r.numTiles) return true;
        if(l.numTiles > r.numTiles) return false;
        return l.maxgcups > r.maxgcups;
    });

    // vec.erase(
    //     std::unique(vec.begin(), vec.end(), [](const auto& l, const auto& r){
    //         return l.tilesize == r.tilesize;
    //     }),
    //     vec.end()
    // );
    
    int deviceId = 0;
    CUDACHECK(cudaGetDevice(&deviceId));
    int smMajor = 0;
    int smMinor = 0;
    CUDACHECK(cudaDeviceGetAttribute(&smMajor, cudaDevAttrComputeCapabilityMajor, deviceId));
    CUDACHECK(cudaDeviceGetAttribute(&smMinor, cudaDevAttrComputeCapabilityMinor, deviceId));

    const int gpuArch = smMajor * 100 + smMinor * 10;

    std::cout 
            << "gpuArch"
            << sep << "alphabetSize"
            << sep << "blocksize"
            << sep << "tilesize"
            << sep << "groupsize"
            << sep << "numItems"
            << sep << "mingcups"
            << sep << "avggcups"
            << sep << "maxgcups"
            << sep << "datatype"
            << sep << "approach"
            << sep << "resultkind"
            << sep << "numTiles"
            << '\n';

    for(const auto& config : vec){
        std::cout 
            << gpuArch
            << sep << alphabetSize
            << sep << config.blocksize
            << sep << config.tilesize
            << sep << config.groupsize
            << sep << config.numItems
            << sep << config.mingcups
            << sep << config.avggcups
            << sep << config.maxgcups
            << sep << config.datatypeString
            << sep << config.approachString
            << sep << config.resultkindString
            << sep << config.numTiles
            << '\n';
    }
}

void printBenchmarkData_bestPerTile(std::vector<BenchmarkData> vec, char sep){
    //sort results for same tile size by average gcups
    std::sort(vec.begin(), vec.end(), [](const auto& l, const auto& r){
        if(l.tilesize < r.tilesize) return true;
        if(l.tilesize > r.tilesize) return false;
        if(l.numTiles < r.numTiles) return true;
        if(l.numTiles > r.numTiles) return false;
        return l.maxgcups > r.maxgcups;
    });

    vec.erase(
        std::unique(vec.begin(), vec.end(), [](const auto& l, const auto& r){
            return l.tilesize == r.tilesize && l.numTiles == r.numTiles;
        }),
        vec.end()
    );

    vec.erase(
        std::remove_if(vec.begin(), vec.end(), [](const auto& l){
            return l.numTiles == 0;
        }),
        vec.end()
    );
    
    int deviceId = 0;
    CUDACHECK(cudaGetDevice(&deviceId));
    int smMajor = 0;
    int smMinor = 0;
    CUDACHECK(cudaDeviceGetAttribute(&smMajor, cudaDevAttrComputeCapabilityMajor, deviceId));
    CUDACHECK(cudaDeviceGetAttribute(&smMinor, cudaDevAttrComputeCapabilityMinor, deviceId));

    const int gpuArch = smMajor * 100 + smMinor * 10;

    std::cout 
            << "gpuArch"
            << sep << "alphabetSize"
            << sep << "blocksize"
            << sep << "tilesize"
            << sep << "groupsize"
            << sep << "numItems"
            << sep << "mingcups"
            << sep << "avggcups"
            << sep << "maxgcups"
            << sep << "datatype"
            << sep << "approach"
            << sep << "resultkind"
            << sep << "numTiles"
            << '\n';

    for(const auto& config : vec){
        std::cout 
            << gpuArch
            << sep << alphabetSize
            << sep << config.blocksize
            << sep << config.tilesize
            << sep << config.groupsize
            << sep << config.numItems
            << sep << config.mingcups
            << sep << config.avggcups
            << sep << config.maxgcups
            << sep << config.datatypeString
            << sep << config.approachString
            << sep << config.resultkindString
            << sep << config.numTiles
            << '\n';
    }
}

void printBestAsConfig(std::vector<BenchmarkData> vec){
    //sort results for same tile size by average gcups
    std::sort(vec.begin(), vec.end(), [](const auto& l, const auto& r){
        if(l.tilesize < r.tilesize) return true;
        if(l.tilesize > r.tilesize) return false;
        if(l.numTiles < r.numTiles) return true;
        if(l.numTiles > r.numTiles) return false;
        return l.maxgcups > r.maxgcups;
    });

    vec.erase(
        std::unique(vec.begin(), vec.end(), [](const auto& l, const auto& r){
            return l.tilesize == r.tilesize && l.numTiles == r.numTiles;
        }),
        vec.end()
    );

    vec.erase(
        std::remove_if(vec.begin(), vec.end(), [](const auto& l){
            return l.numTiles == 0;
        }),
        vec.end()
    );

    std::cout << "using type = std::tuple<\n";
    for(const auto& config : vec){
        std::cout << "\t" << "GaplessConfig<" << alphabetSize 
            << ", Approach::" << config.approachString
            << ", Datatype::" << config.datatypeString
            << ", " << config.blocksize
            << ", " << config.groupsize
            << ", " << config.numItems
            //<< ", " << int(config.maxgcups)
            << ">,\n";
    }
    std::cout << ">;\n";

    std::cout << "std::array<int, " << vec.size() << ">{";
    for(const auto& config : vec){
        std::cout << int(config.maxgcups) << ",";
    }
    std::cout << "};\n";
    
}