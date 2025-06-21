#include <algorithm>
#include <vector>
#include <iostream>

#include "options.hpp"
#include "cudasw4.cuh"

#include "hpc_helpers/all_helpers.cuh"
#include "hpc_helpers/peer_access.cuh"




struct BenchmarkData{
    bool dpx;
    int tilesize;
    int groupsize;
    int numRegs;
    float gcups;
    GaplessKernelConfig::Approach kernelApproach;
};

struct BenchmarkDataSW{
    bool dpx;
    int tilesize;
    int groupsize;
    int numRegs;
    float gcups;
    SmithWatermanKernelConfig::Approach kernelApproach;
};

void writeBenchmarkDataHeader(std::ostream& os){
    os << "tilesize groupsize numRegs dpx kernelApproach gcups" << "\n";
}

std::ostream& operator<<(std::ostream& os, const BenchmarkData& data){

    os << data.tilesize << " " << data.groupsize << " " << data.numRegs 
        << " " << data.dpx << " " << int(data.kernelApproach) << " " << data.gcups;
    return os;
}

std::ostream& operator<<(std::ostream& os, const BenchmarkDataSW& data){

    os << data.tilesize << " " << data.groupsize << " " << data.numRegs 
        << " " << data.dpx << " " << int(data.kernelApproach) << " " << data.gcups;
    return os;
}


void gapless_search(){
    std::cout << "gapless_search\n";

    ProgramOptions options;

    const int numTopOutputs = 0;
    const auto blosumType = cudasw4::BlosumType::BLOSUM62_20;
    const bool verbose = false;

    std::vector<int> deviceIds;
    {
        int num = 0;
        cudaGetDeviceCount(&num); CUERR
        for(int i = 0; i < num; i++){
            deviceIds.push_back(i);
        }
        if(deviceIds.size() > 0){
            if(verbose){
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
 
    using MemoryConfig = cudasw4::MemoryConfig;
    using ScanResult = cudasw4::ScanResult;

    MemoryConfig memoryConfig;
    memoryConfig.maxBatchBytes = options.maxBatchBytes;
    memoryConfig.maxBatchSequences = options.maxBatchSequences;
    memoryConfig.maxTempBytes = options.maxTempBytes;
    memoryConfig.maxGpuMem = options.maxGpuMem;

    const char* letters = "ARNDCQEGHILKMFPSTWYV";

    std::mt19937 gen(42);
    std::uniform_int_distribution<> dist(0,19);





    std::vector<BenchmarkData> allBenchmarkData;

    writeBenchmarkDataHeader(std::cout);

    auto execute = [&](GaplessKernelConfig::Approach kernelApproach, bool useDPX){
        KernelConfigFilenames kernelConfigFilenames;
        cudasw4::CudaSW4 cudaSW4(
            deviceIds, 
            numTopOutputs,
            blosumType, 
            memoryConfig, 
            verbose,
            kernelConfigFilenames
        );

        cudaSW4.setScanType(ScanType::Gapless);
        // cudaSW4.setScanType(ScanType::SW_Endpos);

        // std::ofstream logfile(outputfilename);

        std::vector<BenchmarkData> benchmarkDataVec;

        std::vector<std::tuple<int,int>> validRegConfigs;
        #define X(g,r)\
            validRegConfigs.push_back(std::make_tuple(g,r));
        
        PSSM_GAPLESS_SINGLETILE_FOR_EACH_VALID_CONFIG_DO_X

        #undef X



        for(auto regConfig : validRegConfigs){
            const int groupsize = std::get<0>(regConfig);
            const int numRegs = std::get<1>(regConfig);
        

        // for(int groupsize : {4,8,16}){
        // // for(int groupsize : {4}){
        // // for(int groupsize : {8,16}){
        //     for(int numRegs : {4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64}){
        //     // for(int groupsize : {4,8,16}){
        //     //     for(int numRegs : {32}){



        // // for(int groupsize : {8,16}){
        // //     for(int numRegs : {4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64}){


                const int l = groupsize * numRegs * 2;
                // if(l <= 2048) continue;
 

                BenchmarkData benchmarkData;
                benchmarkData.dpx = useDPX;
                benchmarkData.tilesize = l;
                benchmarkData.groupsize = groupsize;
                benchmarkData.numRegs = numRegs;
                benchmarkData.kernelApproach = kernelApproach;

                std::string sequence(l, ' ');
                for(int i = 0; i < l; i++){
                    sequence[i] = letters[dist(gen)];
                }

                const int pseudoDBLength = l;
                const int pseudoDBSize = 5000000;
                const bool pseudoDBSameSequence = false;

                helpers::CpuTimer timer_read_db("Generate DB");
                auto fullDB_tmp = std::make_shared<cudasw4::PseudoDB>(cudasw4::loadPseudoDB(
                    pseudoDBSize, 
                    pseudoDBLength,
                    pseudoDBSameSequence
                ));
                timer_read_db.stop();
                
                cudaSW4.setDatabase(fullDB_tmp);
                cudaSW4.prefetchDBToGpus();
                
                GaplessKernelConfig config;
                config.dpx = useDPX;
                config.tilesize = l;
                config.groupsize = groupsize;
                config.numRegs = numRegs;
                config.approach =kernelApproach;

                cudaSW4.setCustomKernelConfig_Gapless(config);
                cudasw4::DecodedQueryView queryView(sequence.data(), sequence.size());
                ScanResult scanResult = cudaSW4.scan(queryView, std::nullopt);

                benchmarkData.gcups = scanResult.stats.gcups;

                benchmarkDataVec.push_back(benchmarkData);

                std::cout << benchmarkData << "\n";

        }

        return benchmarkDataVec;
    };

    for(auto kernelApproach : {GaplessKernelConfig::Approach::hardcodedzero, GaplessKernelConfig::Approach::kernelparamzero}){
        auto resultNoDpx = execute(kernelApproach, false);
        allBenchmarkData.insert(allBenchmarkData.end(), resultNoDpx.begin(), resultNoDpx.end());

        int ccMajor = 0;
        cudaDeviceGetAttribute(&ccMajor, cudaDevAttrComputeCapabilityMajor, 0);
        const bool supportsDPX = ccMajor >= 9;
        if(supportsDPX){
            auto resultDpx = execute(kernelApproach, true);
            allBenchmarkData.insert(allBenchmarkData.end(), resultDpx.begin(), resultDpx.end());
        }
    }

    auto bestConfigs = allBenchmarkData;
    std::sort(bestConfigs.begin(), bestConfigs.end(), [](const auto& l, const auto& r){
        if(l.tilesize < r.tilesize) return true;
        if(l.tilesize > r.tilesize) return false;
        return l.gcups > r.gcups;
    });

    std::cout << "sorted\n";
    std::copy(bestConfigs.begin(), bestConfigs.end(), std::ostream_iterator<BenchmarkData>(std::cout, "\n"));

    //only keep best for each tilesize
    bestConfigs.erase(
        std::unique(bestConfigs.begin(), bestConfigs.end(), [](const auto& l, const auto& r){
            return l.tilesize == r.tilesize;
        }),
        bestConfigs.end()
    );

    std::cout << "best\n";
    std::copy(bestConfigs.begin(), bestConfigs.end(), std::ostream_iterator<BenchmarkData>(std::cout, "\n"));
}









void sw_search(){
    std::cout << "sw_search\n";
    ProgramOptions options;

    const int numTopOutputs = 0;
    const auto blosumType = cudasw4::BlosumType::BLOSUM62_20;
    const bool verbose = false;

    std::vector<int> deviceIds;
    {
        int num = 0;
        cudaGetDeviceCount(&num); CUERR
        for(int i = 0; i < num; i++){
            deviceIds.push_back(i);
        }
        if(deviceIds.size() > 0){
            if(verbose){
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

    using MemoryConfig = cudasw4::MemoryConfig;
    using ScanResult = cudasw4::ScanResult;

    MemoryConfig memoryConfig;
    memoryConfig.maxBatchBytes = options.maxBatchBytes;
    memoryConfig.maxBatchSequences = options.maxBatchSequences;
    memoryConfig.maxTempBytes = options.maxTempBytes;
    memoryConfig.maxGpuMem = options.maxGpuMem;

    const char* letters = "ARNDCQEGHILKMFPSTWYV";

    std::mt19937 gen(42);
    std::uniform_int_distribution<> dist(0,19);





    std::vector<BenchmarkDataSW> allBenchmarkData;

    writeBenchmarkDataHeader(std::cout);

    auto execute = [&](bool useDPX){

        KernelConfigFilenames kernelConfigFilenames;
        cudasw4::CudaSW4 cudaSW4(
            deviceIds, 
            numTopOutputs,
            blosumType, 
            memoryConfig, 
            verbose,
            kernelConfigFilenames
        );

        cudaSW4.setScanType(ScanType::SW_Endpos);

        // std::ofstream logfile(outputfilename);

        std::vector<BenchmarkDataSW> benchmarkDataVec;

        std::vector<std::tuple<int,int>> validRegConfigs;
        #define X(g,r)\
            validRegConfigs.push_back(std::make_tuple(g,r));
        
        PSSM_SW_ENDPOS_SINGLETILE_FLOAT_OR_INT_FOR_EACH_VALID_CONFIG_DO_X

        #undef X



        for(auto regConfig : validRegConfigs){
            const int groupsize = std::get<0>(regConfig);
            const int numRegs = std::get<1>(regConfig);

            const int l = groupsize * numRegs;
            // if(l <= 2048) continue;


            BenchmarkDataSW benchmarkData;
            benchmarkData.dpx = useDPX;
            benchmarkData.tilesize = l;
            benchmarkData.groupsize = groupsize;
            benchmarkData.numRegs = numRegs;
            benchmarkData.kernelApproach = SmithWatermanKernelConfig::Approach::Unused;

            std::string sequence(l, ' ');
            for(int i = 0; i < l; i++){
                sequence[i] = letters[dist(gen)];
            }

            const int pseudoDBLength = l;
            const int pseudoDBSize = 5000000;
            const bool pseudoDBSameSequence = false;

            helpers::CpuTimer timer_read_db("Generate DB");
            auto fullDB_tmp = std::make_shared<cudasw4::PseudoDB>(cudasw4::loadPseudoDB(
                pseudoDBSize, 
                pseudoDBLength,
                pseudoDBSameSequence
            ));
            timer_read_db.stop();
            
            cudaSW4.setDatabase(fullDB_tmp);
            cudaSW4.prefetchDBToGpus();
            

            SmithWatermanKernelConfig config;
            config.dpx = useDPX;
            config.tilesize = l;
            config.groupsize = groupsize;
            config.numRegs = numRegs;
            config.approach = SmithWatermanKernelConfig::Approach::Unused;

            cudaSW4.setCustomKernelConfig_SW(config);
            cudasw4::DecodedQueryView queryView(sequence.data(), sequence.size());
            ScanResult scanResult = cudaSW4.scan(queryView, std::nullopt);

            benchmarkData.gcups = scanResult.stats.gcups;

            benchmarkDataVec.push_back(benchmarkData);

            std::cout << benchmarkData << "\n";
        }

        return benchmarkDataVec;
    };

    {
        auto resultNoDpx = execute(false);
        allBenchmarkData.insert(allBenchmarkData.end(), resultNoDpx.begin(), resultNoDpx.end());

        int ccMajor = 0;
        cudaDeviceGetAttribute(&ccMajor, cudaDevAttrComputeCapabilityMajor, 0);
        const bool supportsDPX = ccMajor >= 9;
        if(supportsDPX){
            auto resultDpx = execute(true);
            allBenchmarkData.insert(allBenchmarkData.end(), resultDpx.begin(), resultDpx.end());
        }
    }

    auto bestConfigs = allBenchmarkData;
    std::sort(bestConfigs.begin(), bestConfigs.end(), [](const auto& l, const auto& r){
        if(l.tilesize < r.tilesize) return true;
        if(l.tilesize > r.tilesize) return false;
        return l.gcups > r.gcups;
    });

    std::cout << "sorted\n";
    std::copy(bestConfigs.begin(), bestConfigs.end(), std::ostream_iterator<BenchmarkDataSW>(std::cout, "\n"));

    //only keep best for each tilesize
    bestConfigs.erase(
        std::unique(bestConfigs.begin(), bestConfigs.end(), [](const auto& l, const auto& r){
            return l.tilesize == r.tilesize;
        }),
        bestConfigs.end()
    );

    std::cout << "best\n";
    std::copy(bestConfigs.begin(), bestConfigs.end(), std::ostream_iterator<BenchmarkDataSW>(std::cout, "\n"));
}




int main(int argc, char* argv[]){

    bool gapless = false;
    bool sw = false;
    for(int x = 1; x < argc; x++){
        std::string argstring = argv[x];
        if(argstring == "--gapless"){
            gapless = true;
        }
        if(argstring == "--sw"){
            sw = true;
        }
    }

    if(gapless){
        gapless_search();
    }

    if(sw){
        sw_search();
    }

    return 0;
}

