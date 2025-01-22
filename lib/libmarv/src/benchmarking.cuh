int peakBenchmark(int argc, char* argv[]){
     // {
    //     std::string tempname = "query128.fasta";
    //     std::ofstream tempfile(tempname);

    //     const char* letters = "ARNDCQEGHILKMFPSTWYV";

    //     std::mt19937 gen(42);
    //     std::uniform_int_distribution<> dist(0,19);

    //     //for(int l = 8193; l < 4*4096; l += 64){
    //     {
    //         int l=128;
    //         std::string sequence(l, ' ');
    //         for(int i = 0; i < l; i++){
    //             sequence[i] = letters[dist(gen)];
    //         }
    //         tempfile << ">length " << l << "\n";
    //         tempfile << sequence << "\n";
    //     }
    //     tempfile.flush();

    //     options.queryFiles = {tempname};
    // }
    

    ProgramOptions options;
    bool parseSuccess = parseArgs(argc, argv, options);

    if(!parseSuccess || options.help){
        printHelp(argc, argv);
        return 0;
    }

    options.usePseudoDB = true;
    options.verbose = true;
    options.loadFullDBToGpu = true;
    options.numTopOutputs = 0;
    options.verbose = false;

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
 
    using KernelTypeConfig = cudasw4::KernelTypeConfig;
    using MemoryConfig = cudasw4::MemoryConfig;
    using ScanResult = cudasw4::ScanResult;

    KernelTypeConfig kernelTypeConfig;
    kernelTypeConfig.singlePassType = options.singlePassType;
    kernelTypeConfig.manyPassType_small = options.manyPassType_small;
    kernelTypeConfig.manyPassType_large = options.manyPassType_large;
    kernelTypeConfig.overflowType = options.overflowType;

    MemoryConfig memoryConfig;
    memoryConfig.maxBatchBytes = options.maxBatchBytes;
    memoryConfig.maxBatchSequences = options.maxBatchSequences;
    memoryConfig.maxTempBytes = options.maxTempBytes;
    memoryConfig.maxGpuMem = options.maxGpuMem;

    cudasw4::CudaSW4 cudaSW4(
        deviceIds, 
        options.numTopOutputs,
        options.blosumType, 
        kernelTypeConfig, 
        memoryConfig, 
        options.verbose
    );

    const char* letters = "ARNDCQEGHILKMFPSTWYV";

    std::mt19937 gen(42);
    std::uniform_int_distribution<> dist(0,19);

    for(int l = 1; l <= 65536; l *= 2){
        std::string sequence(l, ' ');
        for(int i = 0; i < l; i++){
            sequence[i] = letters[dist(gen)];
        }

        options.pseudoDBLength = l;
        options.pseudoDBSize = 1024*1024*1024 / l;

        if(options.verbose){
            std::cout << "Generating pseudo db\n";
        }
        helpers::CpuTimer timer_read_db("Generate DB");
        auto fullDB_tmp = std::make_shared<cudasw4::PseudoDB>(cudasw4::loadPseudoDB(
            options.pseudoDBSize, 
            options.pseudoDBLength,
            options.pseudoDBSameSequence
        ));
        if(options.verbose){
            timer_read_db.print();
        }
        
        cudaSW4.setDatabase(fullDB_tmp);


        if(options.verbose){
            cudaSW4.printDBInfo();
            if(options.printLengthPartitions){
                cudaSW4.printDBLengthPartitions();
            }
        }

        if(options.loadFullDBToGpu){
            cudaSW4.prefetchDBToGpus();
        }

        std::cout << "Processing query length " << sequence.size() << "\n";

        cudaSW4.totalTimerStart();

        std::cout.flush();

        cudasw4::DecodedQueryView queryView(sequence.data(), sequence.size());
        ScanResult scanResult = cudaSW4.scan(queryView, std::nullopt);

        std::cout << "Done. Scan time: " << scanResult.stats.seconds << " s, " << scanResult.stats.gcups << " GCUPS\n";

    }

    return 0;
}










int gridsearchPseudo(int argc, char* argv[]){
    ProgramOptions options;
    bool parseSuccess = parseArgs(argc, argv, options);

    if(!parseSuccess || options.help){
        printHelp(argc, argv);
        return 0;
    }

    options.usePseudoDB = true;
    options.pseudoDBSameSequence = false;
    options.verbose = true;
    options.loadFullDBToGpu = true;
    options.numTopOutputs = 0;
    options.verbose = false;

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
 
    using KernelTypeConfig = cudasw4::KernelTypeConfig;
    using MemoryConfig = cudasw4::MemoryConfig;
    using ScanResult = cudasw4::ScanResult;

    KernelTypeConfig kernelTypeConfig;
    kernelTypeConfig.singlePassType = options.singlePassType;
    kernelTypeConfig.manyPassType_small = options.manyPassType_small;
    kernelTypeConfig.manyPassType_large = options.manyPassType_large;
    kernelTypeConfig.overflowType = options.overflowType;

    MemoryConfig memoryConfig;
    memoryConfig.maxBatchBytes = options.maxBatchBytes;
    memoryConfig.maxBatchSequences = options.maxBatchSequences;
    memoryConfig.maxTempBytes = options.maxTempBytes;
    memoryConfig.maxGpuMem = options.maxGpuMem;

    cudasw4::CudaSW4 cudaSW4(
        deviceIds, 
        options.numTopOutputs,
        options.blosumType, 
        kernelTypeConfig, 
        memoryConfig, 
        options.verbose
    );

    const char* letters = "ARNDCQEGHILKMFPSTWYV";

    std::mt19937 gen(42);
    std::uniform_int_distribution<> dist(0,19);

    auto execute = [&](std::string outputfilename, bool usekernelparamzero){

        std::ofstream logfile(outputfilename);

        for(int l = 128; l <= 384; l += 8){
            std::string sequence(l, ' ');
            for(int i = 0; i < l; i++){
                sequence[i] = letters[dist(gen)];
            }

            options.pseudoDBLength = l;
            options.pseudoDBSize = 5000000;

            if(options.verbose){
                std::cout << "Generating pseudo db\n";
            }
            helpers::CpuTimer timer_read_db("Generate DB");
            auto fullDB_tmp = std::make_shared<cudasw4::PseudoDB>(cudasw4::loadPseudoDB(
                options.pseudoDBSize, 
                options.pseudoDBLength,
                options.pseudoDBSameSequence
            ));
            if(options.verbose){
                timer_read_db.print();
            }
            
            cudaSW4.setDatabase(fullDB_tmp);
            if(usekernelparamzero){
                cudaSW4.setGaplessHalf2KernelApproach(KernelApproach::kernelparamzero);
                cudaSW4.setGaplessDPXKernelApproach(KernelApproach::kernelparamzero);
            }else{
                cudaSW4.setGaplessHalf2KernelApproach(KernelApproach::hardcodedzero);
                cudaSW4.setGaplessDPXKernelApproach(KernelApproach::hardcodedzero);
            }

            if(options.verbose){
                cudaSW4.printDBInfo();
                if(options.printLengthPartitions){
                    cudaSW4.printDBLengthPartitions();
                }
            }

            if(options.loadFullDBToGpu){
                cudaSW4.prefetchDBToGpus();
            }

            //std::cout << "Processing query length " << sequence.size() << "\n";
            std::cout << "Processing query length " << sequence.size();

            logfile << sequence.size() << "\n";

            constexpr int smallestValidGroupSize = 4;

            std::vector<std::vector<float>> gcupsMatrix;
            for(int groupsize = smallestValidGroupSize; groupsize <= 16; groupsize *= 2){
                std::cout.flush();

                std::vector<float> gcupsVector;
                for(int numRegs = 4; numRegs <= 64; numRegs += 4){
                    if(groupsize * numRegs * 2 >= int(sequence.size())){
                    //if(groupsize >= 4){
                        cudaSW4.setGroupConfig(groupsize, numRegs);

                        cudasw4::DecodedQueryView queryView(sequence.data(), sequence.size());
                        ScanResult scanResult = cudaSW4.scan(queryView, std::nullopt);
                        gcupsVector.push_back(scanResult.stats.gcups);
                    }else{
                        gcupsVector.push_back(0);
                    }
                }
                gcupsMatrix.push_back(gcupsVector);
            }



            float bestGcups = 0;
            int bestgroupsizeIndex = 0;
            int bestregsIndex = 0;
            for(int r = 0; r < int(gcupsMatrix.size()); r++){
                const auto& vec = gcupsMatrix[r];
                for(int c = 0; c < int(vec.size()); c++){
                    logfile << vec[c] << " ";
                    //std::cout << vec[c] << " ";
                    //printf("%5f ", vec[c]);
                    if(vec[c] > bestGcups){
                        bestgroupsizeIndex = r;
                        bestregsIndex = c;
                        bestGcups = vec[c];
                    }
                }
                logfile << "\n";
                //printf("\n");
                //std::cout << "\n";
            }
            logfile << "Best: " << bestGcups << " " << smallestValidGroupSize * (1u << bestgroupsizeIndex) << " " << (4 * (1+bestregsIndex)) << "\n"; 
            std::cout << ", Best: " << bestGcups << " " << smallestValidGroupSize * (1u << bestgroupsizeIndex) << " " << (4 * (1+bestregsIndex)) << "\n";
            //<< " " << bestgroupsizeIndex << " " << bestregsIndex << "\n";

        }

    };

    execute("gridsearch_128_384_8_alldifferentsubjects_withkernelparamzero.txt", true);
    execute("gridsearch_128_384_8_alldifferentsubjects_withoutkernelparamzero.txt", false);

    return 0;
}


int gridsearchPseudo_SW(int argc, char* argv[]){
    ProgramOptions options;
    bool parseSuccess = parseArgs(argc, argv, options);

    if(!parseSuccess || options.help){
        printHelp(argc, argv);
        return 0;
    }

    options.usePseudoDB = true;
    options.pseudoDBSameSequence = false;
    options.verbose = true;
    options.loadFullDBToGpu = true;
    options.numTopOutputs = 0;
    options.verbose = false;

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
 
    using KernelTypeConfig = cudasw4::KernelTypeConfig;
    using MemoryConfig = cudasw4::MemoryConfig;
    using ScanResult = cudasw4::ScanResult;
    using ScanType = cudasw4::ScanType;

    KernelTypeConfig kernelTypeConfig;
    kernelTypeConfig.singlePassType = options.singlePassType;
    kernelTypeConfig.manyPassType_small = options.manyPassType_small;
    kernelTypeConfig.manyPassType_large = options.manyPassType_large;
    kernelTypeConfig.overflowType = options.overflowType;

    MemoryConfig memoryConfig;
    memoryConfig.maxBatchBytes = options.maxBatchBytes;
    memoryConfig.maxBatchSequences = options.maxBatchSequences;
    memoryConfig.maxTempBytes = options.maxTempBytes;
    memoryConfig.maxGpuMem = options.maxGpuMem;

    cudasw4::CudaSW4 cudaSW4(
        deviceIds, 
        options.numTopOutputs,
        options.blosumType, 
        kernelTypeConfig, 
        memoryConfig, 
        options.verbose
    );

    cudaSW4.setScanType(ScanType::SW_Endpos);

    const char* letters = "ARNDCQEGHILKMFPSTWYV";

    std::mt19937 gen(42);
    std::uniform_int_distribution<> dist(0,19);

    std::ofstream logfile("gridsearch_SW.txt");

    for(int l = 1232+16; l <= 2048; l += 16){
        std::string sequence(l, ' ');
        for(int i = 0; i < l; i++){
            sequence[i] = letters[dist(gen)];
        }

        options.pseudoDBLength = l;
        options.pseudoDBSize = 1024*1024*1024 / l;

        if(options.verbose){
            std::cout << "Generating pseudo db\n";
        }
        helpers::CpuTimer timer_read_db("Generate DB");
        auto fullDB_tmp = std::make_shared<cudasw4::PseudoDB>(cudasw4::loadPseudoDB(
            options.pseudoDBSize, 
            options.pseudoDBLength,
            options.pseudoDBSameSequence
        ));
        if(options.verbose){
            timer_read_db.print();
        }
        
        cudaSW4.setDatabase(fullDB_tmp);


        if(options.verbose){
            cudaSW4.printDBInfo();
            if(options.printLengthPartitions){
                cudaSW4.printDBLengthPartitions();
            }
        }

        if(options.loadFullDBToGpu){
            cudaSW4.prefetchDBToGpus();
        }

        //std::cout << "Processing query length " << sequence.size() << "\n";
        std::cout << "Processing query length " << sequence.size();

        logfile << sequence.size() << "\n";

        std::vector<std::vector<float>> gcupsMatrix;
        for(int groupsize = 4; groupsize <= 32; groupsize *= 2){
            std::cout.flush();

            std::vector<float> gcupsVector;
            for(int numRegs = 4; numRegs <= 64; numRegs += 4){
                if(groupsize * numRegs <= 1024){
                //if(groupsize * numRegs >= int(sequence.size()) && groupsize * numRegs <= 1024){
                //if(groupsize >= 4){
                    cudaSW4.setGroupConfig(groupsize, numRegs);

                    cudasw4::DecodedQueryView queryView(sequence.data(), sequence.size());
                    ScanResult scanResult = cudaSW4.scan(queryView, std::nullopt);
                    gcupsVector.push_back(scanResult.stats.gcups);
                }else{
                    gcupsVector.push_back(0);
                }
            }
            gcupsMatrix.push_back(gcupsVector);
        }



        float bestGcups = 0;
        int bestgroupsizeIndex = 0;
        int bestregsIndex = 0;
        for(int r = 0; r < int(gcupsMatrix.size()); r++){
            const auto& vec = gcupsMatrix[r];
            for(int c = 0; c < int(vec.size()); c++){
                logfile << vec[c] << " ";
                //std::cout << vec[c] << " ";
                //printf("%5f ", vec[c]);
                if(vec[c] > bestGcups){
                    bestgroupsizeIndex = r;
                    bestregsIndex = c;
                    bestGcups = vec[c];
                }
            }
            logfile << "\n";
            //printf("\n");
            //std::cout << "\n";
        }
        logfile << "Best: " << bestGcups << " " << (1u << bestgroupsizeIndex) << " " << (4 * (1+bestregsIndex)) << "\n"; 
        std::cout << ", Best: " << bestGcups << " " << (1u << bestgroupsizeIndex) << " " << (4 * (1+bestregsIndex)) << "\n";
        //<< " " << bestgroupsizeIndex << " " << bestregsIndex << "\n";

    }

    return 0;
}




int gridsearchReal(int argc, char* argv[]){
    ProgramOptions options;
    bool parseSuccess = parseArgs(argc, argv, options);

    if(!parseSuccess || options.help){
        printHelp(argc, argv);
        return 0;
    }

    options.usePseudoDB = false;
    options.pseudoDBSameSequence = false;
    options.verbose = true;
    options.loadFullDBToGpu = true;
    options.numTopOutputs = 0;
    options.verbose = false;

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
 
    using KernelTypeConfig = cudasw4::KernelTypeConfig;
    using MemoryConfig = cudasw4::MemoryConfig;
    using ScanResult = cudasw4::ScanResult;

    KernelTypeConfig kernelTypeConfig;
    kernelTypeConfig.singlePassType = options.singlePassType;
    kernelTypeConfig.manyPassType_small = options.manyPassType_small;
    kernelTypeConfig.manyPassType_large = options.manyPassType_large;
    kernelTypeConfig.overflowType = options.overflowType;

    MemoryConfig memoryConfig;
    memoryConfig.maxBatchBytes = options.maxBatchBytes;
    memoryConfig.maxBatchSequences = options.maxBatchSequences;
    memoryConfig.maxTempBytes = options.maxTempBytes;
    memoryConfig.maxGpuMem = options.maxGpuMem;

    cudasw4::CudaSW4 cudaSW4(
        deviceIds, 
        options.numTopOutputs,
        options.blosumType, 
        kernelTypeConfig, 
        memoryConfig, 
        options.verbose
    );

    if(options.verbose){
        std::cout << "Reading Database: \n";
    }
    try{
        helpers::CpuTimer timer_read_db("Read DB");
        constexpr bool writeAccess = false;
        const bool prefetchSeq = options.prefetchDBFile;
        auto fullDB_tmp = std::make_shared<cudasw4::DB>(cudasw4::loadDB(options.dbPrefix, writeAccess, prefetchSeq));
        if(options.verbose){
            timer_read_db.print();
        }

        cudaSW4.setDatabase(fullDB_tmp);
    }catch(cudasw4::LoadDBException& ex){
        if(options.verbose){
            std::cout << "Failed to map db files. Using fallback db. Error message: " << ex.what() << "\n";
        }
        helpers::CpuTimer timer_read_db("Read DB");
        auto fullDB_tmp = std::make_shared<cudasw4::DBWithVectors>(cudasw4::loadDBWithVectors(options.dbPrefix));
        if(options.verbose){
            timer_read_db.print();
        }

        cudaSW4.setDatabase(fullDB_tmp);
    }

    const char* letters = "ARNDCQEGHILKMFPSTWYV";

    std::mt19937 gen(42);
    std::uniform_int_distribution<> dist(0,19);

    std::ofstream logfile("log1.txt", std::ios::app);

    for(int l = 8192+4096; l <= 65536; l += 4096){
        if(l == 16384 || l == 32768 || l == 65536) continue;

        std::string sequence(l, ' ');
        for(int i = 0; i < l; i++){
            sequence[i] = letters[dist(gen)];
        }

        if(options.verbose){
            cudaSW4.printDBInfo();
            if(options.printLengthPartitions){
                cudaSW4.printDBLengthPartitions();
            }
        }

        if(options.loadFullDBToGpu){
            cudaSW4.prefetchDBToGpus();
        }

        //std::cout << "Processing query length " << sequence.size() << "\n";
        std::cout << "Processing query length " << sequence.size();

        logfile << sequence.size() << "\n";

        std::vector<std::vector<float>> gcupsMatrix;
        for(int groupsize = 1; groupsize <= 16; groupsize *= 2){
            std::cout.flush();

            std::vector<float> gcupsVector;
            for(int numRegs = 4; numRegs <= 64; numRegs += 4){
                //if(groupsize * numRegs * 2 >= int(sequence.size())){
                if(groupsize >= 8){
                    cudaSW4.setGroupConfig(groupsize, numRegs);

                    cudasw4::DecodedQueryView queryView(sequence.data(), sequence.size());
                    ScanResult scanResult = cudaSW4.scan(queryView, std::nullopt);
                    gcupsVector.push_back(scanResult.stats.gcups);
                }else{
                    gcupsVector.push_back(0);
                }
            }
            gcupsMatrix.push_back(gcupsVector);
        }



        float bestGcups = 0;
        int bestgroupsizeIndex = 0;
        int bestregsIndex = 0;
        for(int r = 0; r < int(gcupsMatrix.size()); r++){
            const auto& vec = gcupsMatrix[r];
            for(int c = 0; c < int(vec.size()); c++){
                logfile << vec[c] << " ";
                //std::cout << vec[c] << " ";
                //printf("%5f ", vec[c]);
                if(vec[c] > bestGcups){
                    bestgroupsizeIndex = r;
                    bestregsIndex = c;
                    bestGcups = vec[c];
                }
            }
            logfile << "\n";
            //printf("\n");
            //std::cout << "\n";
        }
        logfile << "Best: " << bestGcups << " " << (1u << bestgroupsizeIndex) << " " << (4 * (1+bestregsIndex)) << "\n"; 
        std::cout << ", Best: " << bestGcups << " " << (1u << bestgroupsizeIndex) << " " << (4 * (1+bestregsIndex)) << "\n";
        //<< " " << bestgroupsizeIndex << " " << bestregsIndex << "\n";

    }

    return 0;
}












int lengthbenchmarkReal(int argc, char* argv[], int firstLength, int lastLength, int stepLength){
    ProgramOptions options;
    bool parseSuccess = parseArgs(argc, argv, options);

    if(!parseSuccess || options.help){
        printHelp(argc, argv);
        return 0;
    }

    options.usePseudoDB = false;
    options.pseudoDBSameSequence = false;
    options.verbose = true;
    options.loadFullDBToGpu = true;
    options.numTopOutputs = 0;
    options.verbose = false;

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
 
    using KernelTypeConfig = cudasw4::KernelTypeConfig;
    using MemoryConfig = cudasw4::MemoryConfig;
    using ScanResult = cudasw4::ScanResult;

    KernelTypeConfig kernelTypeConfig;
    kernelTypeConfig.singlePassType = options.singlePassType;
    kernelTypeConfig.manyPassType_small = options.manyPassType_small;
    kernelTypeConfig.manyPassType_large = options.manyPassType_large;
    kernelTypeConfig.overflowType = options.overflowType;

    MemoryConfig memoryConfig;
    memoryConfig.maxBatchBytes = options.maxBatchBytes;
    memoryConfig.maxBatchSequences = options.maxBatchSequences;
    memoryConfig.maxTempBytes = options.maxTempBytes;
    memoryConfig.maxGpuMem = options.maxGpuMem;

    cudasw4::CudaSW4 cudaSW4(
        deviceIds, 
        options.numTopOutputs,
        options.blosumType, 
        kernelTypeConfig, 
        memoryConfig, 
        options.verbose
    );

    if(options.verbose){
        std::cout << "Reading Database: \n";
    }
    try{
        helpers::CpuTimer timer_read_db("Read DB");
        constexpr bool writeAccess = false;
        const bool prefetchSeq = options.prefetchDBFile;
        auto fullDB_tmp = std::make_shared<cudasw4::DB>(cudasw4::loadDB(options.dbPrefix, writeAccess, prefetchSeq));
        if(options.verbose){
            timer_read_db.print();
        }

        cudaSW4.setDatabase(fullDB_tmp);
    }catch(cudasw4::LoadDBException& ex){
        if(options.verbose){
            std::cout << "Failed to map db files. Using fallback db. Error message: " << ex.what() << "\n";
        }
        helpers::CpuTimer timer_read_db("Read DB");
        auto fullDB_tmp = std::make_shared<cudasw4::DBWithVectors>(cudasw4::loadDBWithVectors(options.dbPrefix));
        if(options.verbose){
            timer_read_db.print();
        }

        cudaSW4.setDatabase(fullDB_tmp);
    }

    const char* letters = "ARNDCQEGHILKMFPSTWYV";

    std::mt19937 gen(42);
    std::uniform_int_distribution<> dist(0,19);

    for(int l = firstLength; l <= lastLength; l += stepLength){

        std::string sequence(l, ' ');
        for(int i = 0; i < l; i++){
            sequence[i] = letters[dist(gen)];
        }

        if(options.verbose){
            cudaSW4.printDBInfo();
            if(options.printLengthPartitions){
                cudaSW4.printDBLengthPartitions();
            }
        }

        if(options.loadFullDBToGpu){
            cudaSW4.prefetchDBToGpus();
        }

        std::cout << "Processing query length " << sequence.size();


        cudasw4::DecodedQueryView queryView(sequence.data(), sequence.size());
        ScanResult scanResult = cudaSW4.scan(queryView, std::nullopt);
        std::cout << scanResult.stats.gcups << " GCUPS\n";
                    
    }

    return 0;
}






struct BenchmarkData{
    bool dpx;
    int tilesize;
    int groupsize;
    int numRegs;
    float gcups;
    KernelApproach kernelApproach;
};

void writeBenchmarkDataHeader(std::ostream& os){
    os << "tilesize groupsize numRegs dpx kernelApproach gcups" << "\n";
}

std::ostream& operator<<(std::ostream& os, const BenchmarkData& data){

    os << data.tilesize << " " << data.groupsize << " " << data.numRegs 
        << " " << data.dpx << " " << int(data.kernelApproach) << " " << data.gcups;
    return os;
}


int peakbenchmarkAllSingleTileConfigs(int argc, char* argv[]){
    ProgramOptions options;
    bool parseSuccess = parseArgs(argc, argv, options);

    if(!parseSuccess || options.help){
        printHelp(argc, argv);
        return 0;
    }

    options.usePseudoDB = true;
    options.pseudoDBSameSequence = false;
    options.verbose = true;
    options.loadFullDBToGpu = true;
    options.numTopOutputs = 0;
    options.verbose = false;

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
 
    using KernelTypeConfig = cudasw4::KernelTypeConfig;
    using MemoryConfig = cudasw4::MemoryConfig;
    using ScanResult = cudasw4::ScanResult;

    KernelTypeConfig kernelTypeConfig;
    kernelTypeConfig.singlePassType = options.singlePassType;
    kernelTypeConfig.manyPassType_small = options.manyPassType_small;
    kernelTypeConfig.manyPassType_large = options.manyPassType_large;
    kernelTypeConfig.overflowType = options.overflowType;

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

    auto execute = [&](std::string /*outputfilename*/, KernelApproach kernelApproach, bool useDPX){
        kernelTypeConfig.singlePassType = useDPX ? cudasw4::KernelType::DPXs16 : cudasw4::KernelType::Half2;

        cudasw4::CudaSW4 cudaSW4(
            deviceIds, 
            options.numTopOutputs,
            options.blosumType, 
            kernelTypeConfig, 
            memoryConfig, 
            options.verbose
        );

        // std::ofstream logfile(outputfilename);

        std::vector<BenchmarkData> benchmarkDataVec;

        for(int groupsize : {4,8,16}){
        // for(int groupsize : {4}){
        // for(int groupsize : {8,16}){
            for(int numRegs : {4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64}){
            // for(int groupsize : {4,8,16}){
            //     for(int numRegs : {32}){
                const int l = groupsize * numRegs * 2;
                if(l > 2048) continue;

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

                options.pseudoDBLength = l;
                options.pseudoDBSize = 5000000;

                if(options.verbose){
                    std::cout << "Generating pseudo db\n";
                }
                helpers::CpuTimer timer_read_db("Generate DB");
                auto fullDB_tmp = std::make_shared<cudasw4::PseudoDB>(cudasw4::loadPseudoDB(
                    options.pseudoDBSize, 
                    options.pseudoDBLength,
                    options.pseudoDBSameSequence
                ));
                if(options.verbose){
                    timer_read_db.print();
                }
                
                cudaSW4.setDatabase(fullDB_tmp);
                cudaSW4.setGaplessHalf2KernelApproach(kernelApproach);
                cudaSW4.setGaplessDPXKernelApproach(kernelApproach);

                if(options.verbose){
                    cudaSW4.printDBInfo();
                    if(options.printLengthPartitions){
                        cudaSW4.printDBLengthPartitions();
                    }
                }

                if(options.loadFullDBToGpu){
                    cudaSW4.prefetchDBToGpus();
                }

                // std::cout << "Processing query length " << sequence.size() << "\n";
                // std::cout << "Processing query length " << sequence.size();

                // logfile << sequence.size() << "\n";

                cudaSW4.setGroupConfig(groupsize, numRegs);
                cudasw4::DecodedQueryView queryView(sequence.data(), sequence.size());
                ScanResult scanResult = cudaSW4.scan(queryView, std::nullopt);

                benchmarkData.gcups = scanResult.stats.gcups;

                benchmarkDataVec.push_back(benchmarkData);

                std::cout << benchmarkData << "\n";

                // std::cout << l << " " << groupsize << " " << numRegs << " " << scanResult.stats.gcups << " GCUPS\n"; 


            }
        }

        return benchmarkDataVec;
    };

    for(auto kernelApproach : {KernelApproach::hardcodedzero, KernelApproach::kernelparamzero}){
        auto resultNoDpx = execute("", kernelApproach, false);
        allBenchmarkData.insert(allBenchmarkData.end(), resultNoDpx.begin(), resultNoDpx.end());

        int ccMajor = 0;
        cudaDeviceGetAttribute(&ccMajor, cudaDevAttrComputeCapabilityMajor, 0);
        const bool supportsDPX = ccMajor >= 9;
        if(supportsDPX){
            auto resultDpx = execute("", kernelApproach, true);
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

    return 0;
}