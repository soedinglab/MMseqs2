


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
#include "cudasw4.cuh"
#include "config.hpp"
#include "target_subject_ids.cuh"

// #include "benchmarking.cuh"

std::vector<std::string> split(const std::string& str, char c){
	std::vector<std::string> result;

	std::stringstream ss(str);
	std::string s;

	while (std::getline(ss, s, c)) {
		result.emplace_back(s);
	}

	return result;
}

void printScanResultPlain(std::ostream& os, const cudasw4::ScanResult& scanResult, const cudasw4::CudaSW4& cudaSW4){
    const int n = scanResult.scores.size();
    for(int i = 0; i < n; i++){
        const auto referenceId = scanResult.referenceIds[i];
        os << "Result " << i << ".";
        os << " Score: " << scanResult.scores[i] << ".";
        os << " Length: " << cudaSW4.getReferenceLength(referenceId) << ".";
        os << " Header " << cudaSW4.getReferenceHeader(referenceId) << ".";
        os << " referenceId " << referenceId << ".";
        os << " Alignment_end_query " << scanResult.endPositions[i].getQueryEndInclusive() << ".";
        os << " Alignment_end_ref " << scanResult.endPositions[i].getSubjectEndInclusive();
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
    const cudasw4::ScanResult& scanResult, 
    const cudasw4::CudaSW4& cudaSW4, 
    int64_t queryId,
    cudasw4::SequenceLengthT queryLength,
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
            << cudaSW4.getReferenceLength(referenceId) << sep
            << cudaSW4.getReferenceHeader(referenceId) << sep
            << referenceId << sep
            << scanResult.endPositions[i].getQueryEndInclusive() << sep
            << scanResult.endPositions[i].getSubjectEndInclusive()
            << "\n";

        //std::cout << " Sequence " << cudaSW4.getReferenceSequence(referenceId) << "\n";
    }
}

struct BatchOfQueries{
    std::vector<char> chars;               
    std::vector<std::size_t> offsets;  
    std::vector<cudasw4::SequenceLengthT> lengths;  
    std::vector<std::string> headers;  
};



int main(int argc, char* argv[])
{
    ProgramOptions options;
    bool parseSuccess = parseArgs(argc, argv, options);

    if(!parseSuccess || options.help){
        printHelp(argc, argv);
        return 0;
    }

    // peakbenchmarkAllSingleTileConfigs(argc, argv);

    // peakBenchmark(argc, argv);
    // gridsearchPseudo(argc, argv);
    //gridsearchReal(argc, argv);

    // lengthbenchmarkReal(argc, argv, 32, 4096, 32);
    // lengthbenchmarkReal(argc, argv, 4096+1024, 16384, 1024);
    // lengthbenchmarkReal(argc, argv, 16384+16384, 65536, 16384);
    // return 0;

    // gridsearchPseudo_SW(argc, argv);
    // return 0;
   



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
 
    using MemoryConfig = cudasw4::MemoryConfig;
    using ScanResult = cudasw4::ScanResult;
    using ScanType = cudasw4::ScanType;

    MemoryConfig memoryConfig;
    memoryConfig.maxBatchBytes = options.maxBatchBytes;
    memoryConfig.maxBatchSequences = options.maxBatchSequences;
    memoryConfig.maxTempBytes = options.maxTempBytes;
    memoryConfig.maxGpuMem = options.maxGpuMem;

    //ScanType scanType = ScanType::Gapless;
    //ScanType scanType = ScanType::SW_Endpos;

    std::shared_ptr<cudasw4::TargetSubjectIds> targetSubjectIds;
    if(options.subjectIdsFilename.has_value()){
        targetSubjectIds = std::make_shared<cudasw4::TargetSubjectIds>(options.subjectIdsFilename.value());
        // for(auto x : targetSubjectIds->subjectIds){
        //     std::cout << x << ", ";
        // }
        // std::cout << "\n";
    }


    std::ofstream outputfile(options.outputfile);
    if(!bool(outputfile)){
        throw std::runtime_error("Cannot open file " + options.outputfile);
    }
    if(options.outputMode == ProgramOptions::OutputMode::TSV){
        printTSVHeader(outputfile);
    }

    int numTopOutputs = options.numTopOutputs;
    if(targetSubjectIds){
        numTopOutputs = targetSubjectIds->subjectIds.size();
    }
    numTopOutputs = std::min(numTopOutputs, cudasw4::MaxNumberOfResults::value());
    // std::cout << "Will output up to " << numTopOutputs << " results\n";

    KernelConfigFilenames kernelConfigFilenames;
    kernelConfigFilenames.gapless = options.kernelConfigsFile_gapless;
    kernelConfigFilenames.sw = options.kernelConfigsFile_sw;

    cudasw4::CudaSW4 cudaSW4(
        deviceIds, 
        numTopOutputs,
        options.blosumType, 
        memoryConfig, 
        options.verbose,
        kernelConfigFilenames
    );

    cudaSW4.setScanType(options.scanType);
    if(targetSubjectIds){
        cudaSW4.setTargetSubjectIds(targetSubjectIds);
    }

    if(!options.usePseudoDB){
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
    }else{
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

    if(!options.interactive){

        for(const auto& queryFile : options.queryFiles){
            std::cout << "Processing query file " << queryFile << "\n";
        // 0 load all queries into memory, then process.
        // 1 load and process queries one after another
        #if 1
            kseqpp::KseqPP reader(queryFile);
            int64_t query_num = 0;

            cudaSW4.totalTimerStart();

            while(reader.next() >= 0){
                std::cout << "Processing query " << query_num << " ... ";
                std::cout.flush();
                const std::string& header = reader.getCurrentHeader();
                const std::string& sequence = reader.getCurrentSequence();

                cudasw4::DecodedQueryView queryView(sequence.data(), sequence.size());

                ScanResult scanResult = cudaSW4.scan(queryView, std::nullopt);
                if(options.verbose){
                    std::cout << "Done. Scan time: " << scanResult.stats.seconds << " s, " << scanResult.stats.gcups << " GCUPS\n";
                }else{
                    std::cout << "Done.\n";
                }

                if(numTopOutputs > 0){
                    if(options.outputMode == ProgramOptions::OutputMode::Plain){
                        outputfile << "Query " << query_num << ", header" <<  header
                            << ", length " << sequence.size()
                            << ", num overflows " << scanResult.stats.numOverflows << "\n";

                        printScanResultPlain(outputfile, scanResult, cudaSW4);
                    }else{
                        printScanResultTSV(outputfile, scanResult, cudaSW4, query_num, sequence.size(), header);
                    }
                    outputfile.flush();
                }

                query_num++;
            }

            auto totalBenchmarkStats = cudaSW4.totalTimerStop();
            if(options.verbose){
                std::cout << "Total time: " << totalBenchmarkStats.seconds << " s, " << totalBenchmarkStats.gcups << " GCUPS\n";
            }

        #else

            BatchOfQueries batchOfQueries;
            {
                
                constexpr int ALIGN = 4;
                kseqpp::KseqPP reader(queryFile);
                batchOfQueries.offsets.push_back(0);
                while(reader.next() >= 0){
                    const std::string& header = reader.getCurrentHeader();
                    const std::string& sequence = reader.getCurrentSequence();
                    //we ignore quality
                    //const std::string& quality = reader.getCurrentQuality();

                    batchOfQueries.chars.insert(batchOfQueries.chars.end(), sequence.begin(), sequence.end());
                    //padding
                    if(batchOfQueries.chars.size() % ALIGN != 0){
                        batchOfQueries.chars.insert(batchOfQueries.chars.end(), ALIGN - batchOfQueries.chars.size() % ALIGN, ' ');
                    }

                    batchOfQueries.offsets.push_back(batchOfQueries.chars.size());
                    batchOfQueries.lengths.push_back(sequence.size());
                    batchOfQueries.headers.push_back(header);
                }
            }

            int64_t numQueries = batchOfQueries.lengths.size();
            const char* maxNumQueriesString = std::getenv("ALIGNER_MAX_NUM_QUERIES");
            if(maxNumQueriesString != nullptr){
                int64_t maxNumQueries = std::atoi(maxNumQueriesString);
                numQueries = std::min(numQueries, maxNumQueries);
            }
        
            std::vector<ScanResult> scanResults(numQueries);

            cudaSW4.totalTimerStart();

            for(int64_t query_num = 0; query_num < numQueries; ++query_num) {
                std::cout << "Processing query " << query_num << " ... ";
                std::cout.flush();
                const size_t offset = batchOfQueries.offsets[query_num];
                const cudasw4::SequenceLengthT length = batchOfQueries.lengths[query_num];
                const char* sequence = batchOfQueries.chars.data() + offset;
                cudasw4::DecodedQueryView queryView(sequence, length);
                ScanResult scanResult = cudaSW4.scan(queryView, std::nullopt);
                scanResults[query_num] = scanResult;
                if(options.verbose){
                    std::cout << "Done. Scan time: " << scanResult.stats.seconds << " s, " << scanResult.stats.gcups << " GCUPS\n";
                }else{
                    std::cout << "Done.\n";
                }
            }

            auto totalBenchmarkStats = cudaSW4.totalTimerStop();

            if(options.verbose){
                std::cout << "Total time: " << totalBenchmarkStats.seconds << " s, " << totalBenchmarkStats.gcups << " GCUPS\n";
            }
            if(numTopOutputs > 0){
                for(int64_t query_num = 0; query_num < numQueries; ++query_num) {
                    const ScanResult& scanResult = scanResults[query_num];

                    if(options.outputMode == ProgramOptions::OutputMode::Plain){
                        outputfile << "Query " << query_num << ", header" <<  batchOfQueries.headers[query_num] 
                            << ", length " << batchOfQueries.lengths[query_num]
                            << ", num overflows " << scanResult.stats.numOverflows << "\n";
                        printScanResultPlain(outputfile, scanResult, cudaSW4);
                    }else{
                        printScanResultTSV(outputfile, scanResult, cudaSW4, query_num, batchOfQueries.lengths[query_num], batchOfQueries.headers[query_num]);
                    }
                }
            }
        #endif

        }
    }else{
        std::cout << "Interactive mode ready\n";
        std::cout << "Use 's inputsequence' to query inputsequence against the database. Press ENTER twice to begin.\n";
        std::cout << "Use 'f inputfile' to query all sequences in inputfile\n";
        std::cout << "Use 'exit' to terminate\n";
        std::cout << "Waiting for command...\n";

        std::string line;
        while(std::getline(std::cin, line)){
            auto tokens = split(line, ' ');
            if(tokens.size() == 0) continue;

            const auto& command = tokens[0];
            if(command == "exit"){
                break;
            }else if(command == "s"){
                if(tokens.size() > 1){
                    auto& sequence = tokens[1];

                    //read the remaining lines to catch multi-line sequence input (for example copy&paste fasta sequence)
                    while(std::getline(std::cin, line)){
                        if(line.empty()) break;
                        sequence += line;
                    }

                    std::cout << "sequence: " << sequence << "\n";
                    std::cout << "Processing query " << 0 << " ... ";
                    std::cout.flush();
                    cudasw4::DecodedQueryView queryView(sequence.data(), sequence.size());
                    ScanResult scanResult = cudaSW4.scan(queryView, std::nullopt);
                    if(options.verbose){
                        std::cout << "Done. Scan time: " << scanResult.stats.seconds << " s, " << scanResult.stats.gcups << " GCUPS\n";
                    }else{
                        std::cout << "Done.\n";
                    }

                    if(options.outputMode == ProgramOptions::OutputMode::Plain){
                        printScanResultPlain(outputfile, scanResult, cudaSW4);
                    }else{
                        printScanResultTSV(outputfile, scanResult, cudaSW4, -1, sequence.size(), "-");
                    }
                }else{
                    std::cout << "Missing argument for command 's'\n";
                }
            }else if(command == "f"){
                if(tokens.size() > 1){
                    const auto& filename = tokens[1];
                    try{
                        kseqpp::KseqPP reader(filename);
                        int64_t query_num = 0;

                        while(reader.next() >= 0){
                            std::cout << "Processing query " << query_num << " ... ";
                            std::cout.flush();
                            const std::string& header = reader.getCurrentHeader();
                            const std::string& sequence = reader.getCurrentSequence();

                            cudasw4::DecodedQueryView queryView(sequence.data(), sequence.size());
                            ScanResult scanResult = cudaSW4.scan(queryView, std::nullopt);
                            if(options.verbose){
                                std::cout << "Done. Scan time: " << scanResult.stats.seconds << " s, " << scanResult.stats.gcups << " GCUPS\n";
                            }else{
                                std::cout << "Done.\n";
                            }

                            if(options.outputMode == ProgramOptions::OutputMode::Plain){
                                std::cout << "Query " << query_num << ", header" <<  header
                                << ", length " << sequence.size()
                                << ", num overflows " << scanResult.stats.numOverflows << "\n";

                                printScanResultPlain(outputfile, scanResult, cudaSW4);
                            }else{
                                printScanResultTSV(outputfile, scanResult, cudaSW4, -1, sequence.size(), "-");
                            }

                            query_num++;
                        }
                    }catch(...){
                        std::cout << "Error\n";
                    }
                }else{
                    std::cout << "Missing argument for command 'f' \n";
                }
            }else{
                std::cout << "Unrecognized command: " << command << "\n";
            }

            std::cout << "Waiting for command...\n";
        }

    }

}
