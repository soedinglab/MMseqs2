#include "options.hpp"
#include "types.hpp"
#include "hpc_helpers/all_helpers.cuh"

#include <string>
#include <iostream>



void printOptions(const ProgramOptions& options){
    std::cout << "Selected options:\n";
    std::cout << "verbose: " << options.verbose << "\n";
    std::cout << "interactive: " << options.interactive << "\n";
    std::cout << "loadFullDBToGpu: " << options.loadFullDBToGpu << "\n";
    std::cout << "prefetchDBFile: " << options.prefetchDBFile << "\n";
    std::cout << "numTopOutputs: " << options.numTopOutputs << "\n";
    std::cout << "gop: " << options.gop << "\n";
    std::cout << "gex: " << options.gex << "\n";
    std::cout << "maxBatchBytes: " << options.maxBatchBytes << "\n";
    std::cout << "maxBatchSequences: " << options.maxBatchSequences << "\n";
    std::cout << "maxTempBytes: " << options.maxTempBytes << "\n";
    for(size_t i = 0; i < options.queryFiles.size(); i++){
        std::cout << "queryFile " << i  << " : " << options.queryFiles[i] << "\n";
    }
    #ifdef CAN_USE_FULL_BLOSUM
    std::cout << "blosum: " << to_string(options.blosumType) << "\n";
    #else
    std::cout << "blosum: " << to_string_nodim(options.blosumType) << "\n";
    #endif
    if(options.usePseudoDB){
        std::cout << "Using built-in pseudo db with " << options.pseudoDBSize 
            << " sequences of length " << options.pseudoDBLength << ". ";
        if(options.pseudoDBSameSequence){
            std::cout << "All sequences are identical\n";
        }else{
            std::cout << "All sequences are different\n";
        }
    }else{
        std::cout << "Using db file: " << options.dbPrefix << "\n";
    }
    std::cout << "memory limit per gpu: " << (options.maxGpuMem == std::numeric_limits<size_t>::max() ? 
        "unlimited" : std::to_string(options.maxGpuMem)) << "\n"; 

    std::cout << "Output mode: " << options.outputModeString() << "\n";
    std::cout << "Output file: " << options.outputfile << "\n";
    std::cout << "Scan type: " << to_string(options.scanType) << "\n";
    std::cout << "File with subject Ids: " << (options.subjectIdsFilename.has_value() ? options.subjectIdsFilename.value() : " unspecified") << "\n";
    std::cout << "kernelConfigsFile_gapless: " << options.kernelConfigsFile_gapless.value_or("unspecified") << "\n";
    std::cout << "kernelConfigsFile_sw: " << options.kernelConfigsFile_sw.value_or("unspecified") << "\n";
}

bool parseArgs(int argc, char** argv, ProgramOptions& options){

    auto parseMemoryString = [](const std::string& string){
        std::size_t result = 0;
        if(string.length() > 0){
            std::size_t factor = 1;
            bool foundSuffix = false;
            switch(string.back()){
                case 'K':{
                    factor = std::size_t(1) << 10; 
                    foundSuffix = true;
                }break;
                case 'M':{
                    factor = std::size_t(1) << 20;
                    foundSuffix = true;
                }break;
                case 'G':{
                    factor = std::size_t(1) << 30;
                    foundSuffix = true;
                }break;
            }
            if(foundSuffix){
                const auto numberString = string.substr(0, string.size()-1);
                result = factor * std::stoull(numberString);
            }else{
                result = std::stoull(string);
            }
        }else{
            result = 0;
        }
        return result;
    };

    auto stringToScanType = [&](const std::string& string){
        if(string == "Gapless") return cudasw4::ScanType::Gapless;
        if(string == "SW_Endpos") return cudasw4::ScanType::SW_Endpos;
        if(string == "Gapless+SW_Endpos") return cudasw4::ScanType::GaplessPlusSW_Endpos;
        std::cout << "Unknown scan type " << string << ". Using Gapless.\n";
        return cudasw4::ScanType::Gapless;
    };

    bool gotQuery = false;
    bool gotDB = false;
    bool gotGex = false;
    bool gotGop = false;

    options.queryFiles.clear();

    for(int i = 1; i < argc; i++){
        const std::string arg = argv[i];
        if(arg == "--help"){
            options.help = true;
        }else if(arg == "--uploadFull"){
            options.loadFullDBToGpu = true;
        }else if(arg == "--verbose"){
            options.verbose = true;            
        }else if(arg == "--interactive"){
            options.interactive = true;            
        }else if(arg == "--printLengthPartitions"){
            options.printLengthPartitions = true;
        }else if(arg == "--prefetchDBFile"){
            options.prefetchDBFile = true;
        }else if(arg == "--top"){
            options.numTopOutputs = std::atoi(argv[++i]);
        }else if(arg == "--gop"){
            options.gop = std::atoi(argv[++i]);
            gotGop = true;
        }else if(arg == "--gex"){
            options.gex = std::atoi(argv[++i]);
            gotGex = true;
        }else if(arg == "--maxBatchBytes"){
            options.maxBatchBytes = parseMemoryString(argv[++i]);
        }else if(arg == "--maxBatchSequences"){
            options.maxBatchSequences = std::atoi(argv[++i]);
        }else if(arg == "--maxTempBytes"){
            options.maxTempBytes = parseMemoryString(argv[++i]);
        }else if(arg == "--maxGpuMem"){
            options.maxGpuMem = parseMemoryString(argv[++i]);
        }else if(arg == "--query"){
            options.queryFiles.push_back(argv[++i]);
            gotQuery = true;
        }else if(arg == "--db"){
            options.dbPrefix = argv[++i];
            gotDB = true;
        }else if(arg == "--mat"){
            const std::string val = argv[++i];
            #ifdef CAN_USE_FULL_BLOSUM
            if(val == "blosum45") options.blosumType = cudasw4::BlosumType::BLOSUM45;
            if(val == "blosum50") options.blosumType = cudasw4::BlosumType::BLOSUM50;
            if(val == "blosum62") options.blosumType = cudasw4::BlosumType::BLOSUM62;
            if(val == "blosum80") options.blosumType = cudasw4::BlosumType::BLOSUM80;
            if(val == "blosum45_20") options.blosumType = cudasw4::BlosumType::BLOSUM45_20;
            if(val == "blosum50_20") options.blosumType = cudasw4::BlosumType::BLOSUM50_20;
            if(val == "blosum62_20") options.blosumType = cudasw4::BlosumType::BLOSUM62_20;
            if(val == "blosum80_20") options.blosumType = cudasw4::BlosumType::BLOSUM80_20;
            #else
            if(val == "blosum45") options.blosumType = cudasw4::BlosumType::BLOSUM45_20;
            if(val == "blosum50") options.blosumType = cudasw4::BlosumType::BLOSUM50_20;
            if(val == "blosum62") options.blosumType = cudasw4::BlosumType::BLOSUM62_20;
            if(val == "blosum80") options.blosumType = cudasw4::BlosumType::BLOSUM80_20;
            if(val == "blosum45_20") options.blosumType = cudasw4::BlosumType::BLOSUM45_20;
            if(val == "blosum50_20") options.blosumType = cudasw4::BlosumType::BLOSUM50_20;
            if(val == "blosum62_20") options.blosumType = cudasw4::BlosumType::BLOSUM62_20;
            if(val == "blosum80_20") options.blosumType = cudasw4::BlosumType::BLOSUM80_20;
            #endif
        }else if(arg == "--pseudodb"){
            options.usePseudoDB = true;
            options.pseudoDBSize = std::atoi(argv[++i]);
            options.pseudoDBLength = std::atoi(argv[++i]);
            int val = std::atoi(argv[++i]);
            options.pseudoDBSameSequence = val != 0;
            gotDB = true;
        }else if(arg == "--tsv"){
            options.outputMode = ProgramOptions::OutputMode::TSV;
        }else if(arg == "--of"){
            options.outputfile = argv[++i];
        }else if(arg == "--scanType"){
            options.scanType = stringToScanType(argv[++i]);
        }else if(arg == "--subjectIdsFile"){
            options.subjectIdsFilename = argv[++i];
        }else if(arg == "--kernelconfigsGapless"){
            options.kernelConfigsFile_gapless = argv[++i];
        }else if(arg == "--kernelconfigsSW"){
            options.kernelConfigsFile_sw = argv[++i];
        }else{
            std::cout << "Unexpected arg " << arg << "\n";
        }
    }

    //set specific gop gex for blosum if no gop gex was set
    if(options.blosumType == cudasw4::BlosumType::BLOSUM45 || options.blosumType == cudasw4::BlosumType::BLOSUM45_20){
        if(!gotGop) options.gop = -13;
        if(!gotGex) options.gex = -2;
    }
    if(options.blosumType == cudasw4::BlosumType::BLOSUM50 || options.blosumType == cudasw4::BlosumType::BLOSUM50_20){
        if(!gotGop) options.gop = -13;
        if(!gotGex) options.gex = -2;
    }
    if(options.blosumType == cudasw4::BlosumType::BLOSUM62 || options.blosumType == cudasw4::BlosumType::BLOSUM62_20){
        if(!gotGop) options.gop = -11;
        if(!gotGex) options.gex = -1;
    }
    if(options.blosumType == cudasw4::BlosumType::BLOSUM80 || options.blosumType == cudasw4::BlosumType::BLOSUM80_20){
        if(!gotGop) options.gop = -10;
        if(!gotGex) options.gex = -1;
    }

    if(!gotQuery){
        std::cout << "Query is missing\n";
        return false;
    }
    if(!gotDB){
        std::cout << "DB prefix is missing\n";
        return false;
    }

    return true;
}

void printHelp(int /*argc*/, char** argv){
    ProgramOptions defaultoptions;

    std::cout << "Usage: " << argv[0] << " [options]\n";
    std::cout << "The GPUs to use are set via CUDA_VISIBLE_DEVICES environment variable.\n";
    std::cout << "Options: \n";

    std::cout << "   Mandatory\n";
    std::cout << "      --query queryfile : Mandatory. Fasta or Fastq. Can be gzip'ed. Repeat this option for multiple query files\n";
    std::cout << "      --db dbPrefix : Mandatory. The DB to query against. The same dbPrefix as used for makedb\n";
    std::cout << "\n";

    std::cout << "   Scoring\n";
    std::cout << "      --top val : Output the val best scores. Default val = " << defaultoptions.numTopOutputs << "\n";
    std::cout << "      --gop val : Gap open score. Overwrites our blosum-dependent default score.\n";
    std::cout << "      --gex val : Gap extend score. Overwrites our blosum-dependent default score.\n";
    #ifdef CAN_USE_FULL_BLOSUM
    std::cout << "      --mat val: Set substitution matrix. Supported values: blosum45, blosum50, blosum62, blosum80, blosum45_20, blosum50_20, blosum62_20, blosum80_20. "
                        "Default: " << "blosum62_20" << "\n";
    #else 
    std::cout << "      --mat val: Set substitution matrix. Supported values: blosum45, blosum50, blosum62, blosum80. "
                        "Default: " << "blosum62" << "\n";
    #endif
    std::cout << "      --scanType val : Set scan type. Supported values = {Gapless, SW_Endpos, Gapless+SW_Endpos}.\n";
    std::cout << "            Gapless: Scan whole DB with gapless alignment. \n";
    std::cout << "            SW_Endpos: Scan whole DB with Smith Waterman Alignment, output score and end position.\n";
    std::cout << "            Gapless+SW_Endpos: Scan whole DB with gapless alignment, then re-scan top results with Smith Waterman. Default val = " << to_string(defaultoptions.scanType) << "\n";
    std::cout << "      --subjectIdsFile val : Only consider database sequences with index specified in file. Must be a text file, one index per line.\n";
    std::cout << "            Do not use together with scanType Gapless+SW_Endpos. When --subjectIdsFile is set, option --top is ignored.\n";
    std::cout << "\n";

    std::cout << "   Memory\n";
    std::cout << "      --maxGpuMem val : Try not to use more than val bytes of gpu memory per gpu. Uses all available gpu memory by default\n";
    std::cout << "      --maxTempBytes val : Size of temp storage in GPU memory. Can use suffix K,M,G. Default val = " << defaultoptions.maxTempBytes << "\n";
    std::cout << "      --maxBatchBytes val : Process DB in batches of at most val bytes. Can use suffix K,M,G. Default val = " << defaultoptions.maxBatchBytes << "\n";
    std::cout << "      --maxBatchSequences val : Process DB in batches of at most val sequences. Default val = " << defaultoptions.maxBatchSequences << "\n";
    std::cout << "\n";
    
    std::cout << "   Misc\n";
    std::cout << "      --of val: Result output file. Parent directory must exist. Default: console output (/dev/stdout)\n";
    std::cout << "      --tsv : Print results as tab-separated values instead of plain text. \n";
    std::cout << "      --verbose : More console output. Shows timings. \n";
    std::cout << "      --printLengthPartitions : Print number of sequences per length partition in db.\n";
    std::cout << "      --interactive : Loads DB, then waits for sequence input by user\n";
    std::cout << "      --help : Print this message\n";
    std::cout << "\n";

    std::cout << "   Performance and benchmarking\n";
    std::cout << "      --prefetchDBFile : Load DB into RAM immediately at program start instead of waiting for the first access.\n";
    std::cout << "      --uploadFull : If enough GPU memory is available to store full db, copy full DB to GPU before processing queries.\n";
    std::cout << "      --pseudodb num length sameSeq: Use a generated DB which contains `num` equal sequences of length `length`."
                        "sameSeq can be 0 or 1. If `sameSeq`!=0, all sequences in DB will be identical\n";
    std::cout << "      --kernelconfigsGapless filename\n";
    std::cout << "      --kernelconfigsSW filename\n";
    std::cout << "\n";

            
}