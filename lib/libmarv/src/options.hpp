#ifndef OPTIONS_HPP
#define OPTIONS_HPP

#include "types.hpp"
#include <string>
#include <iostream>
#include <optional>

struct ProgramOptions{
    enum class OutputMode{
        Plain,
        TSV
    };

    bool help = false;
    bool loadFullDBToGpu = false;
    bool usePseudoDB = false;
    bool printLengthPartitions = false;
    bool interactive = false;
    bool verbose = false;
    bool prefetchDBFile = false;
    bool pseudoDBSameSequence = true;
    int numTopOutputs = 10;
    int gop = -11;
    int gex = -1;
    int pseudoDBLength = 0;
    int pseudoDBSize = 0;
    cudasw4::BlosumType blosumType = cudasw4::BlosumType::BLOSUM62_20;
    OutputMode outputMode = OutputMode::Plain;

    cudasw4::ScanType scanType = cudasw4::ScanType::Gapless;

    size_t maxBatchBytes = 128ull * 1024ull * 1024ull;
    size_t maxBatchSequences = 10'000'000;
    size_t maxTempBytes = 4ull * 1024ull * 1024ull * 1024ull;

    size_t maxGpuMem = std::numeric_limits<size_t>::max();

    std::optional<std::string> subjectIdsFilename;
    std::string outputfile = "/dev/stdout";
    std::string dbPrefix;
    std::vector<std::string> queryFiles;

    std::optional<std::string> kernelConfigsFile_gapless;
    std::optional<std::string> kernelConfigsFile_sw;

    std::string outputModeString() const{
        switch(outputMode){
            case OutputMode::Plain: return "Plain";
            case OutputMode::TSV: return "TSV";
            default: return "Unnamed output mode";
        }
    }
};

void printOptions(const ProgramOptions& options);

bool parseArgs(int argc, char** argv, ProgramOptions& options);

void printHelp(int argc, char** argv);

#endif