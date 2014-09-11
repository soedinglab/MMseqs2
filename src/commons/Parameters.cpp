
#include "Parameters.h"
#include "Sequence.h"
#include "Debug.h"
#include "getoptpp/getopt_pp_standalone.h" // external lib for parsing
#include <iomanip>


void Parameters::printUsageMessage(std::string programUsageHeader,
                                   std::vector<MMseqsParameter> parameters){
    std::stringstream ss;
    ss << programUsageHeader << std::endl;
    for(std::size_t i = 0; i < parameters.size(); i++) {
        ss << std::setw(25) << std::left << parameters[i].name << parameters[i].description << std::endl;
    }
    Debug(Debug::INFO) << ss.str();
}

void Parameters::parseParameters(int argc, char* pargv[],
                                 std::string programUsageHeader,
                                 std::vector<MMseqsParameter> parameters,
                                 size_t requiredParameterCount)
{
    GetOpt::GetOpt_pp ops(argc, pargv);
    ops.exceptions(std::ios::failbit); // throw exceptoin when parsing error
    //ops.exceptions(std::ios::eofbit); // throw exceptoin when parsing error

    try
    {
        ops >> GetOpt::Option('s', sensitivity);
        ops >> GetOpt::Option('k', kmerSize);
        ops >> GetOpt::Option("threads",     threads);
        ops >> GetOpt::Option("max-seq-len", maxSeqLen);
        ops >> GetOpt::Option("alph-size",   alphabetSize);

        if (ops >> GetOpt::OptionPresent("profile")){
            querySeqType  = Sequence::HMM_PROFILE;
            targetSeqType = Sequence::AMINO_ACIDS;
        }
        if (ops >> GetOpt::OptionPresent("nucl")){
            querySeqType  = Sequence::NUCLEOTIDES;
            targetSeqType = Sequence::NUCLEOTIDES;
        }
        
        if(ops >> GetOpt::OptionPresent("z-score") == false){
            // adapt z-score threshold to the sensitivity setting TODO
            // user defined threshold overwrites the automatic setting
            if (1.0 <= sensitivity && sensitivity < 2.0)
                zscoreThr = 500.0;
            else if (2.0 <= sensitivity && sensitivity < 3.0)
                zscoreThr = 300.0;
            else if (3.0 <= sensitivity && sensitivity < 4.0)
                zscoreThr = 100.0;
            else if (4.0 <= sensitivity && sensitivity < 5.0)
                zscoreThr = 50.0;
            else if (5.0 <= sensitivity && sensitivity < 6.0)
                zscoreThr = 40.0;
            else if (6.0 <= sensitivity && sensitivity < 7.0)
                zscoreThr = 30.0;
            else if (7.0 <= sensitivity && sensitivity < 8.0)
                zscoreThr = 20.0;
            else if (8.0 <= sensitivity && sensitivity <= 9.0)
                zscoreThr = 10.0;
        }
        ops >> GetOpt::Option("z-score",  zscoreThr);
        ops >> GetOpt::Option("skip",     skip);
        ops >> GetOpt::Option("max-seqs", maxResListLen);
        ops >> GetOpt::Option("split",    split);
        ops >> GetOpt::Option('m',"sub-mat",  scoringMatrixFile);
        int searchMode = 0;
        ops >> GetOpt::Option("search-mode", searchMode);
        localSearch = (searchMode == 1) ? true : false;

        if (ops >> GetOpt::OptionPresent("no-comp-bias-corr")){
            compBiasCorrection = false;
        }
        if (ops >> GetOpt::OptionPresent("no-spaced-kmer")){
            spacedKmer = false;
        }
    // alignment
        ops >> GetOpt::Option('e', evalThr);
        ops >> GetOpt::Option('c', covThr);
        ops >> GetOpt::Option("max-rejected", maxRejected);
    // clustering
        ops >> GetOpt::Option('g', clusteringMode);
        ops >> GetOpt::Option("min-seq-id", seqIdThr);
        if (ops >> GetOpt::OptionPresent("cascaded")){
            cascaded = true;
        }

    // logging
        ops >> GetOpt::Option('v', verbosity);
        ops.end_of_options();            // I'm done!

        ops.options_remain();
    }
    catch (GetOpt::GetOptEx ex)
    {
        printUsageMessage(programUsageHeader, parameters);
        Debug(Debug::INFO) << "Error in arguments" << "\n";
        EXIT(EXIT_FAILURE);
    }
    
    if (querySeqType == Sequence::NUCLEOTIDES)
        alphabetSize = 5;
    
    // read global parameters
    std::vector<std::string> getFilename;
    ops >> GetOpt::GlobalOption(getFilename);
    std::cout << getFilename.size() << std::endl;
    if(getFilename.size() < requiredParameterCount){
        printUsageMessage(programUsageHeader, parameters);
        Debug(Debug::INFO) << requiredParameterCount << " Database paths are required" << "\n";
        EXIT(EXIT_FAILURE);
    }
    
    switch (getFilename.size()) {
        case 5:
            db5 = getFilename[4];
            db5Index = db5 + ".index";
        case 4:
            db4 = getFilename[3];
            db4Index = db4 + ".index";
        case 3:
            db3 = getFilename[2];
            db3Index = db3 + ".index";
        case 2:
            db2 = getFilename[1];
            db2Index = db2 + ".index";
            db1 = getFilename[0];
            db1Index = db1 + ".index";
            break;
        default:
            printUsageMessage(programUsageHeader, parameters);
            Debug(Debug::INFO) << "Unrecognized parameters!" << "\n";
            EXIT(EXIT_FAILURE);
            break;
    }
}

void Parameters::serialize( std::ostream &stream )  {
    //todo
}

void Parameters::deserialize( std::istream &stream ) {

}

Parameters::Parameters(){
    setDefaults();
}


void Parameters::setDefaults() {
    
    // get environment
    char* mmdir = getenv ("MMDIR");
    if (mmdir == 0){
        std::cerr << "Please set the environment variable $MMDIR to your MMSEQS installation directory.\n";
        exit(1);
    }
    scoringMatrixFile = std::string(mmdir);
    scoringMatrixFile = scoringMatrixFile + "/data/blosum62.out";
    
    kmerSize =  6;
    alphabetSize = 21;
    maxSeqLen = 50000;
    maxResListLen = 300;
    sensitivity = 4.0f;
    split = 1;
    skip = 0;
    querySeqType  = Sequence::AMINO_ACIDS;
    targetSeqType = Sequence::AMINO_ACIDS;
    
    threads = 1;
#ifdef OPENMP
    threads = Util::omp_thread_count();
#endif
    compBiasCorrection = true;
    spacedKmer = true;
    localSearch = true;
    zscoreThr = 50.0f;
    
    evalThr = 0.001;
    covThr = 0.8;
    maxRejected = INT_MAX;
    
    clusteringMode = Parameters::SET_COVER;
    seqIdThr = 0.0;
    validateClustering = 0;
    cascaded = false;

    verbosity = Debug::INFO;

}



