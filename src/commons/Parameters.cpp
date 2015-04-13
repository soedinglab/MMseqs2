
#include "Parameters.h"
#include "Sequence.h"
#include "Debug.h"
#include "getoptpp/getopt_pp_standalone.h" // external lib for parsing
#include "Util.h"
#include <iomanip>



const MMseqsParameter Parameters::PARAM_S=MMseqsParameter(0,"-s",                    "[float]\tSensitivity in the range [1:9]");
const MMseqsParameter Parameters::PARAM_K={1,"-k",                    "[int]\tk-mer size in the range [4:7]"};
const MMseqsParameter Parameters::PARAM_THREADS={2,"--threads",        "[int]\tNumber of cores used for the computation"};
const MMseqsParameter Parameters::PARAM_ALPH_SIZE={3,"--alph-size",    "[int]\tAmino acid alphabet size"};
const MMseqsParameter Parameters::PARAM_MAX_SEQ_LEN={4,"--max-seq-len","[int]\tMaximum sequence length"};
const MMseqsParameter Parameters::PARAM_PROFILE={5,"--profile",        "\tHMM Profile input"};
const MMseqsParameter Parameters::PARAM_NUCL={6,"--nucl",              "\tNucleotide sequences input"};
const MMseqsParameter Parameters::PARAM_Z_SCORE={7,"--z-score",        "[float]\tZ-score threshold "};
const MMseqsParameter Parameters::PARAM_SKIP={8,"--skip",              "[int]\tNumber of skipped k-mers during the index table generation"};
const MMseqsParameter Parameters::PARAM_MAX_SEQS={9,"--max-seqs",      "[int]\tMaximum result sequences per query"};
const MMseqsParameter Parameters::PARAM_SPLIT={10,"--split",            "[int]\tSplits target databases in n equal distrbuted junks"};
const MMseqsParameter Parameters::PARAM_SUB_MAT={11,"--sub-mat",        "[file]\tAmino acid substitution matrix file"};
const MMseqsParameter Parameters::PARAM_SEARCH_MODE={12,"--search-mode","[int]\tSearch mode. Local: 1 Global: 2"};
const MMseqsParameter Parameters::PARAM_NO_COMP_BIAS_CORR={13,"--no-comp-bias-corr","Switch off local amino acid composition bias correction"};
const MMseqsParameter Parameters::PARAM_SPACED_KMER_MODE={14,"--spaced-kmer-mode","Spaced kmers mode (use consecutive pattern). Disable: 0, Enable: 1"};
// alignment
const MMseqsParameter Parameters::PARAM_E={15,"-e",                          "Maximum e-value"};
const MMseqsParameter Parameters::PARAM_C={16,"-c",                          "Minimum alignment coverage"};
const MMseqsParameter Parameters::PARAM_MAX_REJECTED={17,"--max-rejected","Maximum rejected alignments before alignment calculation for a query is aborted"};
// clustering
const MMseqsParameter Parameters::PARAM_G={18,"-g","Greedy clustering by sequence length"};
const MMseqsParameter Parameters::PARAM_A={28,"-a","Affinity clustering"};
const MMseqsParameter Parameters::PARAM_MIN_SEQ_ID={19,"--min-seq-id","Minimum sequence identity of sequences in a cluster"};
const MMseqsParameter Parameters::PARAM_CASCADED={20,"--cascaded", "\tStart the cascaded instead of simple clustering workflow"};
// logging
const MMseqsParameter Parameters::PARAM_V={21,"-v","Verbosity level: 0=NOTHING, 1=ERROR, 2=WARNING, 3=INFO"};
// clustering workflow
const MMseqsParameter Parameters::PARAM_RESTART={22, "--restart","[int]\tRestart the clustering workflow starting with alignment or clustering.\n"
        "\t\tThe value is in the range [1:3]: 1: restart from prefiltering  2: from alignment; 3: from clustering"};
const MMseqsParameter Parameters::PARAM_STEP={23, "--step","[int]\t\tRestart the step of the cascaded clustering. For values in [1:3], the resprective step number, 4 is only the database merging"};

const MMseqsParameter Parameters::PARAM_ORF_MIN_LENGTH={24, "--min-length","[int]\t\tMinimum length of open reading frame to be extracted from fasta file"};
const MMseqsParameter Parameters::PARAM_ORF_MAX_LENGTH={25, "--max-length","[int]\t\tMaximum length of open reading frame to be extracted from fasta file."};
const MMseqsParameter Parameters::PARAM_ORF_MAX_GAP={26, "--max-gaps","[int]\t\tMaximum number of gaps or unknown residues before an open reading frame is rejected"};
const MMseqsParameter Parameters::PARAM_K_SCORE={27,"--k-score","[int]\tSet the K-mer threshold for the K-mer generation"};
const MMseqsParameter Parameters::PARAM_KEEP_TEMP_FILES={28,"--keep-temp-files","\tDo not delete temporary files."};
const MMseqsParameter Parameters::PARAM_ORF_SKIP_INCOMPLETE={29,"--skip-incomplete","\tSkip orfs that have only an end or only a start"};

void Parameters::printUsageMessage(std::string programUsageHeader,
                                   std::vector<MMseqsParameter> parameters){
    std::stringstream ss;
    ss << programUsageHeader << std::endl;
    for(std::size_t i = 0; i < parameters.size(); i++) {
        ss << std::setw(25) << std::left << parameters[i].name << parameters[i].description << std::endl;
    }
    Debug(Debug::INFO) << ss.str();
}

void Parameters::parseParameters(int argc, const char* pargv[],
                                 std::string programUsageHeader,
                                 std::vector<MMseqsParameter> parameters,
                                 size_t requiredParameterCount)
{
    GetOpt::GetOpt_pp ops(argc, pargv);
    ops.exceptions(std::ios::failbit); // throw exception when parsing error
    //ops.exceptions(std::ios::eofbit); // throw exception when parsing error

    try
    {
        ops >> GetOpt::Option('s', sensitivity);
        ops >> GetOpt::Option('k', kmerSize);
        ops >> GetOpt::Option("k-score", kmerScore);
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

        if((ops >> GetOpt::OptionPresent("z-score")) == false){
            // adapt z-score threshold to the sensitivity setting
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
        ops >> GetOpt::Option('l',"max-seqs", maxResListLen);
        ops >> GetOpt::Option("split",    split);
        ops >> GetOpt::Option('m',"sub-mat",  scoringMatrixFile);

        int searchMode = 0;
        if (ops >> GetOpt::OptionPresent("search-mode")){
            ops >> GetOpt::Option("search-mode", searchMode);
            localSearch = (searchMode == 1) ? true : false;
        }

        if (ops >> GetOpt::OptionPresent("no-comp-bias-corr")){
            compBiasCorrection = false;
        }

        int spacedKmerMode = 0;
        if (ops >> GetOpt::OptionPresent("spaced-kmer-mode")){
            ops >> GetOpt::Option("spaced-kmer-mode", spacedKmerMode);
            spacedKmer = (spacedKmerMode == 1) ? true : false;
        }

        if (ops >> GetOpt::OptionPresent("keep-temp-files"))
            keepTempFiles = true;

        // alignment
        ops >> GetOpt::Option('e', evalThr);
        ops >> GetOpt::Option('c', covThr);
        ops >> GetOpt::Option("max-rejected", maxRejected);

        // clustering
        if (ops >> GetOpt::OptionPresent('g')) {
            clusteringMode = Parameters::GREEDY;
        }
        if (ops >> GetOpt::OptionPresent('a')) {
            clusteringMode = Parameters::AFFINITY;
        }
        ops >> GetOpt::Option("min-seq-id", seqIdThr);
        if (ops >> GetOpt::OptionPresent("cascaded")){
            cascaded = true;
        }

        // logging
        ops >> GetOpt::Option('v', verbosity);

        // clustering workflow
        ops >> GetOpt::Option("step", step);
        ops >> GetOpt::Option("restart", restart);

        // extractorf
        ops >> GetOpt::Option("min-length", orfMinLength);
        ops >> GetOpt::Option("max-length", orfMaxLength);
        ops >> GetOpt::Option("max-gaps", orfMaxGaps);
        if (ops >> GetOpt::OptionPresent("skip-incomplete")){
            orfSkipIncomplete = true;
        }


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
    printParameters(argc,pargv,parameters);
}

void Parameters::printParameters(int argc, const char* pargv[],
                                 std::vector<MMseqsParameter> parameters){
    Debug(Debug::WARNING) << "Program call:\n";
    for (int i = 0; i < argc; i++)
        Debug(Debug::WARNING) << pargv[i] << " ";
    Debug(Debug::WARNING) << "\n\n";
    
    for (size_t i = 0; i < parameters.size(); i++) {
        switch (parameters[i].uniqid) {
            case 0:
                Debug(Debug::WARNING) << "Sensitivity:             " << this->sensitivity << "\n";
                break;
            case 1:
                Debug(Debug::WARNING) << "K-mer size:              " << this->kmerSize << "\n";
                break;
            case 2:
                Debug(Debug::WARNING) << "Threads:                 " << this->threads << "\n";
                break;
            case 3:
                Debug(Debug::WARNING) << "Alphabet size:           " << this->alphabetSize << "\n";
                break;
            case 4:
                Debug(Debug::WARNING) << "Max. sequence length:    " << this->maxSeqLen  << "\n";
                break;
            case 5:
                if(this->querySeqType == Sequence::HMM_PROFILE){
                    Debug(Debug::WARNING) << "Query input:              AA Profile\n";
                    Debug(Debug::WARNING) << "DB    input:              AA\n";
                }
                break;
            case 6:
                if(this->querySeqType == Sequence::NUCLEOTIDES){
                    Debug(Debug::WARNING) << "Query input:              Nucleotide\n";
                    Debug(Debug::WARNING) << "DB input:                 Nucleotide\n";
                }
                break;
            case 7:
                Debug(Debug::WARNING) << "Z-Score threshold:       " << this->zscoreThr << "\n";
                break;
            case 8:
                Debug(Debug::WARNING) << "Skip Kmers:              " << this->skip << "\n";
                break;
            case 9:
                Debug(Debug::WARNING) << "Max. results per query:  " << this->maxResListLen  << "\n";
                break;
            case 10:
                Debug(Debug::WARNING) << "Split db:                " << this->split << "\n";
                break;
            case 11:
                Debug(Debug::WARNING) << "Sub Matrix:              " << this->scoringMatrixFile << "\n";
                break;
            case 12:
                if (this->localSearch)
                    Debug(Debug::WARNING) << "Search mode:             local\n";
                else
                    Debug(Debug::WARNING) << "Search mode:             global\n";
                break;
            case 13:
                if (this->compBiasCorrection)
                    Debug(Debug::WARNING) << "Compositional bias:      on\n";
                else
                    Debug(Debug::WARNING) << "Compositional bias:      off\n";
                break;
            case 14:
                if (this->spacedKmer)
                    Debug(Debug::WARNING) << "Spaced kmers:            on\n";
                else
                    Debug(Debug::WARNING) << "Spaced kmers:            off\n";
                break;
            case 15:
                Debug(Debug::WARNING) << "Max. evalue:             " << this->evalThr << "\n";
                break;
            case 16:
                Debug(Debug::WARNING) << "Min. sequence coverage:  " << this->covThr  << "\n";
                break;
            case 17:
                Debug(Debug::WARNING) << "Max. rejected:           ";
                if (this->maxRejected == INT_MAX)
                    Debug(Debug::WARNING) << "off\n";
                else
                    Debug(Debug::WARNING) << this->maxRejected << "\n";
                break;
            case 18:
                if(this->clusteringMode == GREEDY){
                    Debug(Debug::WARNING) << "Cluster type:             " << "greedy" << "\n";
                }else{
                    Debug(Debug::WARNING) << "Cluster type:             " << "simple" << "\n";
                }
                break;
            case 19:
                Debug(Debug::WARNING) << "Min. sequence id:        " << this->seqIdThr  << "\n";
                break;
            case 20:
                if(this->cascaded)
                    Debug(Debug::WARNING) << "Cluster mode:             " << "cascaded" << "\n";
                else
                    Debug(Debug::WARNING) << "Cluster mode:             " << "single" << "\n";
                break;
            case 24:
                Debug(Debug::WARNING) << "Minimum length:      " << this->orfMinLength  << "\n";
                break;
            case 25:
                Debug(Debug::WARNING) << "Maximum length:      " << this->orfMaxLength  << "\n";
                break;
            case 26:
                Debug(Debug::WARNING) << "Maximum gaps in ORF:     " << this->orfMaxGaps  << "\n";
                break;
            case 27:
                if(this->kmerScore != INT_MAX)
                    Debug(Debug::WARNING) << "K-score:             " << this->kmerScore << "\n";
                else
                    Debug(Debug::WARNING) << "K-score:             " << "auto" << "\n";
                break;
            case 28:
                if(this->keepTempFiles)
                    Debug(Debug::WARNING) << "Keep temp files:          yes\n";
                else
                    Debug(Debug::WARNING) << "Keep temp files:          no\n";
                break;
            case 29:
                if(this->orfSkipIncomplete)
                    Debug(Debug::WARNING) << "Skip incomplete ORFs:     yes\n";
                else
                    Debug(Debug::WARNING) << "Skip incomplete ORFs:     no\n";
                break;
            default:
                break;
        }
    }
    Debug(Debug::WARNING) << "\n";
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
        EXIT(1);
    }
    scoringMatrixFile = std::string(mmdir);
    scoringMatrixFile = scoringMatrixFile + "/data/blosum62.out";
    
    kmerSize =  6;
    kmerScore = INT_MAX;
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
    covThr = 0.0;
    maxRejected = INT_MAX;
    
    clusteringMode = Parameters::SET_COVER;
    seqIdThr = 0.0;
    validateClustering = 0;
    cascaded = false;
    restart = 0;
    step = 1;
    keepTempFiles = false;

    verbosity = Debug::INFO;

    //extractorfs
    orfMinLength = 1;
    orfMaxLength = SIZE_MAX;
    orfMaxGaps = SIZE_MAX;
    orfSkipIncomplete = true;
}



