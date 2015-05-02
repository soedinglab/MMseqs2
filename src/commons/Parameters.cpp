
#include "Parameters.h"
#include "Sequence.h"
#include "Debug.h"
#include "getoptpp/getopt_pp_standalone.h" // external lib for parsing
#include "Util.h"
#include <iomanip>



const MMseqsParameter Parameters::PARAM_S=MMseqsParameter(PARAM_S_ID,"-s",                    "[float]\tSensitivity in the range [1:9]");
const MMseqsParameter Parameters::PARAM_K=MMseqsParameter(PARAM_K_ID,"-k",                    "[int]\tk-mer size in the range [4:7]");
const MMseqsParameter Parameters::PARAM_THREADS=MMseqsParameter(PARAM_THREADS_ID,"--threads",        "[int]\tNumber of cores used for the computation");
const MMseqsParameter Parameters::PARAM_ALPH_SIZE=MMseqsParameter(PARAM_ALPH_SIZE_ID,"--alph-size",    "[int]\tAmino acid alphabet size");
const MMseqsParameter Parameters::PARAM_MAX_SEQ_LEN=MMseqsParameter(PARAM_MAX_SEQ_LEN_ID,"--max-seq-len","[int]\tMaximum sequence length");
const MMseqsParameter Parameters::PARAM_PROFILE=MMseqsParameter(PARAM_PROFILE_ID,"--profile",        "\tHMM Profile input");
const MMseqsParameter Parameters::PARAM_NUCL=MMseqsParameter(PARAM_NUCL_ID,"--nucl",              "\tNucleotide sequences input");
const MMseqsParameter Parameters::PARAM_Z_SCORE=MMseqsParameter(PARAM_Z_SCORE_ID,"--z-score",        "[float]\tZ-score threshold ");
const MMseqsParameter Parameters::PARAM_SKIP=MMseqsParameter(PARAM_SKIP_ID,"--skip",              "[int]\tNumber of skipped k-mers during the index table generation");
const MMseqsParameter Parameters::PARAM_MAX_SEQS=MMseqsParameter(PARAM_MAX_SEQS_ID,"--max-seqs",      "[int]\tMaximum result sequences per query");
const MMseqsParameter Parameters::PARAM_SPLIT=MMseqsParameter(PARAM_SPLIT_ID,"--split",            "[int]\tSplits target databases in n equal distrbuted junks");
const MMseqsParameter Parameters::PARAM_SUB_MAT=MMseqsParameter(PARAM_SUB_MAT_ID,"--sub-mat",        "[file]\tAmino acid substitution matrix file");
const MMseqsParameter Parameters::PARAM_SEARCH_MODE=MMseqsParameter(PARAM_SEARCH_MODE_ID,"--search-mode","[int]\tSearch mode. Local: 1 Global: 2");
const MMseqsParameter Parameters::PARAM_NO_COMP_BIAS_CORR=MMseqsParameter(PARAM_NO_COMP_BIAS_CORR_ID,"--no-comp-bias-corr","Switch off local amino acid composition bias correction");
const MMseqsParameter Parameters::PARAM_FAST_MODE=MMseqsParameter(PARAM_FAST_MODE_ID,"--fast-mode","Fast search is using Z-score instead of logP-Value and extracts hits with a score higher than 6");
const MMseqsParameter Parameters::PARAM_SPACED_KMER_MODE=MMseqsParameter(PARAM_SPACED_KMER_MODE_ID,"--spaced-kmer-mode","Spaced kmers mode (use consecutive pattern). Disable: 0, Enable: 1");
// alignment
const MMseqsParameter Parameters::PARAM_E=MMseqsParameter(PARAM_E_ID,"-e",                          "Maximum e-value");
const MMseqsParameter Parameters::PARAM_C=MMseqsParameter(PARAM_C_ID,"-c",                          "Minimum alignment coverage");
const MMseqsParameter Parameters::PARAM_MAX_REJECTED=MMseqsParameter(PARAM_MAX_REJECTED_ID,"--max-rejected","Maximum rejected alignments before alignment calculation for a query is aborted");
// clustering
const MMseqsParameter Parameters::PARAM_G=MMseqsParameter(PARAM_G_ID,"-g","Greedy clustering by sequence length");
const MMseqsParameter Parameters::PARAM_A=MMseqsParameter(PARAM_A_ID,"-a","Affinity clustering");
const MMseqsParameter Parameters::PARAM_MIN_SEQ_ID=MMseqsParameter(PARAM_MIN_SEQ_ID_ID,"--min-seq-id","Minimum sequence identity of sequences in a cluster");
const MMseqsParameter Parameters::PARAM_CASCADED=MMseqsParameter(PARAM_CASCADED_ID,"--cascaded", "\tStart the cascaded instead of simple clustering workflow");
//affinity clustering
const MMseqsParameter Parameters::PARAM_MAXITERATIONS=MMseqsParameter(PARAM_MAXITERATIONS_ID,"--max-iterations","[int]\t Maximum number of iterations in affinity propagation clustering");
const MMseqsParameter Parameters::PARAM_CONVERGENCEITERATIONS=MMseqsParameter(PARAM_CONVERGENCEITERATIONS_ID,"--convergence_iterations","[int]\t Number of iterations the set of representatives has to stay constant");
const MMseqsParameter Parameters::PARAM_DAMPING=MMseqsParameter(PARAM_DAMPING_ID,"--damping","Ratio of previous iteration entering values. Value between [0.5:1).");
const MMseqsParameter Parameters::PARAM_SIMILARITYSCORE=MMseqsParameter(PARAM_SIMILARITYSCORE_ID,"--similarity-type","Type of score used for clustering [1:5]. 1=alignment score. 2=coverage 3=sequence identity 4=E-value 5= Score per Column ");
const MMseqsParameter Parameters::PARAM_PREFERENCE=MMseqsParameter(PARAM_PREFERENCE_ID,"--preference","Preference value influences the number of clusters (default=0). High values lead to more clusters.");
// logging
const MMseqsParameter Parameters::PARAM_V=MMseqsParameter(PARAM_V_ID,"-v","Verbosity level: 0=NOTHING, 1=ERROR, 2=WARNING, 3=INFO");
// clustering workflow
const MMseqsParameter Parameters::PARAM_RESTART=MMseqsParameter(PARAM_RESTART_ID, "--restart","[int]\tRestart the clustering workflow starting with alignment or clustering.\n"
        "\t\tThe value is in the range [1:3]: 1: restart from prefiltering  2: from alignment; 3: from clustering");
const MMseqsParameter Parameters::PARAM_STEP=MMseqsParameter(PARAM_STEP_ID, "--step","[int]\t\tRestart the step of the cascaded clustering. For values in [1:3], the resprective step number, 4 is only the database merging");

const MMseqsParameter Parameters::PARAM_ORF_MIN_LENGTH=MMseqsParameter(PARAM_ORF_MIN_LENGTH_ID, "--min-length","[int]\t\tMinimum length of open reading frame to be extracted from fasta file");
const MMseqsParameter Parameters::PARAM_ORF_MAX_LENGTH=MMseqsParameter(PARAM_ORF_MAX_LENGTH_ID, "--max-length","[int]\t\tMaximum length of open reading frame to be extracted from fasta file.");
const MMseqsParameter Parameters::PARAM_ORF_MAX_GAP=MMseqsParameter(PARAM_ORF_MAX_GAP_ID, "--max-gaps","[int]\t\tMaximum number of gaps or unknown residues before an open reading frame is rejected");
const MMseqsParameter Parameters::PARAM_K_SCORE=MMseqsParameter(PARAM_K_SCORE_ID,"--k-score","[int]\tSet the K-mer threshold for the K-mer generation");
const MMseqsParameter Parameters::PARAM_KEEP_TEMP_FILES=MMseqsParameter(PARAM_KEEP_TEMP_FILES_ID,"--delete-tmp-files","\tDo not delete temporary files.");
const MMseqsParameter Parameters::PARAM_ORF_SKIP_INCOMPLETE=MMseqsParameter(PARAM_ORF_SKIP_INCOMPLETE_ID,"--skip-incomplete","\tSkip orfs that have only an end or only a start");

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
                                 size_t requiredParameterCount,
                                 bool printPar)
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

        if (ops >> GetOpt::OptionPresent("fast-mode")){
            fastMode =  true;
        }

        int spacedKmerMode = 0;
        if (ops >> GetOpt::OptionPresent("spaced-kmer-mode")){
            ops >> GetOpt::Option("spaced-kmer-mode", spacedKmerMode);
            spacedKmer = (spacedKmerMode == 1) ? true : false;
        }

        if (ops >> GetOpt::OptionPresent("delete-tmp-files"))
            keepTempFiles = false;

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
        ops >> GetOpt::Option("max-iterations", maxIteration);
        ops >> GetOpt::Option("convergence_iterations", convergenceIterations);
        ops >> GetOpt::Option("damping", dampingFactor);
        ops >> GetOpt::Option("similarity-type", similarityScoreType);
        ops >> GetOpt::Option("preference", preference);

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
    if(printPar == true)
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
            case PARAM_S_ID:
                Debug(Debug::WARNING) << "Sensitivity:             " << this->sensitivity << "\n";
                break;
            case PARAM_K_ID:
                Debug(Debug::WARNING) << "K-mer size:              " << this->kmerSize << "\n";
                break;
            case PARAM_THREADS_ID:
                Debug(Debug::WARNING) << "Threads:                 " << this->threads << "\n";
                break;
            case PARAM_ALPH_SIZE_ID:
                Debug(Debug::WARNING) << "Alphabet size:           " << this->alphabetSize << "\n";
                break;
            case PARAM_MAX_SEQ_LEN_ID:
                Debug(Debug::WARNING) << "Max. sequence length:    " << this->maxSeqLen  << "\n";
                break;
            case PARAM_PROFILE_ID:
                if(this->querySeqType == Sequence::HMM_PROFILE){
                    Debug(Debug::WARNING) << "Query input:              AA Profile\n";
                    Debug(Debug::WARNING) << "DB    input:              AA\n";
                }
                break;
            case PARAM_NUCL_ID:
                if(this->querySeqType == Sequence::NUCLEOTIDES){
                    Debug(Debug::WARNING) << "Query input:              Nucleotide\n";
                    Debug(Debug::WARNING) << "DB input:                 Nucleotide\n";
                }
                break;
            case PARAM_Z_SCORE_ID:
                Debug(Debug::WARNING) << "Z-Score threshold:       " << this->zscoreThr << "\n";
                break;
            case PARAM_SKIP_ID:
                Debug(Debug::WARNING) << "Skip Kmers:              " << this->skip << "\n";
                break;
            case PARAM_MAX_SEQS_ID:
                Debug(Debug::WARNING) << "Max. results per query:  " << this->maxResListLen  << "\n";
                break;
            case PARAM_SPLIT_ID:
                Debug(Debug::WARNING) << "Split db:                " << this->split << "\n";
                break;
            case PARAM_SUB_MAT_ID:
                Debug(Debug::WARNING) << "Sub Matrix:              " << this->scoringMatrixFile << "\n";
                break;
            case PARAM_SEARCH_MODE_ID:
                if (this->localSearch)
                    Debug(Debug::WARNING) << "Search mode:             local\n";
                else
                    Debug(Debug::WARNING) << "Search mode:             global\n";
                break;
            case PARAM_NO_COMP_BIAS_CORR_ID:
                if (this->compBiasCorrection)
                    Debug(Debug::WARNING) << "Compositional bias:      on\n";
                else
                    Debug(Debug::WARNING) << "Compositional bias:      off\n";
                break;
            case PARAM_FAST_MODE_ID:
                if(this->fastMode)
                    Debug(Debug::WARNING) << "Fastmode:                yes\n";
                else
                    Debug(Debug::WARNING) << "Fastmode:                off\n";
                break;
            case PARAM_SPACED_KMER_MODE_ID:
                if (this->spacedKmer)
                    Debug(Debug::WARNING) << "Spaced kmers:            on\n";
                else
                    Debug(Debug::WARNING) << "Spaced kmers:            off\n";
                break;
            case PARAM_E_ID:
                Debug(Debug::WARNING) << "Max. evalue:             " << this->evalThr << "\n";
                break;
            case PARAM_C_ID:
                Debug(Debug::WARNING) << "Min. sequence coverage:  " << this->covThr  << "\n";
                break;
            case PARAM_MAX_REJECTED_ID:
                Debug(Debug::WARNING) << "Max. rejected:           ";
                if (this->maxRejected == INT_MAX)
                    Debug(Debug::WARNING) << "off\n";
                else
                    Debug(Debug::WARNING) << this->maxRejected << "\n";
                break;
            case PARAM_G_ID:
                if(this->clusteringMode == GREEDY){
                    Debug(Debug::WARNING) << "Cluster type:             " << "greedy" << "\n";
                }else{
                    Debug(Debug::WARNING) << "Cluster type:             " << "simple" << "\n";
                }
                break;
            case PARAM_MIN_SEQ_ID_ID:
                Debug(Debug::WARNING) << "Min. sequence id:        " << this->seqIdThr  << "\n";
                break;
            case PARAM_CASCADED_ID:
                if(this->cascaded)
                    Debug(Debug::WARNING) << "Cluster mode:             " << "cascaded" << "\n";
                else
                    Debug(Debug::WARNING) << "Cluster mode:             " << "single" << "\n";
                break;
            case PARAM_ORF_MIN_LENGTH_ID:
                Debug(Debug::WARNING) << "Minimum length:      " << this->orfMinLength  << "\n";
                break;
            case PARAM_ORF_MAX_LENGTH_ID:
                Debug(Debug::WARNING) << "Maximum length:      " << this->orfMaxLength  << "\n";
                break;
            case PARAM_ORF_MAX_GAP_ID:
                Debug(Debug::WARNING) << "Maximum gaps in ORF:     " << this->orfMaxGaps  << "\n";
                break;
            case PARAM_K_SCORE_ID:
                if(this->kmerScore != INT_MAX)
                    Debug(Debug::WARNING) << "K-score:                 " << this->kmerScore << "\n";
                else
                    Debug(Debug::WARNING) << "K-score:                 " << "auto" << "\n";
                break;
            case PARAM_KEEP_TEMP_FILES_ID:
                if(this->keepTempFiles)
                    Debug(Debug::WARNING) << "Delete tmp files:          yes\n";
                else
                    Debug(Debug::WARNING) << "Delete tmp files:          no\n";
                break;
            case PARAM_ORF_SKIP_INCOMPLETE_ID:
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
    fastMode = false;
    spacedKmer = true;
    localSearch = true;
    zscoreThr = 50.0f;
    
    evalThr = 0.001;
    covThr = 0.0;
    maxRejected = INT_MAX;
    seqIdThr = 0.0;
    
    clusteringMode = Parameters::SET_COVER;
    validateClustering = 0;
    cascaded = false;

    maxIteration=1000;
    convergenceIterations=100;
    dampingFactor=0.6;
    similarityScoreType=APC_BITSCORE;
    preference=0;


    restart = 0;
    step = 1;
    keepTempFiles = true;

    verbosity = Debug::INFO;

    //extractorfs
    orfMinLength = 1;
    orfMaxLength = SIZE_MAX;
    orfMaxGaps = SIZE_MAX;
    orfSkipIncomplete = true;
}



