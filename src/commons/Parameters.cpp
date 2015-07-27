#include "Parameters.h"
#include "Sequence.h"
#include "Debug.h"
#include "getoptpp/getopt_pp_standalone.h" // external lib for parsing
#include "Util.h"
#include <iomanip>


Parameters::Parameters():
        PARAM_S(PARAM_S_ID,"-s", "Sensitivity", "[float]\tSensitivity in the range [1:9]", typeid(float), (void *) &sensitivity),
        PARAM_K(PARAM_K_ID,"-k", "K-mer size", "[int]\tk-mer size in the range [4:7]",typeid(int),  (void *) &kmerSize),
        PARAM_THREADS(PARAM_THREADS_ID,"--threads", "Threads", "[int]\tNumber of cores used for the computation",typeid(int), (void *) &threads),
        PARAM_ALPH_SIZE(PARAM_ALPH_SIZE_ID,"--alph-size", "Alphabet size", "[int]\tAmino acid alphabet size",typeid(int),(void *) &alphabetSize),
        PARAM_MAX_SEQ_LEN(PARAM_MAX_SEQ_LEN_ID,"--max-seq-len","Max. sequence length", "[int]\tMaximum sequence length",typeid(int), (void *) &maxSeqLen),
        PARAM_PROFILE(PARAM_PROFILE_ID,"--profile", "Profile", "\tHMM Profile input",typeid(bool),(void *) &profile),
        PARAM_NUCL(PARAM_NUCL_ID,"--nucl", "Nucleotid","\tNucleotide sequences input",typeid(bool),(void *) &nucl),
        PARAM_Z_SCORE(PARAM_Z_SCORE_ID,"--z-score", "Z-Score threshold", "[float]\tZ-score threshold ",typeid(float),(void *) &zscoreThr ),
        PARAM_K_SCORE(PARAM_K_SCORE_ID,"--k-score", "K-score", "[int]\tSet the K-mer threshold for the K-mer generation",typeid(int),(void *) &kmerScore),
        PARAM_SKIP(PARAM_SKIP_ID,"--skip", "Skip", "[int]\tNumber of skipped k-mers during the index table generation",typeid(int),(void *) &skip),
        PARAM_MAX_SEQS(PARAM_MAX_SEQS_ID,"--max-seqs", "Max. results per query", "[int]\tMaximum result sequences per query",typeid(int),(void *) &maxResListLen),
        PARAM_SPLIT(PARAM_SPLIT_ID,"--split", "Split DB", "[int]\tSplits target databases in n equally distributed chunks",typeid(int),(void *) &split),
        PARAM_SPLIT_AMINOACID(PARAM_SPLIT_AMINOACID_ID,"--split-aa", "Split by amino acid","\tTry to find the best split for the target database by amino acid count instead.",typeid(int), (void *) &splitAA),
        PARAM_SUB_MAT(PARAM_SUB_MAT_ID,"--sub-mat", "Sub Matrix", "[file]\tAmino acid substitution matrix file",typeid(std::string),(void *) &scoringMatrixFile),
        PARAM_SEARCH_MODE(PARAM_SEARCH_MODE_ID,"--search-mode", "Search mode", "[int]\tSearch mode. Global: 0 Local: 1 Local fast: 2",typeid(int), (void *) &searchMode),
        PARAM_NO_COMP_BIAS_CORR(PARAM_NO_COMP_BIAS_CORR_ID,"--no-comp-bias-corr", "Compositional bias","Switch off local amino acid composition bias correction",typeid(bool), (void *) &compBiasCorrection),
        PARAM_FAST_MODE(PARAM_FAST_MODE_ID,"--fast-mode", "Fast Mode", "Fast search is using Z-score instead of logP-Value and extracts hits with a score higher than 6",typeid(bool), (void *) &fastMode),
        PARAM_SPACED_KMER_MODE(PARAM_SPACED_KMER_MODE_ID,"--spaced-kmer-mode", "Spaced Kmer", "Spaced kmers mode (use consecutive pattern). Disable: 0, Enable: 1",typeid(int), (void *) &spacedKmer ),
        PARAM_KEEP_TEMP_FILES(PARAM_KEEP_TEMP_FILES_ID,"--keep-tmp-files", "Keep-tmp-files" ,"\tDo not delete temporary files.",typeid(bool),(void *) &keepTempFiles),
// alignment
        PARAM_E(PARAM_E_ID,"-e", "E-value threshold", "Maximum e-value",typeid(float), (void *) &evalThr),
        PARAM_C(PARAM_C_ID,"-c", "Coverage threshold", "Minimum alignment coverage",typeid(float), (void *) &covThr),
        PARAM_MAX_REJECTED(PARAM_MAX_REJECTED_ID,"--max-rejected", "Max Reject", "Maximum rejected alignments before alignment calculation for a query is aborted",typeid(int),(void *) &maxRejected),
// clustering
        PARAM_MIN_SEQ_ID(PARAM_MIN_SEQ_ID_ID,"--min-seq-id", "Seq. Id Threshold","Minimum sequence identity of sequences in a cluster",typeid(float), (void *) &seqIdThr),
        PARAM_CLUSTER_MODE(PARAM_CLUSTER_MODE_ID,"--cluster-mode", "Cluster mode", "0 Setcover, 1 affinity clustering, 2 Greedy clustering by sequence length",typeid(int), (void *) &clusteringMode),
        PARAM_CASCADED(PARAM_CASCADED_ID,"--cascaded", "Cascaded clustering", "\tStart the cascaded instead of simple clustering workflow",typeid(bool), (void *) &cascaded),
//affinity clustering
        PARAM_MAXITERATIONS(PARAM_MAXITERATIONS_ID,"--max-iterations", "Max iterations affinity", "[int]\tMaximum number of iterations in affinity propagation clustering",typeid(int), (void *) &maxIteration),
        PARAM_CONVERGENCEITERATIONS(PARAM_CONVERGENCEITERATIONS_ID,"--convergence_iterations", "Convergence iterations", "[int]\t Number of iterations the set of representatives has to stay constant",typeid(int), (void *) &convergenceIterations),
        PARAM_DAMPING(PARAM_DAMPING_ID,"--damping", "Damping", "Ratio of previous iteration entering values. Value between [0.5:1).",typeid(float), (void *) &dampingFactor),
        PARAM_SIMILARITYSCORE(PARAM_SIMILARITYSCORE_ID,"--similarity-type", "Similarity type", "Type of score used for clustering [1:5]. 1=alignment score. 2=coverage 3=sequence identity 4=E-value 5= Score per Column ",typeid(int),(void *) &similarityScoreType),
        PARAM_PREFERENCE(PARAM_PREFERENCE_ID,"--preference", "Preference", "Preference value influences the number of clusters (default=0). High values lead to more clusters.",typeid(float), (void *) &preference),
// logging
        PARAM_V(PARAM_V_ID,"-v", "Verbosity","Verbosity level: 0=NOTHING, 1=ERROR, 2=WARNING, 3=INFO",typeid(int), (void *) &verbosity),
// clustering workflow
        PARAM_RESTART(PARAM_RESTART_ID, "--restart", "Restart","[int]\tRestart the clustering workflow starting with alignment or clustering.\n"
                "\t\tThe value is in the range [1:3]: 1: restart from prefiltering  2: from alignment; 3: from clustering",typeid(int),(void *) &restart),
        PARAM_STEP(PARAM_STEP_ID, "--step","Step","[int]\t\tRestart the step of the cascaded clustering. For values in [1:3], the resprective step number, 4 is only the database merging",typeid(int),(void *) &step),
// search workflow
        PARAM_NUM_ITERATIONS(PARAM_NUM_ITERATIONS_ID, "--num-iterations", "Number search iterations","[int]\tSearch iterations",typeid(int),(void *) &numIterations),

        PARAM_ORF_MIN_LENGTH(PARAM_ORF_MIN_LENGTH_ID, "--min-length", "Min orf length", "[int]\t\tMinimum length of open reading frame to be extracted from fasta file",typeid(int),(void *) &orfMinLength),
        PARAM_ORF_MAX_LENGTH(PARAM_ORF_MAX_LENGTH_ID, "--max-length", "Max orf length", "[int]\t\tMaximum length of open reading frame to be extracted from fasta file.",typeid(int),(void *) &orfMaxLength),
        PARAM_ORF_MAX_GAP(PARAM_ORF_MAX_GAP_ID, "--max-gaps", "Max orf gaps", "[int]\t\tMaximum number of gaps or unknown residues before an open reading frame is rejected",typeid(int),(void *) &orfMaxGaps),
        PARAM_ORF_SKIP_INCOMPLETE(PARAM_ORF_SKIP_INCOMPLETE_ID,"--skip-incomplete", "Skip incomplete orfs", "\tSkip orfs that have only an end or only a start",typeid(bool),(void *) &orfSkipIncomplete)
{
    // alignment
    alignment.push_back(PARAM_E);
    alignment.push_back(PARAM_C);
    alignment.push_back(PARAM_MIN_SEQ_ID);
    alignment.push_back(PARAM_MAX_SEQ_LEN);
    alignment.push_back(PARAM_MAX_SEQS);
    alignment.push_back(PARAM_MAX_REJECTED);
    alignment.push_back(PARAM_NUCL);
    alignment.push_back(PARAM_PROFILE);
    alignment.push_back(PARAM_SUB_MAT);
    alignment.push_back(PARAM_THREADS);
    alignment.push_back(PARAM_V);

    // prefilter
    prefilter.push_back(PARAM_S);
    prefilter.push_back(PARAM_K);
    prefilter.push_back(PARAM_K_SCORE);
    prefilter.push_back(PARAM_ALPH_SIZE);
    prefilter.push_back(PARAM_MAX_SEQ_LEN);
    prefilter.push_back(PARAM_PROFILE);
    prefilter.push_back(PARAM_NUCL);
    prefilter.push_back(PARAM_Z_SCORE);
    prefilter.push_back(PARAM_SKIP);
    prefilter.push_back(PARAM_MAX_SEQS);
    prefilter.push_back(PARAM_SPLIT);
    prefilter.push_back(PARAM_SEARCH_MODE);
    prefilter.push_back(PARAM_NO_COMP_BIAS_CORR);
    prefilter.push_back(PARAM_FAST_MODE);
    prefilter.push_back(PARAM_SPACED_KMER_MODE);
    prefilter.push_back(PARAM_SUB_MAT);
    prefilter.push_back(PARAM_THREADS);
    prefilter.push_back(PARAM_V);

    // clustering
    clustering.push_back(PARAM_CLUSTER_MODE);
    clustering.push_back(PARAM_MAX_SEQS);
    clustering.push_back(PARAM_V);
    clustering.push_back(PARAM_MAXITERATIONS);
    clustering.push_back(PARAM_CONVERGENCEITERATIONS);
    clustering.push_back(PARAM_DAMPING);
    clustering.push_back(PARAM_SIMILARITYSCORE);
    clustering.push_back(PARAM_PREFERENCE);
    clustering.push_back(PARAM_MIN_SEQ_ID);

    // find orf
    onlyverbosity.push_back(PARAM_V);

    // create profile db
    createprofiledb.push_back(PARAM_SUB_MAT);
    createprofiledb.push_back(PARAM_V);

    // extract orf
    extractorf.push_back(PARAM_ORF_MIN_LENGTH);
    extractorf.push_back(PARAM_ORF_MAX_LENGTH);
    extractorf.push_back(PARAM_ORF_MAX_GAP);
    extractorf.push_back(PARAM_ORF_SKIP_INCOMPLETE);

    // splitffindex
    splitffindex.push_back(PARAM_SPLIT);
    splitffindex.push_back(PARAM_SPLIT_AMINOACID);

    // create index
    createindex.push_back(PARAM_K);
    createindex.push_back(PARAM_ALPH_SIZE);
    createindex.push_back(PARAM_MAX_SEQ_LEN);
    createindex.push_back(PARAM_SPLIT);
    createindex.push_back(PARAM_SUB_MAT);
    createindex.push_back(PARAM_SEARCH_MODE);
    createindex.push_back(PARAM_SKIP);
    createindex.push_back(PARAM_SPACED_KMER_MODE);
    createindex.push_back(PARAM_V);

    setDefaults();
}


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
                                 std::vector<MMseqsParameter> par,
                                 size_t requiredParameterCount,
                                 bool printPar)
{
    GetOpt::GetOpt_pp ops(argc, pargv);
    ops.exceptions(std::ios::failbit); // throw exception when parsing error
    //ops.exceptions(std::ios::eofbit); // throw exception when parsing error
    try
    {
        for(size_t i = 0; i < par.size(); i++) {
            if (par[i].name[0] == '-' && par[i].name[1] == '-') {
                if (typeid(int) == par[i].type) {
                    ops >> GetOpt::Option(par[i].name + 2, *((int *) par[i].value));
                } else if (typeid(float) == par[i].type) {
                    ops >> GetOpt::Option(par[i].name + 2, *((float *) par[i].value));
                } else if (typeid(std::string) == par[i].type) {
                    std::string val;
                    ops >> GetOpt::Option(par[i].name + 2, val );
                    if(val.length() != 0)
                        par[i].value = (void *) &val;
                } else if (typeid(bool) == par[i].type) {
                    if (ops >> GetOpt::OptionPresent(par[i].name + 2)) {
                        bool * value = (bool *) par[i].value;
                        // toggle Value
                        *value = !*value;
                    }
                } else {
                    Debug(Debug::ERROR) << "Wrong parameter type in parseParameters. Please inform developer\n";
                    EXIT(EXIT_FAILURE);
                }
            } else if (par[i].name[0] == '-') {
                if (typeid(int) == par[i].type) {
                    ops >> GetOpt::Option(par[i].name[1], *((int *) par[i].value));
                } else if (typeid(float) == par[i].type) {
                    ops >> GetOpt::Option(par[i].name[1], *((float *) par[i].value));
                } else if (typeid(std::string) == par[i].type) {
                    std::string val;
                    ops >> GetOpt::Option(par[i].name[1], val);
                    if(val.length() != 0){
                        std::string currVal = *((std::string *)par[i].value);
                        currVal.assign( val );
                    }
                } else if (typeid(bool) == par[i].type) {
                    if (ops >> GetOpt::OptionPresent(par[i].name[1])) {
                        bool * value = (bool *) par[i].value;
                        // toggle Value
                        *value = !*value;
                    }
                } else {
                    Debug(Debug::ERROR) << "Wrong parameter type in parseParameters. Please inform developer\n";
                    EXIT(EXIT_FAILURE);
                }
            }
        }
        ops.end_of_options();            // I'm done!
        ops.options_remain();
    }
    catch (const GetOpt::GetOptEx &ex) {
        printUsageMessage(programUsageHeader, par);
        Debug(Debug::INFO) << "Error in arguments" << "\n";
        EXIT(EXIT_FAILURE);
    }


    if (profile){
        querySeqType  = Sequence::HMM_PROFILE;
        targetSeqType = Sequence::AMINO_ACIDS;
    }
    if (nucl){
        querySeqType  = Sequence::NUCLEOTIDES;
        targetSeqType = Sequence::NUCLEOTIDES;
    }

    if (querySeqType == Sequence::NUCLEOTIDES){
        alphabetSize = 5;
    }
    // read global parameters
    std::vector<std::string> getFilename;
    ops >> GetOpt::GlobalOption(getFilename);

    if(getFilename.size() < requiredParameterCount){
        printUsageMessage(programUsageHeader, par);
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
            printUsageMessage(programUsageHeader, par);
            Debug(Debug::INFO) << "Unrecognized parameters!" << "\n";
            EXIT(EXIT_FAILURE);
            break;
    }
    if(printPar == true)
        printParameters(argc, pargv, par);
}

void Parameters::printParameters(int argc, const char* pargv[],
                                 std::vector<MMseqsParameter> par){
    Debug(Debug::WARNING) << "Program call:\n";
    for (int i = 0; i < argc; i++)
        Debug(Debug::WARNING) << pargv[i] << " ";
    Debug(Debug::WARNING) << "\n\n";
    std::stringstream ss;
    for (size_t i = 0; i < par.size(); i++) {
        ss << std::setw(25) << std::left << par[i].display;
        if(typeid(int) == par[i].type ){
            ss << *((int *)par[i].value);
        } else if(typeid(float) == par[i].type ){
            ss << *((float *)par[i].value);
        }else if(typeid(std::string) == par[i].type ){
            ss << *((std::string *) par[i].value);
        }else if (typeid(bool) == par[i].type){
            bool value = *((bool *)par[i].value);
            std::string valueStr = (value == true) ? "on" : "off";
            ss << valueStr;
        }
        ss << "\n";
    }

    Debug(Debug::WARNING) << ss.str() << "\n";
}

void Parameters::serialize( std::ostream &stream )  {
}

void Parameters::deserialize( std::istream &stream ) {

}

void Parameters::setDefaults() {

    // get environment
    char* mmdirStr = getenv ("MMDIR");
    if (mmdirStr == 0){
        std::cerr << "Please set the environment variable $MMDIR to your MMSEQS installation directory.\n";
        EXIT(EXIT_FAILURE);
    }
    mmdir = std::string(mmdirStr);
    scoringMatrixFile += mmdir + "/data/blosum62.out";

    kmerSize =  6;
    kmerScore = INT_MAX;
    alphabetSize = 21;
    maxSeqLen = 50000;
    maxResListLen = 300;
    sensitivity = 4.0f;
    split = 1;
    splitAA = false;
    skip = 0;
    querySeqType  = Sequence::AMINO_ACIDS;
    targetSeqType = Sequence::AMINO_ACIDS;
    numIterations = 1;
    threads = 1;
#ifdef OPENMP
    threads = Util::omp_thread_count();
#endif
    compBiasCorrection = true;
    fastMode = false;
    spacedKmer = true;
    searchMode = true;
    profile = false;
    nucl = false;
    zscoreThr = 50.0f;

    evalThr = 0.001;
    covThr = 0.0;
    maxRejected = INT_MAX;
    seqIdThr = 0.0;

    clusteringMode = Parameters::SET_COVER;
    validateClustering = 0;
    cascaded = false;

    // affinity clustering
    maxIteration=1000;
    convergenceIterations=100;
    dampingFactor=0.6;
    similarityScoreType=APC_BITSCORE;
    preference=0;

    // Clustering workflow
    restart = 0;
    step = 1;
    keepTempFiles = true;

    // logging
    verbosity = Debug::INFO;

    //extractorfs
    orfMinLength = 1;
    orfMaxLength = SIZE_MAX;
    orfMaxGaps = SIZE_MAX;
    orfSkipIncomplete = true;
}

std::vector<MMseqsParameter> Parameters::combineList(std::vector<MMseqsParameter> par1,
                                                     std::vector<MMseqsParameter> par2) {
    std::vector<MMseqsParameter> retVec;
    std::vector< std::vector<MMseqsParameter>> tmp;
    tmp.push_back(par1);
    tmp.push_back(par2);
    for(size_t i = 0; i < tmp.size(); i++) {
        std::vector<MMseqsParameter> currPar = tmp[i];
        for (size_t i = 0; i < currPar.size(); i++) {
            bool addPar = true;
            for (size_t j = 0; j < retVec.size(); j++) {
                if (currPar[i].uniqid == retVec[j].uniqid)
                    addPar = false;
            }
            if (addPar == true) {
                retVec.push_back(currPar[i]);
            }
        }
    }
    return retVec;
}

std::string Parameters::createParameterString(std::vector<MMseqsParameter> par) {
    std::stringstream ss;
    for (size_t i = 0; i < par.size(); i++) {
        if(typeid(int) == par[i].type ){
            ss << par[i].name << " ";
            ss << *((int *)par[i].value) << " ";
        } else if(typeid(float) == par[i].type ){
            ss << par[i].name << " ";
            ss << *((float *)par[i].value) << " ";
        }else if(typeid(std::string) == par[i].type ){
            ss << par[i].name << " ";
            ss << *((std::string *) par[i].value) << " ";
        }else if (typeid(bool) == par[i].type){
            bool val = *((bool *)(par[i].value));
            if(val == true){
                ss << par[i].name << " ";
            }
        }else
            Debug(Debug::ERROR) << "Wrong parameter type. Please inform developer\n";
    }
    return ss.str();
}
