#include "Parameters.h"
#include "Sequence.h"
#include "Debug.h"
#include "Util.h"
#include <iomanip>

Parameters::Parameters():
        PARAM_S(PARAM_S_ID,"-s", "Sensitivity","[int]\tSensitivity in the range [1:10]. From low (1) to high (10) sensitivity.", typeid(int), (void *) &sensitivity, "^[1-9]\{1\}$"),
        PARAM_K(PARAM_K_ID,"-k", "K-mer size", "[int]\tk-mer size in the range [4:7]",typeid(int),  (void *) &kmerSize, "^[6-7]\{1\}$"),
        PARAM_THREADS(PARAM_THREADS_ID,"--threads", "Threads", "[int]\tNumber of cores used for the computation",typeid(int), (void *) &threads, "^[1-9]\{1\}[0-9]*$"),
        PARAM_ALPH_SIZE(PARAM_ALPH_SIZE_ID,"--alph-size", "Alphabet size", "[int]\tAmino acid alphabet size",typeid(int),(void *) &alphabetSize, "^[1-9]\{1\}[0-9]\{1\}$"),
        PARAM_MAX_SEQ_LEN(PARAM_MAX_SEQ_LEN_ID,"--max-seq-len","Max. sequence length", "[int]\tMaximum sequence length",typeid(int), (void *) &maxSeqLen, "^[1-9]\{1\}[0-9]*$"),
        PARAM_PROFILE(PARAM_PROFILE_ID,"--profile", "Profile", "\tHMM Profile input",typeid(bool),(void *) &profile, ""),
        PARAM_NUCL(PARAM_NUCL_ID,"--nucl", "Nucleotid","\tNucleotide sequences input",typeid(bool),(void *) &nucl , ""),
        PARAM_Z_SCORE(PARAM_Z_SCORE_ID,"--z-score", "Z-Score threshold", "[float]\tZ-score threshold (only implemented in global search mode)",typeid(float),(void *) &zscoreThr, "^[0-9]*(\\.[0-9])?$"),
        PARAM_K_SCORE(PARAM_K_SCORE_ID,"--k-score", "K-score", "[int]\tSet the K-mer threshold for the K-mer generation",typeid(int),(void *) &kmerScore,  "^[1-9]\{1\}[0-9]*$"),
        PARAM_SKIP(PARAM_SKIP_ID,"--skip", "Skip", "[int]\tNumber of skipped k-mers during the index table generation",typeid(int),(void *) &skip,  "^[0-9]\{1\}[0-9]*$"),
        PARAM_MAX_SEQS(PARAM_MAX_SEQS_ID,"--max-seqs", "Max. results per query", "[int]\tMaximum result sequences per query",typeid(int),(void *) &maxResListLen, "^[1-9]\{1\}[0-9]*$"),
        PARAM_SPLIT(PARAM_SPLIT_ID,"--split", "Split DB", "[int]\tSplits target set in n equally distributed chunks",typeid(int),(void *) &split,  "^[1-9]\{1\}[0-9]*$"),
        PARAM_SPLIT_MODE(PARAM_SPLIT_MODE_ID,"--split-mode", "Split mode", "[int]\tMPI Option: Target set: 0 (low memory) or query set: 1 (faster but memory intensiv)",typeid(int),(void *) &splitMode,  "^[0-1]\{1\}$"),
        PARAM_SPLIT_AMINOACID(PARAM_SPLIT_AMINOACID_ID,"--split-aa", "Split by amino acid","\tTry to find the best split for the target database by amino acid count instead",typeid(bool), (void *) &splitAA, "$"),
        PARAM_SUB_MAT(PARAM_SUB_MAT_ID,"--sub-mat", "Sub Matrix", "[file]\tAmino acid substitution matrix file",typeid(std::string),(void *) &scoringMatrixFile, ""),
        PARAM_SEARCH_MODE(PARAM_SEARCH_MODE_ID,"--search-mode", "Search mode", "[int]\tSearch mode. Global: 0 Local: 1 Local fast: 2",typeid(int), (void *) &searchMode, "^[0-2]\{1\}$"),
        PARAM_NO_COMP_BIAS_CORR(PARAM_NO_COMP_BIAS_CORR_ID,"--no-comp-bias-corr", "Compositional bias","Switch off local amino acid composition bias correction",typeid(bool), (void *) &compBiasCorrection, ""),
        PARAM_SPACED_KMER_MODE(PARAM_SPACED_KMER_MODE_ID,"--spaced-kmer-mode", "Spaced Kmer", "[int]\tSpaced kmers mode (use consecutive pattern). Disable: 0, Enable: 1",typeid(int), (void *) &spacedKmer,  "^[0-1]\{1\}" ),
        PARAM_KEEP_TEMP_FILES(PARAM_KEEP_TEMP_FILES_ID,"--keep-tmp-files", "Keep-tmp-files" ,"\tDo not delete temporary files.",typeid(bool),(void *) &keepTempFiles, ""),
// alignment
        PARAM_E(PARAM_E_ID,"-e", "E-value threshold", "Maximum e-value",typeid(float), (void *) &evalThr, "^[0-9]*(\\.[0-9]+)?$"),
        PARAM_C(PARAM_C_ID,"-c", "Coverage threshold", "Minimum alignment coverage",typeid(float), (void *) &covThr, "^0(\\.[0-9]+)?|1\\.0$"),
        PARAM_MAX_REJECTED(PARAM_MAX_REJECTED_ID,"--max-rejected", "Max Reject", "Maximum rejected alignments before alignment calculation for a query is aborted",typeid(int),(void *) &maxRejected, "^[1-9]\{1\}[0-9]*$"),
// clustering
        PARAM_MIN_SEQ_ID(PARAM_MIN_SEQ_ID_ID,"--min-seq-id", "Seq. Id Threshold","Minimum sequence identity of sequences in a cluster",typeid(float), (void *) &seqIdThr, "[0-9]*(\\.[0-9]+)?$"),
        PARAM_CLUSTER_MODE(PARAM_CLUSTER_MODE_ID,"--cluster-mode", "Cluster mode", "0 Setcover, 1 affinity clustering, 2 Greedy clustering by sequence length",typeid(int), (void *) &clusteringMode, "[0-8]\{1\}$"),
        PARAM_CASCADED(PARAM_CASCADED_ID,"--cascaded", "Cascaded clustering", "\tStart the cascaded instead of simple clustering workflow",typeid(bool), (void *) &cascaded, ""),
//affinity clustering
        PARAM_MAXITERATIONS(PARAM_MAXITERATIONS_ID,"--max-iterations", "Max iterations affinity", "[int]\tMaximum number of iterations in affinity propagation clustering",typeid(int), (void *) &maxIteration,  "^[1-9]\{1\}[0-9]*$"),
        PARAM_CONVERGENCEITERATIONS(PARAM_CONVERGENCEITERATIONS_ID,"--convergence_iterations", "Convergence iterations", "[int]\t Number of iterations the set of representatives has to stay constant",typeid(int), (void *) &convergenceIterations,  "^[1-9]\{1\}[0-9]*$"),
        PARAM_DAMPING(PARAM_DAMPING_ID,"--damping", "Damping", "Ratio of previous iteration entering values. Value between [0.5:1).",typeid(float), (void *) &dampingFactor, "^[0-9]*(\\.[0-9]+)?$"),
        PARAM_SIMILARITYSCORE(PARAM_SIMILARITYSCORE_ID,"--similarity-type", "Similarity type", "Type of score used for clustering [1:5]. 1=alignment score. 2=coverage 3=sequence identity 4=E-value 5= Score per Column ",typeid(int),(void *) &similarityScoreType,  "^[1-9]\{1\}[0-9]*$"),
        PARAM_PREFERENCE(PARAM_PREFERENCE_ID,"--preference", "Preference", "Preference value influences the number of clusters (default=0). High values lead to more clusters.",typeid(float), (void *) &preference, "^[0-9]*(\\.[0-9]+)?$"),
// logging
        PARAM_V(PARAM_V_ID,"-v", "Verbosity","Verbosity level: 0=NOTHING, 1=ERROR, 2=WARNING, 3=INFO",typeid(int), (void *) &verbosity, "^[0-3]\{1\}$"),
// create profile (HMM, PSSM)
        PARAM_PROFILE_TYPE(PARAM_PROFILE_TYPE_ID,"--profile-type", "Profile type", "[int]\tMPI Option: HMM 0 or PSSM",typeid(int),(void *) &profileMode,  "^[0-1]\{1\}$"),

// clustering workflow
        PARAM_RESTART(PARAM_RESTART_ID, "--restart", "Restart","[int]\tRestart the clustering workflow starting with alignment or clustering.\n"
                "\t\tThe value is in the range [1:3]: 1: restart from prefiltering  2: from alignment; 3: from clustering",typeid(int),(void *) &restart, "^[0-3]\{1\}$"),
        PARAM_STEP(PARAM_STEP_ID, "--step","Step","[int]\t\tRestart the step of the cascaded clustering. For values in [1:3], the resprective step number, 4 is only the set merging",typeid(int),(void *) &step, "^[0-4]\{1\}$"),
// search workflow
        PARAM_NUM_ITERATIONS(PARAM_NUM_ITERATIONS_ID, "--num-iterations", "Number search iterations","[int]\tSearch iterations",typeid(int),(void *) &numIterations, "^[1-9]\{1\}[0-9]*$"),
// Orfs
        PARAM_ORF_MIN_LENGTH(PARAM_ORF_MIN_LENGTH_ID, "--min-length", "Min orf length", "[int]\t\tMinimum length of open reading frame to be extracted from fasta file",typeid(int),(void *) &orfMinLength, "^[1-9]\{1\}[0-9]*$"),
        PARAM_ORF_MAX_LENGTH(PARAM_ORF_MAX_LENGTH_ID, "--max-length", "Max orf length", "[int]\t\tMaximum length of open reading frame to be extracted from fasta file.",typeid(int),(void *) &orfMaxLength, "^[1-9]\{1\}[0-9]*$"),
        PARAM_ORF_MAX_GAP(PARAM_ORF_MAX_GAP_ID, "--max-gaps", "Max orf gaps", "[int]\t\tMaximum number of gaps or unknown residues before an open reading frame is rejected",typeid(int),(void *) &orfMaxGaps, "^(0|[1-9]{1}[0-9]*)$"),
        PARAM_ORF_SKIP_INCOMPLETE(PARAM_ORF_SKIP_INCOMPLETE_ID,"--skip-incomplete", "Skip incomplete orfs", "\tSkip orfs that have only an end or only a start",typeid(bool),(void *) &orfSkipIncomplete, ""),
        PARAM_ORF_NUMERIC_INDICES(PARAM_ORF_NUMERIC_INDICES_ID,"--numeric-indices", "Use numeric indices", "\tUse numeric indices as the ffindex key instead of trying to parse fasta headers",typeid(bool),(void *) &orfUseNumericIndices, "")
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
    prefilter.push_back(PARAM_SPLIT_MODE);
    prefilter.push_back(PARAM_SEARCH_MODE);
    prefilter.push_back(PARAM_NO_COMP_BIAS_CORR);
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
    clustering.push_back(PARAM_THREADS);


    // find orf
    onlyverbosity.push_back(PARAM_V);

    // create profile db
    createprofiledb.push_back(PARAM_SUB_MAT);
    createprofiledb.push_back(PARAM_PROFILE_TYPE);
    createprofiledb.push_back(PARAM_THREADS);
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

Parameters::~Parameters(){
    createindex.clear();
    splitffindex.clear();
    extractorf.clear();
    onlyverbosity.clear();
    clustering.clear();
    alignment.clear();
    prefilter.clear();
    createprofiledb.clear();
}


int Parameters::compileRegex(regex_t * regex, const char * regexText){
    int status = regcomp(regex, regexText, REG_EXTENDED | REG_NEWLINE);
    if (status != 0 ){
        Debug(Debug::INFO) << "Error in regex " << regexText << "\n";
        EXIT(EXIT_FAILURE);
    }
    return 0;
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
    //ops.exceptions(std::ios::eofbit); // throw exception when parsing error
    std::vector<std::string> getFilename;
    size_t parametersFound = 0;
    for(int argIdx = 0; argIdx < argc; argIdx++ ){
        // it is a parameter if it starts with - or --
        if ((pargv[argIdx][0] == '-' && pargv[argIdx][1] == '-') || (pargv[argIdx][0] == '-')) {
            std::string parameter(pargv[argIdx]);
            for(size_t parIdx = 0; parIdx < par.size(); parIdx++){
                if(parameter.compare(par[parIdx].name) == 0) {
                    if (typeid(int) == par[parIdx].type) {
                        regex_t regex;
                        compileRegex(&regex, par[parIdx].regex);
                        int nomatch = regexec(&regex, pargv[argIdx+1], 0, NULL, 0);
                        regfree(&regex);
                        // if no match found or two matches found (we want exactly one match)
                        if (nomatch){
                            printUsageMessage(programUsageHeader, par);
                            Debug(Debug::INFO) << "Error in argument " << par[parIdx].name << "\n";
                            EXIT(EXIT_FAILURE);
                        }else{
                            *((int *) par[parIdx].value) = atoi(pargv[argIdx+1]);
                        }
                        argIdx++;
                    } else if (typeid(float) == par[parIdx].type) {
                        regex_t regex;
                        compileRegex(&regex, par[parIdx].regex);
                        int nomatch = regexec(&regex, pargv[argIdx+1], 0, NULL, 0);
                        regfree(&regex);
                        if (nomatch){
                            printUsageMessage(programUsageHeader, par);
                            Debug(Debug::INFO) << "Error in argument " << par[parIdx].name << "\n";
                            EXIT(EXIT_FAILURE);
                        }else{
                            *((float *) par[parIdx].value) = atof(pargv[argIdx+1]);
                        }
                        argIdx++;
                    } else if (typeid(std::string) == par[parIdx].type) {
                        std::string val(pargv[argIdx+1]);
                        if(val.length() != 0){
                            std::string * currVal = ((std::string *)par[parIdx].value);
                            currVal->assign( val );
                        }
                        argIdx++;
                    } else if (typeid(bool) == par[parIdx].type) {
                        bool * value = (bool *) par[parIdx].value;
                        // toggle Value
                        *value = !*value;
                    } else {
                        Debug(Debug::ERROR) << "Wrong parameter type in parseParameters. Please inform developer\n";
                        EXIT(EXIT_FAILURE);
                    }
                }
            }
            parametersFound++;
        } else { // it is a filename if its not a parameter
            getFilename.push_back(pargv[argIdx]);
        }
    }
    //TODO make remain parameter logic (
    //if(parametersFound + getFilename.size() < argc)

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
    sensitivity = 4;
    split = 1;
    splitMode = TARGET_DB_SPLIT;
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
    similarityScoreType=APC_SEQID;
    preference=0;

    // Clustering workflow
    restart = 0;
    step = 1;
    keepTempFiles = true;

    // create profile
    profileMode = PROFILE_MODE_HMM;

    // logging
    verbosity = Debug::INFO;

    //extractorfs
    orfMinLength = 1;
    orfMaxLength = SIZE_MAX;
    orfMaxGaps = SIZE_MAX;
    orfUseNumericIndices = false;
    orfSkipIncomplete = false;
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
