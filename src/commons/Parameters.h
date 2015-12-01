// Written by Martin Steinegger Martin.Steinegger@campus.lmu.de
//
// Represents a parameter of MMseqs
//
#ifndef MMSEQS_PARAMETERS
#define MMSEQS_PARAMETERS
#include <string>
#include <vector>
#include <map>
#include <typeinfo>
#include <regex.h>

#define PARAMETER(x) const static int x##_ID = __COUNTER__; \
    				 MMseqsParameter x;


struct MMseqsParameter {
    const int uniqid;
    const char *name;
    const char *display;
    const char *description;
    const std::type_info &type;
    void * value;
    const char * regex;
    MMseqsParameter(int uid,const char * n, const char *display,
                    const char * d, const std::type_info &hash, void * value, const char * regex):
                    uniqid(uid), name(n), display(display), description(d), type(hash), value(value), regex(regex){}
};


class Parameters          // Parameters for gap penalties and pseudocounts
{
public:
    static const int SEARCH_GLOBAL = 0;
    static const int SEARCH_LOCAL = 1;
    static const int SEARCH_LOCAL_FAST = 2;

    static const int PROFILE_MODE_HMM = 0;
    static const int PROFILE_MODE_PSSM = 1;

    static const int SET_COVER = 0;
    static const int CONNECTED_COMPONENT = 1;
    static const int GREEDY = 2;
    static const int AFFINITY = 3;


    static const int APC_ALIGNMENTSCORE=1;
    static const int APC_COVERAGE=2;
    static const int APC_SEQID=3;
    static const int APC_EVAL=4;
    static const int APC_BITSCORE=5;

    static const int TARGET_DB_SPLIT = 0;
    static const int QUERY_DB_SPLIT = 1;
    // COMMON
    const char** argv;            //command line parameters
    char argc;              //dimension of argv
    
    // path to databases
    std::string db1;
    std::string db1Index;
    
    std::string db2;
    std::string db2Index;

    std::string db3;
    std::string db3Index;
    
    std::string db4;
    std::string db4Index;
    
    std::string db5;
    std::string db5Index;

    std::string mmdir;

    std::string scoringMatrixFile;       // path to scoring matrix
    size_t maxSeqLen;                    // sequence length
    size_t maxResListLen;                // Maximal result list length per query
    int    verbosity;                    // log level
    int    querySeqType;                 // Query sequence type (PROFILE, AMINOACIDE, NUCLEOTIDE)
    int    targetSeqType;                // Target sequence type (PROFILE, AMINOACIDE, NUCLEOTIDE)
    int    threads;                      // Amounts of threads
    bool   keepTempFiles;                // Do not delete temp files
    
    // PREFILTER
    int    sensitivity;                  // target sens
    int    kmerSize;                     // kmer size for the prefilter
    int    kmerScore;                    // kmer score for the prefilter
    int    alphabetSize;                 // alphabet size for the prefilter
    float  zscoreThr;                    // z score threshold for global matching
    int    searchMode;                   // Local search type
    bool   profile;                      // using profile information
    bool   nucl;                         // using nucl informatoin
    bool   compBiasCorrection;           // Aminoacid composiont correction
    bool   mask;                         // Mask kmer index or sequence
    int    spacedKmer;                   // Spaced Kmers
    int    split;                        // Split database in n equal chunks
    int    splitMode;                    // Split by query or target DB (MPI only)
    bool   splitAA;                      // Split database by amino acid count instead
    int    skip;                         // Skip amino acid positions
    
    // ALIGNMENT
    std::string ffindexPrefDB;           // prefilter database (input for alignment module)
    float  evalThr;                      // e-value threshold for acceptance
    float  covThr;                       // coverage threshold for acceptance
    int    maxRejected;                  // after n sequences that are above eval stop
    float  seqIdThr;                     // sequence identity threshold for acceptance

    // CLUSTERING
    std::string ffindexAlnDBBase;
    int    clusteringMode;
    int    validateClustering;
    bool   cascaded;

    // SEARCH WORKFLOW
    int numIterations;

    //CLUSTERING
    int maxIteration;                   // Maximum depth of breadth first search in connected component
    int similarityScoreType;            // Type of score to use for reassignment 1=alignment score. 2=coverage 3=sequence identity 4=E-value 5= Score per Column

    //extractorf
    int orfMinLength;
    int orfMaxLength;
    int orfMaxGaps;
    bool   orfSkipIncomplete;

    // CREATE PROFILE
    int profileMode;

    // createdb
    bool useHeader;
    int identifierOffset;

    // rebuildfasta
    bool useHeaderFile;

    // CLUSTERING WORKFLOW
    int restart;
    int step;
    
    void setDefaultPaths();
    void setDefaults();
    void serialize( std::ostream &stream );
    void deserialize( std::istream &stream );
    void parseParameters(int argc, const char* argv[],
                         std::string programUsageHeader,
                         std::vector<MMseqsParameter> par,
                         size_t requiredParameterCount,
                         bool printParameters = true);
    void printUsageMessage(std::string programUsageHeader,
                           std::vector<MMseqsParameter> parameters);
    void printParameters(int argc, const char* pargv[],
                         std::vector<MMseqsParameter> par);
    Parameters();
    ~Parameters();
    // parameter constants
    //    "-s              \t[float]\tSensitivity in the range [1:9] (default=4).\n"
    //    "-k              \t[int]\tk-mer size in the range [4:7] (default=6).\n"
    //    "-cpu            \t[int]\tNumber of cores used for the computation (default=all cores).\n"
    //    "--alph-size     \t[int]\tAmino acid alphabet size (default=21).\n"
    //    "--z-score       \t[float]\tZ-score threshold (default: 50.0).\n"
    //    "--max-seq-len   \t[int]\tMaximum sequence length (default=50000).\n"
    //    "--profile       \t\tHMM Profile input.\n"
    //    "--nucl          \t\tNucleotide sequences input.\n"
    //    "--search-mode   \t[int]\tSearch mode loc: 1 glob: 2 (default=1).\n"
    //    "--no-comp-bias-corr \t\tSwitch off local amino acid composition bias correction.\n"
    //    "--no-spaced-kmer \t\tSwitch off spaced kmers (consecutive pattern).\n"
    //    "--split         \t[int]\tSplits target databases in n equally distributed chunks (default=1)\n"
    //    "--threads       \t[int]\tNumber of threads used to compute. (Default=all cpus)\n"
    //    "--max-seqs      \t[int]\tMaximum result sequences per query (default=300).\n"
    //    "--skip          \t[int]\tNumber of skipped k-mers during the index table generation.\n"
    //    "--sub-mat       \t[file]\tAmino acid substitution matrix file.\n"
    //    "-v              \t[int]\tVerbosity level: 0=NOTHING, 1=ERROR, 2=WARNING, 3=INFO (default=3).\n");
    // alignment
    //    "-e          \t[float]\tMaximum e-value (default=0.01).\n"
    //    "-c          \t[float]\tMinimum alignment coverage (default=0.8).\n"
    //    "--max-rejected\t[int]\tMaximum rejected alignments before alignment calculation for a query is aborted. (default=INT_MAX)\n"
    // clustering
    //    "-g              \t[int]\tgreedy clustering by sequence length (default: set cover clustering algorithm).\n"
    //    "-a              \t[int]\taffinity clustering (default: set cover clustering algorithm).\n"
    //    "--min-seq-id    \t[float]\tMinimum sequence identity of query to target cluster (default = 0.0)\n"
    PARAMETER(PARAM_S);
    PARAMETER(PARAM_K);
    PARAMETER(PARAM_THREADS);
    PARAMETER(PARAM_ALPH_SIZE);
    PARAMETER(PARAM_MAX_SEQ_LEN);
    PARAMETER(PARAM_PROFILE);
    PARAMETER(PARAM_NUCL);
    PARAMETER(PARAM_Z_SCORE);
    PARAMETER(PARAM_K_SCORE);
    PARAMETER(PARAM_SKIP);
    PARAMETER(PARAM_MAX_SEQS);
    PARAMETER(PARAM_SPLIT);
    PARAMETER(PARAM_SPLIT_MODE);
    PARAMETER(PARAM_SPLIT_AMINOACID);
    PARAMETER(PARAM_SUB_MAT);
    PARAMETER(PARAM_SEARCH_MODE);
    PARAMETER(PARAM_NO_COMP_BIAS_CORR);
    PARAMETER(PARAM_MASK);
    PARAMETER(PARAM_SPACED_KMER_MODE);
    PARAMETER(PARAM_KEEP_TEMP_FILES);
    std::vector<MMseqsParameter> prefilter;

    // alignment
    PARAMETER(PARAM_E);
    PARAMETER(PARAM_C);
    PARAMETER(PARAM_MAX_REJECTED);
    PARAMETER(PARAM_MIN_SEQ_ID);
    std::vector<MMseqsParameter> alignment;
    // clustering
    PARAMETER(PARAM_CLUSTER_MODE);

    PARAMETER(PARAM_CASCADED);

    //afinity clustering
    PARAMETER(PARAM_MAXITERATIONS);
    PARAMETER(PARAM_SIMILARITYSCORE);
    std::vector<MMseqsParameter> clustering;
    // create profile (HMM, PSSM)
    PARAMETER(PARAM_PROFILE_TYPE);
    // logging
    PARAMETER(PARAM_V);
    // clustering workflow
    PARAMETER(PARAM_RESTART);
    PARAMETER(PARAM_STEP);
    // search workflow
    PARAMETER(PARAM_NUM_ITERATIONS);
    // extractorfs
    PARAMETER(PARAM_ORF_MIN_LENGTH);
    PARAMETER(PARAM_ORF_MAX_LENGTH);
    PARAMETER(PARAM_ORF_MAX_GAP);
    PARAMETER(PARAM_ORF_SKIP_INCOMPLETE);

    // createdb
    PARAMETER(PARAM_USE_HEADER); // also used by extractorf
    PARAMETER(PARAM_ID_OFFSET);  // same

    // rebuildfasta
    PARAMETER(PARAM_USE_HEADER_FILE)
    
    std::vector<MMseqsParameter> onlyverbosity;
    std::vector<MMseqsParameter> createprofiledb;
    std::vector<MMseqsParameter> extractorf;
    std::vector<MMseqsParameter> splitffindex;
    std::vector<MMseqsParameter> createindex;
    std::vector<MMseqsParameter> formatalignment;
    std::vector<MMseqsParameter> createdb;
    std::vector<MMseqsParameter> rebuildfasta;

    std::vector <MMseqsParameter> combineList(std::vector < MMseqsParameter > par1,
                                              std::vector < MMseqsParameter > par2);

    std::string createParameterString(std::vector < MMseqsParameter > vector);

    int compileRegex(regex_t *regex, const char *regexText);

};

#endif
