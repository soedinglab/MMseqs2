// Written by Martin Steinegger martin.steinegger@mpibpc.mpg.de
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
    bool wasSet;
    MMseqsParameter(int uid, const char * n, const char *display,
                    const char * d, const std::type_info &hash, void * value, const char * regex):
                    uniqid(uid), name(n), display(display), description(d), type(hash), value(value), regex(regex), wasSet(false){}
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
    float  sensitivity;                  // target sens
    int    kmerSize;                     // kmer size for the prefilter
    int    kmerScore;                    // kmer score for the prefilter
    int    alphabetSize;                 // alphabet size for the prefilter
    float  zscoreThr;                    // z score threshold for global matching
    int    searchMode;                   // Local search type
    bool   profile;                      // using profile information
    bool   nucl;                         // using nucl informatoin
    int    compBiasCorrection;           // Aminoacid composiont correction
    int    diagonalScoring;              // switch diagonal scoring
    int    minDiagScoreThr;              // min diagonal score
    int    spacedKmer;                   // Spaced Kmers
    int    split;                        // Split database in n equal chunks
    int    splitMode;                    // Split by query or target DB (MPI only)
    bool   splitAA;                      // Split database by amino acid count instead

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
    std::string forwardFrames;
    std::string reverseFrames;

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

    // gff2ffindex
    std::string gffType;

    // translate nucleotide
    int translationTable;

    // addSequences
    int minSequences;
    
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

    PARAMETER(PARAM_S);
    PARAMETER(PARAM_K);
    PARAMETER(PARAM_THREADS);
    PARAMETER(PARAM_ALPH_SIZE);
    PARAMETER(PARAM_MAX_SEQ_LEN);
    PARAMETER(PARAM_PROFILE);
    PARAMETER(PARAM_NUCL);
    PARAMETER(PARAM_DIAGONAL_SCORING);
    PARAMETER(PARAM_MIN_DIAG_SCORE);
    PARAMETER(PARAM_K_SCORE);
    PARAMETER(PARAM_MAX_SEQS);
    PARAMETER(PARAM_SPLIT);
    PARAMETER(PARAM_SPLIT_MODE);
    PARAMETER(PARAM_SPLIT_AMINOACID);
    PARAMETER(PARAM_SUB_MAT);
    PARAMETER(PARAM_SEARCH_MODE);
    PARAMETER(PARAM_NO_COMP_BIAS_CORR);
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
    // logging
    PARAMETER(PARAM_V);
    std::vector<MMseqsParameter> clustering;

    // create profile (HMM, PSSM)
    PARAMETER(PARAM_PROFILE_TYPE);

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
    PARAMETER(PARAM_ORF_FORWARD_FRAMES);
    PARAMETER(PARAM_ORF_REVERSE_FRAMES);

    // createdb
    PARAMETER(PARAM_USE_HEADER); // also used by extractorf
    PARAMETER(PARAM_ID_OFFSET);  // same

    // rebuildfasta
    PARAMETER(PARAM_USE_HEADER_FILE)

    // gff2ffindex
    PARAMETER(PARAM_GFF_TYPE)

    // translate_nucleotide
    PARAMETER(PARAM_TRANSLATION_TABLE)

    // addsequences
    PARAMETER(PARAM_MIN_SEQUENCES)

    std::vector<MMseqsParameter> empty;

    std::vector<MMseqsParameter> onlyverbosity;
    std::vector<MMseqsParameter> createprofiledb;
    std::vector<MMseqsParameter> extractorf;
    std::vector<MMseqsParameter> splitffindex;
    std::vector<MMseqsParameter> createindex;
    std::vector<MMseqsParameter> formatalignment;
    std::vector<MMseqsParameter> createdb;
    std::vector<MMseqsParameter> rebuildfasta;
    std::vector<MMseqsParameter> gff2ffindex;

    std::vector<MMseqsParameter> searchworkflow;
    std::vector<MMseqsParameter> clusteringWorkflow;
    std::vector<MMseqsParameter> clusterUpdate;
    std::vector<MMseqsParameter> translateNucleotide;
    std::vector<MMseqsParameter> addSequences;

    std::vector <MMseqsParameter> combineList(std::vector < MMseqsParameter > par1,
                                              std::vector < MMseqsParameter > par2);

    std::string createParameterString(std::vector < MMseqsParameter > vector);

    int compileRegex(regex_t *regex, const char *regexText);

};

#endif
