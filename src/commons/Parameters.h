// Written by Martin Steinegger Martin.Steinegger@campus.lmu.de
//
// Represents a parameter of MMseqs
//
#ifndef __MMseqs__Parameters__
#define __MMseqs__Parameters__

#include <iostream>
#include <vector>
#include <map>


#define PARAMETER(x) const static int x##_ID = __COUNTER__; \
                     const static MMseqsParameter x;


struct MMseqsParameter {
    const int uniqid;
    const char *name;
    const char *description;
    MMseqsParameter(int uid,const char * n,const char * d): uniqid(uid), name(n),description(d){}
};


class Parameters          // Parameters for gap penalties and pseudocounts
{
public:
    
    static const int SET_COVER = 0;
    static const int GREEDY = 1;
    static const int AFFINITY = 2;

    static const int APC_ALIGNMENTSCORE=1;
    static const int APC_COVERAGE=2;
    static const int APC_SEQID=3;
    static const int APC_EVAL=4;
    static const int APC_BITSCORE=5;


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
   
    
    std::string scoringMatrixFile;     // path to scoring matrix
    size_t maxSeqLen;                   // sequence length
    size_t maxResListLen;               // Maximal result list length per query
    int    verbosity;                   // log level
    int    querySeqType;                // Query sequence type (PROFILE, AMINOACIDE, NUCLEOTIDE)
    int    targetSeqType;               // Target sequence type (PROFILE, AMINOACIDE, NUCLEOTIDE)
    int    threads;                     // Amounts of threads
    bool   keepTempFiles;               // Do not delete temp files
    
    // PREFILTER
    float  sensitivity;                  // target sens
    int    kmerSize;                     // kmer size for the prefilter
    int    kmerScore;                    // kmer score for the prefilter
    int    alphabetSize;                 // alphabet size for the prefilter
    float  zscoreThr;                    // z score threshold for global matching
    bool   localSearch;                  // Local search type
    bool   compBiasCorrection;           // Aminoacid composiont correction
    bool   fastMode;                     // Search 20.000 times faster than BLAST in (local search only)
    bool   spacedKmer;                   // Spaced Kmers
    int    split;                        // Splite database in n equal Junks
    int    skip;                         // Skip amino acid positions
    
    // ALIGNMENT
    std::string ffindexPrefDB;         // prefilter database (input for alignment module)
    double  evalThr;                     // e-value threshold for acceptance
    double  covThr;                      // coverage threshold for acceptance
    int     maxRejected;                 // after n sequences that are above eval stop
    float  seqIdThr;                     // sequence identity threshold for acceptance

    // CLUSTERING
    std::string ffindexAlnDBBase;
    int    clusteringMode;
    int    validateClustering;
    bool   cascaded;
    //AFFINITYCLUSTERING
    int maxIteration;                   // Maximum number of iterations of affinity clustering.
    int convergenceIterations;          // Number of iterations the representatives have to stay constant.
    float dampingFactor;                  // Reduces oscillation. Value in range of 0.5< <1.
    int similarityScoreType;            // Type of score to use for affinity clustering. (1) alignment score. (2) coverage (3)sequence identity (4)E-value.
    double preference;                  //Preference value influences the number of clusters (default=0). High values lead to more clusters.

    //extractorf
    size_t orfMinLength;
    size_t orfMaxLength;
    size_t orfMaxGaps;
    bool   orfSkipIncomplete;

    // CLUSTERING WORKFLOW
    int restart;
    int step;
    
    void setDefaultPaths();
    void setDefaults();
    void serialize( std::ostream &stream );
    void deserialize( std::istream &stream );
    void parseParameters(int argc, const char* argv[],
                         std::string programUsageHeader,
                         std::vector<MMseqsParameter> parameters,
                         size_t requiredParameterCount,
                         bool printParameters = true);
    void printUsageMessage(std::string programUsageHeader,
                           std::vector<MMseqsParameter> parameters);
    void printParameters(int argc, const char* pargv[],
                         std::vector<MMseqsParameter> parameters);
    Parameters();
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
    //    "--split         \t[int]\tSplits target databases in n equal distrbuted junks (default=1)\n"
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
    PARAMETER(PARAM_SKIP);
    PARAMETER(PARAM_MAX_SEQS);
    PARAMETER(PARAM_SPLIT);
    PARAMETER(PARAM_SUB_MAT);
    PARAMETER(PARAM_SEARCH_MODE);
    PARAMETER(PARAM_NO_COMP_BIAS_CORR);
    PARAMETER(PARAM_FAST_MODE);
    PARAMETER(PARAM_SPACED_KMER_MODE);
    PARAMETER(PARAM_K_SCORE);
    PARAMETER(PARAM_KEEP_TEMP_FILES);
    // alignment
    PARAMETER(PARAM_E);
    PARAMETER(PARAM_C);
    PARAMETER(PARAM_MAX_REJECTED);
    PARAMETER(PARAM_MIN_SEQ_ID);
    // clustering
    PARAMETER(PARAM_G);
    PARAMETER(PARAM_A);
    PARAMETER(PARAM_CASCADED);
    //afinity clustering
    PARAMETER(PARAM_MAXITERATIONS);
    PARAMETER(PARAM_CONVERGENCEITERATIONS);
    PARAMETER(PARAM_DAMPING);
    PARAMETER(PARAM_SIMILARITYSCORE);
    PARAMETER(PARAM_PREFERENCE);

    // logging
    PARAMETER(PARAM_V);
    // clustering workflow
    PARAMETER(PARAM_RESTART);
    PARAMETER(PARAM_STEP);
    // extractorfs
    PARAMETER(PARAM_ORF_MIN_LENGTH);
    PARAMETER(PARAM_ORF_MAX_LENGTH);
    PARAMETER(PARAM_ORF_MAX_GAP);
    PARAMETER(PARAM_ORF_SKIP_INCOMPLETE);

};

#endif /* defined(__MMseqs__Parameters__) */
