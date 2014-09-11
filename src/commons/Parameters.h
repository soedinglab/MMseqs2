// Written by Martin Steinegger Martin.Steinegger@campus.lmu.de
//
// Represents a parameter of MMseqs
//
#ifndef __MMseqs__Parameters__
#define __MMseqs__Parameters__

#include <iostream>
#include <vector>
#include <map>


struct MMseqsParameter {
    const char *name;
    const char *description;
};


class Parameters          // Parameters for gap penalties and pseudocounts
{
public:
    
    static const int SET_COVER = 0;
    static const int GREEDY = 1;
    // COMMON
    char** argv;            //command line parameters
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
    
    // PREFILTER
    float  sensitivity;                  // target sens
    int    kmerSize;                     // kmer size for the prefilter
    int    alphabetSize;                 // alphabet size for the prefilter
    float  zscoreThr;                    // z score threshold for global matching
    bool   localSearch;                  // Local search type
    bool   compBiasCorrection;           // Aminoacid composiont correction
    bool   spacedKmer;                   // Spaced Kmers
    int    split;                        // Splite database in n equal Junks
    int    skip;                         // Skip amino acid positions
    
    // ALIGNMENT
    std::string ffindexPrefDB;         // prefilter database (input for alignment module)
    double  evalThr;                     // e-value threshold for acceptance
    double  covThr;                      // coverage threshold for acceptance
    int     maxRejected;                 // after n sequences that are above eval stop
    
    // CLUSTERING
    std::string ffindexAlnDBBase;
    int    clusteringMode;
    float  seqIdThr;
    int    validateClustering;
    bool   cascaded;
    
    void setDefaultPaths();
    void setDefaults();
    void serialize( std::ostream &stream );
    void deserialize( std::istream &stream );
    void parseParameters(int argc, char* argv[],
                         std::string programUsageHeader,
                         std::vector<MMseqsParameter> parameters,
                         size_t requiredParameterCount);
    void printUsageMessage(std::string programUsageHeader,
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
    //    "--min-seq-id    \t[float]\tMinimum sequence identity of sequences in a cluster (default = 0.0)\n"
    const constexpr static MMseqsParameter PARAM_S={"-s",                    "[float]\tSensitivity in the range [1:9]"};
    const constexpr static MMseqsParameter PARAM_K={"-k",                    "[int]\tk-mer size in the range [4:7]"};
    const constexpr static MMseqsParameter PARAM_THREADS={"--threads",        "[int]\tNumber of cores used for the computation"};
    const constexpr static MMseqsParameter PARAM_ALPH_SIZE={"--alph-size",    "[int]\tAmino acid alphabet size"};
    const constexpr static MMseqsParameter PARAM_MAX_SEQ_LEN={"--max-seq-len","[int]\tMaximum sequence length"};
    const constexpr static MMseqsParameter PARAM_PROFILE={"--profile",        "\tHMM Profile input"};
    const constexpr static MMseqsParameter PARAM_NUCL={"--nucl",              "\tNucleotide sequences input"};
    const constexpr static MMseqsParameter PARAM_Z_SCORE={"--z-score",        "[float]\tZ-score threshold "};
    const constexpr static MMseqsParameter PARAM_SKIP={"--skip",              "[int]\tNumber of skipped k-mers during the index table generation"};
    const constexpr static MMseqsParameter PARAM_MAX_SEQS={"--max-seqs",      "[int]\tMaximum result sequences per query"};
    const constexpr static MMseqsParameter PARAM_SPLIT={"--split",            "[int]\tSplits target databases in n equal distrbuted junks"};
    const constexpr static MMseqsParameter PARAM_SUB_MAT={"--sub-mat",        "[file]\tAmino acid substitution matrix file"};
    const constexpr static MMseqsParameter PARAM_SEARCH_MODE={"--search-mode","[int]\tSearch mode loc: 1 glob: 2"};
    const constexpr static MMseqsParameter PARAM_NO_COMP_BIAS_CORR={"--no-comp-bias-corr","Switch off local amino acid composition bias correction"};
    const constexpr static MMseqsParameter PARAM_NO_SPACED_KMER={"--no-spaced=kmer","Switch off spaced kmers (use consecutive pattern)"};
    // alignment
    const constexpr static MMseqsParameter PARAM_E={"-e",                          "Maximum e-value"};
    const constexpr static MMseqsParameter PARAM_C={"-c",                          "Minimum alignment coverage"};
    const constexpr static MMseqsParameter PARAM_MAX_REJECTED={"--max-rejected","Maximum rejected alignments before alignment calculation for a query is aborted"};
    // clustering
    const constexpr static MMseqsParameter PARAM_G={"-g","Greedy clustering by sequence length"};
    const constexpr static MMseqsParameter PARAM_MIN_SEQ_ID={"--min-seq-id","Minimum sequence identity of sequences in a cluster"};
    const constexpr static MMseqsParameter PARAM_CASCADED={"--cascaded", "\tStart the cascaded instead of simple clustering workflow"};
    // logging
    const constexpr static MMseqsParameter PARAM_V={"-v","Verbosity level: 0=NOTHING, 1=ERROR, 2=WARNING, 3=INFO"};
    
};

#endif /* defined(__MMseqs__Parameters__) */
