// Written by Martin Steinegger martin.steinegger@mpibpc.mpg.de
//
// Represents a parameter of MMseqs
//
#ifndef MMSEQS_PARAMETERS
#define MMSEQS_PARAMETERS
#include <string>
#include <vector>
#include <typeinfo>
#include "Command.h"

#define PARAMETER(x) const static int x##_ID = __COUNTER__; \
    				 MMseqsParameter x;

struct MMseqsParameter {
    const char *name;
    const char *display;
    const char *description;
    const std::type_info &type;
    void * value;
    const char * regex;
    const int uniqid;
    int category;
    bool wasSet;

    static const int COMMAND_PREFILTER = 1;
    static const int COMMAND_ALIGN = 2;
    static const int COMMAND_CLUST = 4;
    static const int COMMAND_COMMON = 8;
    static const int COMMAND_PROFILE = 16;
    static const int COMMAND_MISC = 32;
    static const int COMMAND_CLUSTLINEAR = 64;
    static const int COMMAND_EXPERT = 128;


    MMseqsParameter(int uid, const char * n, const char *display,
                    const char * d, const std::type_info &hash,
                    void * value, const char * regex, int category = COMMAND_MISC):
            name(n), display(display), description(d), type(hash), value(value),
            regex(regex), uniqid(uid), category(category), wasSet(false){}
};


class Parameters {
public:

    static const unsigned int ALIGNMENT_MODE_FAST_AUTO = 0;
    static const unsigned int ALIGNMENT_MODE_SCORE_ONLY = 1;
    static const unsigned int ALIGNMENT_MODE_SCORE_COV = 2;
    static const unsigned int ALIGNMENT_MODE_SCORE_COV_SEQID = 3;

    // format alignment
    static const int FORMAT_ALIGNMENT_BLAST_TAB = 0;
    static const int FORMAT_ALIGNMENT_PAIRWISE  = 1;
    static const int FORMAT_ALIGNMENT_BLAST_WITH_LEN = 2;
    // NOT IMPLEMENTED YET
    static const int FORMAT_ALIGNMENT_SAM       = 99;

    // convertprofiledb
    static const int PROFILE_MODE_HMM = 0;
    static const int PROFILE_MODE_PSSM = 1;
    static const int PROFILE_MODE_HMM3 = 2;

    // clustering
    static const int SET_COVER = 0;
    static const int CONNECTED_COMPONENT = 1;
    static const int GREEDY = 2;
    static const int GREEDY_MEM = 3;

    // clustering
    static const int APC_ALIGNMENTSCORE=1;
    static const int APC_SEQID=2;
    // split mode
    static const int TARGET_DB_SPLIT = 0;
    static const int QUERY_DB_SPLIT = 1;
    static const int DETECT_BEST_DB_SPLIT = 2;

    static const int TAXONOMY_NO_LCA = 0;
    static const int TAXONOMY_SINGLE_SEARCH = 1;
    static const int TAXONOMY_2BLCA = 2;

    static const int PARSE_VARIADIC = 1;
    static const int PARSE_REST = 2;

    static std::string getSplitModeName(int splitMode) {
        switch (splitMode) {
            case 0: return "Target";
            case 1: return "Query";
            case 2: return "Auto";
            default: return "Error";
        }
    };

    // split
    static const int AUTO_SPLIT_DETECTION = 0;

    static const int MAX_SEQ_LEN = 32000;

    // extractalignedregion
    static const int EXTRACT_QUERY  = 1;
    static const int EXTRACT_TARGET = 2;

    static const int CLUST_HASH_DEFAULT_ALPH_SIZE = 3;
    static const int CLUST_LINEAR_DEFAULT_ALPH_SIZE = 13;
    static const int CLUST_LINEAR_DEFAULT_K = 0;
    static const int CLUST_LINEAR_KMER_PER_SEQ = 0;


    // cov mode
    static const int COV_MODE_BIDIRECTIONAL  = 0;
    static const int COV_MODE_TARGET = 1;
    static const int COV_MODE_QUERY = 2;


    // rescorediagonal
    static const int RESCORE_MODE_HAMMING = 0;
    static const int RESCORE_MODE_SUBSTITUTION = 1;
    static const int RESCORE_MODE_ALIGNMENT = 2;

    // header type
    static const int HEADER_TYPE_UNICLUST = 1;
    static const int HEADER_TYPE_METACLUST = 2;

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

    std::string db6;
    std::string db6Index;

    std::vector<std::string> filenames;

    const char** restArgv;
    int restArgc;

    std::string scoringMatrixFile;       // path to scoring matrix
    size_t maxSeqLen;                    // sequence length
    size_t maxResListLen;                // Maximal result list length per query
    int    verbosity;                    // log level
//    int    querySeqType;                 // Query sequence type (PROFILE, AMINOACIDE, NUCLEOTIDE)
//    int    targetSeqType;                // Target sequence type (PROFILE, AMINOACIDE, NUCLEOTIDE)
    int    threads;                      // Amounts of threads
    bool   removeTmpFiles;               // Do not delete temp files
    bool   includeIdentity;              // include identical ids as hit

    // PREFILTER
    float  sensitivity;                  // target sens
    int    kmerSize;                     // kmer size for the prefilter
    int    kmerScore;                    // kmer score for the prefilter
    int    alphabetSize;                 // alphabet size for the prefilter
    //bool   queryProfile;                 // using queryProfile information
    //bool   targetProfile;                // using targetProfile information
    int    compBiasCorrection;           // Aminoacid composiont correction
    int    diagonalScoring;              // switch diagonal scoring
    int    maskMode;                     // mask low complex areas

    int    minDiagScoreThr;              // min diagonal score
    int    spacedKmer;                   // Spaced Kmers
    int    split;                        // Split database in n equal chunks
    int    splitMode;                    // Split by query or target DB
    int    splitMemoryLimit;             // Maximum amount of memory a split can use
    bool   splitAA;                      // Split database by amino acid count instead
    size_t resListOffset;                // Offsets result list
    bool   noPreload;                    // Do not preload database into memory
    bool   earlyExit;                    // Exit immediately after writing the result

    // ALIGNMENT
    int alignmentMode;                   // alignment mode 0=fastest on parameters,
                                         // 1=score only, 2=score, cov, start/end pos, 3=score, cov, start/end pos, seq.id,
    float  evalThr;                      // e-value threshold for acceptance
    float  covThr;                       // coverage query&target threshold for acceptance
    int    covMode;                      // coverage target threshold for acceptance

    int    maxRejected;                  // after n sequences that are above eval stop
    int    maxAccept;                    // after n accepted sequences stop
    float  seqIdThr;                     // sequence identity threshold for acceptance
    bool   addBacktrace;                 // store backtrace string (M=Match, D=deletion, I=insertion)
    bool   realign;                      // realign hit with more conservative score
	
    // workflow
    std::string runner;

    // CLUSTERING
    int    clusteringMode;
    int    clusterSteps;
    bool   cascaded;

    // SEARCH WORKFLOW
    int numIterations;
    float startSens;
    int sensSteps;

    // easysearch
    bool greedyBestHits;

    //CLUSTERING
    int maxIteration;                   // Maximum depth of breadth first search in connected component
    int similarityScoreType;            // Type of score to use for reassignment 1=alignment score. 2=coverage 3=sequence identity 4=E-value 5= Score per Column

    //mergecluster
    std::string DBfile ;

    //extractorfs
    int orfMinLength;
    int orfMaxLength;
    int orfMaxGaps;
    int contigStartMode;
    int contigEndMode;
    int orfStartMode;
    std::string forwardFrames;
    std::string reverseFrames;
    bool useAllTableStarts;

    // convertprofiledb
    int profileMode;

    // format alignment
    int formatAlignmentMode;            // BLAST_TAB, PAIRWISE or SAM
    bool dbOut;

    // rescorediagonal
    int rescoreMode;
    bool filterHits;
    bool globalAlignment;

    // result2msa
    bool allowDeletion;
    bool addInternalId;
    bool compressMSA;
    bool summarizeHeader;
    std::string summaryPrefix;
    bool omitConsensus;
    bool skipQuery;

    // convertmsa
    int identifierField;

    // msa2profile
    int matchMode;
    float matchRatio;

    // result2profile
    float filterMaxSeqId;
    float evalProfile;
    int filterMsa;
    float qsc;
    float qid;
    float cov;
    int Ndiff;
    bool wg;
    float pca;
    float pcb;

    // sequence2profile
    float neff;
    float tau;

    // createtsv
    bool firstSeqRepr;
    
    //result2stats
    std::string stat;

    // linearcluster
    int kmersPerSequence;
    bool includeOnlyExtendable;
    int hashShift;

    // indexdb
    bool includeHeader;

    // createdb
    int identifierOffset;
    bool splitSeqByLen;

    // convert2fasta
    bool useHeaderFile;

    // result2flat
    bool useHeader;

    // gff2db
    std::string gffType;

    // translate nucleotide
    int translationTable;
    bool addOrfStop;

    // createseqfiledb
    int minSequences;
    int maxSequences;
    bool hhFormat;

    // filterDb
    int filterColumn;
    std::string filterColumnRegex;
	std::string filteringFile;
	std::string mappingFile;
	bool positiveFilter;
	bool trimToOneColumn;
    int extractLines;
    float compValue;
    std::string compOperator;
    int sortEntries;
    bool beatsFirst;
    std::string joinDB;

    //aggregate
    std::string mode ;

    // mergedbs
    std::string mergePrefixes;

    // summarizetabs
    float overlap;
    int msaType;

    // extractalignedregion
    int extractMode;

    // convertkb
    std::string kbColumns;
    
    // concatdbs
    bool preserveKeysB;
    
    // diff
    bool useSequenceId;

    //prefixid
    std::string prefix;
    bool tsvOut;

    // clusterUpdate;
    bool recoverDeleted;

    // summarize headers
    int headerType;

    // lca
    std::string lcaRanks;
    std::string blacklist;

    // taxonomy
    int lcaMode;

    static Parameters& getInstance()
    {
        if (instance == NULL) {
            initInstance();
        }
        return *instance;
    }
    static void initInstance() {
        new Parameters;
    }
    
    void setDefaults();
    void parseParameters(int argc, const char* argv[],
                         const Command& command,
                         size_t requiredParameterCount,
                         bool printParameters = true,
                         int parseFlags = 0,
                         int outputFlags = 0);
    void printUsageMessage(const Command& command,
                           int outputFlag);
    void printParameters(int argc, const char* pargv[],
                         const std::vector<MMseqsParameter> &par);
	
	std::vector<MMseqsParameter> removeParameter(const std::vector<MMseqsParameter>& par, const MMseqsParameter& x);

    PARAMETER(PARAM_S)
    PARAMETER(PARAM_K)
    PARAMETER(PARAM_THREADS)
    PARAMETER(PARAM_ALPH_SIZE)
    PARAMETER(PARAM_MAX_SEQ_LEN)
//    PARAMETER(PARAM_QUERY_PROFILE)
//    PARAMETER(PARAM_TARGET_PROFILE)
    //PARAMETER(PARAM_NUCL)
    PARAMETER(PARAM_DIAGONAL_SCORING)
    PARAMETER(PARAM_MASK_RESIDUES)

    PARAMETER(PARAM_MIN_DIAG_SCORE)
    PARAMETER(PARAM_K_SCORE)
    PARAMETER(PARAM_MAX_SEQS)
    PARAMETER(PARAM_SPLIT)
    PARAMETER(PARAM_SPLIT_MODE)
    PARAMETER(PARAM_SPLIT_MEMORY_LIMIT)
    PARAMETER(PARAM_SPLIT_AMINOACID)
    PARAMETER(PARAM_SUB_MAT)
    PARAMETER(PARAM_NO_COMP_BIAS_CORR)
    PARAMETER(PARAM_SPACED_KMER_MODE)
    PARAMETER(PARAM_REMOVE_TMP_FILES)
    PARAMETER(PARAM_INCLUDE_IDENTITY)
    PARAMETER(PARAM_RES_LIST_OFFSET)
    PARAMETER(PARAM_NO_PRELOAD)
    PARAMETER(PARAM_EARLY_EXIT)
    std::vector<MMseqsParameter> prefilter;

    // alignment
    PARAMETER(PARAM_ALIGNMENT_MODE)
    PARAMETER(PARAM_E)
    PARAMETER(PARAM_C)
    PARAMETER(PARAM_COV_MODE)
    PARAMETER(PARAM_MAX_REJECTED)
    PARAMETER(PARAM_MAX_ACCEPT)
    PARAMETER(PARAM_ADD_BACKTRACE)
    PARAMETER(PARAM_REALIGN)
    PARAMETER(PARAM_MIN_SEQ_ID)
    std::vector<MMseqsParameter> align;

    // clustering
    PARAMETER(PARAM_CLUSTER_MODE)
    PARAMETER(PARAM_CLUSTER_STEPS)
    PARAMETER(PARAM_CASCADED)

    // affinity clustering
    PARAMETER(PARAM_MAXITERATIONS)
    PARAMETER(PARAM_SIMILARITYSCORE)

    //Merge Clusters
    PARAMETER(PARAM_BY_DB)

    // logging
    PARAMETER(PARAM_V)
    std::vector<MMseqsParameter> clust;

    // create profile (HMM, PSSM)
    PARAMETER(PARAM_PROFILE_TYPE)

    // format alignment
    PARAMETER(PARAM_FORMAT_MODE)
    PARAMETER(PARAM_DB_OUTPUT)

    // rescoremode
    PARAMETER(PARAM_RESCORE_MODE)
    PARAMETER(PARAM_FILTER_HITS)
    PARAMETER(PARAM_GLOBAL_ALIGNMENT)

    // result2msa
    PARAMETER(PARAM_ALLOW_DELETION)
    PARAMETER(PARAM_ADD_INTERNAL_ID)
    PARAMETER(PARAM_COMPRESS_MSA)
    PARAMETER(PARAM_SUMMARIZE_HEADER)
    PARAMETER(PARAM_SUMMARY_PREFIX)
    PARAMETER(PARAM_OMIT_CONSENSUS)
    PARAMETER(PARAM_SKIP_QUERY)

    // convertmsa
    PARAMETER(PARAM_IDENTIFIER_FIELD)

    // msa2profile
    PARAMETER(PARAM_MATCH_MODE)
    PARAMETER(PARAM_MATCH_RATIO)

    // result2profile
    PARAMETER(PARAM_E_PROFILE)
    PARAMETER(PARAM_FILTER_MSA)
    PARAMETER(PARAM_FILTER_MAX_SEQ_ID)
    PARAMETER(PARAM_FILTER_QSC)
    PARAMETER(PARAM_FILTER_QID)
    PARAMETER(PARAM_FILTER_COV)
    PARAMETER(PARAM_FILTER_NDIFF)
    PARAMETER(PARAM_WG)
    PARAMETER(PARAM_PCA)
    PARAMETER(PARAM_PCB)

    // sequence2profile
    PARAMETER(PARAM_NEFF)
    PARAMETER(PARAM_TAU)

    // createtsv
    PARAMETER(PARAM_FIRST_SEQ_REP_SEQ)
    
    // result2stat
    PARAMETER(PARAM_STAT)

    // linearcluster
    PARAMETER(PARAM_KMER_PER_SEQ)
    PARAMETER(PARAM_INCLUDE_ONLY_EXTENDABLE)
    PARAMETER(PARAM_HASH_SHIFT)

    // workflow
    PARAMETER(PARAM_RUNNER)

    // search workflow
    PARAMETER(PARAM_NUM_ITERATIONS)
    PARAMETER(PARAM_START_SENS)
    PARAMETER(PARAM_SENS_STEPS)

    // easysearch
    PARAMETER(PARAM_GREEDY_BEST_HITS)

    // extractorfs
    PARAMETER(PARAM_ORF_MIN_LENGTH)
    PARAMETER(PARAM_ORF_MAX_LENGTH)
    PARAMETER(PARAM_ORF_MAX_GAP)
    PARAMETER(PARAM_CONTIG_START_MODE)
    PARAMETER(PARAM_CONTIG_END_MODE)
    PARAMETER(PARAM_ORF_START_MODE)
    PARAMETER(PARAM_ORF_FORWARD_FRAMES)
    PARAMETER(PARAM_ORF_REVERSE_FRAMES)
    PARAMETER(PARAM_USE_ALL_TABLE_STARTS)

    // indexdb
    PARAMETER(PARAM_INCLUDE_HEADER)

    // createdb
    PARAMETER(PARAM_USE_HEADER) // also used by extractorfs
    PARAMETER(PARAM_ID_OFFSET)  // same
    PARAMETER(PARAM_DONT_SPLIT_SEQ_BY_LEN)

    // convert2fasta
    PARAMETER(PARAM_USE_HEADER_FILE)

    // gff2db
    PARAMETER(PARAM_GFF_TYPE)

    // translate_nucleotide
    PARAMETER(PARAM_TRANSLATION_TABLE)
    PARAMETER(PARAM_ADD_ORF_STOP)

    // createseqfiledb
    PARAMETER(PARAM_MIN_SEQUENCES)
    PARAMETER(PARAM_MAX_SEQUENCES)
    PARAMETER(PARAM_HH_FORMAT)

    // filterDb
    PARAMETER(PARAM_FILTER_COL)
    PARAMETER(PARAM_FILTER_REGEX)
    PARAMETER(PARAM_FILTER_POS)
    PARAMETER(PARAM_FILTER_FILE)
    PARAMETER(PARAM_MAPPING_FILE)
    PARAMETER(PARAM_TRIM_TO_ONE_COL)
    PARAMETER(PARAM_EXTRACT_LINES)
    PARAMETER(PARAM_COMP_OPERATOR)
    PARAMETER(PARAM_COMP_VALUE)
    PARAMETER(PARAM_SORT_ENTRIES)
    PARAMETER(PARAM_BEATS_FIRST)
    PARAMETER(PARAM_JOIN_DB)

    //aggregate
    PARAMETER(PARAM_MODE)

    // concatdb
    PARAMETER(PARAM_PRESERVEKEYS)

    // diff
    PARAMETER(PARAM_USESEQID)

    // prefixid
    PARAMETER(PARAM_PREFIX)
    PARAMETER(PARAM_TSV)

    // summarize headers
    PARAMETER(PARAM_HEADER_TYPE)

    // mergedbs
    PARAMETER(PARAM_MERGE_PREFIXES)

    // summarizetabs
    PARAMETER(PARAM_OVERLAP)

    // extractdomains
    PARAMETER(PARAM_MSA_TYPE)

    // extract aligned region
    PARAMETER(PARAM_EXTRACT_MODE)

    // convertkb
    PARAMETER(PARAM_KB_COLUMNS)

    // clusterupdate
    PARAMETER(PARAM_RECOVER_DELETED)

    // lca
    PARAMETER(PARAM_LCA_RANKS)
    PARAMETER(PARAM_BLACKLIST)

    // taxonomy
    PARAMETER(PARAM_LCA_MODE)

    std::vector<MMseqsParameter> empty;
    std::vector<MMseqsParameter> rescorediagonal;
    std::vector<MMseqsParameter> onlyverbosity;
    std::vector<MMseqsParameter> createFasta;
    std::vector<MMseqsParameter> convertprofiledb;
    std::vector<MMseqsParameter> sequence2profile;

    std::vector<MMseqsParameter> result2profile;
    std::vector<MMseqsParameter> result2pp;
    std::vector<MMseqsParameter> result2msa;
    std::vector<MMseqsParameter> convertmsa;
    std::vector<MMseqsParameter> msa2profile;
    std::vector<MMseqsParameter> createtsv;
    std::vector<MMseqsParameter> result2stats;
    std::vector<MMseqsParameter> extractorfs;
    std::vector<MMseqsParameter> splitdb;
    std::vector<MMseqsParameter> indexdb;
    std::vector<MMseqsParameter> createindex;
    std::vector<MMseqsParameter> convertalignments;
    std::vector<MMseqsParameter> createdb;
    std::vector<MMseqsParameter> convert2fasta;
    std::vector<MMseqsParameter> result2flat;
    std::vector<MMseqsParameter> gff2ffindex;
    std::vector<MMseqsParameter> clusthash;
    std::vector<MMseqsParameter> kmermatcher;
    std::vector<MMseqsParameter> linclustworkflow;
    std::vector<MMseqsParameter> assemblerworkflow;
    std::vector<MMseqsParameter> easysearchworkflow;
    std::vector<MMseqsParameter> searchworkflow;
    std::vector<MMseqsParameter> clusteringWorkflow;
    std::vector<MMseqsParameter> clusterUpdateSearch;
    std::vector<MMseqsParameter> clusterUpdateClust;
    std::vector<MMseqsParameter> mergeclusters ;
    std::vector<MMseqsParameter> clusterUpdate;
    std::vector<MMseqsParameter> translatenucs;
    std::vector<MMseqsParameter> swapresult;
    std::vector<MMseqsParameter> createseqfiledb;
    std::vector<MMseqsParameter> filterDb;
    std::vector<MMseqsParameter> onlythreads;
    std::vector<MMseqsParameter> subtractdbs;
    std::vector<MMseqsParameter> diff;
    std::vector<MMseqsParameter> concatdbs;
    std::vector<MMseqsParameter> mergedbs;
    std::vector<MMseqsParameter> summarizeheaders;
    std::vector<MMseqsParameter> prefixid;
    std::vector<MMseqsParameter> summarizeresult;
    std::vector<MMseqsParameter> summarizetabs;
    std::vector<MMseqsParameter> extractdomains;
    std::vector<MMseqsParameter> extractalignedregion;
    std::vector<MMseqsParameter> convertkb;
    std::vector<MMseqsParameter> tsv2db;
    std::vector<MMseqsParameter> lca;
    std::vector<MMseqsParameter> taxonomy;
    std::vector<MMseqsParameter> profile2pssm;
    std::vector<MMseqsParameter> profile2cs;
    std::vector<MMseqsParameter> aggregate ;

    std::vector<MMseqsParameter> combineList(const std::vector<MMseqsParameter> &par1,
                                             const std::vector<MMseqsParameter> &par2);

    size_t hashParameter(const std::vector<std::string> &filenames, const std::vector<MMseqsParameter> &par);

    std::string createParameterString(const std::vector<MMseqsParameter> &vector, bool wasSet = false);

    void overrideParameterDescription(Command& command, int uid, const char* description, const char* regex = NULL, int category = 0);

protected:
    Parameters();
    static Parameters* instance;
    virtual ~Parameters() {};

private:
    Parameters(Parameters const&);
    void operator=(Parameters const&);
};

#endif
