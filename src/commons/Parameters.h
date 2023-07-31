// Written by Martin Steinegger martin.steinegger@snu.ac.kr
//
// Represents parameters of MMseqs2
//
#ifndef MMSEQS_PARAMETERS
#define MMSEQS_PARAMETERS
#include <string>
#include <vector>
#include <map>
#include <typeinfo>
#include <cstddef>
#include <utility>
#include <cstdint>

#include "Command.h"
#include "MultiParam.h"

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
    unsigned int category;
    bool wasSet;

    static const unsigned int COMMAND_PREFILTER = 1;
    static const unsigned int COMMAND_ALIGN = 2;
    static const unsigned int COMMAND_CLUST = 4;
    static const unsigned int COMMAND_COMMON = 8;
    static const unsigned int COMMAND_PROFILE = 16;
    static const unsigned int COMMAND_MISC = 32;
    static const unsigned int COMMAND_CLUSTLINEAR = 64;
    static const unsigned int COMMAND_EXPERT = 128;
    static const unsigned int COMMAND_HIDDEN = 256;


    MMseqsParameter(int uid, const char * n, const char *display,
                    const char * d, const std::type_info &hash,
                    void * value, const char * regex, unsigned int category = COMMAND_MISC):
            name(n), display(display), description(d), type(hash), value(value),
            regex(regex), uniqid(uid), category(category), wasSet(false){}

    void addCategory(unsigned int cat) {
        category |= cat;
    }

    void removeCategory(unsigned int cat) {
        category &= ~cat;
    }

    void replaceCategory(unsigned int cat) {
        category = cat;
    }
};


class Parameters {
public:

    static const int DBTYPE_AMINO_ACIDS = 0;
    static const int DBTYPE_NUCLEOTIDES = 1;
    static const int DBTYPE_HMM_PROFILE = 2;
    //static const int DBTYPE_PROFILE_STATE_SEQ = 3;
    //static const int DBTYPE_PROFILE_STATE_PROFILE = 4;
    static const int DBTYPE_ALIGNMENT_RES = 5;
    static const int DBTYPE_CLUSTER_RES = 6;
    static const int DBTYPE_PREFILTER_RES = 7;
    static const int DBTYPE_TAXONOMICAL_RESULT = 8;
    static const int DBTYPE_INDEX_DB = 9;
    static const int DBTYPE_CA3M_DB = 10;
    static const int DBTYPE_MSA_DB = 11;
    static const int DBTYPE_GENERIC_DB = 12;
    static const int DBTYPE_OMIT_FILE = 13;
    static const int DBTYPE_PREFILTER_REV_RES = 14;
    static const int DBTYPE_OFFSETDB = 15;
    static const int DBTYPE_DIRECTORY = 16; // needed for verification
    static const int DBTYPE_FLATFILE = 17; // needed for verification
    static const int DBTYPE_SEQTAXDB = 18; // needed for verification
    static const int DBTYPE_STDIN = 19; // needed for verification
    static const int DBTYPE_URI = 20; // needed for verification

    static const unsigned int DBTYPE_EXTENDED_COMPRESSED = 1;
    static const unsigned int DBTYPE_EXTENDED_INDEX_NEED_SRC = 2;
    static const unsigned int DBTYPE_EXTENDED_CONTEXT_PSEUDO_COUNTS = 4;

    // don't forget to add new database types to DBReader::getDbTypeName and Parameters::PARAM_OUTPUT_DBTYPE

    static const int SEARCH_TYPE_AUTO = 0;
    static const int SEARCH_TYPE_PROTEIN = 1;
    static const int SEARCH_TYPE_TRANSLATED = 2;
    static const int SEARCH_TYPE_NUCLEOTIDES = 3;
    static const int SEARCH_TYPE_TRANS_NUCL_ALN = 4;
    // flag
    static const int SEARCH_MODE_FLAG_QUERY_AMINOACID = 1;
    static const int SEARCH_MODE_FLAG_TARGET_AMINOACID = 2;
    static const int SEARCH_MODE_FLAG_QUERY_TRANSLATED = 4;
    static const int SEARCH_MODE_FLAG_TARGET_TRANSLATED = 8;
    static const int SEARCH_MODE_FLAG_QUERY_PROFILE = 16;
    static const int SEARCH_MODE_FLAG_TARGET_PROFILE = 32;
    static const int SEARCH_MODE_FLAG_QUERY_NUCLEOTIDE = 64;
    static const int SEARCH_MODE_FLAG_TARGET_NUCLEOTIDE = 128;

    static const unsigned int ALIGNMENT_MODE_FAST_AUTO = 0;
    static const unsigned int ALIGNMENT_MODE_SCORE_ONLY = 1;
    static const unsigned int ALIGNMENT_MODE_SCORE_COV = 2;
    static const unsigned int ALIGNMENT_MODE_SCORE_COV_SEQID = 3;
    static const unsigned int ALIGNMENT_MODE_UNGAPPED = 4;

    static const unsigned int ALIGNMENT_OUTPUT_ALIGNMENT = 0;
    static const unsigned int ALIGNMENT_OUTPUT_CLUSTER = 1;

    static const unsigned int EXPAND_TRANSFER_EVALUE = 0;
    static const unsigned int EXPAND_RESCORE_BACKTRACE = 1;

    static const unsigned int PCMODE_SUBSTITUTION_SCORE = 0;
    static const unsigned int PCMODE_CONTEXT_SPECIFIC = 1;


    static const unsigned int WRITER_ASCII_MODE = 0;
    static const unsigned int WRITER_COMPRESSED_MODE = 1;
    static const unsigned int WRITER_LEXICOGRAPHIC_MODE = 2;

    // convertalis alignment
    static const int FORMAT_ALIGNMENT_BLAST_TAB = 0;
    static const int FORMAT_ALIGNMENT_SAM = 1;
    static const int FORMAT_ALIGNMENT_BLAST_WITH_LEN = 2;
    static const int FORMAT_ALIGNMENT_HTML = 3;
    static const int FORMAT_ALIGNMENT_BLAST_TAB_WITH_HEADERS = 4;

    // result2msa
    static const int FORMAT_MSA_CA3M = 0;
    static const int FORMAT_MSA_CA3M_CONSENSUS = 1;
    static const int FORMAT_MSA_FASTADB = 2;
    static const int FORMAT_MSA_FASTADB_SUMMARY = 3;
    static const int FORMAT_MSA_STOCKHOLM_FLAT = 4;
    static const int FORMAT_MSA_A3M = 5;
    static const int FORMAT_MSA_A3M_ALN_INFO = 6;
    // outfmt
    static const int OUTFMT_QUERY = 0;
    static const int OUTFMT_TARGET = 1;
    static const int OUTFMT_EVALUE = 2;
    static const int OUTFMT_GAPOPEN = 3;
    static const int OUTFMT_PIDENT = 4;
    static const int OUTFMT_NIDENT = 5;
    static const int OUTFMT_QSTART = 6;
    static const int OUTFMT_QEND = 7;
    static const int OUTFMT_QLEN = 8;
    static const int OUTFMT_TSTART = 9;
    static const int OUTFMT_TEND = 10;
    static const int OUTFMT_TLEN = 11;
    static const int OUTFMT_ALNLEN = 12;
    static const int OUTFMT_RAW = 13;
    static const int OUTFMT_BITS = 14;
    static const int OUTFMT_CIGAR = 15;
    static const int OUTFMT_QSEQ = 16;
    static const int OUTFMT_TSEQ = 17;
    static const int OUTFMT_QHEADER = 18;
    static const int OUTFMT_THEADER = 19;
    static const int OUTFMT_QALN = 20;
    static const int OUTFMT_TALN = 21;
    static const int OUTFMT_QFRAME = 22;
    static const int OUTFMT_TFRAME = 23;
    static const int OUTFMT_MISMATCH = 24;
    static const int OUTFMT_QCOV = 25;
    static const int OUTFMT_TCOV = 26;
    static const int OUTFMT_EMPTY = 27;
    static const int OUTFMT_QSET = 28;
    static const int OUTFMT_QSETID = 29;
    static const int OUTFMT_TSET = 30;
    static const int OUTFMT_TSETID = 31;
    static const int OUTFMT_TAXID = 32;
    static const int OUTFMT_TAXNAME = 33;
    static const int OUTFMT_TAXLIN = 34;
    static const int OUTFMT_QORFSTART = 35;
    static const int OUTFMT_QORFEND = 36;
    static const int OUTFMT_TORFSTART = 37;
    static const int OUTFMT_TORFEND = 38;
    static const int OUTFMT_FIDENT = 39;

    static const int INDEX_SUBSET_NORMAL = 0;
    static const int INDEX_SUBSET_NO_HEADERS = 1;
    static const int INDEX_SUBSET_NO_PREFILTER = 2;

    static std::vector<int> getOutputFormat(int formatMode, const std::string &outformat, bool &needSequences, bool &needBacktrace, bool &needFullHeaders,
                                            bool &needLookup, bool &needSource, bool &needTaxonomyMapping, bool &needTaxonomy);

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

    // taxonomy output
    static const int TAXONOMY_OUTPUT_LCA = 0;
    static const int TAXONOMY_OUTPUT_ALIGNMENT = 1;
    static const int TAXONOMY_OUTPUT_BOTH = 2;

    // aggregate taxonomy
    static const int AGG_TAX_UNIFORM = 0;
    static const int AGG_TAX_MINUS_LOG_EVAL = 1;
    static const int AGG_TAX_SCORE = 2;

    // pairaln dummy mode
    static const int PAIRALN_DUMMY_MODE_OFF = 0;
    static const int PAIRALN_DUMMY_MODE_ON = 1;

    // pairaln mode
    static const int PAIRALN_MODE_ALL_PER_SPECIES = 0;
    static const int PAIRALN_MODE_COVER_ALL_CHAINS = 1;

    // taxonomy search strategy
    static const int TAXONOMY_SINGLE_SEARCH = 1;
    static const int TAXONOMY_2BLCA = 2;
    static const int TAXONOMY_APPROX_2BLCA = 3;
    static const int TAXONOMY_TOP_HIT = 4;

    static const int PARSE_VARIADIC = 1;
    static const int PARSE_REST = 2;
    static const int PARSE_ALLOW_EMPTY = 4;

    // preload mode
    static const int PRELOAD_MODE_AUTO = 0;
    static const int PRELOAD_MODE_FREAD = 1;
    static const int PRELOAD_MODE_MMAP = 2;
    static const int PRELOAD_MODE_MMAP_TOUCH = 3;

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

    static const int MAX_SEQ_LEN = 65535;

    // extractalignedregion
    static const int EXTRACT_QUERY  = 1;
    static const int EXTRACT_TARGET = 2;

    static const int CLUST_HASH_DEFAULT_ALPH_SIZE = 3;
    static const int CLUST_HASH_DEFAULT_MIN_SEQ_ID = 99;
    static const int CLUST_LINEAR_DEFAULT_ALPH_SIZE = 13;
    static const int CLUST_LINEAR_DEFAULT_K = 0;
    static const int CLUST_LINEAR_KMER_PER_SEQ = 0;

    // cov mode
    static const int COV_MODE_BIDIRECTIONAL  = 0;
    static const int COV_MODE_TARGET = 1;
    static const int COV_MODE_QUERY = 2;
    static const int COV_MODE_LENGTH_QUERY = 3;
    static const int COV_MODE_LENGTH_TARGET = 4;
    static const int COV_MODE_LENGTH_SHORTER = 5;

    // seq. id mode
    static const int SEQ_ID_ALN_LEN  = 0;
    static const int SEQ_ID_SHORT = 1;
    static const int SEQ_ID_LONG = 2;

    // seq. split mode
    static const int SEQUENCE_SPLIT_MODE_HARD = 0;
    static const int SEQUENCE_SPLIT_MODE_SOFT = 1;

    // rescorediagonal
    static const int RESCORE_MODE_HAMMING = 0;
    static const int RESCORE_MODE_SUBSTITUTION = 1;
    static const int RESCORE_MODE_ALIGNMENT = 2;
    static const int RESCORE_MODE_END_TO_END_ALIGNMENT = 3;
    static const int RESCORE_MODE_WINDOW_QUALITY_ALIGNMENT = 4;

    // combinepvalperset
    static const int AGGREGATION_MODE_MULTIHIT = 0;
    static const int AGGREGATION_MODE_MIN_PVAL = 1;
    static const int AGGREGATION_MODE_PRODUCT = 2;
    static const int AGGREGATION_MODE_TRUNCATED_PRODUCT = 3;

    // header type
    static const int HEADER_TYPE_UNICLUST = 1;
    static const int HEADER_TYPE_METACLUST = 2;

    // createsubdb, filtertaxseqdb type
    static const int SUBDB_MODE_HARD = 0;
    static const int SUBDB_MODE_SOFT = 1;

    static const int ID_MODE_KEYS = 0;
    static const int ID_MODE_LOOKUP = 1;

    // unpackdb
    static const int UNPACK_NAME_KEY = 0;
    static const int UNPACK_NAME_ACCESSION = 1;

    // result direction
    static const int PARAM_RESULT_DIRECTION_QUERY  = 0;
    static const int PARAM_RESULT_DIRECTION_TARGET = 1;

    // path to databases
    std::string db1;
    std::string db1Index;
    std::string db1dbtype;


    std::string hdr1;
    std::string hdr1Index;
    std::string hdr1dbtype;


    std::string db2;
    std::string db2Index;
    std::string db2dbtype;


    std::string hdr2;
    std::string hdr2Index;
    std::string hdr2dbtype;


    std::string db3;
    std::string db3Index;
    std::string db3dbtype;


    std::string hdr3;
    std::string hdr3Index;
    std::string hdr3dbtype;

    std::string db4;
    std::string db4Index;
    std::string db4dbtype;

    std::string hdr4;
    std::string hdr4Index;
    std::string hdr4dbtype;

    std::string db5;
    std::string db5Index;
    std::string db5dbtype;

    std::string hdr5;
    std::string hdr5Index;
    std::string hdr5dbtype;

    std::string db6;
    std::string db6Index;
    std::string db6dbtype;

    std::string hdr6;
    std::string hdr6Index;
    std::string hdr6dbtype;


    std::vector<std::string> filenames;

    const char** restArgv;
    int restArgc;

    MultiParam<NuclAA<std::string>> scoringMatrixFile;       // path to scoring matrix
    MultiParam<NuclAA<std::string>> seedScoringMatrixFile;   // seed sub. matrix
    size_t maxSeqLen;                    // sequence length
    size_t maxResListLen;                // Maximal result list length per query
    int    verbosity;                    // log level
    int    threads;                      // Amounts of threads
    int    compressed;                   // compressed writer
    bool   removeTmpFiles;               // Do not delete temp files
    bool   includeIdentity;              // include identical ids as hit

    // PREFILTER
    float  sensitivity;                  // target sens
    int    kmerSize;                     // kmer size for the prefilter
    int targetSearchMode;                // target search mode
    MultiParam<SeqProf<int>> kmerScore;   // kmer score for the prefilter
    MultiParam<NuclAA<int>> alphabetSize; // alphabet size for the prefilter
    int    compBiasCorrection;           // Aminoacid composiont correction
    float    compBiasCorrectionScale;    // Aminoacid composiont correction scale factor

    bool   diagonalScoring;              // switch diagonal scoring
    int    exactKmerMatching;            // only exact k-mer matching
    int    maskMode;                     // mask low complex areas
    float  maskProb;                     // mask probability
    int    maskLowerCaseMode;            // mask lowercase letters in prefilter and kmermatchers

    int    minDiagScoreThr;              // min diagonal score
    int    spacedKmer;                   // Spaced Kmers
    int    split;                        // Split database in n equal chunks
    int    splitMode;                    // Split by query or target DB
    size_t splitMemoryLimit;             // Maximum memory in bytes a split can use
    size_t diskSpaceLimit;               // Maximum disk space in bytes for sliced reverse profile search
    bool   splitAA;                      // Split database by amino acid count instead
    int    preloadMode;                  // Preload mode of database
    float  scoreBias;                    // Add this bias to the score when computing the alignements
    float  realignScoreBias;             // Add this bias additionally when realigning
    int    realignMaxSeqs;               // Max alignments to realign
    std::string spacedKmerPattern;       // User-specified kmer pattern
    std::string localTmp;                // Local temporary path


    // ALIGNMENT
    int alignmentMode;                   // alignment mode 0=fastest on parameters,
                                         // 1=score only, 2=score, cov, start/end pos, 3=score, cov, start/end pos, seq.id,
    int alignmentOutputMode;             // alignment output mode 0=alignment, 1=cluster
    double evalThr;                      // e-value threshold for acceptance
    float  covThr;                       // coverage query&target threshold for acceptance
    int    covMode;                      // coverage target threshold for acceptance
    int    seqIdMode;                    // seq. id. normalize mode

    int    maxRejected;                  // after n sequences that are above eval stop
    int    maxAccept;                    // after n accepted sequences stop
    int    altAlignment;                 // show up to this many alternative alignments
    float  seqIdThr;                     // sequence identity threshold for acceptance
    int    alnLenThr;                    // min. alignment length
    bool   addBacktrace;                 // store backtrace string (M=Match, D=deletion, I=insertion)
    bool   realign;                      // realign hit with more conservative score
    MultiParam<NuclAA<int>> gapOpen;             // gap open cost
    MultiParam<NuclAA<int>> gapExtend;           // gap extension cost
    float correlationScoreWeight; // correlation score weight
#ifdef GAP_POS_SCORING
    int    gapPseudoCount;               // for calculation of position-specific gap opening penalties
#endif
    int    zdrop;                        // zdrop

    // workflow
    std::string runner;
    bool reuseLatest;

    // CLUSTERING
    int    clusteringMode;
    int    clusterSteps;
    bool   singleStepClustering;
    int    clusterReassignment;

    // SEARCH WORKFLOW
    int numIterations;
    float startSens;
    int sensSteps;
    bool exhaustiveSearch;
    int exhaustiveFilterMsa;
    int strand;
    int orfFilter;
    float orfFilterSens;
    double orfFilterEval;
    bool lcaSearch;

    // easysearch
    bool greedyBestHits;

    //CLUSTERING
    int maxIteration;                   // Maximum depth of breadth first search in connected component
    int similarityScoreType;            // Type of score to use for reassignment 1=alignment score. 2=coverage 3=sequence identity 4=E-value 5= Score per Column

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
    int translate;
    int createLookup;

    // convertalis
    int formatAlignmentMode;
    std::string outfmt;
    bool dbOut;

    // rescorediagonal
    int rescoreMode;
    bool wrappedScoring;
    bool filterHits;
    bool globalAlignment;
    int sortResults;

    // result2msa
    int msaFormatMode;
    bool allowDeletion;
    std::string summaryPrefix;
    bool skipQuery;

    // convertmsa
    int identifierField;

    // msa2profile
    int matchMode;
    float matchRatio;

    // result2profile
    int maskProfile;
    float filterMaxSeqId;
    double evalProfile;
    int filterMsa;
    float qsc;
    std::string qid;
    float covMSAThr;
    int Ndiff;
    int filterMinEnable;
    bool wg;
    int pcmode;
    MultiParam<PseudoCounts> pca;
    MultiParam<PseudoCounts> pcb;

    // sequence2profile
    float neff;
    float tau;

    // createtsv
    bool firstSeqRepr;
    int idxSeqSrc;
    bool fullHeader;
    size_t targetTsvColumn;

    //result2stats
    std::string stat;

    // linearcluster
    int kmersPerSequence;
    MultiParam<NuclAA<float>> kmersPerSequenceScale;
    bool includeOnlyExtendable;
    bool ignoreMultiKmer;
    int hashShift;
    int pickNbest;
    int adjustKmerLength;
    int resultDirection;
    float weightThr;
    std::string weightFile;

    // indexdb
    int checkCompatible;
    int searchType;
    int indexSubset;

    // createdb
    int identifierOffset;
    int dbType;
    int createdbMode;
    bool shuffleDatabase;

    // splitsequence
    int sequenceOverlap;
    int sequenceSplitMode;
    int headerSplitMode;

    // convert2fasta
    bool useHeaderFile;
    int writeLookup;

    // result2flat
    bool useHeader;

    // createclusearchdb
    std::string dbSuffixList;

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
    int columnToTake;
    std::string filterColumnRegex;
    std::string filteringFile;
    std::string mappingFile;
    std::string filterExpression;
    bool positiveFilter;
    bool trimToOneColumn;
    int extractLines;
    double compValue;
    std::string compOperator;
    int sortEntries;
    bool beatsFirst;
    std::string joinDB;

    // besthitperset
    bool simpleBestHit;
    float alpha;
    bool shortOutput;
    int aggregationMode;

    // mergedbs
    std::string mergePrefixes;
    bool mergeStopEmpty;

    // summarizetabs
    float overlap;
    int msaType;

    // extractalignedregion
    int extractMode;

    // convertkb
    std::string kbColumns;

    // concatdbs
    bool preserveKeysB;
    bool takeLargerEntry;

    // offsetalignments
    int chainAlignment;
    int mergeQuery;

    // tsv2db
    int outputDbType;

    // diff
    bool useSequenceId;

    // prefixid
    std::string prefix;
    bool tsvOut;

    // clusterUpdate;
    bool recoverDeleted;

    // summarize headers
    int headerType;

    // filtertaxdb, filtertaxseqdb
    std::string taxonList;

    // view
    std::string idList;
    int idxEntryType;

    // lca
    int pickIdFrom;
    std::string lcaRanks;
    int showTaxLineage;
    std::string blacklist;

    // aggregatetax
    float majorityThr;
    int voteMode;

    // pairaln
    int pairdummymode;
    int pairmode;

    // taxonomyreport
    int reportMode;

    // createtaxdb
    std::string ncbiTaxDump;
    std::string taxMappingFile;
    int taxMappingMode;
    int taxDbMode;

    // exapandaln
    int expansionMode;
    int expandFilterClusters;

    // taxonomy
    int taxonomySearchMode;
    int taxonomyOutputMode;

    // createsubdb
    int subDbMode;
    int dbIdMode;

    // tar2db
    std::string tarInclude;
    std::string tarExclude;

    // unpackdb
    std::string unpackSuffix;
    int unpackNameMode;

    // for modules that should handle -h themselves
    bool help;

    // tool citations
    std::map<unsigned int, const char*> citations;

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
    void initMatrices();
    void parseParameters(int argc, const char *pargv[], const Command &command, bool printPar, int parseFlags,
                         int outputFlags);
    void printUsageMessage(const Command& command, unsigned int outputFlag, const char* extraText = NULL);
    void printParameters(const std::string &module, int argc, const char* pargv[],
                         const std::vector<MMseqsParameter*> &par);

    void checkIfDatabaseIsValid(const Command& command, int argc, const char** argv, bool isStartVar, bool isMiddleVar, bool isEndVar);

    std::vector<MMseqsParameter*> removeParameter(const std::vector<MMseqsParameter*>& par, const MMseqsParameter& x);

    PARAMETER(PARAM_S)
    PARAMETER(PARAM_K)
    PARAMETER(PARAM_TARGET_SEARCH_MODE)
    PARAMETER(PARAM_THREADS)
    PARAMETER(PARAM_COMPRESSED)
    PARAMETER(PARAM_ALPH_SIZE)
    PARAMETER(PARAM_MAX_SEQ_LEN)
    PARAMETER(PARAM_DIAGONAL_SCORING)
    PARAMETER(PARAM_EXACT_KMER_MATCHING)
    PARAMETER(PARAM_MASK_RESIDUES)
    PARAMETER(PARAM_MASK_PROBABILTY)
    PARAMETER(PARAM_MASK_LOWER_CASE)

    PARAMETER(PARAM_MIN_DIAG_SCORE)
    PARAMETER(PARAM_K_SCORE)
    PARAMETER(PARAM_MAX_SEQS)
    PARAMETER(PARAM_SPLIT)
    PARAMETER(PARAM_SPLIT_MODE)
    PARAMETER(PARAM_SPLIT_MEMORY_LIMIT)
    PARAMETER(PARAM_DISK_SPACE_LIMIT)
    PARAMETER(PARAM_SPLIT_AMINOACID)
    PARAMETER(PARAM_SUB_MAT)
    PARAMETER(PARAM_SEED_SUB_MAT)
    PARAMETER(PARAM_NO_COMP_BIAS_CORR)
    PARAMETER(PARAM_NO_COMP_BIAS_CORR_SCALE)
    PARAMETER(PARAM_SPACED_KMER_MODE)
    PARAMETER(PARAM_REMOVE_TMP_FILES)
    PARAMETER(PARAM_INCLUDE_IDENTITY)
    PARAMETER(PARAM_PRELOAD_MODE)
    PARAMETER(PARAM_SPACED_KMER_PATTERN)
    PARAMETER(PARAM_LOCAL_TMP)
    std::vector<MMseqsParameter*> prefilter;
    std::vector<MMseqsParameter*> ungappedprefilter;

    // alignment
    PARAMETER(PARAM_ALIGNMENT_MODE)
    PARAMETER(PARAM_ALIGNMENT_OUTPUT_MODE)
    PARAMETER(PARAM_E)
    PARAMETER(PARAM_C)
    PARAMETER(PARAM_COV_MODE)
    PARAMETER(PARAM_SEQ_ID_MODE)
    PARAMETER(PARAM_MAX_REJECTED)
    PARAMETER(PARAM_MAX_ACCEPT)
    PARAMETER(PARAM_ADD_BACKTRACE)
    PARAMETER(PARAM_REALIGN)
    PARAMETER(PARAM_MIN_SEQ_ID)
    PARAMETER(PARAM_MIN_ALN_LEN)
    PARAMETER(PARAM_SCORE_BIAS)
    PARAMETER(PARAM_REALIGN_SCORE_BIAS)
    PARAMETER(PARAM_REALIGN_MAX_SEQS)
    PARAMETER(PARAM_CORR_SCORE_WEIGHT)
    PARAMETER(PARAM_ALT_ALIGNMENT)
    PARAMETER(PARAM_GAP_OPEN)
    PARAMETER(PARAM_GAP_EXTEND)
#ifdef GAP_POS_SCORING
    PARAMETER(PARAM_GAP_PSEUDOCOUNT)
#endif
    PARAMETER(PARAM_ZDROP)

    // clustering
    PARAMETER(PARAM_CLUSTER_MODE)
    PARAMETER(PARAM_CLUSTER_STEPS)
    PARAMETER(PARAM_CASCADED)
    PARAMETER(PARAM_CLUSTER_REASSIGN)

    // affinity clustering
    PARAMETER(PARAM_MAXITERATIONS)
    PARAMETER(PARAM_SIMILARITYSCORE)

    // logging
    PARAMETER(PARAM_V)
    std::vector<MMseqsParameter*> clust;

    // format alignment
    PARAMETER(PARAM_FORMAT_MODE)
    PARAMETER(PARAM_FORMAT_OUTPUT)
    PARAMETER(PARAM_DB_OUTPUT)

    // rescoremode
    PARAMETER(PARAM_RESCORE_MODE)
    PARAMETER(PARAM_WRAPPED_SCORING)
    PARAMETER(PARAM_FILTER_HITS)
    PARAMETER(PARAM_SORT_RESULTS)

    // result2msa
    PARAMETER(PARAM_MSA_FORMAT_MODE)
    PARAMETER(PARAM_ALLOW_DELETION)
    PARAMETER(PARAM_SUMMARY_PREFIX)
    PARAMETER(PARAM_SKIP_QUERY)

    // convertmsa
    PARAMETER(PARAM_IDENTIFIER_FIELD)

    // msa2profile
    PARAMETER(PARAM_MATCH_MODE)
    PARAMETER(PARAM_MATCH_RATIO)

    // result2profile
    PARAMETER(PARAM_MASK_PROFILE)
    PARAMETER(PARAM_E_PROFILE)
    PARAMETER(PARAM_FILTER_MSA)
    PARAMETER(PARAM_FILTER_MAX_SEQ_ID)
    PARAMETER(PARAM_FILTER_QSC)
    PARAMETER(PARAM_FILTER_QID)
    PARAMETER(PARAM_FILTER_COV)
    PARAMETER(PARAM_FILTER_NDIFF)
    PARAMETER(PARAM_FILTER_MIN_ENABLE)
    PARAMETER(PARAM_WG)
    PARAMETER(PARAM_PC_MODE)
    PARAMETER(PARAM_PCA)
    PARAMETER(PARAM_PCB)

    // sequence2profile
    PARAMETER(PARAM_NEFF)
    PARAMETER(PARAM_TAU)

    // createtsv
    PARAMETER(PARAM_TARGET_COLUMN)
    PARAMETER(PARAM_FIRST_SEQ_REP_SEQ)
    PARAMETER(PARAM_FULL_HEADER)
    PARAMETER(PARAM_IDX_SEQ_SRC)

    // result2stat
    PARAMETER(PARAM_STAT)

    // linearcluster
    PARAMETER(PARAM_KMER_PER_SEQ)
    PARAMETER(PARAM_KMER_PER_SEQ_SCALE)
    PARAMETER(PARAM_INCLUDE_ONLY_EXTENDABLE)
    PARAMETER(PARAM_IGNORE_MULTI_KMER)
    PARAMETER(PARAM_HASH_SHIFT)
    PARAMETER(PARAM_PICK_N_SIMILAR)
    PARAMETER(PARAM_ADJUST_KMER_LEN)
    PARAMETER(PARAM_RESULT_DIRECTION)
    PARAMETER(PARAM_WEIGHT_FILE)
    PARAMETER(PARAM_WEIGHT_THR)

    // workflow
    PARAMETER(PARAM_RUNNER)
    PARAMETER(PARAM_REUSELATEST)

    // search workflow
    PARAMETER(PARAM_NUM_ITERATIONS)
    PARAMETER(PARAM_START_SENS)
    PARAMETER(PARAM_SENS_STEPS)
    PARAMETER(PARAM_EXHAUSTIVE_SEARCH)
    PARAMETER(PARAM_EXHAUSTIVE_SEARCH_FILTER)
    PARAMETER(PARAM_STRAND)
    PARAMETER(PARAM_ORF_FILTER)
    PARAMETER(PARAM_ORF_FILTER_S)
    PARAMETER(PARAM_ORF_FILTER_E)
    PARAMETER(PARAM_LCA_SEARCH)

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
    PARAMETER(PARAM_TRANSLATE)
    PARAMETER(PARAM_CREATE_LOOKUP)

    // indexdb
    PARAMETER(PARAM_CHECK_COMPATIBLE)
    PARAMETER(PARAM_SEARCH_TYPE)
    PARAMETER(PARAM_INDEX_SUBSET)

    // createdb
    PARAMETER(PARAM_USE_HEADER) // also used by extractorfs
    PARAMETER(PARAM_ID_OFFSET)  // same
    PARAMETER(PARAM_DB_TYPE)
    PARAMETER(PARAM_CREATEDB_MODE)
    PARAMETER(PARAM_SHUFFLE)
    PARAMETER(PARAM_WRITE_LOOKUP)

    // convert2fasta
    PARAMETER(PARAM_USE_HEADER_FILE)

    // split sequence
    PARAMETER(PARAM_SEQUENCE_OVERLAP)
    PARAMETER(PARAM_SEQUENCE_SPLIT_MODE)
    PARAMETER(PARAM_HEADER_SPLIT_MODE)

    // createclusearchdb
    PARAMETER(PARAM_DB_SUFFIX_LIST)

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
    PARAMETER(PARAM_COLUMN_TO_TAKE)
    PARAMETER(PARAM_FILTER_REGEX)
    PARAMETER(PARAM_FILTER_POS)
    PARAMETER(PARAM_FILTER_FILE)
    PARAMETER(PARAM_FILTER_EXPRESSION)
    PARAMETER(PARAM_MAPPING_FILE)
    PARAMETER(PARAM_TRIM_TO_ONE_COL)
    PARAMETER(PARAM_EXTRACT_LINES)
    PARAMETER(PARAM_COMP_OPERATOR)
    PARAMETER(PARAM_COMP_VALUE)
    PARAMETER(PARAM_SORT_ENTRIES)
    PARAMETER(PARAM_BEATS_FIRST)
    PARAMETER(PARAM_JOIN_DB)

    //besthitperset
    PARAMETER(PARAM_SIMPLE_BEST_HIT)
    PARAMETER(PARAM_ALPHA)
    PARAMETER(PARAM_SHORT_OUTPUT)
    PARAMETER(PARAM_AGGREGATION_MODE)

    // concatdb
    PARAMETER(PARAM_PRESERVEKEYS)
    PARAMETER(PARAM_TAKE_LARGER_ENTRY)

    // offsetalignment
    PARAMETER(PARAM_CHAIN_ALIGNMENT)
    PARAMETER(PARAM_MERGE_QUERY)


    // tsv2db
    PARAMETER(PARAM_OUTPUT_DBTYPE)

    // diff
    PARAMETER(PARAM_USESEQID)

    // prefixid
    PARAMETER(PARAM_PREFIX)
    PARAMETER(PARAM_TSV)

    // summarize headers
    PARAMETER(PARAM_HEADER_TYPE)

    // mergedbs
    PARAMETER(PARAM_MERGE_PREFIXES)
    PARAMETER(PARAM_MERGE_STOP_EMPTY)

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

    // filtertaxdb, filtertaxseqdb
    PARAMETER(PARAM_TAXON_LIST)

    // view
    PARAMETER(PARAM_ID_LIST)
    PARAMETER(PARAM_IDX_ENTRY_TYPE)

    // lca, addtaxonomy and aggregatetax
    PARAMETER(PARAM_PICK_ID_FROM)
    PARAMETER(PARAM_LCA_RANKS)
    PARAMETER(PARAM_BLACKLIST)
    PARAMETER(PARAM_TAXON_ADD_LINEAGE)

    // aggregatetax
    PARAMETER(PARAM_MAJORITY)
    PARAMETER(PARAM_VOTE_MODE)

    // pairaln
    PARAMETER(PARAM_PAIRING_DUMMY_MODE)
    PARAMETER(PARAM_PAIRING_MODE)
    
    // taxonomyreport
    PARAMETER(PARAM_REPORT_MODE)

    // createtaxdb
    PARAMETER(PARAM_NCBI_TAX_DUMP)
    PARAMETER(PARAM_TAX_MAPPING_FILE)
    PARAMETER(PARAM_TAX_MAPPING_MODE)
    PARAMETER(PARAM_TAX_DB_MODE)

    // exapandaln
    PARAMETER(PARAM_EXPANSION_MODE)
    PARAMETER(PARAM_EXPAND_FILTER_CLUSTERS)

    // taxonomy
    PARAMETER(PARAM_LCA_MODE)
    PARAMETER(PARAM_TAX_OUTPUT_MODE)

    // createsubdb
    PARAMETER(PARAM_SUBDB_MODE)
    PARAMETER(PARAM_ID_MODE)

    // tar2db
    PARAMETER(PARAM_TAR_INCLUDE)
    PARAMETER(PARAM_TAR_EXCLUDE)

    // unpackdb
    PARAMETER(PARAM_UNPACK_SUFFIX)
    PARAMETER(PARAM_UNPACK_NAME_MODE)

    // for modules that should handle -h themselves
    PARAMETER(PARAM_HELP)
    PARAMETER(PARAM_HELP_LONG)

    struct PredefinedSubstitutionMatrix{
        std::string name;
        const unsigned char * subMatData;
        unsigned int subMatDataLen;
        PredefinedSubstitutionMatrix(const char * name, const unsigned char * subMatData, const unsigned int subMatDataLen)
                : name(name), subMatData(subMatData), subMatDataLen(subMatDataLen) {}

    };
    std::vector<PredefinedSubstitutionMatrix> substitutionMatrices;

    std::vector<MMseqsParameter*> empty;
    std::vector<MMseqsParameter*> onlyverbosity;
    std::vector<MMseqsParameter*> view;
    std::vector<MMseqsParameter*> verbandcompression;
    std::vector<MMseqsParameter*> onlythreads;
    std::vector<MMseqsParameter*> threadsandcompression;

    std::vector<MMseqsParameter*> alignall;
    std::vector<MMseqsParameter*> align;
    std::vector<MMseqsParameter*> rescorediagonal;
    std::vector<MMseqsParameter*> alignbykmer;
    std::vector<MMseqsParameter*> createFasta;
    std::vector<MMseqsParameter*> convertprofiledb;
    std::vector<MMseqsParameter*> sequence2profile;
    std::vector<MMseqsParameter*> result2profile;
    std::vector<MMseqsParameter*> result2msa;
    std::vector<MMseqsParameter*> result2dnamsa;
    std::vector<MMseqsParameter*> filtera3m;
    std::vector<MMseqsParameter*> filterresult;
    std::vector<MMseqsParameter*> convertmsa;
    std::vector<MMseqsParameter*> msa2profile;
    std::vector<MMseqsParameter*> createtsv;
    std::vector<MMseqsParameter*> result2stats;
    std::vector<MMseqsParameter*> extractorfs;
    std::vector<MMseqsParameter*> extractframes;
    std::vector<MMseqsParameter*> orftocontig;
    std::vector<MMseqsParameter*> reverseseq;
    std::vector<MMseqsParameter*> splitdb;
    std::vector<MMseqsParameter*> splitsequence;
    std::vector<MMseqsParameter*> masksequence;
    std::vector<MMseqsParameter*> indexdb;
    std::vector<MMseqsParameter*> kmerindexdb;
    std::vector<MMseqsParameter*> createindex;
    std::vector<MMseqsParameter*> createlinindex;
    std::vector<MMseqsParameter*> convertalignments;
    std::vector<MMseqsParameter*> createdb;
    std::vector<MMseqsParameter*> convert2fasta;
    std::vector<MMseqsParameter*> result2flat;
    std::vector<MMseqsParameter*> result2repseq;
    std::vector<MMseqsParameter*> gff2db;
    std::vector<MMseqsParameter*> clusthash;
    std::vector<MMseqsParameter*> kmermatcher;
    std::vector<MMseqsParameter*> kmersearch;
    std::vector<MMseqsParameter*> countkmer;
    std::vector<MMseqsParameter*> easylinclustworkflow;
    std::vector<MMseqsParameter*> linclustworkflow;
    std::vector<MMseqsParameter*> easysearchworkflow;
    std::vector<MMseqsParameter*> searchworkflow;
    std::vector<MMseqsParameter*> linsearchworkflow;
    std::vector<MMseqsParameter*> easylinsearchworkflow;
    std::vector<MMseqsParameter*> mapworkflow;
    std::vector<MMseqsParameter*> easyclusterworkflow;
    std::vector<MMseqsParameter*> clusterworkflow;
    std::vector<MMseqsParameter*> clusterUpdateSearch;
    std::vector<MMseqsParameter*> clusterUpdateClust;
    std::vector<MMseqsParameter*> mergeclusters;
    std::vector<MMseqsParameter*> clusterUpdate;
    std::vector<MMseqsParameter*> translatenucs;
    std::vector<MMseqsParameter*> swapresult;
    std::vector<MMseqsParameter*> swapdb;
    std::vector<MMseqsParameter*> createseqfiledb;
    std::vector<MMseqsParameter*> filterDb;
    std::vector<MMseqsParameter*> offsetalignment;
    std::vector<MMseqsParameter*> proteinaln2nucl;
    std::vector<MMseqsParameter*> subtractdbs;
    std::vector<MMseqsParameter*> diff;
    std::vector<MMseqsParameter*> concatdbs;
    std::vector<MMseqsParameter*> mergedbs;
    std::vector<MMseqsParameter*> summarizeheaders;
    std::vector<MMseqsParameter*> prefixid;
    std::vector<MMseqsParameter*> summarizeresult;
    std::vector<MMseqsParameter*> summarizetabs;
    std::vector<MMseqsParameter*> extractdomains;
    std::vector<MMseqsParameter*> createclusearchdb;
    std::vector<MMseqsParameter*> extractalignedregion;
    std::vector<MMseqsParameter*> convertkb;
    std::vector<MMseqsParameter*> tsv2db;
    std::vector<MMseqsParameter*> lca;
    std::vector<MMseqsParameter*> majoritylca;
    std::vector<MMseqsParameter*> addtaxonomy;
    std::vector<MMseqsParameter*> taxonomyreport;
    std::vector<MMseqsParameter*> filtertaxdb;
    std::vector<MMseqsParameter*> filtertaxseqdb;
    std::vector<MMseqsParameter*> aggregatetax;
    std::vector<MMseqsParameter*> aggregatetaxweights;
    std::vector<MMseqsParameter*> taxonomy;
    std::vector<MMseqsParameter*> easytaxonomy;
    std::vector<MMseqsParameter*> createsubdb;
    std::vector<MMseqsParameter*> renamedbkeys;
    std::vector<MMseqsParameter*> createtaxdb;
    std::vector<MMseqsParameter*> profile2pssm;
    std::vector<MMseqsParameter*> profile2neff;
    std::vector<MMseqsParameter*> profile2seq;
    std::vector<MMseqsParameter*> besthitbyset;
    std::vector<MMseqsParameter*> combinepvalbyset;
    std::vector<MMseqsParameter*> multihitdb;
    std::vector<MMseqsParameter*> multihitsearch;
    std::vector<MMseqsParameter*> expandaln;
    std::vector<MMseqsParameter*> expand2profile;
    std::vector<MMseqsParameter*> pairaln;
    std::vector<MMseqsParameter*> sortresult;
    std::vector<MMseqsParameter*> enrichworkflow;
    std::vector<MMseqsParameter*> databases;
    std::vector<MMseqsParameter*> tar2db;
    std::vector<MMseqsParameter*> unpackdbs;
    std::vector<MMseqsParameter*> appenddbtoindex;

    std::vector<MMseqsParameter*> combineList(const std::vector<MMseqsParameter*> &par1,
                                             const std::vector<MMseqsParameter*> &par2);

    size_t hashParameter(const std::vector<DbType> &dbtypes, const std::vector<std::string> &filenames, const std::vector<MMseqsParameter*> &par);

    std::string createParameterString(const std::vector<MMseqsParameter*> &vector, bool wasSet = false);

    void overrideParameterDescription(MMseqsParameter& par, const char *description, const char *regex = NULL, int category = 0);

    static std::vector<std::string> findMissingTaxDbFiles(const std::string &filename);
    static void printTaxDbError(const std::string &filename, const std::vector<std::string>& missingFiles);

    static const uint32_t DBTYPE_MASK = 0x0000FFFF;

    static bool isEqualDbtype(const int type1, const int type2) {
        return ((type1 & DBTYPE_MASK) == (type2 & DBTYPE_MASK));
    }

    static const char* getDbTypeName(int dbtype) {
        switch (dbtype & DBTYPE_MASK) {
            case DBTYPE_AMINO_ACIDS: return "Aminoacid";
            case DBTYPE_NUCLEOTIDES: return "Nucleotide";
            case DBTYPE_HMM_PROFILE: return "Profile";
            case DBTYPE_ALIGNMENT_RES: return "Alignment";
            case DBTYPE_CLUSTER_RES: return "Clustering";
            case DBTYPE_PREFILTER_RES: return "Prefilter";
            case DBTYPE_TAXONOMICAL_RESULT: return "Taxonomy";
            case DBTYPE_INDEX_DB: return "Index";
            case DBTYPE_CA3M_DB: return "CA3M";
            case DBTYPE_MSA_DB: return "MSA";
            case DBTYPE_GENERIC_DB: return "Generic";
            case DBTYPE_PREFILTER_REV_RES: return "Bi-directional prefilter";
            case DBTYPE_OFFSETDB: return "Offsetted headers";
            case DBTYPE_DIRECTORY: return "Directory";
            case DBTYPE_FLATFILE: return "Flatfile";
            case DBTYPE_STDIN: return "stdin";
            case DBTYPE_URI: return "uri";

            default: return "Unknown";
        }
    }

protected:
    Parameters();
    static Parameters* instance;

private:
    Parameters(Parameters const&);
    void operator=(Parameters const&);

};

#endif
