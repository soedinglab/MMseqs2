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

#define USE_ONLY_SET_PARAMETERS true // when createParameterString, generates
// flags only for the param set by the user

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

    static const int DBTYPE_AMINO_ACIDS = 0;
    static const int DBTYPE_NUCLEOTIDES = 1;
    static const int DBTYPE_HMM_PROFILE = 2;
    static const int DBTYPE_PROFILE_STATE_SEQ = 3;
    static const int DBTYPE_PROFILE_STATE_PROFILE = 4;
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
    // don't forget to add new database types to DBReader::getDbTypeName and Parameters::PARAM_OUTPUT_DBTYPE

    static const int SEARCH_TYPE_AUTO = 0;
    static const int SEARCH_TYPE_PROTEIN = 1;
    static const int SEARCH_TYPE_TRANSLATED = 2;
    static const int SEARCH_TYPE_NUCLEOTIDES = 3;


    static const unsigned int ALIGNMENT_MODE_FAST_AUTO = 0;
    static const unsigned int ALIGNMENT_MODE_SCORE_ONLY = 1;
    static const unsigned int ALIGNMENT_MODE_SCORE_COV = 2;
    static const unsigned int ALIGNMENT_MODE_SCORE_COV_SEQID = 3;
    static const unsigned int ALIGNMENT_MODE_UNGAPPED = 4;


    static const unsigned int WRITER_ASCII_MODE = 0;
    static const unsigned int WRITER_COMPRESSED_MODE = 1;
    static const unsigned int WRITER_LEXICOGRAPHIC_MODE = 2;

    // convertalis alignment
    static const int FORMAT_ALIGNMENT_BLAST_TAB = 0;
    static const int FORMAT_ALIGNMENT_SAM = 1;
    static const int FORMAT_ALIGNMENT_BLAST_WITH_LEN = 2;

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
    static std::vector<int> getOutputFormat(const std::string &outformat, bool &needSequences, bool &needBacktrace, bool &needFullHeaders);

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
    static const int TAXONOMY_2BLCA_APPROX = 3;


    static const int PARSE_VARIADIC = 1;
    static const int PARSE_REST = 2;

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
    static const int CLUST_LINEAR_DEFAULT_ALPH_SIZE = 13;
    static const int CLUST_LINEAR_DEFAULT_K = 0;
    static const int CLUST_LINEAR_KMER_PER_SEQ = 0;


    // cov mode
    static const int COV_MODE_BIDIRECTIONAL  = 0;
    static const int COV_MODE_TARGET = 1;
    static const int COV_MODE_QUERY = 2;
    static const int COV_MODE_LENGTH_QUERY = 3;
    static const int COV_MODE_LENGTH_TARGET = 4;

    // seq. id mode
    static const int SEQ_ID_ALN_LEN  = 0;
    static const int SEQ_ID_SHORT = 1;
    static const int SEQ_ID_LONG = 2;


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

    std::string scoringMatrixFile;       // path to scoring matrix
    std::string seedScoringMatrixFile;   // seed sub. matrix
    size_t maxSeqLen;                    // sequence length
    size_t maxResListLen;                // Maximal result list length per query
    int    verbosity;                    // log level
//    int    querySeqType;                 // Query sequence type (PROFILE, AMINOACIDE, NUCLEOTIDE)
//    int    targetSeqType;                // Target sequence type (PROFILE, AMINOACIDE, NUCLEOTIDE)
    int    threads;                      // Amounts of threads
    int    compressed;                   // compressed writer
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
    int    exactKmerMatching;            // only exact k-mer matching
    int    maskMode;                     // mask low complex areas
    int    maskLowerCaseMode;            // maske lowercase letters in prefilter and kmermatchers

    int    minDiagScoreThr;              // min diagonal score
    int    spacedKmer;                   // Spaced Kmers
    int    split;                        // Split database in n equal chunks
    int    splitMode;                    // Split by query or target DB
    int    splitMemoryLimit;             // Maximum amount of memory a split can use
    int    diskSpaceLimit;               // Disk space max usage for sliced reverse profile search
    bool   splitAA;                      // Split database by amino acid count instead
    size_t resListOffset;                // Offsets result list
    int    preloadMode;                  // Preload mode of database
    float  scoreBias;                    // Add this bias to the score when computing the alignements
    std::string spacedKmerPattern;       // User-specified kmer pattern
    std::string localTmp;                // Local temporary path

    // ALIGNMENT
    int alignmentMode;                   // alignment mode 0=fastest on parameters,
                                         // 1=score only, 2=score, cov, start/end pos, 3=score, cov, start/end pos, seq.id,
    float  evalThr;                      // e-value threshold for acceptance
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
    int    gapOpen;                      // gap open
    int    gapExtend;                    // gap extend

    // workflow
    std::string runner;
    bool reuseLatest;

    // CLUSTERING
    int    clusteringMode;
    int    clusterSteps;
    bool   cascaded;

    // SEARCH WORKFLOW
    int numIterations;
    float startSens;
    int sensSteps;
    bool sliceSearch;
    int strand;

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

    // convertprofiledb
    int profileMode;

    // convertalis
    int formatAlignmentMode;            // BLAST_TAB, PAIRWISE or SAM
    std::string outfmt;
    bool dbOut;

    // rescorediagonal
    int rescoreMode;
    bool filterHits;
    bool globalAlignment;
    int sortResults;

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
    int maskProfile;
    float filterMaxSeqId;
    float evalProfile;
    int filterMsa;
    float qsc;
    float qid;
    float covMSAThr;
    int Ndiff;
    bool wg;
    float pca;
    float pcb;

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
    bool includeOnlyExtendable;
    int skipNRepeatKmer;
    int hashShift;
    int pickNbest;
    int adjustKmerLength;

    // indexdb
    bool checkCompatible;
    int searchType;

    // createdb
    int identifierOffset;
    int dbType;
    bool splitSeqByLen;
    bool shuffleDatabase;

    // splitsequence
    int sequenceOverlap;

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
    int columnToTake;
    std::string filterColumnRegex;
    std::string filteringFile;
    std::string mappingFile;
    std::string filterExpression;
    bool positiveFilter;
    bool trimToOneColumn;
    int extractLines;
    float compValue;
    std::string compOperator;
    int sortEntries;
    bool beatsFirst;
    std::string joinDB;
    std::string compPos ;
    std::string clusterFile ;

    // besthitperset
    bool simpleBestHit;
    float alpha;
    bool shortOutput;

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
    bool takeLargerEntry;

    // offsetalignments
    int chainAlignment;
    int mergeQuery;

    // tsv2db
    int outputDbType;

    // diff
    bool useSequenceId;

    //prefixid
    std::string prefix;
    bool tsvOut;

    // clusterUpdate;
    bool recoverDeleted;

    // summarize headers
    int headerType;

    // filtertaxdb
    std::string taxonList;
    bool invertSelection;

    // view
    std::string idList;
    int idxEntryType;
    // lca
    std::string lcaRanks;
    std::string blacklist;

    // exapandaln
    int expansionMode;

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
                           unsigned int outputFlag);
    void printParameters(const std::string &module, int argc, const char* pargv[],
                         const std::vector<MMseqsParameter*> &par);

    std::vector<MMseqsParameter*> removeParameter(const std::vector<MMseqsParameter*>& par, const MMseqsParameter& x);

    PARAMETER(PARAM_S)
    PARAMETER(PARAM_K)
    PARAMETER(PARAM_THREADS)
    PARAMETER(PARAM_COMPRESSED)
    PARAMETER(PARAM_ALPH_SIZE)
    PARAMETER(PARAM_MAX_SEQ_LEN)
//    PARAMETER(PARAM_QUERY_PROFILE)
//    PARAMETER(PARAM_TARGET_PROFILE)
    //PARAMETER(PARAM_NUCL)
    PARAMETER(PARAM_DIAGONAL_SCORING)
    PARAMETER(PARAM_EXACT_KMER_MATCHING)
    PARAMETER(PARAM_MASK_RESIDUES)
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
    PARAMETER(PARAM_SPACED_KMER_MODE)
    PARAMETER(PARAM_REMOVE_TMP_FILES)
    PARAMETER(PARAM_INCLUDE_IDENTITY)
    PARAMETER(PARAM_RES_LIST_OFFSET)
    PARAMETER(PARAM_PRELOAD_MODE)
    PARAMETER(PARAM_SPACED_KMER_PATTERN)
    PARAMETER(PARAM_LOCAL_TMP)
    std::vector<MMseqsParameter*> prefilter;
    std::vector<MMseqsParameter*> ungappedprefilter;

    // alignment
    PARAMETER(PARAM_ALIGNMENT_MODE)
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
    PARAMETER(PARAM_ALT_ALIGNMENT)
    PARAMETER(PARAM_GAP_OPEN)
    PARAMETER(PARAM_GAP_EXTEND)
    std::vector<MMseqsParameter*> align;

    // clustering
    PARAMETER(PARAM_CLUSTER_MODE)
    PARAMETER(PARAM_CLUSTER_STEPS)
    PARAMETER(PARAM_CASCADED)

    // affinity clustering
    PARAMETER(PARAM_MAXITERATIONS)
    PARAMETER(PARAM_SIMILARITYSCORE)

    // logging
    PARAMETER(PARAM_V)
    std::vector<MMseqsParameter*> clust;

    // create profile (HMM, PSSM)
    PARAMETER(PARAM_PROFILE_TYPE)

    // format alignment
    PARAMETER(PARAM_FORMAT_MODE)
    PARAMETER(PARAM_FORMAT_OUTPUT)
    PARAMETER(PARAM_DB_OUTPUT)

    // rescoremode
    PARAMETER(PARAM_RESCORE_MODE)
    PARAMETER(PARAM_FILTER_HITS)
    PARAMETER(PARAM_GLOBAL_ALIGNMENT)
    PARAMETER(PARAM_SORT_RESULTS)

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
    PARAMETER(PARAM_MASK_PROFILE)
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
    PARAMETER(PARAM_TARGET_COLUMN)
    PARAMETER(PARAM_FIRST_SEQ_REP_SEQ)
    PARAMETER(PARAM_FULL_HEADER)
    PARAMETER(PARAM_IDX_SEQ_SRC)

    // result2stat
    PARAMETER(PARAM_STAT)

    // linearcluster
    PARAMETER(PARAM_KMER_PER_SEQ)
    PARAMETER(PARAM_INCLUDE_ONLY_EXTENDABLE)
    PARAMETER(PARAM_SKIP_N_REPEAT_KMER)
    PARAMETER(PARAM_HASH_SHIFT)
    PARAMETER(PARAM_PICK_N_SIMILAR)
    PARAMETER(PARAM_ADJUST_KMER_LEN)

    // workflow
    PARAMETER(PARAM_RUNNER)
    PARAMETER(PARAM_REUSELATEST)

    // search workflow
    PARAMETER(PARAM_NUM_ITERATIONS)
    PARAMETER(PARAM_START_SENS)
    PARAMETER(PARAM_SENS_STEPS)
    PARAMETER(PARAM_SLICE_SEARCH)
    PARAMETER(PARAM_STRAND)


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
    PARAMETER(PARAM_CHECK_COMPATIBLE)
    PARAMETER(PARAM_SEARCH_TYPE)

    // createdb
    PARAMETER(PARAM_USE_HEADER) // also used by extractorfs
    PARAMETER(PARAM_ID_OFFSET)  // same
    PARAMETER(PARAM_DB_TYPE)
    PARAMETER(PARAM_DONT_SPLIT_SEQ_BY_LEN)
    PARAMETER(PARAM_DONT_SHUFFLE)

    // convert2fasta
    PARAMETER(PARAM_USE_HEADER_FILE)

    // split sequence
    PARAMETER(PARAM_SEQUENCE_OVERLAP)

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
    PARAMETER(PARAM_COMPUTE_POSITIONS)
    PARAMETER(PARAM_TRANSITIVE_REPLACE)

    //besthitperset
    PARAMETER(PARAM_SIMPLE_BEST_HIT)
    PARAMETER(PARAM_ALPHA)
    PARAMETER(PARAM_SHORT_OUTPUT)

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

    // filtertaxdb
    PARAMETER(PARAM_TAXON_LIST)
    PARAMETER(PARAM_INVERT_SELECTION)

    // view
    PARAMETER(PARAM_ID_LIST)
    PARAMETER(PARAM_IDX_ENTRY_TYPE)

    // lca
    PARAMETER(PARAM_LCA_RANKS)
    PARAMETER(PARAM_BLACKLIST)

    // exapandaln
    PARAMETER(PARAM_EXPANSION_MODE)

    // taxonomy
    PARAMETER(PARAM_LCA_MODE)

    std::vector<MMseqsParameter*> empty;
    std::vector<MMseqsParameter*> onlyverbosity;
    std::vector<MMseqsParameter*> view;
    std::vector<MMseqsParameter*> verbandcompression;
    std::vector<MMseqsParameter*> onlythreads;
    std::vector<MMseqsParameter*> threadsandcompression;

    std::vector<MMseqsParameter*> rescorediagonal;
    std::vector<MMseqsParameter*> alignbykmer;
    std::vector<MMseqsParameter*> createFasta;
    std::vector<MMseqsParameter*> convertprofiledb;
    std::vector<MMseqsParameter*> sequence2profile;
    std::vector<MMseqsParameter*> result2profile;
    std::vector<MMseqsParameter*> result2pp;
    std::vector<MMseqsParameter*> result2msa;
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
    std::vector<MMseqsParameter*> indexdb;
    std::vector<MMseqsParameter*> kmerindexdb;
    std::vector<MMseqsParameter*> createindex;
    std::vector<MMseqsParameter*> createlinindex;
    std::vector<MMseqsParameter*> convertalignments;
    std::vector<MMseqsParameter*> createdb;
    std::vector<MMseqsParameter*> convert2fasta;
    std::vector<MMseqsParameter*> result2flat;
    std::vector<MMseqsParameter*> gff2ffindex;
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
    std::vector<MMseqsParameter*> subtractdbs;
    std::vector<MMseqsParameter*> diff;
    std::vector<MMseqsParameter*> concatdbs;
    std::vector<MMseqsParameter*> mergedbs;
    std::vector<MMseqsParameter*> summarizeheaders;
    std::vector<MMseqsParameter*> prefixid;
    std::vector<MMseqsParameter*> summarizeresult;
    std::vector<MMseqsParameter*> summarizetabs;
    std::vector<MMseqsParameter*> extractdomains;
    std::vector<MMseqsParameter*> extractalignedregion;
    std::vector<MMseqsParameter*> convertkb;
    std::vector<MMseqsParameter*> tsv2db;
    std::vector<MMseqsParameter*> lca;
    std::vector<MMseqsParameter*> filtertaxdb;
    std::vector<MMseqsParameter*> taxonomy;
    std::vector<MMseqsParameter*> profile2pssm;
    std::vector<MMseqsParameter*> profile2cs;
    std::vector<MMseqsParameter*> besthitbyset;
    std::vector<MMseqsParameter*> combinepvalbyset;
    std::vector<MMseqsParameter*> summerizeresultsbyset;
    std::vector<MMseqsParameter*> multihitdb;
    std::vector<MMseqsParameter*> multihitsearch;
    std::vector<MMseqsParameter*> expandaln;
    std::vector<MMseqsParameter*> sortresult;
    std::vector<MMseqsParameter*> enrichworkflow;

    std::vector<MMseqsParameter*> combineList(const std::vector<MMseqsParameter*> &par1,
                                             const std::vector<MMseqsParameter*> &par2);

    size_t hashParameter(const std::vector<std::string> &filenames, const std::vector<MMseqsParameter*> &par);

    std::string createParameterString(const std::vector<MMseqsParameter*> &vector, bool wasSet = false);

    void overrideParameterDescription(Command& command, int uid, const char* description, const char* regex = NULL, int category = 0);

    static bool isEqualDbtype(const int type1, const int type2) {
        return ((type1 & 0x7FFFFFFF) == (type2 & 0x7FFFFFFF));
    }

protected:
    Parameters();
    static Parameters* instance;
    virtual ~Parameters() {};

private:
    Parameters(Parameters const&);
    void operator=(Parameters const&);
};

#endif
