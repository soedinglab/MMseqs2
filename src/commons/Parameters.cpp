#include "Parameters.h"
#include "Sequence.h"
#include "Debug.h"
#include "Util.h"
#include <unistd.h>
#include "DistanceCalculator.h"

#include <iomanip>
#include <regex.h>
#include <string>
#include <sstream>
#include <climits>

#ifdef OPENMP
#include <omp.h>
#endif

Parameters* Parameters::instance = NULL;

extern const char* binary_name;
extern const char* version;

Parameters::Parameters():
        PARAM_S(PARAM_S_ID,"-s", "Sensitivity","sensitivity: 1.0 faster; 4.0 fast default; 7.5 sensitive [1.0,7.5]", typeid(float), (void *) &sensitivity, "^[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_PREFILTER),
        PARAM_K(PARAM_K_ID,"-k", "K-mer size", "k-mer size in the range [6,7] (0: set automatically to optimum)",typeid(int),  (void *) &kmerSize, "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PREFILTER|MMseqsParameter::COMMAND_CLUSTLINEAR|MMseqsParameter::COMMAND_EXPERT),
        PARAM_THREADS(PARAM_THREADS_ID,"--threads", "Threads", "number of cores used for the computation (uses all cores by default)",typeid(int), (void *) &threads, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_COMMON),
        PARAM_ALPH_SIZE(PARAM_ALPH_SIZE_ID,"--alph-size", "Alphabet size", "alphabet size [2,21]",typeid(int),(void *) &alphabetSize, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PREFILTER|MMseqsParameter::COMMAND_CLUSTLINEAR|MMseqsParameter::COMMAND_EXPERT),
        // Regex for Range 1-32768
        // Please do not change manually, use a tool to regenerate
        // e.g.: http://gamon.webfactional.com/regexnumericrangegenerator/
        PARAM_MAX_SEQ_LEN(PARAM_MAX_SEQ_LEN_ID,"--max-seq-len","Max. sequence length", "Maximum sequence length [1,32768]",typeid(int), (void *) &maxSeqLen, "^0*([1-9]|[1-8][0-9]|9[0-9]|[1-8][0-9]{2}|9[0-8][0-9]|99[0-9]|[1-8][0-9]{3}|9[0-8][0-9]{2}|99[0-8][0-9]|999[0-9]|[12][0-9]{4}|3[01][0-9]{3}|32[0-6][0-9]{2}|327[0-5][0-9]|3276[0-8])$", MMseqsParameter::COMMAND_COMMON|MMseqsParameter::COMMAND_EXPERT),
        PARAM_DIAGONAL_SCORING(PARAM_DIAGONAL_SCORING_ID,"--diag-score", "Diagonal Scoring", "use diagonal score for sorting the prefilter results [0,1]", typeid(int),(void *) &diagonalScoring, "^[0-1]{1}$", MMseqsParameter::COMMAND_PREFILTER|MMseqsParameter::COMMAND_EXPERT),
        PARAM_MASK_RESIDUES(PARAM_MASK_RESIDUES_ID,"--mask", "Mask Residues", "0: w/o low complexity masking, 1: with low complexity masking, (createindex only) 2: add both masked and unmasked sequences to index", typeid(int),(void *) &maskMode, "^[0-2]{1}", MMseqsParameter::COMMAND_PREFILTER|MMseqsParameter::COMMAND_EXPERT),
        PARAM_MIN_DIAG_SCORE(PARAM_MIN_DIAG_SCORE_ID,"--min-ungapped-score", "Minimum Diagonal score", "accept only matches with ungapped alignment score above this threshold", typeid(int),(void *) &minDiagScoreThr, "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PREFILTER|MMseqsParameter::COMMAND_EXPERT),
        PARAM_K_SCORE(PARAM_K_SCORE_ID,"--k-score", "K-score", "k-mer threshold for generating similar-k-mer lists",typeid(int),(void *) &kmerScore,  "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PREFILTER|MMseqsParameter::COMMAND_EXPERT),
        PARAM_MAX_SEQS(PARAM_MAX_SEQS_ID,"--max-seqs", "Max. results per query", "maximum result sequences per query (this parameter affects the sensitivity)",typeid(int),(void *) &maxResListLen, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_COMMON|MMseqsParameter::COMMAND_EXPERT),
        PARAM_SPLIT(PARAM_SPLIT_ID,"--split", "Split DB", "Splits input sets into N equally distributed chunks. The default value sets the best split automatically. createindex can only be used with split 1.",typeid(int),(void *) &split,  "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PREFILTER|MMseqsParameter::COMMAND_EXPERT),
        PARAM_SPLIT_MODE(PARAM_SPLIT_MODE_ID,"--split-mode", "Split mode", "0: split target db; 1: split query db;  2: auto, depending on main memory",typeid(int),(void *) &splitMode,  "^[0-2]{1}$", MMseqsParameter::COMMAND_PREFILTER|MMseqsParameter::COMMAND_EXPERT),
        PARAM_SPLIT_MEMORY_LIMIT(PARAM_SPLIT_MEMORY_LIMIT_ID, "--split-memory-limit", "Split Memory Limit", "Maximum system memory in megabyte that one split may use. Defaults (0) to all available system memory.", typeid(int), (void*) &splitMemoryLimit, "^(0|[1-9]{1}[0-9]*)$", MMseqsParameter::COMMAND_PREFILTER|MMseqsParameter::COMMAND_EXPERT),
        PARAM_SPLIT_AMINOACID(PARAM_SPLIT_AMINOACID_ID,"--split-aa", "Split by amino acid","Try to find the best split for the target database by amino acid count instead",typeid(bool), (void *) &splitAA, "$", MMseqsParameter::COMMAND_EXPERT),
        PARAM_SUB_MAT(PARAM_SUB_MAT_ID,"--sub-mat", "Sub Matrix", "amino acid substitution matrix file",typeid(std::string),(void *) &scoringMatrixFile, "", MMseqsParameter::COMMAND_COMMON|MMseqsParameter::COMMAND_EXPERT),
        PARAM_NO_COMP_BIAS_CORR(PARAM_NO_COMP_BIAS_CORR_ID,"--comp-bias-corr", "Compositional bias","correct for locally biased amino acid composition [0,1]",typeid(int), (void *) &compBiasCorrection, "^[0-1]{1}$", MMseqsParameter::COMMAND_PREFILTER|MMseqsParameter::COMMAND_ALIGN|MMseqsParameter::COMMAND_PROFILE|MMseqsParameter::COMMAND_EXPERT),
        PARAM_SPACED_KMER_MODE(PARAM_SPACED_KMER_MODE_ID,"--spaced-kmer-mode", "Spaced Kmer", "0: use consecutive positions a k-mers; 1: use spaced k-mers",typeid(int), (void *) &spacedKmer,  "^[0-1]{1}", MMseqsParameter::COMMAND_PREFILTER|MMseqsParameter::COMMAND_EXPERT),
        PARAM_REMOVE_TMP_FILES(PARAM_REMOVE_TMP_FILES_ID, "--remove-tmp-files", "Remove Temporary Files" , "Delete temporary files", typeid(bool), (void *) &removeTmpFiles, "",MMseqsParameter::COMMAND_EXPERT),
        PARAM_INCLUDE_IDENTITY(PARAM_INCLUDE_IDENTITY_ID,"--add-self-matches", "Include identical Seq. Id.","artificially add entries of queries with themselves (for clustering)",typeid(bool), (void *) &includeIdentity, "", MMseqsParameter::COMMAND_PREFILTER|MMseqsParameter::COMMAND_ALIGN|MMseqsParameter::COMMAND_EXPERT),
        PARAM_RES_LIST_OFFSET(PARAM_RES_LIST_OFFSET_ID,"--offset-result", "Offset result","Offset result list",typeid(int), (void *) &resListOffset, "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PREFILTER|MMseqsParameter::COMMAND_EXPERT),
        PARAM_NO_PRELOAD(PARAM_NO_PRELOAD_ID, "--no-preload", "No preload", "Do not preload database", typeid(bool), (void*) &noPreload, "", MMseqsParameter::COMMAND_MISC|MMseqsParameter::COMMAND_EXPERT),
        PARAM_EARLY_EXIT(PARAM_EARLY_EXIT_ID, "--early-exit", "Early exit", "Exit immediately after writing the result", typeid(bool), (void*) &earlyExit, "", MMseqsParameter::COMMAND_MISC|MMseqsParameter::COMMAND_EXPERT),
        // alignment
        PARAM_ALIGNMENT_MODE(PARAM_ALIGNMENT_MODE_ID,"--alignment-mode", "Alignment mode", "What to compute: 0: automatic; 1: score+end_pos; 2:+start_pos+cov; 3: +seq.id",typeid(int), (void *) &alignmentMode, "^[0-4]{1}$", MMseqsParameter::COMMAND_ALIGN|MMseqsParameter::COMMAND_EXPERT),
        PARAM_E(PARAM_E_ID,"-e", "E-value threshold", "list matches below this E-value [0.0, inf]",typeid(float), (void *) &evalThr, "^([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)|[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_ALIGN),
        PARAM_C(PARAM_C_ID,"-c", "Coverage threshold", "list matches above this fraction of aligned (covered) residues (see --cov-mode)",typeid(float), (void *) &covThr, "^0(\\.[0-9]+)?|^1(\\.0+)?$", MMseqsParameter::COMMAND_ALIGN| MMseqsParameter::COMMAND_CLUSTLINEAR),
        PARAM_COV_MODE(PARAM_COV_MODE_ID, "--cov-mode", "Coverage Mode", "0: coverage of query and target, 1: coverage of target, 2: coverage of query", typeid(int), (void *) &covMode, "^[0-2]{1}$", MMseqsParameter::COMMAND_ALIGN),
        PARAM_MAX_REJECTED(PARAM_MAX_REJECTED_ID,"--max-rejected", "Max Reject", "maximum rejected alignments before alignment calculation for a query is aborted",typeid(int),(void *) &maxRejected, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_ALIGN),
        PARAM_MAX_ACCEPT(PARAM_MAX_ACCEPT_ID,"--max-accept", "Max Accept", "maximum accepted alignments before alignment calculation for a query is stopped",typeid(int),(void *) &maxAccept, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_ALIGN),
        PARAM_ADD_BACKTRACE(PARAM_ADD_BACKTRACE_ID, "-a", "Add backtrace", "add backtrace string (convert to alignments with mmseqs convertalis utility)", typeid(bool), (void *) &addBacktrace, "", MMseqsParameter::COMMAND_ALIGN),
        PARAM_REALIGN(PARAM_REALIGN_ID, "--realign", "Realign hit", "compute more conservative, shorter alignments (scores and E-values not changed)", typeid(bool), (void *) &realign, "", MMseqsParameter::COMMAND_ALIGN|MMseqsParameter::COMMAND_EXPERT),
        PARAM_MIN_SEQ_ID(PARAM_MIN_SEQ_ID_ID,"--min-seq-id", "Seq. Id Threshold","list matches above this sequence identity (for clustering) [0.0,1.0]",typeid(float), (void *) &seqIdThr, "^0(\\.[0-9]+)?|1(\\.0+)?$", MMseqsParameter::COMMAND_ALIGN),
        // clustering
        PARAM_CLUSTER_MODE(PARAM_CLUSTER_MODE_ID,"--cluster-mode", "Cluster mode", "0: Setcover, 1: connected component, 2: Greedy clustering by sequence length",typeid(int), (void *) &clusteringMode, "[0-2]{1}$", MMseqsParameter::COMMAND_CLUST),
        PARAM_CLUSTER_STEPS(PARAM_CLUSTER_STEPS_ID,"--cluster-steps", "Cascaded clustering steps", "cascaded clustering steps from 1 to -s",typeid(int), (void *) &clusterSteps, "^[1-9]{1}$", MMseqsParameter::COMMAND_CLUST|MMseqsParameter::COMMAND_EXPERT),
        PARAM_CASCADED(PARAM_CASCADED_ID,"--single-step-clustering", "Single step clustering", "switches from cascaded to simple clustering workflow",typeid(bool), (void *) &cascaded, "", MMseqsParameter::COMMAND_CLUST),
        //affinity clustering
        PARAM_MAXITERATIONS(PARAM_MAXITERATIONS_ID,"--max-iterations", "Max depth connected component", "maximum depth of breadth first search in connected component",typeid(int), (void *) &maxIteration,  "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_CLUST|MMseqsParameter::COMMAND_EXPERT),
        PARAM_SIMILARITYSCORE(PARAM_SIMILARITYSCORE_ID,"--similarity-type", "Similarity type", "type of score used for clustering [1:2]. 1=alignment score. 2=sequence identity ",typeid(int),(void *) &similarityScoreType,  "^[1-2]{1}$", MMseqsParameter::COMMAND_CLUST|MMseqsParameter::COMMAND_EXPERT),
        // logging
        PARAM_V(PARAM_V_ID,"-v", "Verbosity","verbosity level: 0=nothing, 1: +errors, 2: +warnings, 3: +info",typeid(int), (void *) &verbosity, "^[0-3]{1}$", MMseqsParameter::COMMAND_COMMON),
        // create profile (HMM)
        PARAM_PROFILE_TYPE(PARAM_PROFILE_TYPE_ID,"--profile-type", "Profile type", "0: HMM (HHsuite) 1: PSSM or 2: HMMER3",typeid(int),(void *) &profileMode,  "^[0-2]{1}$"),
        // convertalignments
        PARAM_FORMAT_MODE(PARAM_FORMAT_MODE_ID,"--format-mode", "Alignment Format", "output format 0: BLAST-TAB, 1: PAIRWISE, 2: BLAST-TAB + query/db length", typeid(int), (void*) &formatAlignmentMode, "^[0-2]{1}$"),
        PARAM_DB_OUTPUT(PARAM_DB_OUTPUT_ID, "--db-output", "Database Output", "Output a result db instead of a text file", typeid(bool), (void*) &dbOut, ""),
        // rescorediagonal
        PARAM_RESCORE_MODE(PARAM_RESCORE_MODE_ID,"--rescore-mode", "Rescore mode", "rescore diagonal by: 0 hamming distance, 1 local alignment (score only) or 2 local alignment", typeid(int), (void *) &rescoreMode, "^[0-2]{1}$"),
        PARAM_FILTER_HITS(PARAM_FILTER_HITS_ID,"--filter-hits", "Remove hits by seq.id. and coverage", "filter hits by seq.id. and coverage", typeid(bool), (void *) &filterHits, "", MMseqsParameter::COMMAND_EXPERT),
        PARAM_GLOBAL_ALIGNMENT(PARAM_GLOBAL_ALIGNMENT_ID,"--global-alignment", "In substitution scoring mode, performs global alignment along the diagonal", "Rescore the complete diagonal", typeid(bool), (void *) &globalAlignment, "", MMseqsParameter::COMMAND_EXPERT),
        // result2msa
        PARAM_ALLOW_DELETION(PARAM_ALLOW_DELETION_ID,"--allow-deletion", "Allow Deletion", "allow deletions in a MSA", typeid(bool), (void*) &allowDeletion, ""),
        PARAM_ADD_INTERNAL_ID(PARAM_ADD_INTERNAL_ID_ID,"--add-iternal-id", "Add internal id", "add internal id as comment to MSA", typeid(bool), (void*) &addInternalId, ""),
        PARAM_COMPRESS_MSA(PARAM_COMPRESS_MSA_ID,"--compress", "Compress MSA", "create MSA in ca3m format", typeid(bool), (void*) &compressMSA, ""),
        PARAM_SUMMARIZE_HEADER(PARAM_SUMMARIZE_HEADER_ID,"--summarize", "Summarize headers", "summarize cluster headers into a single header description", typeid(bool), (void*) &summarizeHeader, ""),
        PARAM_SUMMARY_PREFIX(PARAM_SUMMARY_PREFIX_ID, "--summary-prefix", "Summary prefix","sets the cluster summary prefix",typeid(std::string),(void *) &summaryPrefix, ""),
        PARAM_OMIT_CONSENSUS(PARAM_OMIT_CONSENSUS_ID, "--omit-consensus", "Omit Consensus", "Omit consensus sequence in alignment", typeid(bool), (void*) &omitConsensus, "", MMseqsParameter::COMMAND_EXPERT),
        PARAM_SKIP_QUERY(PARAM_SKIP_QUERY_ID, "--skip-query", "Skip Query", "Skip the query sequence", typeid(bool), (void*) &skipQuery, ""),
        // result2profile
        PARAM_E_PROFILE(PARAM_E_PROFILE_ID,"--e-profile", "Profile e-value threshold", "includes sequences matches with < e-value thr. into the profile [>=0.0]", typeid(float), (void *) &evalProfile, "^([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)|([0-9]*(\\.[0-9]+)?)$", MMseqsParameter::COMMAND_PROFILE),
        PARAM_FILTER_MSA(PARAM_FILTER_MSA_ID,"--filter-msa", "Filter MSA", "filter msa: 0: do not filter, 1: filter", typeid(int), (void*) &filterMsa, "^[0-1]{1}$", MMseqsParameter::COMMAND_PROFILE|MMseqsParameter::COMMAND_EXPERT),
        PARAM_FILTER_MAX_SEQ_ID(PARAM_FILTER_MAX_SEQ_ID_ID,"--max-seq-id", "Maximum sequence identity threshold", "reduce redundancy of output MSA using max. pairwise sequence identity [0.0,1.0]", typeid(float), (void*) &filterMaxSeqId, "^0(\\.[0-9]+)?|1(\\.0+)?$", MMseqsParameter::COMMAND_PROFILE|MMseqsParameter::COMMAND_EXPERT),
        PARAM_FILTER_QSC(PARAM_FILTER_QSC_ID, "--qsc", "Minimum score per column", "reduce diversity of output MSAs using min. score per aligned residue with query sequences [-50.0,100.0]", typeid(float), (void*) &qsc, "^\\-*[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_PROFILE|MMseqsParameter::COMMAND_EXPERT),
        PARAM_FILTER_QID(PARAM_FILTER_QID_ID, "--qid", "Minimum seq. id.", "reduce diversity of output MSAs using min.seq. identity with query sequences [0.0,1.0]", typeid(float), (void*) &qid, "^0(\\.[0-9]+)?|1(\\.0+)?$", MMseqsParameter::COMMAND_PROFILE|MMseqsParameter::COMMAND_EXPERT),
        PARAM_FILTER_COV(PARAM_FILTER_COV_ID, "--cov", "Minimum coverage", "filter output MSAs using min. fraction of query residues covered by matched sequences [0.0,1.0]", typeid(float), (void*) &cov, "^0(\\.[0-9]+)?|1(\\.0+)?$", MMseqsParameter::COMMAND_PROFILE|MMseqsParameter::COMMAND_EXPERT),
        PARAM_FILTER_NDIFF(PARAM_FILTER_NDIFF_ID, "--diff", "Select n most diverse seqs", "filter MSAs by selecting most diverse set of sequences, keeping at least this many seqs in each MSA block of length 50", typeid(int), (void*) &Ndiff, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PROFILE|MMseqsParameter::COMMAND_EXPERT),
        PARAM_WG(PARAM_WG_ID, "--wg", "Use global sequence weighting", "use global sequence weighting for profile calculation", typeid(bool), (void*) &wg, "", MMseqsParameter::COMMAND_PROFILE|MMseqsParameter::COMMAND_EXPERT),
        PARAM_PCA(PARAM_PCA_ID, "--pca", "Pseudo count a", "pseudo count admixture strength", typeid(float), (void*) &pca, "^[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_PROFILE|MMseqsParameter::COMMAND_EXPERT),
        PARAM_PCB(PARAM_PCB_ID, "--pcb", "Pseudo count b", "pseudo counts: Neff at half of maximum admixture (0.0,infinity)", typeid(float), (void*) &pcb, "^[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_PROFILE|MMseqsParameter::COMMAND_EXPERT),
        // msa2profile
        PARAM_MATCH_MODE(PARAM_MATCH_MODE_ID, "--match-mode", "Match mode", "0: Columns that have a residue in the first sequence are kept, 1: columns that have a residue in --match-ratio of all sequences are kept.", typeid(int), (void*)&matchMode, "^(0|1)$", MMseqsParameter::COMMAND_PROFILE),
        PARAM_MATCH_RATIO(PARAM_MATCH_RATIO_ID, "--match-ratio", "Match ratio", "columns that have a residue in this ratio of all sequences are kept", typeid(float), (void*)&matchRatio, "^0(\\.[0-9]+)?|1(\\.0+)?$", MMseqsParameter::COMMAND_PROFILE),
        // sequence2profile
        PARAM_NEFF(PARAM_NEFF_ID, "--neff", "Neff", "Neff included into context state profile (1.0,20.0)", typeid(float), (void*) &neff, "^[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_PROFILE),
        PARAM_TAU(PARAM_TAU_ID, "--tau", "Tau", "Tau: context state pseudo count mixture (0.0,1.0)", typeid(float), (void*) &tau, "[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_PROFILE),
        //createtesv
        PARAM_FIRST_SEQ_REP_SEQ(PARAM_FIRST_SEQ_REP_SEQ_ID, "--first-seq-as-repr", "first sequence as respresentative", "Use the first sequence of the clustering result as representative sequence", typeid(bool), (void*) &firstSeqRepr, "", MMseqsParameter::COMMAND_MISC),
        // result2stats
        PARAM_STAT(PARAM_STAT_ID, "--stat", "Statistics to be computed", "can be one of: linecount, mean, doolittle, charges, seqlen, firstline.", typeid(std::string), (void*) &stat, ""),
        // linearcluster
        PARAM_KMER_PER_SEQ(PARAM_KMER_PER_SEQ_ID, "--kmer-per-seq", "Kmer per sequence", "kmer per sequence", typeid(int), (void*) &kmersPerSequence, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_CLUSTLINEAR),
        PARAM_INCLUDE_ONLY_EXTENDABLE(PARAM_INCLUDE_ONLY_EXTENDABLE_ID, "--inlude-only-extendable", "Include only extendable", "Include only extendable", typeid(bool), (void*) &includeOnlyExtendable, "", MMseqsParameter::COMMAND_CLUSTLINEAR),
        PARAM_HASH_SHIFT(PARAM_HASH_SHIFT_ID, "--hash-shift", "Shift hash", "Shift k-mer hash", typeid(int), (void*) &hashShift, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_CLUSTLINEAR|MMseqsParameter::COMMAND_EXPERT),
        // workflow
        PARAM_RUNNER(PARAM_RUNNER_ID, "--mpi-runner", "Sets the MPI runner","use MPI on compute grid with this MPI command (e.g. \"mpirun -np 42\")",typeid(std::string),(void *) &runner, "", MMseqsParameter::COMMAND_EXPERT),
        // search workflow
        PARAM_NUM_ITERATIONS(PARAM_NUM_ITERATIONS_ID, "--num-iterations", "Number search iterations","Search iterations",typeid(int),(void *) &numIterations, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PROFILE),
        PARAM_START_SENS(PARAM_START_SENS_ID, "--start-sens", "Start sensitivity","start sensitivity",typeid(float),(void *) &startSens, "^[0-9]*(\\.[0-9]+)?$"),
        PARAM_SENS_STEPS(PARAM_SENS_STEPS_ID, "--sens-steps", "Search steps","Search steps performed from --start-sense and -s.",typeid(int),(void *) &sensSteps, "^[1-9]{1}$"),
        // Orfs
        PARAM_ORF_MIN_LENGTH(PARAM_ORF_MIN_LENGTH_ID, "--min-length", "Min codons in orf", "minimum codon number in open reading frames",typeid(int),(void *) &orfMinLength, "^[1-9]{1}[0-9]*$"),
        PARAM_ORF_MAX_LENGTH(PARAM_ORF_MAX_LENGTH_ID, "--max-length", "Max codons in length", "maximum codon number in open reading frames",typeid(int),(void *) &orfMaxLength, "^[1-9]{1}[0-9]*$"),
        PARAM_ORF_MAX_GAP(PARAM_ORF_MAX_GAP_ID, "--max-gaps", "Max orf gaps", "maximum number of codons with gaps or unknown residues before an open reading frame is rejected",typeid(int),(void *) &orfMaxGaps, "^(0|[1-9]{1}[0-9]*)$"),
        PARAM_ORF_START_STATE(PARAM_ORF_START_STATE_ID,"--orf-start-state", "Orf start state", "Orf start can be 0: incomplete, 1: complete, 2: both",typeid(int),(void *) &orfStartState, "^[0-2]{1}"),
        PARAM_ORF_END_STATE(PARAM_ORF_END_STATE_ID,"--orf-end-state", "Orf end state", "Orf end can be 0: incomplete, 1: complete, 2: both ",typeid(int),(void *) &orfEndState, "^[0-2]{1}"),
        PARAM_ORF_LONGEST(PARAM_ORF_LONGEST_ID,"--longest-orf", "Find longest orf", "does the first found start codon start an orf (results in the longst possible orf)",typeid(bool),(void *) &orfLongest, ""),
        PARAM_ORF_EXTENDMIN(PARAM_ORF_EXTENDMIN_ID,"--extend-min", "Extend short orfs", "if an orf would be rejected because of the min length threshold, allow it to be extended to the next stop codon",typeid(bool),(void *) &orfExtendMin, ""),
        PARAM_ORF_FORWARD_FRAMES(PARAM_ORF_FORWARD_FRAMES_ID, "--forward-frames", "Forward Frames", "comma-seperated list of ORF frames on the forward strand to be extracted", typeid(std::string), (void *) &forwardFrames, ""),
        PARAM_ORF_REVERSE_FRAMES(PARAM_ORF_REVERSE_FRAMES_ID, "--reverse-frames", "Reverse Frames", "comma-seperated list of ORF frames on the reverse strand to be extracted", typeid(std::string), (void *) &reverseFrames, ""),
        // indexdb
        PARAM_INCLUDE_HEADER(PARAM_INCLUDE_HEADER_ID, "--include-headers", "Include Header", "Include the header index into the index", typeid(bool), (void *) &includeHeader, ""),
        // createdb
        PARAM_USE_HEADER(PARAM_USE_HEADER_ID,"--use-fasta-header", "Use fasta header", "use the id parsed from the fasta header as the index key instead of using incrementing numeric identifiers",typeid(bool),(void *) &useHeader, ""),
        PARAM_ID_OFFSET(PARAM_ID_OFFSET_ID, "--id-offset", "Offset of numeric ids", "numeric ids in index file are offset by this value ",typeid(int),(void *) &identifierOffset, "^(0|[1-9]{1}[0-9]*)$"),
        PARAM_DONT_SPLIT_SEQ_BY_LEN(PARAM_DONT_SPLIT_SEQ_BY_LEN_ID,"--dont-split-seq-by-len", "Split Seq. by len", "Dont split sequences by --max-seq-len",typeid(bool),(void *) &splitSeqByLen, ""),
        PARAM_USE_HEADER_FILE(PARAM_USE_HEADER_FILE_ID, "--use-header-file", "Use ffindex header", "use the ffindex header file instead of the body to map the entry keys",typeid(bool),(void *) &useHeaderFile, ""),
        PARAM_GFF_TYPE(PARAM_GFF_TYPE_ID,"--gff-type", "GFF Type", "type in the GFF file to filter by",typeid(std::string),(void *) &gffType, ""),
        PARAM_TRANSLATION_TABLE(PARAM_TRANSLATION_TABLE_ID,"--translation-table", "Translation Table", "1) CANONICAL, 2) VERT_MITOCHONDRIAL, 3) YEAST_MITOCHONDRIAL, 4) MOLD_MITOCHONDRIAL, 5) INVERT_MITOCHONDRIAL, 6) CILIATE, 9) FLATWORM_MITOCHONDRIAL, 10) EUPLOTID, 11) PROKARYOTE, 12) ALT_YEAST, 13) ASCIDIAN_MITOCHONDRIAL, 14) ALT_FLATWORM_MITOCHONDRIAL, 15) BLEPHARISMA, 16) CHLOROPHYCEAN_MITOCHONDRIAL, 21) TREMATODE_MITOCHONDRIAL, 22) SCENEDESMUS_MITOCHONDRIAL, 23) THRAUSTOCHYTRIUM_MITOCHONDRIAL, 24) PTEROBRANCHIA_MITOCHONDRIAL, 25) GRACILIBACTERIA, 26) PACHYSOLEN, 27) KARYORELICT, 28) CONDYLOSTOMA, 29) MESODINIUM, 30) PERTRICH, 31) BLASTOCRITHIDIA", typeid(int),(void *) &translationTable, "^[1-9]{1}[0-9]*$"),
        PARAM_ADD_ORF_STOP(PARAM_ADD_ORF_STOP_ID,"--add-orf-stop", "Add Orf Stop", "add * at complete start and end", typeid(bool),(void *) &addOrfStop, ""),
        PARAM_MIN_SEQUENCES(PARAM_MIN_SEQUENCES_ID,"--min-sequences", "Min Sequences", "minimum number of sequences a cluster may contain", typeid(int),(void *) &minSequences,"^[1-9]{1}[0-9]*$"),
        PARAM_MAX_SEQUENCES(PARAM_MAX_SEQUENCES_ID,"--max-sequences", "Max Sequences", "maximum number of sequences a cluster may contain", typeid(int),(void *) &maxSequences,"^[1-9]{1}[0-9]*$"),
        PARAM_HH_FORMAT(PARAM_HH_FORMAT_ID,"--hh-format", "HH format", "format entries to use with hhsuite (for singleton clusters)", typeid(bool), (void *) &hhFormat, ""),
        // filterdb
        PARAM_FILTER_COL(PARAM_FILTER_COL_ID,"--filter-column", "Filter column", "column", typeid(int),(void *) &filterColumn,"^[1-9]{1}[0-9]*$"),
        PARAM_FILTER_REGEX(PARAM_FILTER_REGEX_ID,"--filter-regex", "Filter regex", "regex to select column (example float: [0-9]*(.[0-9]+)? int:[1-9]{1}[0-9])", typeid(std::string),(void *) &filterColumnRegex,"^.*$"),
        PARAM_FILTER_POS(PARAM_FILTER_POS_ID,"--positive-filter", "Positive filter", "used in conjunction with --filter-file. If true, out  = in \\intersect filter ; if false, out = in - filter", typeid(bool),(void *) &positiveFilter,""),
        PARAM_FILTER_FILE(PARAM_FILTER_FILE_ID,"--filter-file", "Filter file", "specify a file that contains the filtering elements", typeid(std::string),(void *) &filteringFile,""),
        PARAM_MAPPING_FILE(PARAM_MAPPING_FILE_ID,"--mapping-file", "Mapping file", "specify a file that translates the keys of a DB to new keys", typeid(std::string),(void *) &mappingFile,""),
        PARAM_TRIM_TO_ONE_COL(PARAM_TRIM_TO_ONE_COL_ID,"--trim-to-one-column", "trim the results to one column","Output only the column specified by --filter-column.",typeid(bool), (void *) &trimToOneColumn, ""),
        PARAM_EXTRACT_LINES(PARAM_EXTRACT_LINES_ID,"--extract-lines", "Extract n lines", "extract n lines of each entry.",typeid(int), (void *) &extractLines, "^[1-9]{1}[0-9]*$"),
        PARAM_COMP_OPERATOR(PARAM_COMP_OPERATOR_ID, "--comparison-operator", "Numerical comparison operator", "Filter by comparing each entry row numerically by using the le) less-than-equal, ge) greater-than-equal or e) equal operator.", typeid(std::string), (void *) &compOperator, ""),
        PARAM_COMP_VALUE(PARAM_COMP_VALUE_ID, "--comparison-value", "Numerical comparison value", "Filter by comparing each entry to this value.", typeid(float), (void *) &compValue, ""),
        PARAM_SORT_ENTRIES(PARAM_SORT_ENTRIES_ID, "--sort-entries", "Sort entries", "Sort column set by --filter-column, by 0) no sorting, 1) increasing,  2) decreasing or 3) random shuffle.", typeid(int), (void *) &sortEntries, "^[1-9]{1}[0-9]*$"),
        PARAM_BEATS_FIRST(PARAM_BEATS_FIRST_ID, "--beats-first", "Beats first", "Filter by comparing each entry to the first entry.", typeid(bool), (void*) &beatsFirst, ""),
        // concatdb
        PARAM_PRESERVEKEYS(PARAM_PRESERVEKEYS_ID,"--preserve-keys", "Preserve the keys", "the keys of the two DB should be distinct, and they will be preserved in the concatenation.",typeid(bool), (void *) &preserveKeysB, ""),
        //diff
        PARAM_USESEQID(PARAM_USESEQID_ID,"--use-seq-id", "Match sequences by their ID", "Sequence ID (Uniprot, GenBank, ...) is used for identifying matches between the old and the new DB.",typeid(bool), (void *) &useSequenceId, ""),
        // prefixid
        PARAM_PREFIX(PARAM_PREFIX_ID, "--prefix", "Prefix", "Use this prefix for all entries", typeid(std::string),(void *) &prefix,""),
        // summarize headers
        PARAM_HEADER_TYPE(PARAM_HEADER_TYPE_ID,"--header-type", "Header type", "Header Type: 1 Uniclust, 2 Metaclust",typeid(int), (void *) &headerType, "[1-2]{1}"),
        // mergedbs
        PARAM_MERGE_PREFIXES(PARAM_MERGE_PREFIXES_ID, "--prefixes", "Merge prefixes", "Comma separated list of prefixes for each entry", typeid(std::string),(void *) &mergePrefixes,""),
        // summarizeresult
        PARAM_OVERLAP(PARAM_OVERLAP_ID, "--overlap", "Overlap", "maximum overlap", typeid(float), (void*) &overlap, "^[0-9]*(\\.[0-9]+)?$"),
        // msa2profile
        PARAM_MSA_TYPE(PARAM_MSA_TYPE_ID,"--msa-type", "MSA type", "MSA Type: cA3M 0, A3M 1, FASTA 2", typeid(int), (void *) &msaType, "^[0-2]{1}$"),
        // extractalignedregion
        PARAM_EXTRACT_MODE(PARAM_EXTRACT_MODE_ID,"--extract-mode", "Extract mode", "Query 1, Target 2", typeid(int), (void *) &extractMode, "^[1-2]{1}$"),
        // convertkb
        PARAM_KB_COLUMNS(PARAM_KB_COLUMNS_ID, "--kb-columns", "UniprotKB Columns", "list of indices of UniprotKB columns to be extracted", typeid(std::string), (void *) &kbColumns, ""),
        PARAM_RECOVER_DELETED(PARAM_RECOVER_DELETED_ID, "--recover-deleted", "Recover Deleted", "Indicates if sequences are allowed to be be removed during updating", typeid(bool), (void*) &recoverDeleted, ""),
        PARAM_LCA_RANKS(PARAM_LCA_RANKS_ID, "--lca-ranks", "LCA Ranks", "Ranks to return in LCA computation", typeid(std::string), (void*) &lcaRanks, ""),
        PARAM_BLACKLIST(PARAM_BLACKLIST_ID, "--blacklist", "Blacklisted Taxa", "Comma separted list of ignored taxa in LCA computation", typeid(std::string), (void*)&blacklist, "([0-9]+,)?[0-9]+"),
        PARAM_LCA_MODE(PARAM_LCA_MODE_ID, "--lca-mode", "LCA Mode", "LCA Mode: No LCA 0, Single Search LCA 1, 2bLCA 2", typeid(int), (void*) &lcaMode, "^[0-2]{1}$")
{
    if (instance) {
        Debug(Debug::ERROR) << "Parameter instance already exists!\n";
        abort();
    }
    instance = this;

    // alignment
    align.push_back(PARAM_SUB_MAT);
    align.push_back(PARAM_ADD_BACKTRACE);
    align.push_back(PARAM_ALIGNMENT_MODE);
    align.push_back(PARAM_E);
    align.push_back(PARAM_MIN_SEQ_ID);
    align.push_back(PARAM_C);
    align.push_back(PARAM_COV_MODE);
    align.push_back(PARAM_MAX_SEQ_LEN);
    align.push_back(PARAM_MAX_SEQS);
    align.push_back(PARAM_NO_COMP_BIAS_CORR);
    align.push_back(PARAM_REALIGN);
    align.push_back(PARAM_MAX_REJECTED);
    align.push_back(PARAM_MAX_ACCEPT);
    align.push_back(PARAM_INCLUDE_IDENTITY);
    align.push_back(PARAM_NO_PRELOAD);
    align.push_back(PARAM_EARLY_EXIT);
    align.push_back(PARAM_PCA);
    align.push_back(PARAM_PCB);
    align.push_back(PARAM_THREADS);
    align.push_back(PARAM_V);

    // prefilter
    prefilter.push_back(PARAM_SUB_MAT);
    prefilter.push_back(PARAM_S);
    prefilter.push_back(PARAM_K);
    prefilter.push_back(PARAM_K_SCORE);
    prefilter.push_back(PARAM_ALPH_SIZE);
    prefilter.push_back(PARAM_MAX_SEQ_LEN);
    prefilter.push_back(PARAM_MAX_SEQS);
    prefilter.push_back(PARAM_RES_LIST_OFFSET);
    prefilter.push_back(PARAM_SPLIT);
    prefilter.push_back(PARAM_SPLIT_MODE);
    prefilter.push_back(PARAM_SPLIT_MEMORY_LIMIT);
    prefilter.push_back(PARAM_C);
    prefilter.push_back(PARAM_COV_MODE);
    prefilter.push_back(PARAM_NO_COMP_BIAS_CORR);
    prefilter.push_back(PARAM_DIAGONAL_SCORING);
    prefilter.push_back(PARAM_MASK_RESIDUES);
    prefilter.push_back(PARAM_MIN_DIAG_SCORE);
    prefilter.push_back(PARAM_INCLUDE_IDENTITY);
    prefilter.push_back(PARAM_SPACED_KMER_MODE);
    prefilter.push_back(PARAM_NO_PRELOAD);
    prefilter.push_back(PARAM_EARLY_EXIT);
    prefilter.push_back(PARAM_PCA);
    prefilter.push_back(PARAM_PCB);
    prefilter.push_back(PARAM_THREADS);
    prefilter.push_back(PARAM_V);

    // clustering
    clust.push_back(PARAM_CLUSTER_MODE);
    clust.push_back(PARAM_V);
    clust.push_back(PARAM_MAXITERATIONS);
    clust.push_back(PARAM_SIMILARITYSCORE);
    clust.push_back(PARAM_THREADS);

    // find orf
    onlyverbosity.push_back(PARAM_V);

    // rescorediagonal
    rescorediagonal.push_back(PARAM_RESCORE_MODE);
    rescorediagonal.push_back(PARAM_SUB_MAT);
    rescorediagonal.push_back(PARAM_FILTER_HITS);
    rescorediagonal.push_back(PARAM_GLOBAL_ALIGNMENT);
    rescorediagonal.push_back(PARAM_C);
    rescorediagonal.push_back(PARAM_E);
    rescorediagonal.push_back(PARAM_COV_MODE);
    rescorediagonal.push_back(PARAM_MIN_SEQ_ID);
    rescorediagonal.push_back(PARAM_INCLUDE_IDENTITY);
    rescorediagonal.push_back(PARAM_THREADS);
    rescorediagonal.push_back(PARAM_V);

    // convertprofiledb
    convertprofiledb.push_back(PARAM_SUB_MAT);
    convertprofiledb.push_back(PARAM_PROFILE_TYPE);
    convertprofiledb.push_back(PARAM_THREADS);
    convertprofiledb.push_back(PARAM_V);


    // sequence2profile
    sequence2profile.push_back(PARAM_PCA);
    sequence2profile.push_back(PARAM_PCB);
    sequence2profile.push_back(PARAM_NEFF);
    sequence2profile.push_back(PARAM_TAU);
    sequence2profile.push_back(PARAM_THREADS);
    sequence2profile.push_back(PARAM_SUB_MAT);
    sequence2profile.push_back(PARAM_V);

    // create fasta
    createFasta.push_back(PARAM_V);

    // result2profile
    result2profile.push_back(PARAM_SUB_MAT);
//    result2profile.push_back(PARAM_QUERY_PROFILE);
    result2profile.push_back(PARAM_E_PROFILE);
    result2profile.push_back(PARAM_NO_COMP_BIAS_CORR);
    result2profile.push_back(PARAM_WG);
    result2profile.push_back(PARAM_FILTER_MSA);
    result2profile.push_back(PARAM_FILTER_MAX_SEQ_ID);
    result2profile.push_back(PARAM_FILTER_QID);
    result2profile.push_back(PARAM_FILTER_QSC);
    result2profile.push_back(PARAM_FILTER_COV);
    result2profile.push_back(PARAM_FILTER_NDIFF);
    result2profile.push_back(PARAM_PCA);
    result2profile.push_back(PARAM_PCB);
    result2profile.push_back(PARAM_OMIT_CONSENSUS);
    result2profile.push_back(PARAM_NO_PRELOAD);
    result2profile.push_back(PARAM_EARLY_EXIT);
    result2profile.push_back(PARAM_THREADS);
    result2profile.push_back(PARAM_V);

    // createtsv
    createtsv.push_back(PARAM_FIRST_SEQ_REP_SEQ);

    //result2stats
    result2stats.push_back(PARAM_STAT);
    result2stats.push_back(PARAM_THREADS);
    result2stats.push_back(PARAM_V);

    // format alignment
    convertalignments.push_back(PARAM_FORMAT_MODE);
    convertalignments.push_back(PARAM_NO_PRELOAD);
    convertalignments.push_back(PARAM_EARLY_EXIT);
    convertalignments.push_back(PARAM_DB_OUTPUT);
    convertalignments.push_back(PARAM_THREADS);
    convertalignments.push_back(PARAM_V);

    // result2msa
    result2msa.push_back(PARAM_SUB_MAT);
    result2msa.push_back(PARAM_E_PROFILE);
    result2msa.push_back(PARAM_ALLOW_DELETION);
    result2msa.push_back(PARAM_ADD_INTERNAL_ID);
    result2msa.push_back(PARAM_NO_COMP_BIAS_CORR);
    result2msa.push_back(PARAM_FILTER_MSA);
    result2msa.push_back(PARAM_FILTER_MAX_SEQ_ID);
    result2msa.push_back(PARAM_FILTER_QID);
    result2msa.push_back(PARAM_FILTER_QSC);
    result2msa.push_back(PARAM_FILTER_COV);
    result2msa.push_back(PARAM_FILTER_NDIFF);
    result2msa.push_back(PARAM_THREADS);
    result2msa.push_back(PARAM_V);
    result2msa.push_back(PARAM_COMPRESS_MSA);
    result2msa.push_back(PARAM_SUMMARIZE_HEADER);
    result2msa.push_back(PARAM_SUMMARY_PREFIX);
    result2msa.push_back(PARAM_OMIT_CONSENSUS);
    result2msa.push_back(PARAM_SKIP_QUERY);
    //result2msa.push_back(PARAM_FIRST_SEQ_REP_SEQ);

    // msa2profile
    msa2profile.push_back(PARAM_MSA_TYPE);
    msa2profile.push_back(PARAM_SUB_MAT);
    msa2profile.push_back(PARAM_MATCH_MODE);
    msa2profile.push_back(PARAM_MATCH_RATIO);
    msa2profile.push_back(PARAM_PCA);
    msa2profile.push_back(PARAM_PCB);
    msa2profile.push_back(PARAM_NO_COMP_BIAS_CORR);
    msa2profile.push_back(PARAM_WG);
    msa2profile.push_back(PARAM_FILTER_MSA);
    msa2profile.push_back(PARAM_FILTER_COV);
    msa2profile.push_back(PARAM_FILTER_QID);
    msa2profile.push_back(PARAM_FILTER_QSC);
    msa2profile.push_back(PARAM_FILTER_MAX_SEQ_ID);
    msa2profile.push_back(PARAM_FILTER_NDIFF);
    msa2profile.push_back(PARAM_THREADS);
    msa2profile.push_back(PARAM_V);

    // profile2pssm
    profile2pssm.push_back(PARAM_SUB_MAT);
    profile2pssm.push_back(PARAM_MAX_SEQ_LEN);
    profile2pssm.push_back(PARAM_NO_COMP_BIAS_CORR);
    profile2pssm.push_back(PARAM_DB_OUTPUT);
    profile2pssm.push_back(PARAM_THREADS);
    profile2pssm.push_back(PARAM_V);

    // profile2cs
    profile2cs.push_back(PARAM_SUB_MAT);
    profile2cs.push_back(PARAM_THREADS);
    profile2cs.push_back(PARAM_V);

    // extract orf
    extractorfs.push_back(PARAM_ORF_MIN_LENGTH);
    extractorfs.push_back(PARAM_ORF_MAX_LENGTH);
    extractorfs.push_back(PARAM_ORF_MAX_GAP);
    extractorfs.push_back(PARAM_ORF_START_STATE);
    extractorfs.push_back(PARAM_ORF_END_STATE);
    extractorfs.push_back(PARAM_ORF_LONGEST);
    extractorfs.push_back(PARAM_ORF_EXTENDMIN);
    extractorfs.push_back(PARAM_ORF_FORWARD_FRAMES);
    extractorfs.push_back(PARAM_ORF_REVERSE_FRAMES);
    extractorfs.push_back(PARAM_ID_OFFSET);
    extractorfs.push_back(PARAM_THREADS);

    // splitdb
    splitdb.push_back(PARAM_SPLIT);
    splitdb.push_back(PARAM_SPLIT_AMINOACID);

    // create index
    indexdb.push_back(PARAM_SUB_MAT);
    indexdb.push_back(PARAM_K);
    indexdb.push_back(PARAM_ALPH_SIZE);
    indexdb.push_back(PARAM_MAX_SEQS);
    indexdb.push_back(PARAM_MAX_SEQ_LEN);
    indexdb.push_back(PARAM_MASK_RESIDUES);
    indexdb.push_back(PARAM_SPACED_KMER_MODE);
    indexdb.push_back(PARAM_S);
    indexdb.push_back(PARAM_K_SCORE);
    indexdb.push_back(PARAM_INCLUDE_HEADER);
    indexdb.push_back(PARAM_SPLIT);
    indexdb.push_back(PARAM_SPLIT_MEMORY_LIMIT);
    indexdb.push_back(PARAM_THREADS);
    indexdb.push_back(PARAM_V);

    // create db
    createdb.push_back(PARAM_MAX_SEQ_LEN);
    createdb.push_back(PARAM_DONT_SPLIT_SEQ_BY_LEN);
    createdb.push_back(PARAM_ID_OFFSET);
    createdb.push_back(PARAM_V);

    // convert2fasta
    convert2fasta.push_back(PARAM_USE_HEADER_FILE);
    convert2fasta.push_back(PARAM_V);

    // result2flat
    result2flat.push_back(PARAM_USE_HEADER);
    result2flat.push_back(PARAM_V);

    // gff2db
    gff2ffindex.push_back(PARAM_GFF_TYPE);
    gff2ffindex.push_back(PARAM_ID_OFFSET);
    gff2ffindex.push_back(PARAM_V);


    // translate nucleotide
    translatenucs.push_back(PARAM_TRANSLATION_TABLE);
    translatenucs.push_back(PARAM_ADD_ORF_STOP);
    translatenucs.push_back(PARAM_V);

    // createseqfiledb
    createseqfiledb.push_back(PARAM_MIN_SEQUENCES);
    createseqfiledb.push_back(PARAM_MAX_SEQUENCES);
    createseqfiledb.push_back(PARAM_HH_FORMAT);
    createseqfiledb.push_back(PARAM_THREADS);
    createseqfiledb.push_back(PARAM_V);

    // filterDb
    filterDb.push_back(PARAM_FILTER_COL);
    filterDb.push_back(PARAM_FILTER_REGEX);
    filterDb.push_back(PARAM_FILTER_POS);
    filterDb.push_back(PARAM_FILTER_FILE);
    filterDb.push_back(PARAM_BEATS_FIRST);
    filterDb.push_back(PARAM_MAPPING_FILE);
    filterDb.push_back(PARAM_THREADS);
    filterDb.push_back(PARAM_V);
    filterDb.push_back(PARAM_TRIM_TO_ONE_COL);
    filterDb.push_back(PARAM_EXTRACT_LINES);
    filterDb.push_back(PARAM_COMP_OPERATOR);
    filterDb.push_back(PARAM_COMP_VALUE);
    filterDb.push_back(PARAM_SORT_ENTRIES);
    filterDb.push_back(PARAM_INCLUDE_IDENTITY);


    onlythreads.push_back(PARAM_THREADS);
    onlythreads.push_back(PARAM_V);

    // swap results
    swapresult.push_back(PARAM_SUB_MAT);
    swapresult.push_back(PARAM_E);
    swapresult.push_back(PARAM_SPLIT_MEMORY_LIMIT);
    swapresult.push_back(PARAM_THREADS);

    // subtractdbs
    subtractdbs.push_back(PARAM_THREADS);
    subtractdbs.push_back(PARAM_E_PROFILE);
    subtractdbs.push_back(PARAM_V);

    // clusthash
    clusthash.push_back(PARAM_SUB_MAT);
    clusthash.push_back(PARAM_ALPH_SIZE);
    clusthash.push_back(PARAM_MIN_SEQ_ID);
    clusthash.push_back(PARAM_MAX_SEQ_LEN);
    clusthash.push_back(PARAM_THREADS);
    clusthash.push_back(PARAM_V);

    // kmermatcher
    kmermatcher.push_back(PARAM_SUB_MAT);
    kmermatcher.push_back(PARAM_ALPH_SIZE);
    kmermatcher.push_back(PARAM_MIN_SEQ_ID);
    kmermatcher.push_back(PARAM_KMER_PER_SEQ);
    kmermatcher.push_back(PARAM_MASK_RESIDUES);
    kmermatcher.push_back(PARAM_COV_MODE);
    kmermatcher.push_back(PARAM_K);
    kmermatcher.push_back(PARAM_C);
    kmermatcher.push_back(PARAM_MAX_SEQ_LEN);
    kmermatcher.push_back(PARAM_HASH_SHIFT);
    kmermatcher.push_back(PARAM_SPLIT_MEMORY_LIMIT);
    kmermatcher.push_back(PARAM_THREADS);
    kmermatcher.push_back(PARAM_V);


    // mergedbs
    mergedbs.push_back(PARAM_MERGE_PREFIXES);
    mergedbs.push_back(PARAM_V);

    // summarize
    summarizeheaders.push_back(PARAM_SUMMARY_PREFIX);
    summarizeheaders.push_back(PARAM_HEADER_TYPE);
    summarizeheaders.push_back(PARAM_THREADS);
    summarizeheaders.push_back(PARAM_V);

    // diff
    diff.push_back(PARAM_USESEQID);
    diff.push_back(PARAM_THREADS);
    diff.push_back(PARAM_V);

    // prefixid
    prefixid.push_back(PARAM_PREFIX);
    prefixid.push_back(PARAM_MAPPING_FILE);
    prefixid.push_back(PARAM_THREADS);
    prefixid.push_back(PARAM_V);

    // summarizeresult
    summarizeresult.push_back(PARAM_ADD_BACKTRACE);
    summarizeresult.push_back(PARAM_OVERLAP);
    summarizeresult.push_back(PARAM_E);
    summarizeresult.push_back(PARAM_C);
    summarizeresult.push_back(PARAM_THREADS);
    summarizeresult.push_back(PARAM_V);

    // summarizetabs
    summarizetabs.push_back(PARAM_OVERLAP);
    summarizetabs.push_back(PARAM_E);
    summarizetabs.push_back(PARAM_C);
    summarizetabs.push_back(PARAM_THREADS);
    summarizetabs.push_back(PARAM_V);

    // annoate
    extractdomains.push_back(PARAM_SUB_MAT);
    extractdomains.push_back(PARAM_MSA_TYPE);
    extractdomains.push_back(PARAM_E);
    extractdomains.push_back(PARAM_C);
    extractdomains.push_back(PARAM_THREADS);
    extractdomains.push_back(PARAM_V);

    // concatdbs
    concatdbs.push_back(PARAM_PRESERVEKEYS);
    concatdbs.push_back(PARAM_THREADS);
    concatdbs.push_back(PARAM_V);

    // extractalignedregion
    extractalignedregion.push_back(PARAM_EXTRACT_MODE);
    extractalignedregion.push_back(PARAM_THREADS);
    extractalignedregion.push_back(PARAM_V);

    // convertkb
    convertkb.push_back(PARAM_MAPPING_FILE);
    convertkb.push_back(PARAM_KB_COLUMNS);
    convertkb.push_back(PARAM_V);

    // lca
    lca.push_back(PARAM_LCA_RANKS);
    lca.push_back(PARAM_BLACKLIST);
    lca.push_back(PARAM_V);
    lca.push_back(PARAM_THREADS);

    // WORKFLOWS
    searchworkflow = combineList(align, prefilter);
    searchworkflow = combineList(searchworkflow, result2profile);
    searchworkflow = combineList(searchworkflow, extractorfs);
    searchworkflow = combineList(searchworkflow, translatenucs);
    searchworkflow.push_back(PARAM_NUM_ITERATIONS);
    searchworkflow.push_back(PARAM_START_SENS);
    searchworkflow.push_back(PARAM_SENS_STEPS);
    searchworkflow.push_back(PARAM_RUNNER);
    searchworkflow.push_back(PARAM_REMOVE_TMP_FILES);


    // createindex workflow
    createindex = combineList(indexdb, extractorfs);
    createindex =  combineList(createindex, translatenucs);
    createindex.push_back(PARAM_REMOVE_TMP_FILES);

    // linclust workflow
    linclustworkflow = combineList(clust, align);
    linclustworkflow = combineList(linclustworkflow, kmermatcher);
    linclustworkflow = combineList(linclustworkflow, rescorediagonal);
    linclustworkflow.push_back(PARAM_REMOVE_TMP_FILES);
    linclustworkflow.push_back(PARAM_RUNNER);


    // assembler workflow
    assemblerworkflow = combineList(rescorediagonal, kmermatcher);
    assemblerworkflow.push_back(PARAM_NUM_ITERATIONS);
    assemblerworkflow.push_back(PARAM_REMOVE_TMP_FILES);
    assemblerworkflow.push_back(PARAM_RUNNER);

    // clustering workflow
    clusteringWorkflow = combineList(prefilter, align);
    clusteringWorkflow = combineList(clusteringWorkflow, clust);
    clusteringWorkflow.push_back(PARAM_CASCADED);
    clusteringWorkflow.push_back(PARAM_CLUSTER_STEPS);
    clusteringWorkflow.push_back(PARAM_REMOVE_TMP_FILES);
    clusteringWorkflow.push_back(PARAM_RUNNER);
    clusteringWorkflow = combineList(clusteringWorkflow, linclustworkflow);

    // taxonomy
    taxonomy = combineList(searchworkflow, lca);
    taxonomy.push_back(PARAM_LCA_MODE);
    taxonomy.push_back(PARAM_REMOVE_TMP_FILES);
    taxonomy.push_back(PARAM_RUNNER);


    clusterUpdateSearch = removeParameter(searchworkflow,PARAM_MAX_SEQS);
    clusterUpdateClust = removeParameter(clusteringWorkflow,PARAM_MAX_SEQS);
    clusterUpdate = combineList(clusterUpdateSearch, clusterUpdateClust);
    clusterUpdate.push_back(PARAM_USESEQID);
    clusterUpdate.push_back(PARAM_RECOVER_DELETED);

    //checkSaneEnvironment();
    setDefaults();
}

void Parameters::printUsageMessage(const Command& command,
                                   const int outputFlag){
    const std::vector<MMseqsParameter>& parameters = *command.params;

    std::ostringstream ss;
    ss << "mmseqs " << command.cmd << ":\n";
    ss << (command.longDescription != NULL ? command.longDescription : command.shortDescription) << "\n\n";

    if(command.citations > 0) {
        ss << "Please cite:\n";
        if(command.citations & CITATION_LINCLUST) {
            ss << "Steinegger, M. & Soding, J. Linclust: clustering billions of protein sequences per day on a single server. bioRxiv 104034 (2017)\n\n";
        }
        if(command.citations & CITATION_MMSEQS1) {
            ss << "Hauser, M., Steinegger, M. & Soding, J. MMseqs software suite for fast and deep clustering and searching of large protein sequence sets. Bioinformatics, 32(9), 1323-1330 (2016). \n\n";
        }
        if(command.citations & CITATION_UNICLUST) {
            ss << "Mirdita, M., von den Driesch, L., Galiez, C., Martin M., Soding J. & Steinegger M. Uniclust databases of clustered and deeply annotated protein sequences and alignments. Nucleic Acids Res (2017), D170-D176 (2016).\n\n";
        }
        if(command.citations & CITATION_MMSEQS2) {
            ss << "Steinegger, M. & Soding, J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, doi:10.1038/nbt.3988 (2017)\n\n";
        }
    }
    ss << "© " << command.author << "\n\n";
    ss << "Usage: " << command.usage << (parameters.size() > 0 ? " [options]" : "") << "\n\n";

    struct {
        const char* title;
        int category;
    } categories[] = {
            {"prefilter",MMseqsParameter::COMMAND_PREFILTER},
            {"align",    MMseqsParameter::COMMAND_ALIGN},
            {"clust",    MMseqsParameter::COMMAND_CLUST},
            {"kmermatcher", MMseqsParameter::COMMAND_CLUSTLINEAR},
            {"profile",  MMseqsParameter::COMMAND_PROFILE},
            {"misc",     MMseqsParameter::COMMAND_MISC},
            {"common",   MMseqsParameter::COMMAND_COMMON},
    };

    size_t maxWidth = 0;
    for(size_t i = 0; i < parameters.size(); i++) {
        maxWidth = std::max(strlen(parameters[i].name), maxWidth);
    }
    maxWidth+=2; // space in front of options
    std::map<int, bool> alreadyPrintMap;
    // header
    ss << std::setprecision(3) << std::fixed;
    bool printExpert = (MMseqsParameter::COMMAND_EXPERT & outputFlag) ;

    for(size_t i = 0; i < ARRAY_SIZE(categories); ++i) {
        bool categoryFound = false;
        for (size_t j = 0; j < parameters.size(); j++) {
            const MMseqsParameter &par = parameters[j];
            bool isExpert =  (par.category & MMseqsParameter::COMMAND_EXPERT);
            bool alreadyPrint = alreadyPrintMap[par.uniqid];
            if (par.category & categories[i].category && (printExpert || isExpert == false ) && alreadyPrint == false ) {
                int others = (par.category ^ categories[i].category);
                if(others & outputFlag  )
                    continue;
                categoryFound = true;
                break;
            }
        }
        if(categoryFound){
            std::string title(categories[i].title);
            title.append(" options");
            ss << std::left << std::setw(maxWidth) << title << "\t";
            ss << std::left << std::setw(10) << "default" << "\t";
            ss << "description [value range]" << std::endl;

            // body
            for (size_t j = 0; j < parameters.size(); j++) {
                const MMseqsParameter &par = parameters[j];
                bool isExpert =  (par.category & MMseqsParameter::COMMAND_EXPERT);
                bool alreadyPrint = alreadyPrintMap[par.uniqid];
                if(par.category & categories[i].category && (printExpert || isExpert == false ) && alreadyPrint == false ){
                    ss << std::left << std::setw(maxWidth) << "  " + std::string(par.name) << "\t";
                    ss << std::boolalpha << std::left << std::setw(6);
                    if (par.type == typeid(int)) {
                        ss << std::left << std::setw(10) << *((int *) par.value);
                    } else if (par.type == typeid(float)) {
                        ss << std::left << std::setw(10) << *((float *) par.value);
                    } else if (par.type == typeid(bool)) {
                        bool flag = *((bool *) par.value);
                        std::string flagOutput = (flag) ? "true" : "false";
                        ss << std::left << std::setw(10) << flagOutput;
                    } else if (par.type == typeid(std::string)) {
                        std::string &out = *((std::string *) par.value);
                        ss << std::left << std::setw(10) << out;
                        //            for(std::string::const_iterator j = out.begin(); j != out.end(); ++j) {
                        //                if(*j == '\n') {
                        //                    ss << "\\n";
                        //                } else if (*j == '\t') {
                        //                    ss << "\\t";
                        //                } else {
                        //                    ss << *j;
                        //                }
                        //            }
                    }
                    ss << "\t";
                    ss << std::left << std::setw(60) << par.description << std::endl;
                    alreadyPrintMap[par.uniqid] = true;
                }
            }
            ss << "\n";
        }
    }
    if (printExpert==false) {
        ss << "An extended list of options can be obtained by calling '" << binary_name << " " << command.cmd << " -h'.\n";
    }
    Debug(Debug::INFO) << ss.str();
}

int compileRegex(regex_t * regex, const char * regexText){
    int status = regcomp(regex, regexText, REG_EXTENDED | REG_NEWLINE);
    if (status != 0 ){
        Debug(Debug::ERROR) << "Error in regex " << regexText << "\n";
        EXIT(EXIT_FAILURE);
    }
    return 0;
}

void Parameters::parseParameters(int argc, const char* pargv[],
                                 const Command& command,
                                 size_t requiredParameterCount,
                                 bool printPar,
                                 int parseFlags,
                                 int outputFlags) {
    filenames.clear();
    std::vector<MMseqsParameter>& par = *command.params;
    size_t parametersFound = 0;
    for(int argIdx = 0; argIdx < argc; argIdx++ ){
        // it is a parameter if it starts with - or --
        const bool longParameter = (pargv[argIdx][0] == '-' && pargv[argIdx][1] == '-');
        if (longParameter || (pargv[argIdx][0] == '-')) {
            if ((parseFlags & PARSE_REST) && longParameter && pargv[argIdx][2] == '\0') {
                restArgv = pargv + argIdx + 1;
                restArgc = argc - (argIdx + 1);
                break;
            }
            std::string parameter(pargv[argIdx]);
            bool hasUnrecognizedParameter = true;
            for(size_t parIdx = 0; parIdx < par.size(); parIdx++){
                if (parameter.compare("-h") == 0) {
                    printUsageMessage(command, 0xFFFFFFFF);
                    EXIT(EXIT_SUCCESS);
                }

                if(parameter.compare(par[parIdx].name) == 0) {
                    if (typeid(bool) != par[parIdx].type && argIdx + 1 == argc) {
                        printUsageMessage(command, outputFlags);
                        Debug(Debug::ERROR) << "Missing argument " << par[parIdx].name << "\n";
                        EXIT(EXIT_FAILURE);
                    }

                    if (par[parIdx].wasSet) {
                        printUsageMessage(command, outputFlags);
                        Debug(Debug::ERROR) << "Duplicate parameter " << par[parIdx].name << "\n";
                        EXIT(EXIT_FAILURE);
                    }

                    if (typeid(int) == par[parIdx].type) {
                        regex_t regex;
                        compileRegex(&regex, par[parIdx].regex);
                        int nomatch = regexec(&regex, pargv[argIdx+1], 0, NULL, 0);
                        regfree(&regex);
                        // if no match found or two matches found (we want exactly one match)
                        if (nomatch){
                            printUsageMessage(command, outputFlags);
                            Debug(Debug::ERROR) << "Error in argument " << par[parIdx].name << "\n";
                            EXIT(EXIT_FAILURE);
                        }else{
                            *((int *) par[parIdx].value) = atoi(pargv[argIdx+1]);
                            par[parIdx].wasSet = true;
                        }
                        argIdx++;
                    } else if (typeid(float) == par[parIdx].type) {
                        regex_t regex;
                        compileRegex(&regex, par[parIdx].regex);
                        int nomatch = regexec(&regex, pargv[argIdx+1], 0, NULL, 0);
                        regfree(&regex);
                        if (nomatch){
                            printUsageMessage(command, outputFlags);
                            Debug(Debug::ERROR) << "Error in argument " << par[parIdx].name << "\n";
                            EXIT(EXIT_FAILURE);
                        }else{
                            double input = strtod(pargv[argIdx+1], NULL);
                            *((float *) par[parIdx].value) = static_cast<float>(input);
                            par[parIdx].wasSet = true;
                        }
                        argIdx++;
                    } else if (typeid(std::string) == par[parIdx].type) {
                        std::string val(pargv[argIdx+1]);
                        if(val.length() != 0){
                            std::string * currVal = ((std::string *)par[parIdx].value);
                            currVal->assign( val );
                            par[parIdx].wasSet = true;
                        }
                        argIdx++;
                    } else if (typeid(bool) == par[parIdx].type) {
                        bool * value = (bool *) par[parIdx].value;
                        par[parIdx].wasSet = true;
                        // toggle Value
                        *value = !*value;
                    } else {
                        Debug(Debug::ERROR) << "Wrong parameter type in parseParameters. Please inform the developers\n";
                        EXIT(EXIT_FAILURE);
                    }

                    hasUnrecognizedParameter = false;
                    continue;
                }
            }

            if (hasUnrecognizedParameter) {
                printUsageMessage(command, outputFlags);
                Debug(Debug::ERROR) << "Unrecognized parameter " << parameter << "\n";

                // Suggest some parameter that the user might have meant
                std::vector<MMseqsParameter>::const_iterator index = par.end();
                size_t minDistance = SIZE_MAX;
                for (std::vector<MMseqsParameter>::const_iterator it = par.begin(); it != par.end(); ++it) {
                    size_t distance = DistanceCalculator::levenshteinDistance(parameter, (*it).name);
                    if(distance < minDistance) {
                        minDistance = distance;
                        index = it;
                    }
                }

                if(index != par.end()) {
                    Debug(Debug::ERROR) << "Did you mean \"" << (*index).name << "\"?\n";
                }

                EXIT(EXIT_FAILURE);
            }

            parametersFound++;
        } else { // it is a filename if its not a parameter
            filenames.emplace_back(pargv[argIdx]);
        }
    }

    if (MMseqsMPI::isMaster()) {
        Debug::setDebugLevel(verbosity);
    }

#ifdef OPENMP
    omp_set_num_threads(threads);
#endif

    const size_t MAX_DB_PARAMETER = 6;

    if (requiredParameterCount > MAX_DB_PARAMETER) {
        Debug(Debug::ERROR) << "Use argv if you need more than " << MAX_DB_PARAMETER << " db parameters" << "\n";
        EXIT(EXIT_FAILURE);
    }

    if (filenames.size() < requiredParameterCount){
        printUsageMessage(command, outputFlags);
        Debug(Debug::ERROR) << requiredParameterCount << " Database paths are required" << "\n";
        EXIT(EXIT_FAILURE);
    }

    switch (std::min(filenames.size(), MAX_DB_PARAMETER)) {
        case 6:
            db6 = filenames[5];
            db6Index = db6;
            db6Index.append(".index");
            // FALLTHROUGH
        case 5:
            db5 = filenames[4];
            db5Index = db5;
            db5Index.append(".index");
            // FALLTHROUGH
        case 4:
            db4 = filenames[3];
            db4Index = db4;
            db4Index.append(".index");
            // FALLTHROUGH
        case 3:
            db3 = filenames[2];
            db3Index = db3;
            db3Index.append(".index");
            // FALLTHROUGH
        case 2:
            db2 = filenames[1];
            db2Index = db2;
            db2Index.append(".index");
            // FALLTHROUGH
        case 1:
            db1 = filenames[0];
            db1Index = db1;
            db1Index.append(".index");
            break;
        default:
            // Do not abort execution if we expect a variable amount of parameters
            if (parseFlags & PARSE_VARIADIC)
                break;
            // FALLTHROUGH
        case 0:
            printUsageMessage(command, outputFlags);
            Debug(Debug::ERROR) << "Unrecognized parameters!" << "\n";
            printParameters(argc, pargv, par);
            EXIT(EXIT_FAILURE);
    }
    if(printPar == true) {
        printParameters(argc, pargv, par);
    }
}

void Parameters::printParameters(int argc, const char* pargv[],
                                 const std::vector<MMseqsParameter> &par){
    if (Debug::debugLevel < 3) {
        return;
    }

    Debug(Debug::INFO) << "Program call:\n";
    for (int i = 0; i < argc; i++)
        Debug(Debug::INFO) << pargv[i] << " ";
    Debug(Debug::INFO) << "\n\n";

    size_t maxWidth = 0;
    for(size_t i = 0; i < par.size(); i++) {
        maxWidth = std::max(strlen(par[i].display), maxWidth);
    }

    std::stringstream ss;
    ss << std::boolalpha;

    ss << std::setw(maxWidth) << std::left  << "MMseqs Version:" << "\t" << version << "\n";


    for (size_t i = 0; i < par.size(); i++) {
        ss << std::setw(maxWidth) << std::left << par[i].display << "\t";
        if(typeid(int) == par[i].type ){
            ss << *((int *)par[i].value);
        } else if(typeid(float) == par[i].type ){
            ss << *((float *)par[i].value);
        }else if(typeid(std::string) == par[i].type ){
            ss << *((std::string *) par[i].value);
        }else if (typeid(bool) == par[i].type){
            ss << *((bool *)par[i].value);
        }
        ss << "\n";
    }

    Debug(Debug::INFO) << ss.str() << "\n";
}

void Parameters::setDefaults() {
    restArgv = NULL;
    restArgc = 0;

    scoringMatrixFile = "blosum62.out";

    kmerSize =  0;
    kmerScore = INT_MAX;
    alphabetSize = 21;
    maxSeqLen = MAX_SEQ_LEN; // 2^16
    maxResListLen = 300;
    sensitivity = 4;
    split = AUTO_SPLIT_DETECTION;
    splitMode = DETECT_BEST_DB_SPLIT;
    splitMemoryLimit = 0;
    splitAA = false;

    // search workflow
    numIterations = 1;
    startSens = 4;
    sensSteps = 1;

    threads = 1;
#ifdef OPENMP
    #ifdef _SC_NPROCESSORS_ONLN
    threads = sysconf(_SC_NPROCESSORS_ONLN);
#endif
    if(threads  <= 1){
        threads = Util::omp_thread_count();
    }
#endif
    compBiasCorrection = 1;
    diagonalScoring = 1;
    maskMode = 1;
    minDiagScoreThr = 15;
    spacedKmer = true;
    includeIdentity = false;
    alignmentMode = ALIGNMENT_MODE_FAST_AUTO;
    evalThr = 0.001;
    covThr = 0.0;
    covMode = COV_MODE_BIDIRECTIONAL;
    maxRejected = INT_MAX;
    maxAccept   = INT_MAX;
    seqIdThr = 0.0;
    addBacktrace = false;
    realign = false;
    clusteringMode = SET_COVER;
    cascaded = true;
    clusterSteps = 3;
    resListOffset = 0;
    noPreload = false;
    earlyExit = false;

    // affinity clustering
    maxIteration=1000;
    similarityScoreType=APC_SEQID;

    // workflow
    const char *runnerEnv = getenv("RUNNER");
    if (runnerEnv != NULL) {
        runner = runnerEnv;
    } else {
        runner = "";
    }

    // Clustering workflow
    removeTmpFiles = false;

    // convertprofiledb
    profileMode = PROFILE_MODE_HMM;

    // createdb
    splitSeqByLen = true;


    // format alignment
    formatAlignmentMode = FORMAT_ALIGNMENT_BLAST_TAB;
    dbOut = false;

    // rescore diagonal
    globalAlignment = false;

    // result2msa
    allowDeletion = false;
    addInternalId = false;
    compressMSA = false;
    summarizeHeader = false;
    summaryPrefix = "cl";
    compressMSA = false;
    omitConsensus = false;
    skipQuery = false;

    // msa2profile
    matchMode = 0;
    matchRatio = 0.5;

    // result2profile
    evalProfile = evalThr;
    filterMsa = 1;
    filterMaxSeqId = 0.9;
    qid = 0.0;           // default for minimum sequence identity with query
    qsc = -20.0f;        // default for minimum score per column with query
    cov = 0.0;           // default for minimum coverage threshold
    Ndiff = 1000;        // pick Ndiff most different sequences from alignment
    wg = false;
    pca = 1.0;
    pcb = 1.5;

    // sequence2profile
    neff = 1.0;
    tau = 0.9;

    // logging
    verbosity = Debug::INFO;

    //extractorfs
    orfMinLength = 1;
    orfMaxLength = INT_MAX;
    orfMaxGaps = INT_MAX;
    orfStartState = 2;
    orfEndState = 2;
    orfLongest = false;
    orfExtendMin = false;
    forwardFrames = "1,2,3";
    reverseFrames = "1,2,3";

    // createdb
    identifierOffset = 0;

    // convert2fasta
    useHeaderFile = false;

    // result2flat
    useHeader = false;

    // translate nucleotide
    addOrfStop = false;
    translationTable = 1;

    // createseqfiledb
    minSequences = 1;
    maxSequences = INT_MAX;
    hhFormat = false;

    // rescorediagonal
    rescoreMode = Parameters::RESCORE_MODE_HAMMING;
    filterHits = false;
    // filterDb
    filterColumn = 1;
    filterColumnRegex = "^.*$";
    positiveFilter = true;
    filteringFile = "";
    trimToOneColumn = false;
    extractLines = 0;
    compOperator = "";
    compValue = 0;
    sortEntries = 0;

    // concatdbs
    preserveKeysB = false;

    // diff
    useSequenceId = false;

    // prefixid
    prefix = "";

    // mergedbs
    mergePrefixes = "";

    // summarizetabs
    overlap = 0.0f;
    msaType = 0;

    // summarize header
    headerType = Parameters::HEADER_TYPE_UNICLUST;

    // extractalignedregion
    extractMode = EXTRACT_TARGET;

    // convertkb
    kbColumns = "";

    // linearcluster
    kmersPerSequence = 20;
    includeOnlyExtendable = false;
    hashShift = 5;

    // result2stats
    stat = "";

    // createtsv
    firstSeqRepr = false;

    // lca
    lcaRanks = "";
    // bin for all unclassified sequences
    // https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=12908
    // other sequences (plasmids, etc)
    // https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=28384
    blacklist = "12908,28384";

    // taxonomy
    lcaMode = 2;
}

std::vector<MMseqsParameter> Parameters::combineList(const std::vector<MMseqsParameter> &par1,
                                                     const std::vector<MMseqsParameter> &par2) {
    std::vector<MMseqsParameter> retVec;
    std::vector< std::vector<MMseqsParameter>> tmp;
    tmp.push_back(par1);
    tmp.push_back(par2);
    for(size_t z = 0; z < tmp.size(); z++) {
        std::vector<MMseqsParameter> currPar = tmp[z];
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

size_t Parameters::hashParameter(const std::vector<std::string> &filenames, const std::vector<MMseqsParameter> &par){
    std::string hashString;
    hashString.reserve(1024);
    for (size_t i = 0; i < filenames.size(); ++i){
        hashString.append(filenames[i]);
        hashString.append(" ");
    }
    hashString.append(createParameterString(par));
    hashString.append(version);
    for (int i = 0; i < restArgc; ++i) {
        hashString.append(restArgv[i]);
    }
    return Util::hash(hashString.c_str(), hashString.size());
}

std::string Parameters::createParameterString(const std::vector<MMseqsParameter> &par) {
    std::ostringstream ss;
    for (size_t i = 0; i < par.size(); ++i) {
        // Never pass the MPI parameters along, they are passed by the environment
        if (par[i].uniqid == PARAM_RUNNER_ID) {
            continue;
        }

        if (typeid(int) == par[i].type){
            ss << par[i].name << " ";
            ss << *((int *)par[i].value) << " ";
        } else if (typeid(float) == par[i].type){
            ss << par[i].name << " ";
            ss << *((float *)par[i].value) << " ";
        } else if (typeid(std::string) == par[i].type){
            if (*((std::string *) par[i].value) != "") {
                ss << par[i].name << " ";
                ss << *((std::string *) par[i].value) << " ";
            }
        } else if (typeid(bool) == par[i].type){
            bool val = *((bool *)(par[i].value));
            if (val == true){
                ss << par[i].name << " ";
            }
        } else {
            Debug(Debug::ERROR) << "Wrong parameter type. Please inform the developers!\n";
            EXIT(EXIT_FAILURE);
        }
    }

    return ss.str();
}

std::vector<MMseqsParameter> Parameters::removeParameter(const std::vector<MMseqsParameter> &par, const MMseqsParameter &x){
    std::vector<MMseqsParameter> newParamList;
    for (std::vector<MMseqsParameter>::const_iterator i = par.begin();i!=par.end();i++) {
        if (i->name != x.name)
            newParamList.push_back(*i);
    }
    return newParamList;
}

void Parameters::overrideParameterDescription(Command &command, const int uid,
                                              const char *description, const char *regex, const int category) {
    for (std::vector<MMseqsParameter>::iterator i = command.params->begin(); i != command.params->end(); i++) {
        if (i->uniqid == uid) {
            if(description != NULL){
                i->description = description;
            }
            if (regex != NULL) {
                i->regex = regex;
            }
            if (category != 0) {
                i->category = category;
            }
            break;
        }
    }
}
