#include "Parameters.h"
#include "Util.h"
#include "DistanceCalculator.h"
#include "Debug.h"
#include "CommandCaller.h"
#include "ByteParser.h"
#include "FileUtil.h"

#include <map>
#include <iomanip>
#include <regex.h>
#include <unistd.h>

#include "blosum62.out.h"
#include "PAM30.out.h"
#include "VTML80.out.h"
#include "VTML40.out.h"
#include "nucleotide.out.h"
#include "base64/base64.h"

#ifdef __CYGWIN__
#include <sys/cygwin.h>
#endif

#ifdef OPENMP
#include <omp.h>
#endif

Parameters* Parameters::instance = NULL;

extern const char* binary_name;
extern const char* version;

Parameters::Parameters():
        scoringMatrixFile(NuclAA<std::string>("INVALID", "INVALID")),
        seedScoringMatrixFile(NuclAA<std::string>("INVALID", "INVALID")),
        alphabetSize(NuclAA<int>(INT_MAX,INT_MAX)),
        PARAM_S(PARAM_S_ID, "-s", "Sensitivity", "Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive", typeid(float), (void *) &sensitivity, "^[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_PREFILTER),
        PARAM_K(PARAM_K_ID, "-k", "k-mer length", "k-mer length (0: automatically set to optimum)", typeid(int), (void *) &kmerSize, "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_CLUSTLINEAR | MMseqsParameter::COMMAND_EXPERT),
        PARAM_TARGET_SEARCH_MODE(PARAM_TARGET_SEARCH_MODE_ID, "--target-search-mode", "Target search mode", "target search mode (0: regular k-mer, 1: similar k-mer)", typeid(int), (void *) &targetSearchMode, "^[0-1]{1}$", MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_CLUSTLINEAR | MMseqsParameter::COMMAND_EXPERT),
        PARAM_THREADS(PARAM_THREADS_ID, "--threads", "Threads", "Number of CPU-cores used (all by default)", typeid(int), (void *) &threads, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_COMMON),
        PARAM_COMPRESSED(PARAM_COMPRESSED_ID, "--compressed", "Compressed", "Write compressed output", typeid(int), (void *) &compressed, "^[0-1]{1}$", MMseqsParameter::COMMAND_COMMON),
        PARAM_ALPH_SIZE(PARAM_ALPH_SIZE_ID, "--alph-size", "Alphabet size", "Alphabet size (range 2-21)", typeid(MultiParam<NuclAA<int>>), (void *) &alphabetSize, "", MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_CLUSTLINEAR | MMseqsParameter::COMMAND_EXPERT),
        PARAM_MAX_SEQ_LEN(PARAM_MAX_SEQ_LEN_ID, "--max-seq-len", "Max sequence length", "Maximum sequence length", typeid(size_t), (void *) &maxSeqLen, "^[0-9]{1}[0-9]*", MMseqsParameter::COMMAND_COMMON | MMseqsParameter::COMMAND_EXPERT),
        PARAM_DIAGONAL_SCORING(PARAM_DIAGONAL_SCORING_ID, "--diag-score", "Diagonal scoring", "Use ungapped diagonal scoring during prefilter", typeid(bool), (void *) &diagonalScoring, "", MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_EXPERT),
        PARAM_EXACT_KMER_MATCHING(PARAM_EXACT_KMER_MATCHING_ID, "--exact-kmer-matching", "Exact k-mer matching", "Extract only exact k-mers for matching (range 0-1)", typeid(int), (void *) &exactKmerMatching, "^[0-1]{1}$", MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_EXPERT),
        PARAM_MASK_RESIDUES(PARAM_MASK_RESIDUES_ID, "--mask", "Mask residues", "Mask sequences in k-mer stage: 0: w/o low complexity masking, 1: with low complexity masking", typeid(int), (void *) &maskMode, "^[0-1]{1}", MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_EXPERT),
        PARAM_MASK_PROBABILTY(PARAM_MASK_PROBABILTY_ID, "--mask-prob", "Mask residues probability", "Mask sequences is probablity is above threshold", typeid(float), (void *) &maskProb, "^0(\\.[0-9]+)?|^1(\\.0+)?$", MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_EXPERT),
        PARAM_MASK_LOWER_CASE(PARAM_MASK_LOWER_CASE_ID, "--mask-lower-case", "Mask lower case residues", "Lowercase letters will be excluded from k-mer search 0: include region, 1: exclude region", typeid(int), (void *) &maskLowerCaseMode, "^[0-1]{1}", MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_EXPERT),
        PARAM_MIN_DIAG_SCORE(PARAM_MIN_DIAG_SCORE_ID, "--min-ungapped-score", "Minimum diagonal score", "Accept only matches with ungapped alignment score above threshold", typeid(int), (void *) &minDiagScoreThr, "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_EXPERT),
        PARAM_K_SCORE(PARAM_K_SCORE_ID, "--k-score", "k-score", "k-mer threshold for generating similar k-mer lists", typeid(MultiParam<SeqProf<int>>), (void *) &kmerScore, "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_EXPERT),
        PARAM_MAX_SEQS(PARAM_MAX_SEQS_ID, "--max-seqs", "Max results per query", "Maximum results per query sequence allowed to pass the prefilter (affects sensitivity)", typeid(size_t), (void *) &maxResListLen, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PREFILTER),
        PARAM_SPLIT(PARAM_SPLIT_ID, "--split", "Split database", "Split input into N equally distributed chunks. 0: set the best split automatically", typeid(int), (void *) &split, "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_EXPERT),
        PARAM_SPLIT_MODE(PARAM_SPLIT_MODE_ID, "--split-mode", "Split mode", "0: split target db; 1: split query db; 2: auto, depending on main memory", typeid(int), (void *) &splitMode, "^[0-2]{1}$", MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_EXPERT),
        PARAM_SPLIT_MEMORY_LIMIT(PARAM_SPLIT_MEMORY_LIMIT_ID, "--split-memory-limit", "Split memory limit", "Set max memory per split. E.g. 800B, 5K, 10M, 1G. Default (0) to all available system memory", typeid(ByteParser), (void *) &splitMemoryLimit, "^(0|[1-9]{1}[0-9]*(B|K|M|G|T)?)$", MMseqsParameter::COMMAND_COMMON | MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_EXPERT),
        PARAM_DISK_SPACE_LIMIT(PARAM_DISK_SPACE_LIMIT_ID, "--disk-space-limit", "Disk space limit", "Set max disk space to use for reverse profile searches. E.g. 800B, 5K, 10M, 1G. Default (0) to all available disk space in the temp folder", typeid(ByteParser), (void *) &diskSpaceLimit, "^(0|[1-9]{1}[0-9]*(B|K|M|G|T)?)$", MMseqsParameter::COMMAND_COMMON | MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_EXPERT),
        PARAM_SPLIT_AMINOACID(PARAM_SPLIT_AMINOACID_ID, "--split-aa", "Split by amino acid", "Try to find the best split boundaries by entry lengths", typeid(bool), (void *) &splitAA, "$", MMseqsParameter::COMMAND_EXPERT),
        PARAM_SUB_MAT(PARAM_SUB_MAT_ID, "--sub-mat", "Substitution matrix", "Substitution matrix file", typeid(MultiParam<NuclAA<std::string>>), (void *) &scoringMatrixFile, "", MMseqsParameter::COMMAND_COMMON | MMseqsParameter::COMMAND_EXPERT),
        PARAM_SEED_SUB_MAT(PARAM_SEED_SUB_MAT_ID, "--seed-sub-mat", "Seed substitution matrix", "Substitution matrix file for k-mer generation", typeid(MultiParam<NuclAA<std::string>>), (void *) &seedScoringMatrixFile, "", MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_EXPERT),
        PARAM_NO_COMP_BIAS_CORR(PARAM_NO_COMP_BIAS_CORR_ID, "--comp-bias-corr", "Compositional bias", "Correct for locally biased amino acid composition (range 0-1)", typeid(int), (void *) &compBiasCorrection, "^[0-1]{1}$", MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PROFILE | MMseqsParameter::COMMAND_EXPERT),
        PARAM_NO_COMP_BIAS_CORR_SCALE(PARAM_NO_COMP_BIAS_CORR_SCALE_ID, "--comp-bias-corr-scale", "Compositional bias", "Correct for locally biased amino acid composition (range 0-1)", typeid(float), (void *) &compBiasCorrectionScale,  "^0(\\.[0-9]+)?|^1(\\.0+)?$", MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PROFILE | MMseqsParameter::COMMAND_EXPERT),

        PARAM_SPACED_KMER_MODE(PARAM_SPACED_KMER_MODE_ID, "--spaced-kmer-mode", "Spaced k-mers", "0: use consecutive positions in k-mers; 1: use spaced k-mers", typeid(int), (void *) &spacedKmer, "^[0-1]{1}", MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_EXPERT),
        PARAM_REMOVE_TMP_FILES(PARAM_REMOVE_TMP_FILES_ID, "--remove-tmp-files", "Remove temporary files", "Delete temporary files", typeid(bool), (void *) &removeTmpFiles, "", MMseqsParameter::COMMAND_COMMON | MMseqsParameter::COMMAND_EXPERT),
        PARAM_INCLUDE_IDENTITY(PARAM_INCLUDE_IDENTITY_ID, "--add-self-matches", "Include identical seq. id.", "Artificially add entries of queries with themselves (for clustering)", typeid(bool), (void *) &includeIdentity, "", MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
        PARAM_PRELOAD_MODE(PARAM_PRELOAD_MODE_ID, "--db-load-mode", "Preload mode", "Database preload mode 0: auto, 1: fread, 2: mmap, 3: mmap+touch", typeid(int), (void *) &preloadMode, "[0-3]{1}", MMseqsParameter::COMMAND_COMMON | MMseqsParameter::COMMAND_EXPERT),
        PARAM_SPACED_KMER_PATTERN(PARAM_SPACED_KMER_PATTERN_ID, "--spaced-kmer-pattern", "Spaced k-mer pattern", "User-specified spaced k-mer pattern", typeid(std::string), (void *) &spacedKmerPattern, "^1[01]*1$", MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_EXPERT),
        PARAM_LOCAL_TMP(PARAM_LOCAL_TMP_ID, "--local-tmp", "Local temporary path", "Path where some of the temporary files will be created", typeid(std::string), (void *) &localTmp, "", MMseqsParameter::COMMAND_PREFILTER | MMseqsParameter::COMMAND_EXPERT),
        // alignment
        PARAM_ALIGNMENT_MODE(PARAM_ALIGNMENT_MODE_ID, "--alignment-mode", "Alignment mode", "How to compute the alignment:\n0: automatic\n1: only score and end_pos\n2: also start_pos and cov\n3: also seq.id\n4: only ungapped alignment", typeid(int), (void *) &alignmentMode, "^[0-5]{1}$", MMseqsParameter::COMMAND_ALIGN),
        PARAM_ALIGNMENT_OUTPUT_MODE(PARAM_ALIGNMENT_OUTPUT_MODE_ID, "--alignment-output-mode", "Alignment mode", "How to compute the alignment:\n0: automatic\n1: only score and end_pos\n2: also start_pos and cov\n3: also seq.id\n4: only ungapped alignment\n5: score only (output) cluster format", typeid(int), (void *) &alignmentOutputMode, "^[0-1]{1}$", MMseqsParameter::COMMAND_ALIGN),
        PARAM_E(PARAM_E_ID, "-e", "E-value threshold", "List matches below this E-value (range 0.0-inf)", typeid(double), (void *) &evalThr, "^([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)|[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_ALIGN),
        PARAM_C(PARAM_C_ID, "-c", "Coverage threshold", "List matches above this fraction of aligned (covered) residues (see --cov-mode)", typeid(float), (void *) &covThr, "^0(\\.[0-9]+)?|^1(\\.0+)?$", MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_CLUSTLINEAR),
        PARAM_COV_MODE(PARAM_COV_MODE_ID, "--cov-mode", "Coverage mode", "0: coverage of query and target\n1: coverage of target\n2: coverage of query\n3: target seq. length has to be at least x% of query length\n4: query seq. length has to be at least x% of target length\n5: short seq. needs to be at least x% of the other seq. length", typeid(int), (void *) &covMode, "^[0-5]{1}$", MMseqsParameter::COMMAND_ALIGN),
        PARAM_SEQ_ID_MODE(PARAM_SEQ_ID_MODE_ID, "--seq-id-mode", "Seq. id. mode", "0: alignment length 1: shorter, 2: longer sequence", typeid(int), (void *) &seqIdMode, "^[0-2]{1}$", MMseqsParameter::COMMAND_ALIGN),
        PARAM_MAX_REJECTED(PARAM_MAX_REJECTED_ID, "--max-rejected", "Max reject", "Maximum rejected alignments before alignment calculation for a query is stopped", typeid(int), (void *) &maxRejected, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_ALIGN),
        PARAM_MAX_ACCEPT(PARAM_MAX_ACCEPT_ID, "--max-accept", "Max accept", "Maximum accepted alignments before alignment calculation for a query is stopped", typeid(int), (void *) &maxAccept, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_ALIGN),
        PARAM_ADD_BACKTRACE(PARAM_ADD_BACKTRACE_ID, "-a", "Add backtrace", "Add backtrace string (convert to alignments with mmseqs convertalis module)", typeid(bool), (void *) &addBacktrace, "", MMseqsParameter::COMMAND_ALIGN),
        PARAM_REALIGN(PARAM_REALIGN_ID, "--realign", "Realign hits", "Compute more conservative, shorter alignments (scores and E-values not changed)", typeid(bool), (void *) &realign, "", MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
        PARAM_MIN_SEQ_ID(PARAM_MIN_SEQ_ID_ID, "--min-seq-id", "Seq. id. threshold", "List matches above this sequence identity (for clustering) (range 0.0-1.0)", typeid(float), (void *) &seqIdThr, "^0(\\.[0-9]+)?|1(\\.0+)?$", MMseqsParameter::COMMAND_ALIGN),
        PARAM_MIN_ALN_LEN(PARAM_MIN_ALN_LEN_ID, "--min-aln-len", "Min alignment length", "Minimum alignment length (range 0-INT_MAX)", typeid(int), (void *) &alnLenThr, "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_ALIGN),
        PARAM_SCORE_BIAS(PARAM_SCORE_BIAS_ID, "--score-bias", "Score bias", "Score bias when computing SW alignment (in bits)", typeid(float), (void *) &scoreBias, "^-?[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
        PARAM_REALIGN_SCORE_BIAS(PARAM_REALIGN_SCORE_BIAS_ID, "--realign-score-bias", "Realign score bias", "Additional bias when computing realignment", typeid(float), (void *) &realignScoreBias, "^-?[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
        PARAM_REALIGN_MAX_SEQS(PARAM_REALIGN_MAX_SEQS_ID, "--realign-max-seqs", "Realign max seqs", "Maximum number of results to return in realignment", typeid(int), (void *) &realignMaxSeqs, "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
        PARAM_CORR_SCORE_WEIGHT(PARAM_CORR_SCORE_WEIGHT_ID, "--corr-score-weight", "Correlation score weight", "Weight of backtrace correlation score that is added to the alignment score", typeid(float), (void *) &correlationScoreWeight, "^-?[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
        PARAM_ALT_ALIGNMENT(PARAM_ALT_ALIGNMENT_ID, "--alt-ali", "Alternative alignments", "Show up to this many alternative alignments", typeid(int), (void *) &altAlignment, "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_ALIGN),
        PARAM_GAP_OPEN(PARAM_GAP_OPEN_ID, "--gap-open", "Gap open cost", "Gap open cost", typeid(MultiParam<NuclAA<int>>), (void *) &gapOpen, "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
        PARAM_GAP_EXTEND(PARAM_GAP_EXTEND_ID, "--gap-extend", "Gap extension cost", "Gap extension cost", typeid(MultiParam<NuclAA<int>>), (void *) &gapExtend, "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
#ifdef GAP_POS_SCORING
        PARAM_GAP_PSEUDOCOUNT(PARAM_GAP_PSEUDOCOUNT_ID, "--gap-pc", "Gap pseudo count", "Pseudo count for calculating position-specific gap opening penalties", typeid(int), &gapPseudoCount, "^[0-9]+$", MMseqsParameter::COMMAND_ALIGN|MMseqsParameter::COMMAND_EXPERT),
#endif
        PARAM_ZDROP(PARAM_ZDROP_ID, "--zdrop", "Zdrop", "Maximal allowed difference between score values before alignment is truncated  (nucleotide alignment only)", typeid(int), (void*) &zdrop, "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
        // clustering
        PARAM_CLUSTER_MODE(PARAM_CLUSTER_MODE_ID, "--cluster-mode", "Cluster mode", "0: Set-Cover (greedy)\n1: Connected component (BLASTclust)\n2,3: Greedy clustering by sequence length (CDHIT)", typeid(int), (void *) &clusteringMode, "[0-3]{1}$", MMseqsParameter::COMMAND_CLUST),
        PARAM_CLUSTER_STEPS(PARAM_CLUSTER_STEPS_ID, "--cluster-steps", "Cascaded clustering steps", "Cascaded clustering steps from 1 to -s", typeid(int), (void *) &clusterSteps, "^[1-9]{1}$", MMseqsParameter::COMMAND_CLUST | MMseqsParameter::COMMAND_EXPERT),
        PARAM_CASCADED(PARAM_CASCADED_ID, "--single-step-clustering", "Single step clustering", "Switch from cascaded to simple clustering workflow", typeid(bool), (void *) &singleStepClustering, "", MMseqsParameter::COMMAND_CLUST),
        PARAM_CLUSTER_REASSIGN(PARAM_CLUSTER_REASSIGN_ID, "--cluster-reassign", "Cluster reassign", "Cascaded clustering can cluster sequence that do not fulfill the clustering criteria.\nCluster reassignment corrects these errors", typeid(bool), (void *) &clusterReassignment, "", MMseqsParameter::COMMAND_CLUST),
        // affinity clustering
        PARAM_MAXITERATIONS(PARAM_MAXITERATIONS_ID, "--max-iterations", "Max connected component depth", "Maximum depth of breadth first search in connected component clustering", typeid(int), (void *) &maxIteration, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_CLUST | MMseqsParameter::COMMAND_EXPERT),
        PARAM_SIMILARITYSCORE(PARAM_SIMILARITYSCORE_ID, "--similarity-type", "Similarity type", "Type of score used for clustering. 1: alignment score 2: sequence identity", typeid(int), (void *) &similarityScoreType, "^[1-2]{1}$", MMseqsParameter::COMMAND_CLUST | MMseqsParameter::COMMAND_EXPERT),
        // logging
        PARAM_V(PARAM_V_ID, "-v", "Verbosity", "Verbosity level: 0: quiet, 1: +errors, 2: +warnings, 3: +info", typeid(int), (void *) &verbosity, "^[0-3]{1}$", MMseqsParameter::COMMAND_COMMON),
        // convertalignments
        PARAM_FORMAT_MODE(PARAM_FORMAT_MODE_ID, "--format-mode", "Alignment format", "Output format:\n0: BLAST-TAB\n1: SAM\n2: BLAST-TAB + query/db length\n3: Pretty HTML\n4: BLAST-TAB + column headers\nBLAST-TAB (0) and BLAST-TAB + column headers (4) support custom output formats (--format-output)", typeid(int), (void *) &formatAlignmentMode, "^[0-4]{1}$"),
        PARAM_FORMAT_OUTPUT(PARAM_FORMAT_OUTPUT_ID, "--format-output", "Format alignment output", "Choose comma separated list of output columns from: query,target,evalue,gapopen,pident,fident,nident,qstart,qend,qlen\ntstart,tend,tlen,alnlen,raw,bits,cigar,qseq,tseq,qheader,theader,qaln,taln,qframe,tframe,mismatch,qcov,tcov\nqset,qsetid,tset,tsetid,taxid,taxname,taxlineage,qorfstart,qorfend,torfstart,torfend", typeid(std::string), (void *) &outfmt, ""),
        PARAM_DB_OUTPUT(PARAM_DB_OUTPUT_ID, "--db-output", "Database output", "Return a result DB instead of a text file", typeid(bool), (void *) &dbOut, "", MMseqsParameter::COMMAND_EXPERT),
        // --include-only-extendablediagonal
        PARAM_RESCORE_MODE(PARAM_RESCORE_MODE_ID, "--rescore-mode", "Rescore mode", "Rescore diagonals with:\n0: Hamming distance\n1: local alignment (score only)\n2: local alignment\n3: global alignment\n4: longest alignment fulfilling window quality criterion", typeid(int), (void *) &rescoreMode, "^[0-4]{1}$"),
        PARAM_WRAPPED_SCORING(PARAM_WRAPPED_SCORING_ID, "--wrapped-scoring", "Allow wrapped scoring", "Double the (nucleotide) query sequence during the scoring process to allow wrapped diagonal scoring around end and start", typeid(bool), (void *) &wrappedScoring, "", MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
        PARAM_FILTER_HITS(PARAM_FILTER_HITS_ID, "--filter-hits", "Remove hits by seq. id. and coverage", "Filter hits by seq.id. and coverage", typeid(bool), (void *) &filterHits, "", MMseqsParameter::COMMAND_EXPERT),
        PARAM_SORT_RESULTS(PARAM_SORT_RESULTS_ID, "--sort-results", "Sort results", "Sort results: 0: no sorting, 1: sort by E-value (Alignment) or seq.id. (Hamming)", typeid(int), (void *) &sortResults, "^[0-1]{1}$", MMseqsParameter::COMMAND_EXPERT),
        // result2msa
        PARAM_MSA_FORMAT_MODE(PARAM_MSA_FORMAT_MODE_ID, "--msa-format-mode", "MSA format mode", "Format MSA as: 0: binary cA3M DB\n1: binary ca3m w. consensus DB\n2: aligned FASTA DB\n3: aligned FASTA w. header summary\n4: STOCKHOLM flat file\n5: A3M format\n6: A3M format w. alignment info", typeid(int), (void *) &msaFormatMode, "^[0-6]{1}$"),
        PARAM_ALLOW_DELETION(PARAM_ALLOW_DELETION_ID, "--allow-deletion", "Allow deletions", "Allow deletions in a MSA", typeid(bool), (void *) &allowDeletion, ""),
        PARAM_SUMMARY_PREFIX(PARAM_SUMMARY_PREFIX_ID, "--summary-prefix", "Summary prefix", "Set the cluster summary prefix", typeid(std::string), (void *) &summaryPrefix, "", MMseqsParameter::COMMAND_EXPERT),
        PARAM_SKIP_QUERY(PARAM_SKIP_QUERY_ID, "--skip-query", "Skip query", "Skip the query sequence", typeid(bool), (void *) &skipQuery, "", MMseqsParameter::COMMAND_EXPERT),
        // convertmsa
        PARAM_IDENTIFIER_FIELD(PARAM_IDENTIFIER_FIELD_ID, "--identifier-field", "Identifier field", "Field from STOCKHOLM comments for choosing the MSA identifier: 0: ID, 1: AC. If the respective comment does not exist, the name of the first sequence will become the identifier", typeid(int), (void *) &identifierField, "^[0-1]{1}$", MMseqsParameter::COMMAND_COMMON),
        // msa2profile
        PARAM_MATCH_MODE(PARAM_MATCH_MODE_ID, "--match-mode", "Match mode", "0: Columns that have a residue in the first sequence are kept, 1: columns that have a residue in --match-ratio of all sequences are kept", typeid(int), (void *) &matchMode, "^(0|1)$", MMseqsParameter::COMMAND_PROFILE),
        PARAM_MATCH_RATIO(PARAM_MATCH_RATIO_ID, "--match-ratio", "Match ratio", "Columns that have a residue in this ratio of all sequences are kept", typeid(float), (void *) &matchRatio, "^0(\\.[0-9]+)?|1(\\.0+)?$", MMseqsParameter::COMMAND_PROFILE),
        // result2profile
        PARAM_MASK_PROFILE(PARAM_MASK_PROFILE_ID, "--mask-profile", "Mask profile", "Mask query sequence of profile using tantan [0,1]", typeid(int), (void *) &maskProfile, "^[0-1]{1}$", MMseqsParameter::COMMAND_PROFILE | MMseqsParameter::COMMAND_EXPERT),
        PARAM_E_PROFILE(PARAM_E_PROFILE_ID, "--e-profile", "Profile E-value threshold", "Include sequences matches with < E-value thr. into the profile (>=0.0)", typeid(double), (void *) &evalProfile, "^([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)|([0-9]*(\\.[0-9]+)?)$", MMseqsParameter::COMMAND_PROFILE),
        PARAM_FILTER_MSA(PARAM_FILTER_MSA_ID, "--filter-msa", "Filter MSA", "Filter msa: 0: do not filter, 1: filter", typeid(int), (void *) &filterMsa, "^[0-1]{1}$", MMseqsParameter::COMMAND_PROFILE | MMseqsParameter::COMMAND_EXPERT),
        PARAM_FILTER_MAX_SEQ_ID(PARAM_FILTER_MAX_SEQ_ID_ID, "--max-seq-id", "Maximum seq. id. threshold", "Reduce redundancy of output MSA using max. pairwise sequence identity [0.0,1.0]", typeid(float), (void *) &filterMaxSeqId, "^0(\\.[0-9]+)?|1(\\.0+)?$", MMseqsParameter::COMMAND_PROFILE | MMseqsParameter::COMMAND_EXPERT),
        PARAM_FILTER_QSC(PARAM_FILTER_QSC_ID, "--qsc", "Minimum score per column", "Reduce diversity of output MSAs using min. score per aligned residue with query sequences [-50.0,100.0]", typeid(float), (void *) &qsc, "^\\-*[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_PROFILE | MMseqsParameter::COMMAND_EXPERT),
        PARAM_FILTER_QID(PARAM_FILTER_QID_ID, "--qid", "Minimum seq. id.", "Reduce diversity of output MSAs using min.seq. identity with query sequences [0.0,1.0]\nAlternatively, can be a list of multiple thresholds:\nE.g.: 0.15,0.30,0.50 to defines filter buckets of ]0.15-0.30] and ]0.30-0.50]", typeid(std::string), (void *) &qid, "", MMseqsParameter::COMMAND_PROFILE | MMseqsParameter::COMMAND_EXPERT),
        PARAM_FILTER_COV(PARAM_FILTER_COV_ID, "--cov", "Minimum coverage", "Filter output MSAs using min. fraction of query residues covered by matched sequences [0.0,1.0]", typeid(float), (void *) &covMSAThr, "^0(\\.[0-9]+)?|1(\\.0+)?$", MMseqsParameter::COMMAND_PROFILE | MMseqsParameter::COMMAND_EXPERT),
        PARAM_FILTER_NDIFF(PARAM_FILTER_NDIFF_ID, "--diff", "Select N most diverse seqs", "Filter MSAs by selecting most diverse set of sequences, keeping at least this many seqs in each MSA block of length 50", typeid(int), (void *) &Ndiff, "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PROFILE | MMseqsParameter::COMMAND_EXPERT),
        PARAM_FILTER_MIN_ENABLE(PARAM_FILTER_MIN_ENABLE_ID, "--filter-min-enable", "Use filter only at N seqs", "Only filter MSAs with more than N sequences, 0 always filters", typeid(int), (void *) &filterMinEnable, "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PROFILE | MMseqsParameter::COMMAND_EXPERT),

        PARAM_WG(PARAM_WG_ID, "--wg", "Global sequence weighting", "Use global sequence weighting for profile calculation", typeid(bool), (void *) &wg, "", MMseqsParameter::COMMAND_PROFILE | MMseqsParameter::COMMAND_EXPERT),
        PARAM_PC_MODE(PARAM_PC_MODE_ID, "--pseudo-cnt-mode", "Pseudo count mode", "use 0: substitution-matrix or 1: context-specific pseudocounts", typeid(int), (void *) &pcmode, "^[0-1]{1}$", MMseqsParameter::COMMAND_PROFILE | MMseqsParameter::COMMAND_EXPERT),
        PARAM_PCA(PARAM_PCA_ID, "--pca", "Pseudo count a", "Pseudo count admixture strength", typeid(MultiParam<PseudoCounts>), (void *) &pca, "^[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_PROFILE | MMseqsParameter::COMMAND_EXPERT),
        PARAM_PCB(PARAM_PCB_ID, "--pcb", "Pseudo count b", "Pseudo counts: Neff at half of maximum admixture (range 0.0-inf)", typeid(MultiParam<PseudoCounts>), (void *) &pcb, "^[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_PROFILE | MMseqsParameter::COMMAND_EXPERT),
        // sequence2profile
        PARAM_NEFF(PARAM_NEFF_ID, "--neff", "Neff", "Neff included into context state profile (1.0,20.0)", typeid(float), (void *) &neff, "^[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_PROFILE),
        PARAM_TAU(PARAM_TAU_ID, "--tau", "Tau", "Tau: context state pseudo count mixture (0.0,1.0)", typeid(float), (void *) &tau, "[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_PROFILE),
        //createtsv
        PARAM_TARGET_COLUMN(PARAM_TARGET_COLUMN_ID, "--target-column", "Target column", "Select a target column (default 1), 0 if no target id exists", typeid(int), (void *) &targetTsvColumn, "^[0-9]*$"),
        PARAM_FIRST_SEQ_REP_SEQ(PARAM_FIRST_SEQ_REP_SEQ_ID, "--first-seq-as-repr", "First sequence as representative", "Use the first sequence of the clustering result as representative sequence", typeid(bool), (void *) &firstSeqRepr, "", MMseqsParameter::COMMAND_MISC),
        PARAM_FULL_HEADER(PARAM_FULL_HEADER_ID, "--full-header", "Add full header", "Replace DB ID by its corresponding Full Header", typeid(bool), (void *) &fullHeader, ""),
        PARAM_IDX_SEQ_SRC(PARAM_IDX_SEQ_SRC_ID, "--idx-seq-src", "Sequence source", "0: auto, 1: split/translated sequences, 2: input sequences", typeid(int), (void *) &idxSeqSrc, "^[0-2]{1}$", MMseqsParameter::COMMAND_MISC),

        // result2stats
        PARAM_STAT(PARAM_STAT_ID, "--stat", "Statistics to be computed", "One of: linecount, mean, min, max, doolittle, charges, seqlen, firstline", typeid(std::string), (void *) &stat, ""),
        // linearcluster
        PARAM_KMER_PER_SEQ(PARAM_KMER_PER_SEQ_ID, "--kmer-per-seq", "k-mers per sequence", "k-mers per sequence", typeid(int), (void *) &kmersPerSequence, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_CLUSTLINEAR),
        PARAM_KMER_PER_SEQ_SCALE(PARAM_KMER_PER_SEQ_SCALE_ID, "--kmer-per-seq-scale", "Scale k-mers per sequence", "Scale k-mer per sequence based on sequence length as kmer-per-seq val + scale x seqlen", typeid(MultiParam<NuclAA<float>>), (void *) &kmersPerSequenceScale, "^0(\\.[0-9]+)?|1(\\.0+)?$", MMseqsParameter::COMMAND_CLUSTLINEAR),
        PARAM_INCLUDE_ONLY_EXTENDABLE(PARAM_INCLUDE_ONLY_EXTENDABLE_ID, "--include-only-extendable", "Include only extendable", "Include only extendable", typeid(bool), (void *) &includeOnlyExtendable, "", MMseqsParameter::COMMAND_CLUSTLINEAR),
        PARAM_IGNORE_MULTI_KMER(PARAM_IGNORE_MULTI_KMER_ID, "--ignore-multi-kmer", "Skip repeating k-mers", "Skip k-mers occurring multiple times (>=2)", typeid(bool), (void *) &ignoreMultiKmer, "", MMseqsParameter::COMMAND_CLUSTLINEAR | MMseqsParameter::COMMAND_EXPERT),
        PARAM_HASH_SHIFT(PARAM_HASH_SHIFT_ID, "--hash-shift", "Shift hash", "Shift k-mer hash initialization", typeid(int), (void *) &hashShift, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_CLUSTLINEAR | MMseqsParameter::COMMAND_EXPERT),
        PARAM_PICK_N_SIMILAR(PARAM_PICK_N_SIMILAR_ID, "--pick-n-sim-kmer", "Add N similar to search", "Add N similar k-mers to search", typeid(int), (void *) &pickNbest, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_CLUSTLINEAR | MMseqsParameter::COMMAND_EXPERT),
        PARAM_ADJUST_KMER_LEN(PARAM_ADJUST_KMER_LEN_ID, "--adjust-kmer-len", "Adjust k-mer length", "Adjust k-mer length based on specificity (only for nucleotides)", typeid(bool), (void *) &adjustKmerLength, "", MMseqsParameter::COMMAND_CLUSTLINEAR | MMseqsParameter::COMMAND_EXPERT),
        PARAM_RESULT_DIRECTION(PARAM_RESULT_DIRECTION_ID, "--result-direction", "Result direction", "result is 0: query, 1: target centric", typeid(int), (void *) &resultDirection, "^[0-1]{1}$", MMseqsParameter::COMMAND_CLUSTLINEAR | MMseqsParameter::COMMAND_EXPERT),
        PARAM_WEIGHT_FILE(PARAM_WEIGHT_FILE_ID, "--weights", "Weight file name", "Weights used for cluster priorization", typeid(std::string), (void*) &weightFile, "", MMseqsParameter::COMMAND_CLUSTLINEAR | MMseqsParameter::COMMAND_EXPERT ),
        PARAM_WEIGHT_THR(PARAM_WEIGHT_THR_ID, "--cluster-weight-threshold", "Cluster Weight threshold", "Weight threshold used for cluster priorization", typeid(float), (void*) &weightThr, "^[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_CLUSTLINEAR | MMseqsParameter::COMMAND_EXPERT ),
        // workflow
        PARAM_RUNNER(PARAM_RUNNER_ID, "--mpi-runner", "MPI runner", "Use MPI on compute cluster with this MPI command (e.g. \"mpirun -np 42\")", typeid(std::string), (void *) &runner, "", MMseqsParameter::COMMAND_COMMON | MMseqsParameter::COMMAND_EXPERT),
        PARAM_REUSELATEST(PARAM_REUSELATEST_ID, "--force-reuse", "Force restart with latest tmp", "Reuse tmp filse in tmp/latest folder ignoring parameters and version changes", typeid(bool), (void *) &reuseLatest, "", MMseqsParameter::COMMAND_COMMON | MMseqsParameter::COMMAND_EXPERT),
        // search workflow
        PARAM_NUM_ITERATIONS(PARAM_NUM_ITERATIONS_ID, "--num-iterations", "Search iterations", "Number of iterative profile search iterations", typeid(int), (void *) &numIterations, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PROFILE),
        PARAM_START_SENS(PARAM_START_SENS_ID, "--start-sens", "Start sensitivity", "Start sensitivity", typeid(float), (void *) &startSens, "^[0-9]*(\\.[0-9]+)?$"),
        PARAM_SENS_STEPS(PARAM_SENS_STEPS_ID, "--sens-steps", "Search steps", "Number of search steps performed from --start-sens to -s", typeid(int), (void *) &sensSteps, "^[1-9]{1}$"),
        PARAM_EXHAUSTIVE_SEARCH(PARAM_EXHAUSTIVE_SEARCH_ID, "--exhaustive-search", "Exhaustive search mode", "For bigger profile DB, run iteratively the search by greedily swapping the search results", typeid(bool), (void *) &exhaustiveSearch, "", MMseqsParameter::COMMAND_PROFILE | MMseqsParameter::COMMAND_EXPERT),
        PARAM_EXHAUSTIVE_SEARCH_FILTER(PARAM_EXHAUSTIVE_SEARCH_FILTER_ID, "--exhaustive-search-filter", "Filter results during exhaustive search", "Filter result during search: 0: do not filter, 1: filter", typeid(int), (void *) &exhaustiveFilterMsa, "^[0-1]{1}$", MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),

        PARAM_STRAND(PARAM_STRAND_ID, "--strand", "Strand selection", "Strand selection only works for DNA/DNA search 0: reverse, 1: forward, 2: both", typeid(int), (void *) &strand, "^[0-2]{1}$", MMseqsParameter::COMMAND_EXPERT),
        PARAM_ORF_FILTER(PARAM_ORF_FILTER_ID, "--orf-filter", "ORF filter", "Prefilter query ORFs with non-selective search\nOnly used during nucleotide-vs-protein classification\nNOTE: Consider disabling when classifying short reads", typeid(int), (void *) &orfFilter, "^[0-1]{1}$"),
        PARAM_ORF_FILTER_S(PARAM_ORF_FILTER_S_ID, "--orf-filter-s", "ORF filter sensitivity", "Sensitivity used for query ORF prefiltering", typeid(float), (void *) &orfFilterSens, "^[0-9]*(\\.[0-9]+)?$"),
        PARAM_ORF_FILTER_E(PARAM_ORF_FILTER_E_ID, "--orf-filter-e", "ORF filter e-value", "E-value threshold used for query ORF prefiltering", typeid(double), (void *) &orfFilterEval, "^([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)|[0-9]*(\\.[0-9]+)?$"),
        PARAM_LCA_SEARCH(PARAM_LCA_SEARCH_ID, "--lca-search", "LCA search mode", "Efficient search for LCA candidates", typeid(bool), (void *) &lcaSearch, "", MMseqsParameter::COMMAND_PROFILE | MMseqsParameter::COMMAND_EXPERT),
        // easysearch
        PARAM_GREEDY_BEST_HITS(PARAM_GREEDY_BEST_HITS_ID, "--greedy-best-hits", "Greedy best hits", "Choose the best hits greedily to cover the query", typeid(bool), (void *) &greedyBestHits, ""),
        // extractorfs
        PARAM_ORF_MIN_LENGTH(PARAM_ORF_MIN_LENGTH_ID, "--min-length", "Min codons in orf", "Minimum codon number in open reading frames", typeid(int), (void *) &orfMinLength, "^[1-9]{1}[0-9]*$"),
        PARAM_ORF_MAX_LENGTH(PARAM_ORF_MAX_LENGTH_ID, "--max-length", "Max codons in length", "Maximum codon number in open reading frames", typeid(int), (void *) &orfMaxLength, "^[1-9]{1}[0-9]*$"),
        PARAM_ORF_MAX_GAP(PARAM_ORF_MAX_GAP_ID, "--max-gaps", "Max orf gaps", "Maximum number of codons with gaps or unknown residues before an open reading frame is rejected", typeid(int), (void *) &orfMaxGaps, "^(0|[1-9]{1}[0-9]*)$"),
        PARAM_CONTIG_START_MODE(PARAM_CONTIG_START_MODE_ID, "--contig-start-mode", "Contig start mode", "Contig start can be 0: incomplete, 1: complete, 2: both", typeid(int), (void *) &contigStartMode, "^[0-2]{1}"),
        PARAM_CONTIG_END_MODE(PARAM_CONTIG_END_MODE_ID, "--contig-end-mode", "Contig end mode", "Contig end can be 0: incomplete, 1: complete, 2: both", typeid(int), (void *) &contigEndMode, "^[0-2]{1}"),
        PARAM_ORF_START_MODE(PARAM_ORF_START_MODE_ID, "--orf-start-mode", "Orf start mode", "Orf fragment can be 0: from start to stop, 1: from any to stop, 2: from last encountered start to stop (no start in the middle)", typeid(int), (void *) &orfStartMode, "^[0-2]{1}"),
        PARAM_ORF_FORWARD_FRAMES(PARAM_ORF_FORWARD_FRAMES_ID, "--forward-frames", "Forward frames", "Comma-separated list of frames on the forward strand to be extracted", typeid(std::string), (void *) &forwardFrames, ""),
        PARAM_ORF_REVERSE_FRAMES(PARAM_ORF_REVERSE_FRAMES_ID, "--reverse-frames", "Reverse frames", "Comma-separated list of frames on the reverse strand to be extracted", typeid(std::string), (void *) &reverseFrames, ""),
        PARAM_USE_ALL_TABLE_STARTS(PARAM_USE_ALL_TABLE_STARTS_ID, "--use-all-table-starts", "Use all table starts", "Use all alternatives for a start codon in the genetic table, if false - only ATG (AUG)", typeid(bool), (void *) &useAllTableStarts, ""),
        PARAM_TRANSLATE(PARAM_TRANSLATE_ID, "--translate", "Translate orf", "Translate ORF to amino acid", typeid(int), (void *) &translate, "^[0-1]{1}"),
        PARAM_CREATE_LOOKUP(PARAM_CREATE_LOOKUP_ID, "--create-lookup", "Create lookup", "Create database lookup file (can be very large)", typeid(int), (void *) &createLookup, "^[0-1]{1}", MMseqsParameter::COMMAND_EXPERT),
        // indexdb
        PARAM_CHECK_COMPATIBLE(PARAM_CHECK_COMPATIBLE_ID, "--check-compatible", "Check compatible", "0: Always recreate index, 1: Check if recreating index is needed, 2: Fail if index is incompatible", typeid(int), (void *) &checkCompatible, "^[0-2]{1}$", MMseqsParameter::COMMAND_MISC),
        PARAM_SEARCH_TYPE(PARAM_SEARCH_TYPE_ID, "--search-type", "Search type", "Search type 0: auto 1: amino acid, 2: translated, 3: nucleotide, 4: translated nucleotide alignment", typeid(int), (void *) &searchType, "^[0-4]{1}"),
        PARAM_INDEX_SUBSET(PARAM_INDEX_SUBSET_ID, "--index-subset", "Index subset", "Create specialized index with subset of entries\n0: normal index\n1: index without headers\n2: index without prefiltering data\n4: index without aln (for cluster db)\nFlags can be combined bit wise", typeid(int), (void *) &indexSubset, "^[0-7]{1}", MMseqsParameter::COMMAND_EXPERT),
        PARAM_INDEX_DBSUFFIX(PARAM_INDEX_DBSUFFIX_ID, "--index-dbsuffix", "Index dbsuffix", "A suffix of the db (used for cluster dbs)", typeid(std::string), (void *) &indexDbsuffix, "", MMseqsParameter::COMMAND_HIDDEN),
        // createdb
        PARAM_USE_HEADER(PARAM_USE_HEADER_ID, "--use-fasta-header", "Use fasta header", "Use the id parsed from the fasta header as the index key instead of using incrementing numeric identifiers", typeid(bool), (void *) &useHeader, ""),
        PARAM_ID_OFFSET(PARAM_ID_OFFSET_ID, "--id-offset", "Offset of numeric ids", "Numeric ids in index file are offset by this value", typeid(int), (void *) &identifierOffset, "^(0|[1-9]{1}[0-9]*)$"),
        PARAM_DB_TYPE(PARAM_DB_TYPE_ID, "--dbtype", "Database type", "Database type 0: auto, 1: amino acid 2: nucleotides", typeid(int), (void *) &dbType, "[0-2]{1}"),
        PARAM_CREATEDB_MODE(PARAM_CREATEDB_MODE_ID, "--createdb-mode", "Createdb mode", "Createdb mode 0: copy data, 1: soft link data and write new index (works only with single line fasta/q)", typeid(int), (void *) &createdbMode, "^[0-1]{1}$"),
        PARAM_SHUFFLE(PARAM_SHUFFLE_ID, "--shuffle", "Shuffle input database", "Shuffle input database", typeid(bool), (void *) &shuffleDatabase, ""),
        PARAM_WRITE_LOOKUP(PARAM_WRITE_LOOKUP_ID, "--write-lookup", "Write lookup file", "write .lookup file containing mapping from internal id, fasta id and file number", typeid(int), (void *) &writeLookup, "^[0-1]{1}", MMseqsParameter::COMMAND_EXPERT),
        PARAM_USE_HEADER_FILE(PARAM_USE_HEADER_FILE_ID, "--use-header-file", "Use header DB", "use the sequence header DB instead of the body to map the entry keys", typeid(bool), (void *) &useHeaderFile, ""),
        // setextendeddbtype
        PARAM_EXTENDED_DBTYPE(PARAM_EXTENDED_DBTYPE_ID, "--extended-dbtype", "Extended dbtype", "Set extended dbtype 1: compressed, 2: need src, 4: context pseudoe cnts", typeid(int), (void *) &extendedDbtype, "^[0-4]{1}"),
        // splitsequence
        PARAM_SEQUENCE_OVERLAP(PARAM_SEQUENCE_OVERLAP_ID, "--sequence-overlap", "Overlap between sequences", "Overlap between sequences", typeid(int), (void *) &sequenceOverlap, "^(0|[1-9]{1}[0-9]*)$"),
        PARAM_SEQUENCE_SPLIT_MODE(PARAM_SEQUENCE_SPLIT_MODE_ID, "--sequence-split-mode", "Sequence split mode", "Sequence split mode 0: copy data, 1: soft link data and write new index,", typeid(int), (void *) &sequenceSplitMode, "^[0-1]{1}$"),
        PARAM_HEADER_SPLIT_MODE(PARAM_HEADER_SPLIT_MODE_ID, "--headers-split-mode", "Header split mode", "Header split mode: 0: split position, 1: original header", typeid(int), (void *) &headerSplitMode, "^[0-1]{1}$"),
        // createclusearchdb
        PARAM_DB_SUFFIX_LIST(PARAM_DB_SUFFIX_LIST_ID, "--db-suffix-list", "Database suffixes", "Suffixes for database to be split in rep/seq", typeid(std::string), (void *) &dbSuffixList, ""),
        // gff2db
        PARAM_GFF_TYPE(PARAM_GFF_TYPE_ID, "--gff-type", "GFF type", "Comma separated list of feature types in the GFF file to select", typeid(std::string), (void *) &gffType, ""),
        // translatenucs
        PARAM_TRANSLATION_TABLE(PARAM_TRANSLATION_TABLE_ID, "--translation-table", "Translation table", "1) CANONICAL, 2) VERT_MITOCHONDRIAL, 3) YEAST_MITOCHONDRIAL, 4) MOLD_MITOCHONDRIAL, 5) INVERT_MITOCHONDRIAL, 6) CILIATE\n9) FLATWORM_MITOCHONDRIAL, 10) EUPLOTID, 11) PROKARYOTE, 12) ALT_YEAST, 13) ASCIDIAN_MITOCHONDRIAL, 14) ALT_FLATWORM_MITOCHONDRIAL\n15) BLEPHARISMA, 16) CHLOROPHYCEAN_MITOCHONDRIAL, 21) TREMATODE_MITOCHONDRIAL, 22) SCENEDESMUS_MITOCHONDRIAL\n23) THRAUSTOCHYTRIUM_MITOCHONDRIAL, 24) PTEROBRANCHIA_MITOCHONDRIAL, 25) GRACILIBACTERIA, 26) PACHYSOLEN, 27) KARYORELICT, 28) CONDYLOSTOMA\n 29) MESODINIUM, 30) PERTRICH, 31) BLASTOCRITHIDIA", typeid(int), (void *) &translationTable, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_MISC | MMseqsParameter::COMMAND_EXPERT),
        // createseqfiledb
        PARAM_ADD_ORF_STOP(PARAM_ADD_ORF_STOP_ID, "--add-orf-stop", "Add orf stop", "Add stop codon '*' at complete start and end", typeid(bool), (void *) &addOrfStop, ""),
        // createseqfiledb
        PARAM_MIN_SEQUENCES(PARAM_MIN_SEQUENCES_ID, "--min-sequences", "Min sequences", "Minimum number of sequences a cluster may contain", typeid(int), (void *) &minSequences, "^[1-9]{1}[0-9]*$"),
        PARAM_MAX_SEQUENCES(PARAM_MAX_SEQUENCES_ID, "--max-sequences", "Max sequences", "Maximum number of sequences a cluster may contain", typeid(int), (void *) &maxSequences, "^[1-9]{1}[0-9]*$"),
        PARAM_HH_FORMAT(PARAM_HH_FORMAT_ID, "--hh-format", "HH format", "Format entries to use with hhsuite (for singleton clusters)", typeid(bool), (void *) &hhFormat, ""),
        // filterdb
        PARAM_FILTER_COL(PARAM_FILTER_COL_ID, "--filter-column", "Filter column", "column", typeid(int), (void *) &filterColumn, "^[1-9]{1}[0-9]*$"),
        PARAM_COLUMN_TO_TAKE(PARAM_COLUMN_TO_TAKE_ID, "--column-to-take", "Column to take", "column to take in join mode. If -1, the whole line is taken", typeid(int), (void *) &columnToTake, "^(-1|0|[1-9]{1}[0-9]*)$"),
        PARAM_FILTER_REGEX(PARAM_FILTER_REGEX_ID, "--filter-regex", "Filter regex", "Regex to select column (example float: [0-9]*(.[0-9]+)? int:[1-9]{1}[0-9])", typeid(std::string), (void *) &filterColumnRegex, "^.*$"),
        PARAM_FILTER_POS(PARAM_FILTER_POS_ID, "--positive-filter", "Positive filter", "Used in conjunction with --filter-file. If true, out  = in \\intersect filter ; if false, out = in - filter", typeid(bool), (void *) &positiveFilter, ""),
        PARAM_FILTER_FILE(PARAM_FILTER_FILE_ID, "--filter-file", "Filter file", "Specify a file that contains the filtering elements", typeid(std::string), (void *) &filteringFile, ""),
        PARAM_FILTER_EXPRESSION(PARAM_FILTER_EXPRESSION_ID, "--filter-expression", "Filter expression", "Specify a mathematical expression to filter lines", typeid(std::string), (void *) &filterExpression, ""),
        PARAM_MAPPING_FILE(PARAM_MAPPING_FILE_ID, "--mapping-file", "Mapping file", "Specify a file that translates the keys of a DB to new keys, TSV format", typeid(std::string), (void *) &mappingFile, ""),
        PARAM_TRIM_TO_ONE_COL(PARAM_TRIM_TO_ONE_COL_ID, "--trim-to-one-column", "Trim to one column", "Output only the column specified by --filter-column", typeid(bool), (void *) &trimToOneColumn, ""),
        PARAM_EXTRACT_LINES(PARAM_EXTRACT_LINES_ID, "--extract-lines", "Extract N lines", "Extract n lines of each entry", typeid(int), (void *) &extractLines, "^[1-9]{1}[0-9]*$"),
        PARAM_COMP_OPERATOR(PARAM_COMP_OPERATOR_ID, "--comparison-operator", "Numerical comparison operator", "Filter by comparing each entry row numerically by using the le) less-than-equal, ge) greater-than-equal or e) equal operator", typeid(std::string), (void *) &compOperator, ""),
        PARAM_COMP_VALUE(PARAM_COMP_VALUE_ID, "--comparison-value", "Numerical comparison value", "Filter by comparing each entry to this value", typeid(double), (void *) &compValue, "^.*$"),
        PARAM_SORT_ENTRIES(PARAM_SORT_ENTRIES_ID, "--sort-entries", "Sort entries", "Sort column set by --filter-column, by 0: no sorting, 1: increasing, 2: decreasing, 3: random shuffle", typeid(int), (void *) &sortEntries, "^[1-9]{1}[0-9]*$"),
        PARAM_BEATS_FIRST(PARAM_BEATS_FIRST_ID, "--beats-first", "Beats first", "Filter by comparing each entry to the first entry", typeid(bool), (void *) &beatsFirst, ""),
        PARAM_JOIN_DB(PARAM_JOIN_DB_ID, "--join-db", "join to DB", "Join another database entry with respect to the database identifier in the chosen column", typeid(std::string), (void *) &joinDB, ""),
        // besthitperset
        PARAM_SIMPLE_BEST_HIT(PARAM_SIMPLE_BEST_HIT_ID, "--simple-best-hit", "Use simple best hit", "Update the p-value by a single best hit, or by best and second best hits", typeid(bool), (void *) &simpleBestHit, ""),
        PARAM_ALPHA(PARAM_ALPHA_ID, "--alpha", "Alpha", "Set alpha for combining p-values during aggregation", typeid(float), (void *) &alpha, "^0(\\.[0-9]+)?|^1(\\.0+)?$"),
        PARAM_SHORT_OUTPUT(PARAM_SHORT_OUTPUT_ID, "--short-output", "Short output", "The output database will contain only the spread p-value", typeid(bool), (void *) &shortOutput, ""),
        PARAM_AGGREGATION_MODE(PARAM_AGGREGATION_MODE_ID, "--aggregation-mode", "Aggregation mode", "Combined P-values computed from 0: multi-hit, 1: minimum of all P-values, 2: product-of-P-values, 3: truncated product", typeid(int), (void *) &aggregationMode, "^[0-4]{1}$"),
        // concatdb
        PARAM_PRESERVEKEYS(PARAM_PRESERVEKEYS_ID, "--preserve-keys", "Preserve the keys", "The keys of the two DB should be distinct, and they will be preserved in the concatenation", typeid(bool), (void *) &preserveKeysB, ""),
        PARAM_TAKE_LARGER_ENTRY(PARAM_TAKE_LARGER_ENTRY_ID, "--take-larger-entry", "Take the larger entry", "Only keep the larger entry (dataSize >) in the concatenation, both databases need the same keys in the index", typeid(bool), (void *) &takeLargerEntry, ""),
        // offsetalignment
        PARAM_CHAIN_ALIGNMENT(PARAM_CHAIN_ALIGNMENT_ID, "--chain-alignments", "Chain overlapping alignments", "Chain overlapping alignments", typeid(int), (void *) &chainAlignment, "^[0-1]{1}", MMseqsParameter::COMMAND_EXPERT),
        PARAM_MERGE_QUERY(PARAM_MERGE_QUERY_ID, "--merge-query", "Merge query", "Combine ORFs/split sequences to a single entry", typeid(int), (void *) &mergeQuery, "^[0-1]{1}", MMseqsParameter::COMMAND_EXPERT),
        // tsv2db
        PARAM_OUTPUT_DBTYPE(PARAM_OUTPUT_DBTYPE_ID, "--output-dbtype", "Output database type", "Set database type for resulting database: Amino acid sequences 0, Nucl. seq. 1, Profiles 2, Alignment result 5, Clustering result 6, Prefiltering result 7, Taxonomy result 8, Indexed database 9, cA3M MSAs 10, FASTA or A3M MSAs 11, Generic database 12, Omit dbtype file 13, Bi-directional prefiltering result 14, Offsetted headers 15", typeid(int), (void *) &outputDbType, "^(0|[1-9]{1}[0-9]*)$"),
        //diff
        PARAM_USESEQID(PARAM_USESEQID_ID, "--use-seq-id", "Match sequences by their ID", "Sequence ID (Uniprot, GenBank, ...) is used for identifying matches between the old and the new DB", typeid(bool), (void *) &useSequenceId, ""),
        // prefixid
        PARAM_PREFIX(PARAM_PREFIX_ID, "--prefix", "Prefix", "Use this prefix for all entries", typeid(std::string), (void *) &prefix, ""),
        PARAM_TSV(PARAM_TSV_ID, "--tsv", "Tsv", "Return output in TSV format", typeid(bool), (void *) &tsvOut, ""),
        // summarize headers
        PARAM_HEADER_TYPE(PARAM_HEADER_TYPE_ID, "--header-type", "Header type", "Header Type: 1: Uniclust, 2: Metaclust", typeid(int), (void *) &headerType, "[1-2]{1}"),
        // mergedbs
        PARAM_MERGE_PREFIXES(PARAM_MERGE_PREFIXES_ID, "--prefixes", "Merge prefixes", "Comma separated list of prefixes for each entry", typeid(std::string), (void *) &mergePrefixes, "", MMseqsParameter::COMMAND_EXPERT),
        PARAM_MERGE_STOP_EMPTY(PARAM_MERGE_STOP_EMPTY_ID, "--merge-stop-empty", "Stop merge after empty", "Don't continue merging entries after an empty entry", typeid(bool), (void*)&mergeStopEmpty, "", MMseqsParameter::COMMAND_EXPERT),
        // summarizeresult
        PARAM_OVERLAP(PARAM_OVERLAP_ID, "--overlap", "Overlap threshold", "Maximum overlap of covered regions", typeid(float), (void *) &overlap, "^[0-9]*(\\.[0-9]+)?$"),
        // msa2profile
        PARAM_MSA_TYPE(PARAM_MSA_TYPE_ID, "--msa-type", "MSA type", "MSA Type: 0: cA3M, 1: A3M, 2: FASTA", typeid(int), (void *) &msaType, "^[0-2]{1}$"),
        // extractalignedregion
        PARAM_EXTRACT_MODE(PARAM_EXTRACT_MODE_ID, "--extract-mode", "Extract mode", "Extract from 1: Query, 2: Target", typeid(int), (void *) &extractMode, "^[1-2]{1}$"),
        // convertkb
        PARAM_KB_COLUMNS(PARAM_KB_COLUMNS_ID, "--kb-columns", "UniprotKB columns", "list of indices of UniprotKB columns to be extracted", typeid(std::string), (void *) &kbColumns, ""),
        PARAM_RECOVER_DELETED(PARAM_RECOVER_DELETED_ID, "--recover-deleted", "Recover deleted", "Find and recover deleted sequences during updating of clustering", typeid(bool), (void *) &recoverDeleted, ""),
        // filtertaxdb
        PARAM_TAXON_LIST(PARAM_TAXON_LIST_ID, "--taxon-list", "Selected taxa", "Taxonomy ID, possibly multiple values separated by ','", typeid(std::string), (void *) &taxonList, ""),
        // view
        PARAM_ID_LIST(PARAM_ID_LIST_ID, "--id-list", "Selected entries with key", "Entries to be printed separated by ','", typeid(std::string), (void *) &idList, ""),
        PARAM_IDX_ENTRY_TYPE(PARAM_IDX_ENTRY_TYPE_ID, "--idx-entry-type", "Index entry type", "0: sequence, 1: src sequence, 2: header, 3: src header", typeid(int), (void *) &idxEntryType, "^[0-3]{1}$"),
        // lca and addtaxonomy
        PARAM_PICK_ID_FROM(PARAM_PICK_ID_FROM_ID, "--pick-id-from", "Extract mode", "Query 1, Target 2", typeid(int), (void *) &pickIdFrom, "^[1-2]{1}$"),
        PARAM_LCA_RANKS(PARAM_LCA_RANKS_ID, "--lca-ranks", "LCA ranks", "Add column with specified ranks (',' separated)", typeid(std::string), (void *) &lcaRanks, ""),
        PARAM_BLACKLIST(PARAM_BLACKLIST_ID, "--blacklist", "Taxon blacklist", "Comma separated list of ignored taxa in LCA computation", typeid(std::string), (void *) &blacklist, "([0-9]+,)?[0-9]+"),
        PARAM_TAXON_ADD_LINEAGE(PARAM_TAXON_ADD_LINEAGE_ID, "--tax-lineage", "Column with taxonomic lineage", "0: don't show, 1: add all lineage names, 2: add all lineage taxids", typeid(int), (void *) &showTaxLineage, "^[0-2]{1}$"),
        // aggregatetax
        PARAM_MAJORITY(PARAM_MAJORITY_ID, "--majority", "Majority threshold", "minimal fraction of agreement among taxonomically assigned sequences of a set", typeid(float), (void *) &majorityThr, "^0(\\.[0-9]+)?|^1(\\.0+)?$"),
        PARAM_VOTE_MODE(PARAM_VOTE_MODE_ID, "--vote-mode", "Vote mode", "Mode of assigning weights to compute majority. 0: uniform, 1: minus log E-value, 2: score", typeid(int), (void *) &voteMode, "^[0-2]{1}$"),
        // pairaln
        PARAM_PAIRING_DUMMY_MODE(PARAM_PAIRING_DUMMY_MODE_ID, "--pairing-dummy-mode", "Include dummy pairing", "0: dont include, 1: include - an entry that will cause result2msa to write a gap only line", typeid(int), (void *) &pairdummymode, "^[0-1]{1}$", MMseqsParameter::COMMAND_EXPERT),
        PARAM_PAIRING_MODE(PARAM_PAIRING_MODE_ID, "--pairing-mode", "Pairing mode", "0: pair maximal per species, 1: pair only if all chains are covered per species", typeid(int), (void *) &pairmode, "^[0-1]{1}$", MMseqsParameter::COMMAND_EXPERT),
        // taxonomyreport
        PARAM_REPORT_MODE(PARAM_REPORT_MODE_ID, "--report-mode", "Report mode", "Taxonomy report mode 0: Kraken 1: Krona", typeid(int), (void *) &reportMode, "^[0-1]{1}$"),
        // createtaxdb
        PARAM_NCBI_TAX_DUMP(PARAM_NCBI_TAX_DUMP_ID, "--ncbi-tax-dump", "NCBI tax dump directory", "NCBI tax dump directory. The tax dump can be downloaded here \"ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz\"", typeid(std::string), (void *) &ncbiTaxDump, ""),
        PARAM_TAX_MAPPING_FILE(PARAM_TAX_MAPPING_FILE_ID, "--tax-mapping-file", "Taxonomy mapping file", "File to map sequence identifier to taxonomical identifier", typeid(std::string), (void *) &taxMappingFile, ""),
        PARAM_TAX_MAPPING_MODE(PARAM_TAX_MAPPING_MODE_ID, "--tax-mapping-mode", "Taxonomy mapping mode", "Map taxonomy based on sequence database 0: .lookup file 1: .source file", typeid(int), (void *) &taxMappingMode, "^[0-1]{1}$"),
        PARAM_TAX_DB_MODE(PARAM_TAX_DB_MODE_ID, "--tax-db-mode", "Taxonomy db mode", "Create taxonomy database as: 0: .dmp flat files (human readable) 1: binary dump (faster readin)", typeid(int), (void *) &taxDbMode, "^[0-1]{1}$"),
        // expandaln
        PARAM_EXPANSION_MODE(PARAM_EXPANSION_MODE_ID, "--expansion-mode", "Expansion mode", "Update score, E-value, and sequence identity by 0: input alignment 1: rescoring the inferred backtrace", typeid(int), (void *) &expansionMode, "^[0-2]{1}$"),
        PARAM_EXPAND_FILTER_CLUSTERS(PARAM_EXPAND_FILTER_CLUSTERS_ID, "--expand-filter-clusters", "Expand filter clusters", "Filter each target cluster during expansion 0: no filter 1: filter", typeid(int), (void *) &expandFilterClusters, "^[0-1]{1}$"),
        // taxonomy
        PARAM_LCA_MODE(PARAM_LCA_MODE_ID, "--lca-mode", "LCA mode", "LCA Mode 1: single search LCA , 2/3: approximate 2bLCA, 4: top hit", typeid(int), (void *) &taxonomySearchMode, "^[1-4]{1}$"),
        PARAM_TAX_OUTPUT_MODE(PARAM_TAX_OUTPUT_MODE_ID, "--tax-output-mode", "Taxonomy output mode", "0: output LCA, 1: output alignment 2: output both", typeid(int), (void *) &taxonomyOutputMode, "^[0-2]{1}$"),
        // createsubdb, filtertaxseqdb
        PARAM_SUBDB_MODE(PARAM_SUBDB_MODE_ID, "--subdb-mode", "Subdb mode", "Subdb mode 0: copy data 1: soft link data and write index", typeid(int), (void *) &subDbMode, "^[0-1]{1}$"),
        PARAM_ID_MODE(PARAM_ID_MODE_ID, "--id-mode", "Database ID mode", "Select DB entries based on 0: database keys, 1: FASTA identifiers (.lookup)", typeid(int), (void *) &dbIdMode, "^[0-1]{1}$"),
        PARAM_TAR_INCLUDE(PARAM_TAR_INCLUDE_ID, "--tar-include", "Tar Inclusion Regex", "Include file names based on this regex", typeid(std::string), (void *) &tarInclude, "^.*$"),
        PARAM_TAR_EXCLUDE(PARAM_TAR_EXCLUDE_ID, "--tar-exclude", "Tar Exclusion Regex", "Exclude file names based on this regex", typeid(std::string), (void *) &tarExclude, "^.*$"),
        // unpackdb
        PARAM_UNPACK_SUFFIX(PARAM_UNPACK_SUFFIX_ID, "--unpack-suffix", "Unpack suffix", "File suffix for unpacked files.\nAdd .gz suffix to write compressed files.", typeid(std::string), (void *) &unpackSuffix, "^.*$"),
        PARAM_UNPACK_NAME_MODE(PARAM_UNPACK_NAME_MODE_ID, "--unpack-name-mode", "Unpack name mode", "Name unpacked files by 0: DB key, 1: accession (through .lookup)", typeid(int), (void *) &unpackNameMode, "^[0-1]{1}$"),
        // for modules that should handle -h themselves
        PARAM_HELP(PARAM_HELP_ID, "-h", "Help", "Help", typeid(bool), (void *) &help, "", MMseqsParameter::COMMAND_HIDDEN),
        PARAM_HELP_LONG(PARAM_HELP_LONG_ID, "--help", "Help", "Help", typeid(bool), (void *) &help, "", MMseqsParameter::COMMAND_HIDDEN)
{
    if (instance) {
        Debug(Debug::ERROR) << "Parameter instance already exists!\n";
        abort();
    }
    instance = this;


        // onlyverbosity
    onlyverbosity.push_back(&PARAM_V);

    // verbandcompression
    verbandcompression.push_back(&PARAM_COMPRESSED);
    verbandcompression.push_back(&PARAM_V);

    // onlythreads
    onlythreads.push_back(&PARAM_THREADS);
    onlythreads.push_back(&PARAM_V);

    // threadsandcompression
    threadsandcompression.push_back(&PARAM_THREADS);
    threadsandcompression.push_back(&PARAM_COMPRESSED);
    threadsandcompression.push_back(&PARAM_V);

    // createclusearchdb
    createclusearchdb.push_back(&PARAM_THREADS);
    createclusearchdb.push_back(&PARAM_COMPRESSED);
    createclusearchdb.push_back(&PARAM_V);
    createclusearchdb.push_back(&PARAM_DB_SUFFIX_LIST);

    // alignall
    alignall.push_back(&PARAM_SUB_MAT);
    alignall.push_back(&PARAM_ADD_BACKTRACE);
    alignall.push_back(&PARAM_ALIGNMENT_MODE);
//    alignall.push_back(&PARAM_WRAPPED_SCORING);
    alignall.push_back(&PARAM_E);
    alignall.push_back(&PARAM_MIN_SEQ_ID);
    alignall.push_back(&PARAM_MIN_ALN_LEN);
    alignall.push_back(&PARAM_SEQ_ID_MODE);
//    alignall.push_back(&PARAM_ALT_ALIGNMENT);
    alignall.push_back(&PARAM_C);
    alignall.push_back(&PARAM_COV_MODE);
    alignall.push_back(&PARAM_MAX_SEQ_LEN);
    alignall.push_back(&PARAM_NO_COMP_BIAS_CORR);
    alignall.push_back(&PARAM_NO_COMP_BIAS_CORR_SCALE);

//    alignall.push_back(&PARAM_REALIGN);
//    alignall.push_back(&PARAM_MAX_REJECTED);
//    alignall.push_back(&PARAM_MAX_ACCEPT);
    alignall.push_back(&PARAM_INCLUDE_IDENTITY);
    alignall.push_back(&PARAM_PRELOAD_MODE);
    alignall.push_back(&PARAM_PCA);
    alignall.push_back(&PARAM_PCB);
    alignall.push_back(&PARAM_SCORE_BIAS);
    alignall.push_back(&PARAM_GAP_OPEN);
    alignall.push_back(&PARAM_GAP_EXTEND);
    alignall.push_back(&PARAM_ZDROP);
    alignall.push_back(&PARAM_THREADS);
    alignall.push_back(&PARAM_COMPRESSED);
    alignall.push_back(&PARAM_V);

    // alignment
    align.push_back(&PARAM_SUB_MAT);
    align.push_back(&PARAM_ADD_BACKTRACE);
    align.push_back(&PARAM_ALIGNMENT_MODE);
    align.push_back(&PARAM_ALIGNMENT_OUTPUT_MODE);
    align.push_back(&PARAM_WRAPPED_SCORING);
    align.push_back(&PARAM_E);
    align.push_back(&PARAM_MIN_SEQ_ID);
    align.push_back(&PARAM_MIN_ALN_LEN);
    align.push_back(&PARAM_SEQ_ID_MODE);
    align.push_back(&PARAM_ALT_ALIGNMENT);
    align.push_back(&PARAM_C);
    align.push_back(&PARAM_COV_MODE);
    align.push_back(&PARAM_MAX_SEQ_LEN);
    align.push_back(&PARAM_NO_COMP_BIAS_CORR);
    align.push_back(&PARAM_NO_COMP_BIAS_CORR_SCALE);

    align.push_back(&PARAM_MAX_REJECTED);
    align.push_back(&PARAM_MAX_ACCEPT);
    align.push_back(&PARAM_INCLUDE_IDENTITY);
    align.push_back(&PARAM_PRELOAD_MODE);
    align.push_back(&PARAM_PCA);
    align.push_back(&PARAM_PCB);
    align.push_back(&PARAM_SCORE_BIAS);
    align.push_back(&PARAM_REALIGN);
    align.push_back(&PARAM_REALIGN_SCORE_BIAS);
    align.push_back(&PARAM_REALIGN_MAX_SEQS);
    align.push_back(&PARAM_CORR_SCORE_WEIGHT);
    align.push_back(&PARAM_GAP_OPEN);
    align.push_back(&PARAM_GAP_EXTEND);
    align.push_back(&PARAM_ZDROP);
    align.push_back(&PARAM_THREADS);
    align.push_back(&PARAM_COMPRESSED);
    align.push_back(&PARAM_V);

    // prefilter
    prefilter.push_back(&PARAM_SUB_MAT);
    prefilter.push_back(&PARAM_SEED_SUB_MAT);
    prefilter.push_back(&PARAM_S);
    prefilter.push_back(&PARAM_K);
    prefilter.push_back(&PARAM_TARGET_SEARCH_MODE);
    prefilter.push_back(&PARAM_K_SCORE);
    prefilter.push_back(&PARAM_ALPH_SIZE);
    prefilter.push_back(&PARAM_MAX_SEQ_LEN);
    prefilter.push_back(&PARAM_MAX_SEQS);
    prefilter.push_back(&PARAM_SPLIT);
    prefilter.push_back(&PARAM_SPLIT_MODE);
    prefilter.push_back(&PARAM_SPLIT_MEMORY_LIMIT);
    prefilter.push_back(&PARAM_C);
    prefilter.push_back(&PARAM_COV_MODE);
    prefilter.push_back(&PARAM_NO_COMP_BIAS_CORR);
    prefilter.push_back(&PARAM_NO_COMP_BIAS_CORR_SCALE);
    prefilter.push_back(&PARAM_DIAGONAL_SCORING);
    prefilter.push_back(&PARAM_EXACT_KMER_MATCHING);
    prefilter.push_back(&PARAM_MASK_RESIDUES);
    prefilter.push_back(&PARAM_MASK_PROBABILTY);
    prefilter.push_back(&PARAM_MASK_LOWER_CASE);
    prefilter.push_back(&PARAM_MIN_DIAG_SCORE);
    prefilter.push_back(&PARAM_TAXON_LIST);
    prefilter.push_back(&PARAM_INCLUDE_IDENTITY);
    prefilter.push_back(&PARAM_SPACED_KMER_MODE);
    prefilter.push_back(&PARAM_PRELOAD_MODE);
    prefilter.push_back(&PARAM_PCA);
    prefilter.push_back(&PARAM_PCB);
    prefilter.push_back(&PARAM_SPACED_KMER_PATTERN);
    prefilter.push_back(&PARAM_LOCAL_TMP);
    prefilter.push_back(&PARAM_THREADS);
    prefilter.push_back(&PARAM_COMPRESSED);
    prefilter.push_back(&PARAM_V);

    // ungappedprefilter
    ungappedprefilter.push_back(&PARAM_SUB_MAT);
    ungappedprefilter.push_back(&PARAM_C);
    ungappedprefilter.push_back(&PARAM_E);
    ungappedprefilter.push_back(&PARAM_COV_MODE);
    ungappedprefilter.push_back(&PARAM_NO_COMP_BIAS_CORR);
    ungappedprefilter.push_back(&PARAM_NO_COMP_BIAS_CORR_SCALE);
    ungappedprefilter.push_back(&PARAM_MIN_DIAG_SCORE);
    ungappedprefilter.push_back(&PARAM_MAX_SEQS);
    ungappedprefilter.push_back(&PARAM_TAXON_LIST);
    ungappedprefilter.push_back(&PARAM_PRELOAD_MODE);
    ungappedprefilter.push_back(&PARAM_THREADS);
    ungappedprefilter.push_back(&PARAM_COMPRESSED);
    ungappedprefilter.push_back(&PARAM_V);

    // gappedprefilter
    gappedprefilter.push_back(&PARAM_SUB_MAT);
    gappedprefilter.push_back(&PARAM_GAP_OPEN);
    gappedprefilter.push_back(&PARAM_GAP_EXTEND);
    gappedprefilter.push_back(&PARAM_E);
    gappedprefilter.push_back(&PARAM_C);
    gappedprefilter.push_back(&PARAM_COV_MODE);
    gappedprefilter.push_back(&PARAM_NO_COMP_BIAS_CORR);
    gappedprefilter.push_back(&PARAM_NO_COMP_BIAS_CORR_SCALE);
    gappedprefilter.push_back(&PARAM_MIN_DIAG_SCORE);
    gappedprefilter.push_back(&PARAM_MAX_SEQS);
    gappedprefilter.push_back(&PARAM_TAXON_LIST);
    gappedprefilter.push_back(&PARAM_PRELOAD_MODE);
    gappedprefilter.push_back(&PARAM_THREADS);
    gappedprefilter.push_back(&PARAM_COMPRESSED);
    gappedprefilter.push_back(&PARAM_V);

    // clustering
    clust.push_back(&PARAM_CLUSTER_MODE);
    clust.push_back(&PARAM_MAXITERATIONS);
    clust.push_back(&PARAM_SIMILARITYSCORE);
    clust.push_back(&PARAM_THREADS);
    clust.push_back(&PARAM_COMPRESSED);
    clust.push_back(&PARAM_V);
    clust.push_back(&PARAM_WEIGHT_FILE);
    clust.push_back(&PARAM_WEIGHT_THR);

    // rescorediagonal
    rescorediagonal.push_back(&PARAM_SUB_MAT);
    rescorediagonal.push_back(&PARAM_RESCORE_MODE);
    rescorediagonal.push_back(&PARAM_WRAPPED_SCORING);
    rescorediagonal.push_back(&PARAM_FILTER_HITS);
    rescorediagonal.push_back(&PARAM_E);
    rescorediagonal.push_back(&PARAM_C);
    rescorediagonal.push_back(&PARAM_ADD_BACKTRACE);
    rescorediagonal.push_back(&PARAM_COV_MODE);
    rescorediagonal.push_back(&PARAM_MIN_SEQ_ID);
    rescorediagonal.push_back(&PARAM_MIN_ALN_LEN);
    rescorediagonal.push_back(&PARAM_SEQ_ID_MODE);
    rescorediagonal.push_back(&PARAM_INCLUDE_IDENTITY);
    rescorediagonal.push_back(&PARAM_SORT_RESULTS);
    rescorediagonal.push_back(&PARAM_PRELOAD_MODE);
    rescorediagonal.push_back(&PARAM_THREADS);
    rescorediagonal.push_back(&PARAM_COMPRESSED);
    rescorediagonal.push_back(&PARAM_V);

    // alignbykmer
    alignbykmer.push_back(&PARAM_SUB_MAT);
    alignbykmer.push_back(&PARAM_K);
    alignbykmer.push_back(&PARAM_SPACED_KMER_MODE);
    alignbykmer.push_back(&PARAM_SPACED_KMER_PATTERN);
    alignbykmer.push_back(&PARAM_ALPH_SIZE);
    alignbykmer.push_back(&PARAM_FILTER_HITS);
    alignbykmer.push_back(&PARAM_C);
    alignbykmer.push_back(&PARAM_E);
    alignbykmer.push_back(&PARAM_COV_MODE);
    alignbykmer.push_back(&PARAM_MIN_SEQ_ID);
    alignbykmer.push_back(&PARAM_MIN_ALN_LEN);
    alignbykmer.push_back(&PARAM_INCLUDE_IDENTITY);
    alignbykmer.push_back(&PARAM_GAP_OPEN);
    alignbykmer.push_back(&PARAM_GAP_EXTEND);
    alignbykmer.push_back(&PARAM_THREADS);
    alignbykmer.push_back(&PARAM_COMPRESSED);
    alignbykmer.push_back(&PARAM_V);

    // convertprofiledb
    convertprofiledb.push_back(&PARAM_SUB_MAT);
    convertprofiledb.push_back(&PARAM_THREADS);
    convertprofiledb.push_back(&PARAM_COMPRESSED);
    convertprofiledb.push_back(&PARAM_V);


    // sequence2profile
    sequence2profile.push_back(&PARAM_PCA);
    sequence2profile.push_back(&PARAM_PCB);
    sequence2profile.push_back(&PARAM_NEFF);
    sequence2profile.push_back(&PARAM_TAU);
    sequence2profile.push_back(&PARAM_THREADS);
    sequence2profile.push_back(&PARAM_SUB_MAT);
    sequence2profile.push_back(&PARAM_COMPRESSED);
    sequence2profile.push_back(&PARAM_V);

    // create fasta
    createFasta.push_back(&PARAM_V);

    // result2profile
    result2profile.push_back(&PARAM_SUB_MAT);
    result2profile.push_back(&PARAM_E);
    result2profile.push_back(&PARAM_MASK_PROFILE);
    result2profile.push_back(&PARAM_E_PROFILE);
    result2profile.push_back(&PARAM_NO_COMP_BIAS_CORR);
    result2profile.push_back(&PARAM_NO_COMP_BIAS_CORR_SCALE);
    result2profile.push_back(&PARAM_WG);
    result2profile.push_back(&PARAM_ALLOW_DELETION);
    result2profile.push_back(&PARAM_FILTER_MSA);
    result2profile.push_back(&PARAM_FILTER_MIN_ENABLE);
    result2profile.push_back(&PARAM_FILTER_MAX_SEQ_ID);
    result2profile.push_back(&PARAM_FILTER_QID);
    result2profile.push_back(&PARAM_FILTER_QSC);
    result2profile.push_back(&PARAM_FILTER_COV);
    result2profile.push_back(&PARAM_FILTER_NDIFF);
    result2profile.push_back(&PARAM_PC_MODE);
    result2profile.push_back(&PARAM_PCA);
    result2profile.push_back(&PARAM_PCB);
    result2profile.push_back(&PARAM_PRELOAD_MODE);
    result2profile.push_back(&PARAM_GAP_OPEN);
    result2profile.push_back(&PARAM_GAP_EXTEND);
#ifdef GAP_POS_SCORING
    result2profile.push_back(&PARAM_GAP_PSEUDOCOUNT);
#endif
    result2profile.push_back(&PARAM_THREADS);
    result2profile.push_back(&PARAM_COMPRESSED);
    result2profile.push_back(&PARAM_V);

    // createtsv
    createtsv.push_back(&PARAM_FIRST_SEQ_REP_SEQ);
    createtsv.push_back(&PARAM_TARGET_COLUMN);
    createtsv.push_back(&PARAM_FULL_HEADER);
    createtsv.push_back(&PARAM_IDX_SEQ_SRC);
    createtsv.push_back(&PARAM_DB_OUTPUT);
    createtsv.push_back(&PARAM_THREADS);
    createtsv.push_back(&PARAM_COMPRESSED);
    createtsv.push_back(&PARAM_V);

    //result2stats
    result2stats.push_back(&PARAM_STAT);
    result2stats.push_back(&PARAM_TSV);
    result2stats.push_back(&PARAM_COMPRESSED);
    result2stats.push_back(&PARAM_THREADS);
    result2stats.push_back(&PARAM_V);

    // format alignment
    convertalignments.push_back(&PARAM_SUB_MAT);
    convertalignments.push_back(&PARAM_FORMAT_MODE);
    convertalignments.push_back(&PARAM_FORMAT_OUTPUT);
    convertalignments.push_back(&PARAM_TRANSLATION_TABLE);
    convertalignments.push_back(&PARAM_GAP_OPEN);
    convertalignments.push_back(&PARAM_GAP_EXTEND);
    convertalignments.push_back(&PARAM_DB_OUTPUT);
    convertalignments.push_back(&PARAM_PRELOAD_MODE);
    convertalignments.push_back(&PARAM_SEARCH_TYPE);
    convertalignments.push_back(&PARAM_THREADS);
    convertalignments.push_back(&PARAM_COMPRESSED);
    convertalignments.push_back(&PARAM_V);

    // result2msa
    result2msa.push_back(&PARAM_SUB_MAT);
    result2msa.push_back(&PARAM_GAP_OPEN);
    result2msa.push_back(&PARAM_GAP_EXTEND);
    result2msa.push_back(&PARAM_ALLOW_DELETION);
    result2msa.push_back(&PARAM_NO_COMP_BIAS_CORR);
    result2msa.push_back(&PARAM_NO_COMP_BIAS_CORR_SCALE);
    result2msa.push_back(&PARAM_MSA_FORMAT_MODE);
    result2msa.push_back(&PARAM_SUMMARY_PREFIX);
    result2msa.push_back(&PARAM_SKIP_QUERY);
    result2msa.push_back(&PARAM_FILTER_MSA);
    result2msa.push_back(&PARAM_FILTER_MIN_ENABLE);
    result2msa.push_back(&PARAM_FILTER_MAX_SEQ_ID);
    result2msa.push_back(&PARAM_FILTER_QID);
    result2msa.push_back(&PARAM_FILTER_QSC);
    result2msa.push_back(&PARAM_FILTER_COV);
    result2msa.push_back(&PARAM_FILTER_NDIFF);
    result2msa.push_back(&PARAM_PRELOAD_MODE);
    result2msa.push_back(&PARAM_THREADS);
    result2msa.push_back(&PARAM_COMPRESSED);
    result2msa.push_back(&PARAM_V);

    // result2dnamsa
    result2dnamsa.push_back(&PARAM_SKIP_QUERY);
    result2dnamsa.push_back(&PARAM_THREADS);
    result2dnamsa.push_back(&PARAM_COMPRESSED);
    //result2msa.push_back(&PARAM_FIRST_SEQ_REP_SEQ);
    result2dnamsa.push_back(&PARAM_V);


    // filtera3m
    filtera3m.push_back(&PARAM_SUB_MAT);
    filtera3m.push_back(&PARAM_GAP_OPEN);
    filtera3m.push_back(&PARAM_GAP_EXTEND);
    filtera3m.push_back(&PARAM_FILTER_MIN_ENABLE);
    filtera3m.push_back(&PARAM_FILTER_MAX_SEQ_ID);
    filtera3m.push_back(&PARAM_FILTER_QID);
    filtera3m.push_back(&PARAM_FILTER_QSC);
    filtera3m.push_back(&PARAM_FILTER_COV);
    filtera3m.push_back(&PARAM_FILTER_NDIFF);
    filtera3m.push_back(&PARAM_V);

    // filterresult
    filterresult.push_back(&PARAM_SUB_MAT);
    filterresult.push_back(&PARAM_GAP_OPEN);
    filterresult.push_back(&PARAM_GAP_EXTEND);
    filterresult.push_back(&PARAM_NO_COMP_BIAS_CORR);
    filterresult.push_back(&PARAM_NO_COMP_BIAS_CORR_SCALE);
    filterresult.push_back(&PARAM_ALLOW_DELETION);
    filterresult.push_back(&PARAM_FILTER_MIN_ENABLE);
    filterresult.push_back(&PARAM_FILTER_MAX_SEQ_ID);
    filterresult.push_back(&PARAM_FILTER_QID);
    filterresult.push_back(&PARAM_FILTER_QSC);
    filterresult.push_back(&PARAM_FILTER_COV);
    filterresult.push_back(&PARAM_FILTER_NDIFF);
    filterresult.push_back(&PARAM_PRELOAD_MODE);
    filterresult.push_back(&PARAM_THREADS);
    filterresult.push_back(&PARAM_COMPRESSED);
    filterresult.push_back(&PARAM_INCLUDE_IDENTITY);
    filterresult.push_back(&PARAM_V);

    // convertmsa
    convertmsa.push_back(&PARAM_IDENTIFIER_FIELD);
    convertmsa.push_back(&PARAM_COMPRESSED);
    convertmsa.push_back(&PARAM_V);

    // msa2profile
    msa2profile.push_back(&PARAM_MSA_TYPE);
    msa2profile.push_back(&PARAM_SUB_MAT);
    msa2profile.push_back(&PARAM_MATCH_MODE);
    msa2profile.push_back(&PARAM_MATCH_RATIO);
    msa2profile.push_back(&PARAM_PC_MODE);
    msa2profile.push_back(&PARAM_PCA);
    msa2profile.push_back(&PARAM_PCB);
    msa2profile.push_back(&PARAM_NO_COMP_BIAS_CORR);
    msa2profile.push_back(&PARAM_NO_COMP_BIAS_CORR_SCALE);
    msa2profile.push_back(&PARAM_WG);
    msa2profile.push_back(&PARAM_FILTER_MSA);
    msa2profile.push_back(&PARAM_FILTER_MIN_ENABLE);
    msa2profile.push_back(&PARAM_FILTER_COV);
    msa2profile.push_back(&PARAM_FILTER_QID);
    msa2profile.push_back(&PARAM_FILTER_QSC);
    msa2profile.push_back(&PARAM_FILTER_MAX_SEQ_ID);
    msa2profile.push_back(&PARAM_FILTER_NDIFF);
    msa2profile.push_back(&PARAM_GAP_OPEN);
    msa2profile.push_back(&PARAM_GAP_EXTEND);
    msa2profile.push_back(&PARAM_SKIP_QUERY);
#ifdef GAP_POS_SCORING
    msa2profile.push_back(&PARAM_GAP_PSEUDOCOUNT);
#endif
    msa2profile.push_back(&PARAM_THREADS);
    msa2profile.push_back(&PARAM_COMPRESSED);
    msa2profile.push_back(&PARAM_V);

    // profile2pssm
    profile2pssm.push_back(&PARAM_SUB_MAT);
    profile2pssm.push_back(&PARAM_MAX_SEQ_LEN);
    profile2pssm.push_back(&PARAM_NO_COMP_BIAS_CORR);
    profile2pssm.push_back(&PARAM_NO_COMP_BIAS_CORR_SCALE);
    profile2pssm.push_back(&PARAM_DB_OUTPUT);
    profile2pssm.push_back(&PARAM_THREADS);
    profile2pssm.push_back(&PARAM_COMPRESSED);
    profile2pssm.push_back(&PARAM_V);

    // profile2seq (profile2consensus + profile2repseq)
    profile2seq.push_back(&PARAM_SUB_MAT);
    profile2seq.push_back(&PARAM_MAX_SEQ_LEN);
    profile2seq.push_back(&PARAM_THREADS);
    profile2seq.push_back(&PARAM_COMPRESSED);
    profile2seq.push_back(&PARAM_V);


    // extract orf
    extractorfs.push_back(&PARAM_ORF_MIN_LENGTH);
    extractorfs.push_back(&PARAM_ORF_MAX_LENGTH);
    extractorfs.push_back(&PARAM_ORF_MAX_GAP);
    extractorfs.push_back(&PARAM_CONTIG_START_MODE);
    extractorfs.push_back(&PARAM_CONTIG_END_MODE);
    extractorfs.push_back(&PARAM_ORF_START_MODE);
    extractorfs.push_back(&PARAM_ORF_FORWARD_FRAMES);
    extractorfs.push_back(&PARAM_ORF_REVERSE_FRAMES);
    extractorfs.push_back(&PARAM_TRANSLATION_TABLE);
    extractorfs.push_back(&PARAM_TRANSLATE);
    extractorfs.push_back(&PARAM_USE_ALL_TABLE_STARTS);
    extractorfs.push_back(&PARAM_ID_OFFSET);
    extractorfs.push_back(&PARAM_CREATE_LOOKUP);
    extractorfs.push_back(&PARAM_THREADS);
    extractorfs.push_back(&PARAM_COMPRESSED);
    extractorfs.push_back(&PARAM_V);

    // extract frames
    extractframes.push_back(&PARAM_ORF_FORWARD_FRAMES);
    extractframes.push_back(&PARAM_ORF_REVERSE_FRAMES);
    extractframes.push_back(&PARAM_CREATE_LOOKUP);
    extractframes.push_back(&PARAM_THREADS);
    extractframes.push_back(&PARAM_COMPRESSED);
    extractframes.push_back(&PARAM_V);

    // orf to contig
    orftocontig.push_back(&PARAM_THREADS);
    orftocontig.push_back(&PARAM_COMPRESSED);
    orftocontig.push_back(&PARAM_V);

    // orf to contig
    reverseseq.push_back(&PARAM_THREADS);
    reverseseq.push_back(&PARAM_COMPRESSED);
    reverseseq.push_back(&PARAM_V);

    // splitsequence
    splitsequence.push_back(&PARAM_MAX_SEQ_LEN);
    splitsequence.push_back(&PARAM_SEQUENCE_OVERLAP);
    splitsequence.push_back(&PARAM_SEQUENCE_SPLIT_MODE);
    splitsequence.push_back(&PARAM_HEADER_SPLIT_MODE);
    splitsequence.push_back(&PARAM_CREATE_LOOKUP);
    splitsequence.push_back(&PARAM_THREADS);
    splitsequence.push_back(&PARAM_COMPRESSED);
    splitsequence.push_back(&PARAM_V);

    // mask sequence
    masksequence.push_back(&PARAM_MASK_PROBABILTY);
    masksequence.push_back(&PARAM_THREADS);
    masksequence.push_back(&PARAM_COMPRESSED);
    masksequence.push_back(&PARAM_V);
    
    // splitdb
    splitdb.push_back(&PARAM_SPLIT);
    splitdb.push_back(&PARAM_SPLIT_AMINOACID);
    splitdb.push_back(&PARAM_COMPRESSED);
    splitdb.push_back(&PARAM_V);

    // create index
    indexdb.push_back(&PARAM_SEED_SUB_MAT);
    indexdb.push_back(&PARAM_K);
    indexdb.push_back(&PARAM_ALPH_SIZE);
    indexdb.push_back(&PARAM_NO_COMP_BIAS_CORR);
    indexdb.push_back(&PARAM_NO_COMP_BIAS_CORR_SCALE);
    indexdb.push_back(&PARAM_MAX_SEQ_LEN);
    indexdb.push_back(&PARAM_MAX_SEQS);
    indexdb.push_back(&PARAM_INDEX_DBSUFFIX);
    indexdb.push_back(&PARAM_MASK_RESIDUES);
    indexdb.push_back(&PARAM_MASK_PROBABILTY);
    indexdb.push_back(&PARAM_MASK_LOWER_CASE);
    indexdb.push_back(&PARAM_SPACED_KMER_MODE);
    indexdb.push_back(&PARAM_SPACED_KMER_PATTERN);
    indexdb.push_back(&PARAM_S);
    indexdb.push_back(&PARAM_K_SCORE);
    indexdb.push_back(&PARAM_CHECK_COMPATIBLE);
    indexdb.push_back(&PARAM_SEARCH_TYPE);
    indexdb.push_back(&PARAM_SPLIT);
    indexdb.push_back(&PARAM_SPLIT_MEMORY_LIMIT);
    indexdb.push_back(&PARAM_INDEX_SUBSET);
    indexdb.push_back(&PARAM_V);
    indexdb.push_back(&PARAM_THREADS);

    // create kmer index
    kmerindexdb.push_back(&PARAM_SEED_SUB_MAT);
    kmerindexdb.push_back(&PARAM_K);
    kmerindexdb.push_back(&PARAM_HASH_SHIFT);
    kmerindexdb.push_back(&PARAM_KMER_PER_SEQ);
    kmerindexdb.push_back(&PARAM_KMER_PER_SEQ_SCALE);
    kmerindexdb.push_back(&PARAM_MIN_SEQ_ID);
    kmerindexdb.push_back(&PARAM_ADJUST_KMER_LEN);
    kmerindexdb.push_back(&PARAM_SPLIT_MEMORY_LIMIT);
    kmerindexdb.push_back(&PARAM_IGNORE_MULTI_KMER);
    kmerindexdb.push_back(&PARAM_ALPH_SIZE);
    kmerindexdb.push_back(&PARAM_MAX_SEQ_LEN);
    kmerindexdb.push_back(&PARAM_MASK_RESIDUES);
    kmerindexdb.push_back(&PARAM_MASK_PROBABILTY);
    kmerindexdb.push_back(&PARAM_MASK_LOWER_CASE);
    kmerindexdb.push_back(&PARAM_CHECK_COMPATIBLE);
    kmerindexdb.push_back(&PARAM_SEARCH_TYPE);
    kmerindexdb.push_back(&PARAM_SPACED_KMER_MODE);
    kmerindexdb.push_back(&PARAM_SPACED_KMER_PATTERN);
    kmerindexdb.push_back(&PARAM_V);
    kmerindexdb.push_back(&PARAM_THREADS);

    // create db
    createdb.push_back(&PARAM_DB_TYPE);
    createdb.push_back(&PARAM_SHUFFLE);
    createdb.push_back(&PARAM_CREATEDB_MODE);
    createdb.push_back(&PARAM_WRITE_LOOKUP);
    createdb.push_back(&PARAM_ID_OFFSET);
    createdb.push_back(&PARAM_COMPRESSED);
    createdb.push_back(&PARAM_V);

    // convert2fasta
    convert2fasta.push_back(&PARAM_USE_HEADER_FILE);
    convert2fasta.push_back(&PARAM_V);

    // result2flat
    result2flat.push_back(&PARAM_USE_HEADER);
    result2flat.push_back(&PARAM_V);

    // result2repseq
    result2repseq.push_back(&PARAM_PRELOAD_MODE);
    result2repseq.push_back(&PARAM_COMPRESSED);
    result2repseq.push_back(&PARAM_THREADS);
    result2repseq.push_back(&PARAM_V);

    // gff2db
    gff2db.push_back(&PARAM_GFF_TYPE);
    gff2db.push_back(&PARAM_ID_OFFSET);
    gff2db.push_back(&PARAM_THREADS);
    gff2db.push_back(&PARAM_V);


    // translate nucleotide
    translatenucs.push_back(&PARAM_TRANSLATION_TABLE);
    translatenucs.push_back(&PARAM_ADD_ORF_STOP);
    translatenucs.push_back(&PARAM_V);
    translatenucs.push_back(&PARAM_COMPRESSED);
    translatenucs.push_back(&PARAM_THREADS);

    // createseqfiledb
    createseqfiledb.push_back(&PARAM_MIN_SEQUENCES);
    createseqfiledb.push_back(&PARAM_MAX_SEQUENCES);
    createseqfiledb.push_back(&PARAM_HH_FORMAT);
    createseqfiledb.push_back(&PARAM_PRELOAD_MODE);
    createseqfiledb.push_back(&PARAM_THREADS);
    createseqfiledb.push_back(&PARAM_COMPRESSED);
    createseqfiledb.push_back(&PARAM_V);

    // filterDb
    filterDb.push_back(&PARAM_FILTER_EXPRESSION);
    filterDb.push_back(&PARAM_FILTER_COL);
    filterDb.push_back(&PARAM_COLUMN_TO_TAKE);
    filterDb.push_back(&PARAM_FILTER_REGEX);
    filterDb.push_back(&PARAM_FILTER_POS);
    filterDb.push_back(&PARAM_FILTER_FILE);
    filterDb.push_back(&PARAM_BEATS_FIRST);
    filterDb.push_back(&PARAM_MAPPING_FILE);
    filterDb.push_back(&PARAM_TRIM_TO_ONE_COL);
    filterDb.push_back(&PARAM_EXTRACT_LINES);
    filterDb.push_back(&PARAM_COMP_OPERATOR);
    filterDb.push_back(&PARAM_COMP_VALUE);
    filterDb.push_back(&PARAM_SORT_ENTRIES);
    filterDb.push_back(&PARAM_INCLUDE_IDENTITY);
    filterDb.push_back(&PARAM_JOIN_DB);
    filterDb.push_back(&PARAM_THREADS);
    filterDb.push_back(&PARAM_COMPRESSED);
    filterDb.push_back(&PARAM_V);

    // besthitperset
    besthitbyset.push_back(&PARAM_SIMPLE_BEST_HIT);
    besthitbyset.push_back(&PARAM_THREADS);
    besthitbyset.push_back(&PARAM_COMPRESSED);
    besthitbyset.push_back(&PARAM_V);


    // combinepvalperset
    combinepvalbyset.push_back(&PARAM_ALPHA);
    combinepvalbyset.push_back(&PARAM_AGGREGATION_MODE);
//    combinepvalperset.push_back(&PARAM_SHORT_OUTPUT);
    combinepvalbyset.push_back(&PARAM_THREADS);
    combinepvalbyset.push_back(&PARAM_COMPRESSED);
    combinepvalbyset.push_back(&PARAM_V);


    // offsetalignment
    offsetalignment.push_back(&PARAM_CHAIN_ALIGNMENT);
    offsetalignment.push_back(&PARAM_MERGE_QUERY);
    offsetalignment.push_back(&PARAM_SEARCH_TYPE);
    offsetalignment.push_back(&PARAM_THREADS);
    offsetalignment.push_back(&PARAM_COMPRESSED);
    offsetalignment.push_back(&PARAM_PRELOAD_MODE);
    offsetalignment.push_back(&PARAM_V);

    // proteinaln2nucl
    proteinaln2nucl.push_back(&PARAM_SUB_MAT);
    proteinaln2nucl.push_back(&PARAM_GAP_OPEN);
    proteinaln2nucl.push_back(&PARAM_GAP_EXTEND);
    proteinaln2nucl.push_back(&PARAM_THREADS);
    proteinaln2nucl.push_back(&PARAM_COMPRESSED);
    proteinaln2nucl.push_back(&PARAM_V);

    // tsv2db
    tsv2db.push_back(&PARAM_INCLUDE_IDENTITY);
    tsv2db.push_back(&PARAM_OUTPUT_DBTYPE);
    tsv2db.push_back(&PARAM_COMPRESSED);
    tsv2db.push_back(&PARAM_V);

    // swap results
    swapresult.push_back(&PARAM_SUB_MAT);
    swapresult.push_back(&PARAM_E);
    swapresult.push_back(&PARAM_SPLIT_MEMORY_LIMIT);
    swapresult.push_back(&PARAM_GAP_OPEN);
    swapresult.push_back(&PARAM_GAP_EXTEND);
    swapresult.push_back(&PARAM_THREADS);
    swapresult.push_back(&PARAM_COMPRESSED);
    swapresult.push_back(&PARAM_PRELOAD_MODE);
    swapresult.push_back(&PARAM_V);

    // swap results
    swapdb.push_back(&PARAM_SPLIT_MEMORY_LIMIT);
    swapdb.push_back(&PARAM_THREADS);
    swapdb.push_back(&PARAM_COMPRESSED);
    swapdb.push_back(&PARAM_V);

    // subtractdbs
    subtractdbs.push_back(&PARAM_THREADS);
    subtractdbs.push_back(&PARAM_E_PROFILE);
    subtractdbs.push_back(&PARAM_E);
    subtractdbs.push_back(&PARAM_COMPRESSED);
    subtractdbs.push_back(&PARAM_V);

    // setextendeddbtype
    extendeddbtype.push_back(&PARAM_EXTENDED_DBTYPE);
    extendeddbtype.push_back(&PARAM_V);

    // clusthash
    clusthash.push_back(&PARAM_SUB_MAT);
    clusthash.push_back(&PARAM_ALPH_SIZE);
    clusthash.push_back(&PARAM_MIN_SEQ_ID);
    clusthash.push_back(&PARAM_MAX_SEQ_LEN);
    clusthash.push_back(&PARAM_PRELOAD_MODE);
    clusthash.push_back(&PARAM_THREADS);
    clusthash.push_back(&PARAM_COMPRESSED);
    clusthash.push_back(&PARAM_V);

    // kmermatcher
    kmermatcher.push_back(&PARAM_SUB_MAT);
    kmermatcher.push_back(&PARAM_ALPH_SIZE);
    kmermatcher.push_back(&PARAM_MIN_SEQ_ID);
    kmermatcher.push_back(&PARAM_KMER_PER_SEQ);
    kmermatcher.push_back(&PARAM_SPACED_KMER_MODE);
    kmermatcher.push_back(&PARAM_SPACED_KMER_PATTERN);
    kmermatcher.push_back(&PARAM_KMER_PER_SEQ_SCALE);
    kmermatcher.push_back(&PARAM_ADJUST_KMER_LEN);
    kmermatcher.push_back(&PARAM_MASK_RESIDUES);
    kmermatcher.push_back(&PARAM_MASK_PROBABILTY);
    kmermatcher.push_back(&PARAM_MASK_LOWER_CASE);
    kmermatcher.push_back(&PARAM_COV_MODE);
    kmermatcher.push_back(&PARAM_K);
    kmermatcher.push_back(&PARAM_C);
    kmermatcher.push_back(&PARAM_MAX_SEQ_LEN);
    kmermatcher.push_back(&PARAM_HASH_SHIFT);
    kmermatcher.push_back(&PARAM_SPLIT_MEMORY_LIMIT);
    kmermatcher.push_back(&PARAM_INCLUDE_ONLY_EXTENDABLE);
    kmermatcher.push_back(&PARAM_IGNORE_MULTI_KMER);
    kmermatcher.push_back(&PARAM_THREADS);
    kmermatcher.push_back(&PARAM_COMPRESSED);
    kmermatcher.push_back(&PARAM_V);
    kmermatcher.push_back(&PARAM_WEIGHT_FILE);
    kmermatcher.push_back(&PARAM_WEIGHT_THR);

    // kmermatcher
    kmersearch.push_back(&PARAM_SEED_SUB_MAT);
    kmersearch.push_back(&PARAM_KMER_PER_SEQ);
    kmersearch.push_back(&PARAM_KMER_PER_SEQ_SCALE);
    kmersearch.push_back(&PARAM_MASK_RESIDUES);
    kmersearch.push_back(&PARAM_MASK_PROBABILTY);
    kmersearch.push_back(&PARAM_MASK_LOWER_CASE);
    kmersearch.push_back(&PARAM_COV_MODE);
    kmersearch.push_back(&PARAM_C);
    kmersearch.push_back(&PARAM_MAX_SEQ_LEN);
    kmersearch.push_back(&PARAM_PICK_N_SIMILAR);
    kmersearch.push_back(&PARAM_RESULT_DIRECTION);
    kmersearch.push_back(&PARAM_SPLIT_MEMORY_LIMIT);
    kmersearch.push_back(&PARAM_THREADS);
    kmersearch.push_back(&PARAM_COMPRESSED);
    kmersearch.push_back(&PARAM_V);

    // countkmer
    countkmer.push_back(&PARAM_K);
    countkmer.push_back(&PARAM_SPACED_KMER_MODE);
    countkmer.push_back(&PARAM_SPACED_KMER_PATTERN);
    countkmer.push_back(&PARAM_THREADS);


    // mergedbs
    mergedbs.push_back(&PARAM_MERGE_PREFIXES);
    mergedbs.push_back(&PARAM_MERGE_STOP_EMPTY);
    mergedbs.push_back(&PARAM_COMPRESSED);
    mergedbs.push_back(&PARAM_V);

    // summarize
    summarizeheaders.push_back(&PARAM_SUMMARY_PREFIX);
    summarizeheaders.push_back(&PARAM_HEADER_TYPE);
    summarizeheaders.push_back(&PARAM_THREADS);
    summarizeheaders.push_back(&PARAM_COMPRESSED);
    summarizeheaders.push_back(&PARAM_V);

    // diff
    diff.push_back(&PARAM_USESEQID);
    diff.push_back(&PARAM_THREADS);
    diff.push_back(&PARAM_COMPRESSED);
    diff.push_back(&PARAM_V);

    // prefixid
    prefixid.push_back(&PARAM_PREFIX);
    prefixid.push_back(&PARAM_MAPPING_FILE);
    prefixid.push_back(&PARAM_TSV);
    prefixid.push_back(&PARAM_THREADS);
    prefixid.push_back(&PARAM_COMPRESSED);
    prefixid.push_back(&PARAM_V);

    // summarizeresult
    summarizeresult.push_back(&PARAM_ADD_BACKTRACE);
    summarizeresult.push_back(&PARAM_OVERLAP);
    summarizeresult.push_back(&PARAM_C);
    summarizeresult.push_back(&PARAM_THREADS);
    summarizeresult.push_back(&PARAM_COMPRESSED);
    summarizeresult.push_back(&PARAM_V);

    // summarizetabs
    summarizetabs.push_back(&PARAM_OVERLAP);
    summarizetabs.push_back(&PARAM_E);
    summarizetabs.push_back(&PARAM_C);
    summarizetabs.push_back(&PARAM_THREADS);
    summarizetabs.push_back(&PARAM_COMPRESSED);
    summarizetabs.push_back(&PARAM_V);

    // annoate
    extractdomains.push_back(&PARAM_SUB_MAT);
    extractdomains.push_back(&PARAM_MSA_TYPE);
    extractdomains.push_back(&PARAM_E);
    extractdomains.push_back(&PARAM_C);
    extractdomains.push_back(&PARAM_THREADS);
    extractdomains.push_back(&PARAM_COMPRESSED);
    extractdomains.push_back(&PARAM_V);

    // concatdbs
    concatdbs.push_back(&PARAM_COMPRESSED);
    concatdbs.push_back(&PARAM_PRESERVEKEYS);
    concatdbs.push_back(&PARAM_TAKE_LARGER_ENTRY);
    concatdbs.push_back(&PARAM_THREADS);
    concatdbs.push_back(&PARAM_V);

    // extractalignedregion
    extractalignedregion.push_back(&PARAM_COMPRESSED);
    extractalignedregion.push_back(&PARAM_EXTRACT_MODE);
    extractalignedregion.push_back(&PARAM_PRELOAD_MODE);
    extractalignedregion.push_back(&PARAM_THREADS);
    extractalignedregion.push_back(&PARAM_V);

    // convertkb
    convertkb.push_back(&PARAM_COMPRESSED);
    convertkb.push_back(&PARAM_MAPPING_FILE);
    convertkb.push_back(&PARAM_KB_COLUMNS);
    convertkb.push_back(&PARAM_V);

    // filtertaxdb
    filtertaxdb.push_back(&PARAM_COMPRESSED);
    filtertaxdb.push_back(&PARAM_TAXON_LIST);
    filtertaxdb.push_back(&PARAM_THREADS);
    filtertaxdb.push_back(&PARAM_V);

    // filtertaxseqdb
    filtertaxseqdb.push_back(&PARAM_COMPRESSED);
    filtertaxseqdb.push_back(&PARAM_TAXON_LIST);
    filtertaxseqdb.push_back(&PARAM_SUBDB_MODE);
    filtertaxseqdb.push_back(&PARAM_THREADS);
    filtertaxseqdb.push_back(&PARAM_V);

    // aggregatetaxweights
    aggregatetaxweights.push_back(&PARAM_MAJORITY);
    aggregatetaxweights.push_back(&PARAM_VOTE_MODE);
    aggregatetaxweights.push_back(&PARAM_LCA_RANKS);
    aggregatetaxweights.push_back(&PARAM_TAXON_ADD_LINEAGE);
    aggregatetaxweights.push_back(&PARAM_COMPRESSED);
    aggregatetaxweights.push_back(&PARAM_THREADS);
    aggregatetaxweights.push_back(&PARAM_V);

    // aggregatetax
    aggregatetax.push_back(&PARAM_LCA_RANKS);
    aggregatetax.push_back(&PARAM_TAXON_ADD_LINEAGE);
    aggregatetax.push_back(&PARAM_COMPRESSED);
    aggregatetax.push_back(&PARAM_THREADS);
    aggregatetax.push_back(&PARAM_V);

    // TODO should we add this in the future?
    //aggregatetax.push_back(&PARAM_BLACKLIST);

    // lca
    lca.push_back(&PARAM_LCA_RANKS);
    lca.push_back(&PARAM_BLACKLIST);
    lca.push_back(&PARAM_TAXON_ADD_LINEAGE);
    lca.push_back(&PARAM_COMPRESSED);
    lca.push_back(&PARAM_THREADS);
    lca.push_back(&PARAM_V);

    // majoritylca
    majoritylca.push_back(&PARAM_MAJORITY);
    majoritylca.push_back(&PARAM_VOTE_MODE);
    majoritylca.push_back(&PARAM_LCA_RANKS);
    majoritylca.push_back(&PARAM_BLACKLIST);
    majoritylca.push_back(&PARAM_TAXON_ADD_LINEAGE);
    majoritylca.push_back(&PARAM_COMPRESSED);
    majoritylca.push_back(&PARAM_THREADS);
    majoritylca.push_back(&PARAM_V);

    // createsubdb
    createsubdb.push_back(&PARAM_SUBDB_MODE);
    createsubdb.push_back(&PARAM_ID_MODE);
    createsubdb.push_back(&PARAM_V);

    // renamedbkeys
    renamedbkeys.push_back(&PARAM_SUBDB_MODE);
    renamedbkeys.push_back(&PARAM_THREADS);
    renamedbkeys.push_back(&PARAM_V);

    // createtaxdb
    createtaxdb.push_back(&PARAM_NCBI_TAX_DUMP);
    createtaxdb.push_back(&PARAM_TAX_MAPPING_FILE);
    createtaxdb.push_back(&PARAM_TAX_MAPPING_MODE);
    createtaxdb.push_back(&PARAM_TAX_DB_MODE);
    createtaxdb.push_back(&PARAM_THREADS);
    createtaxdb.push_back(&PARAM_V);

    // addtaxonomy
    addtaxonomy.push_back(&PARAM_TAXON_ADD_LINEAGE);
    addtaxonomy.push_back(&PARAM_LCA_RANKS);
    addtaxonomy.push_back(&PARAM_PICK_ID_FROM);
    addtaxonomy.push_back(&PARAM_COMPRESSED);
    addtaxonomy.push_back(&PARAM_THREADS);
    addtaxonomy.push_back(&PARAM_V);

    // taxonomyreport
    taxonomyreport.push_back(&PARAM_REPORT_MODE);
    taxonomyreport.push_back(&PARAM_THREADS);
    taxonomyreport.push_back(&PARAM_V);

    // view
    view.push_back(&PARAM_ID_LIST);
    view.push_back(&PARAM_ID_MODE);
    view.push_back(&PARAM_IDX_ENTRY_TYPE);
    view.push_back(&PARAM_V);

    // exapandaln
    expandaln.push_back(&PARAM_EXPANSION_MODE);
    expandaln.push_back(&PARAM_SUB_MAT);
    expandaln.push_back(&PARAM_GAP_OPEN);
    expandaln.push_back(&PARAM_GAP_EXTEND);
    expandaln.push_back(&PARAM_MAX_SEQ_LEN);
    expandaln.push_back(&PARAM_SCORE_BIAS);
    expandaln.push_back(&PARAM_NO_COMP_BIAS_CORR);
    expandaln.push_back(&PARAM_NO_COMP_BIAS_CORR_SCALE);
    expandaln.push_back(&PARAM_E);
    expandaln.push_back(&PARAM_MIN_SEQ_ID);
//    expandaln.push_back(&PARAM_MIN_SEQ_ID);
//    expandaln.push_back(&PARAM_SEQ_ID_MODE);
    expandaln.push_back(&PARAM_C);
    expandaln.push_back(&PARAM_COV_MODE);
    expandaln.push_back(&PARAM_PC_MODE);
    expandaln.push_back(&PARAM_PCA);
    expandaln.push_back(&PARAM_PCB);
    expandaln.push_back(&PARAM_EXPAND_FILTER_CLUSTERS);
    expandaln.push_back(&PARAM_FILTER_MIN_ENABLE);
    expandaln.push_back(&PARAM_FILTER_MAX_SEQ_ID);
    expandaln.push_back(&PARAM_FILTER_QID);
    expandaln.push_back(&PARAM_FILTER_QSC);
    expandaln.push_back(&PARAM_FILTER_COV);
    expandaln.push_back(&PARAM_FILTER_NDIFF);
    expandaln.push_back(&PARAM_PRELOAD_MODE);
    expandaln.push_back(&PARAM_COMPRESSED);
    expandaln.push_back(&PARAM_THREADS);
    expandaln.push_back(&PARAM_V);

    // expand2profile
    expand2profile.push_back(&PARAM_EXPANSION_MODE);
    expand2profile.push_back(&PARAM_SUB_MAT);
    expand2profile.push_back(&PARAM_GAP_OPEN);
    expand2profile.push_back(&PARAM_GAP_EXTEND);
#ifdef GAP_POS_SCORING
    expand2profile.push_back(&PARAM_GAP_PSEUDOCOUNT);
#endif
    expand2profile.push_back(&PARAM_MAX_SEQ_LEN);
    expand2profile.push_back(&PARAM_SCORE_BIAS);
    expand2profile.push_back(&PARAM_NO_COMP_BIAS_CORR);
    expand2profile.push_back(&PARAM_NO_COMP_BIAS_CORR_SCALE);
    expand2profile.push_back(&PARAM_E_PROFILE);
//    expand2profile.push_back(&PARAM_E);
//    expand2profile.push_back(&PARAM_MIN_SEQ_ID);
//    expand2profile.push_back(&PARAM_SEQ_ID_MODE);
    expand2profile.push_back(&PARAM_C);
    expand2profile.push_back(&PARAM_COV_MODE);
    expand2profile.push_back(&PARAM_MASK_PROFILE);
    expand2profile.push_back(&PARAM_WG);
    expand2profile.push_back(&PARAM_ALLOW_DELETION);
    expand2profile.push_back(&PARAM_EXPAND_FILTER_CLUSTERS);
    expand2profile.push_back(&PARAM_FILTER_MSA);
    expand2profile.push_back(&PARAM_FILTER_MIN_ENABLE);
    expand2profile.push_back(&PARAM_FILTER_MAX_SEQ_ID);
    expand2profile.push_back(&PARAM_FILTER_QID);
    expand2profile.push_back(&PARAM_FILTER_QSC);
    expand2profile.push_back(&PARAM_FILTER_COV);
    expand2profile.push_back(&PARAM_FILTER_NDIFF);
    expand2profile.push_back(&PARAM_PC_MODE);
    expand2profile.push_back(&PARAM_PCA);
    expand2profile.push_back(&PARAM_PCB);
    expand2profile.push_back(&PARAM_PRELOAD_MODE);
    expand2profile.push_back(&PARAM_COMPRESSED);
    expand2profile.push_back(&PARAM_THREADS);
    expand2profile.push_back(&PARAM_V);

    pairaln.push_back(&PARAM_PRELOAD_MODE);
    pairaln.push_back(&PARAM_PAIRING_DUMMY_MODE);
    pairaln.push_back(&PARAM_PAIRING_MODE);
    pairaln.push_back(&PARAM_COMPRESSED);
    pairaln.push_back(&PARAM_THREADS);
    pairaln.push_back(&PARAM_V);

    sortresult.push_back(&PARAM_COMPRESSED);
    sortresult.push_back(&PARAM_THREADS);
    sortresult.push_back(&PARAM_V);

    // WORKFLOWS
    searchworkflow = combineList(align, prefilter);
    searchworkflow = combineList(searchworkflow, rescorediagonal);
    searchworkflow = combineList(searchworkflow, result2profile);
    searchworkflow = combineList(searchworkflow, extractorfs);
    searchworkflow = combineList(searchworkflow, translatenucs);
    searchworkflow = combineList(searchworkflow, splitsequence);
    searchworkflow = combineList(searchworkflow, offsetalignment);
    // needed for slice search, however all its parameters are already present in searchworkflow
    // searchworkflow = combineList(searchworkflow, sortresult);
    searchworkflow.push_back(&PARAM_NUM_ITERATIONS);
    searchworkflow.push_back(&PARAM_START_SENS);
    searchworkflow.push_back(&PARAM_SENS_STEPS);
    searchworkflow.push_back(&PARAM_EXHAUSTIVE_SEARCH);
    searchworkflow.push_back(&PARAM_EXHAUSTIVE_SEARCH_FILTER);
    searchworkflow.push_back(&PARAM_STRAND);
    searchworkflow.push_back(&PARAM_LCA_SEARCH);
    searchworkflow.push_back(&PARAM_DISK_SPACE_LIMIT);
    searchworkflow.push_back(&PARAM_RUNNER);
    searchworkflow.push_back(&PARAM_REUSELATEST);
    searchworkflow.push_back(&PARAM_REMOVE_TMP_FILES);

    linsearchworkflow = combineList(align, kmersearch);
    linsearchworkflow = combineList(linsearchworkflow, swapresult);
    linsearchworkflow = combineList(linsearchworkflow, extractorfs);
    linsearchworkflow = combineList(linsearchworkflow, translatenucs);
    linsearchworkflow = combineList(linsearchworkflow, offsetalignment);
    linsearchworkflow.push_back(&PARAM_RUNNER);
    linsearchworkflow.push_back(&PARAM_REUSELATEST);
    linsearchworkflow.push_back(&PARAM_REMOVE_TMP_FILES);

    // easyslinsearch
    easylinsearchworkflow = combineList(createlinindex, linsearchworkflow);
    easylinsearchworkflow = combineList(easylinsearchworkflow, convertalignments);

    // easysearch
    easysearchworkflow = combineList(searchworkflow, convertalignments);
    easysearchworkflow = combineList(easysearchworkflow, summarizeresult);
    easysearchworkflow = combineList(easysearchworkflow, createdb);
    easysearchworkflow.push_back(&PARAM_GREEDY_BEST_HITS);

    // createindex workflow
    createindex = combineList(indexdb, extractorfs);
    createindex = combineList(createindex, translatenucs);
    createindex = combineList(createindex, splitsequence);
    createindex.push_back(&PARAM_STRAND);
    createindex.push_back(&PARAM_REMOVE_TMP_FILES);

    // createindex workflow
    createlinindex = combineList(kmerindexdb, extractorfs);
    createlinindex = combineList(createlinindex, translatenucs);
    createlinindex.push_back(&PARAM_REMOVE_TMP_FILES);

    // linclust workflow
    linclustworkflow = combineList(clust, align);
    linclustworkflow = combineList(linclustworkflow, kmermatcher);
    linclustworkflow = combineList(linclustworkflow, rescorediagonal);
    linclustworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
    linclustworkflow.push_back(&PARAM_REUSELATEST);
    linclustworkflow.push_back(&PARAM_RUNNER);

    // easylinclustworkflow
    easylinclustworkflow = combineList(linclustworkflow, createdb);

    // clustering workflow
    clusterworkflow = combineList(prefilter, align);
    clusterworkflow = combineList(clusterworkflow, rescorediagonal);
    clusterworkflow = combineList(clusterworkflow, clust);
    clusterworkflow.push_back(&PARAM_CASCADED);
    clusterworkflow.push_back(&PARAM_CLUSTER_STEPS);
    clusterworkflow.push_back(&PARAM_CLUSTER_REASSIGN);
    clusterworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
    clusterworkflow.push_back(&PARAM_REUSELATEST);
    clusterworkflow.push_back(&PARAM_RUNNER);
    clusterworkflow = combineList(clusterworkflow, linclustworkflow);

    // easyclusterworkflow
    easyclusterworkflow = combineList(clusterworkflow, createdb);

    // taxonomy
    taxonomy.push_back(&PARAM_ORF_FILTER);
    taxonomy.push_back(&PARAM_ORF_FILTER_E);
    taxonomy.push_back(&PARAM_ORF_FILTER_S);
    taxonomy.push_back(&PARAM_LCA_MODE);
    taxonomy.push_back(&PARAM_TAX_OUTPUT_MODE);
    taxonomy = combineList(taxonomy, aggregatetaxweights);
    taxonomy = combineList(taxonomy, lca);
    taxonomy = combineList(taxonomy, searchworkflow);
    taxonomy = removeParameter(taxonomy, PARAM_NUM_ITERATIONS);
    taxonomy = removeParameter(taxonomy, PARAM_START_SENS);
    taxonomy = removeParameter(taxonomy, PARAM_SENS_STEPS);

    // easy taxonomy
    easytaxonomy = combineList(taxonomy, addtaxonomy);
    easytaxonomy = combineList(easytaxonomy, taxonomyreport);
    easytaxonomy = combineList(easytaxonomy, convertalignments);
    easytaxonomy = combineList(easytaxonomy, createtsv);
    easytaxonomy = combineList(easytaxonomy, createdb);
    easytaxonomy = removeParameter(easytaxonomy, PARAM_TAX_OUTPUT_MODE);
    easytaxonomy = removeParameter(easytaxonomy, PARAM_PICK_ID_FROM);

    // multi hit db
    multihitdb = combineList(createdb, extractorfs);
    multihitdb = combineList(multihitdb, extractorfs);
    multihitdb = combineList(multihitdb, translatenucs);
    multihitdb = combineList(multihitdb, result2stats);

    // multi hit search
    multihitsearch = combineList(searchworkflow, besthitbyset);

    clusterUpdateSearch = removeParameter(searchworkflow, PARAM_MAX_SEQS);
    clusterUpdateClust = removeParameter(clusterworkflow, PARAM_MAX_SEQS);
    clusterUpdate = combineList(clusterUpdateSearch, clusterUpdateClust);
    clusterUpdate.push_back(&PARAM_REUSELATEST);
    clusterUpdate.push_back(&PARAM_USESEQID);
    clusterUpdate.push_back(&PARAM_RECOVER_DELETED);

    mapworkflow = combineList(prefilter, rescorediagonal);
    mapworkflow = combineList(mapworkflow, extractorfs);
    mapworkflow = combineList(mapworkflow, translatenucs);
    mapworkflow.push_back(&PARAM_START_SENS);
    mapworkflow.push_back(&PARAM_SENS_STEPS);
    mapworkflow.push_back(&PARAM_RUNNER);
    mapworkflow.push_back(&PARAM_REUSELATEST);
    mapworkflow.push_back(&PARAM_REMOVE_TMP_FILES);

    enrichworkflow = combineList(searchworkflow, prefilter);
    enrichworkflow = combineList(enrichworkflow, subtractdbs);
    enrichworkflow = combineList(enrichworkflow, align);
    enrichworkflow = combineList(enrichworkflow, expandaln);
    enrichworkflow = combineList(enrichworkflow, result2profile);

    databases.push_back(&PARAM_HELP);
    databases.push_back(&PARAM_HELP_LONG);
    databases.push_back(&PARAM_TSV);
    databases.push_back(&PARAM_REUSELATEST);
    databases.push_back(&PARAM_REMOVE_TMP_FILES);
    databases.push_back(&PARAM_COMPRESSED);
    databases.push_back(&PARAM_THREADS);
    databases.push_back(&PARAM_V);

    // tar2db
    tar2db.push_back(&PARAM_OUTPUT_DBTYPE);
    tar2db.push_back(&PARAM_TAR_INCLUDE);
    tar2db.push_back(&PARAM_TAR_EXCLUDE);
    tar2db.push_back(&PARAM_COMPRESSED);
    tar2db.push_back(&PARAM_THREADS);
    tar2db.push_back(&PARAM_V);

    // unpackdb
    unpackdbs.push_back(&PARAM_UNPACK_NAME_MODE);
    unpackdbs.push_back(&PARAM_UNPACK_SUFFIX);
    unpackdbs.push_back(&PARAM_THREADS);
    unpackdbs.push_back(&PARAM_V);

    // appenddbtoindex
    appenddbtoindex.push_back(&PARAM_ID_LIST);
    appenddbtoindex.push_back(&PARAM_V);

    //checkSaneEnvironment();
    setDefaults();
}


void Parameters::printUsageMessage(const Command& command, const unsigned int outputFlag, const char* extraText) {
    const std::vector<MMseqsParameter*>& parameters = *command.params;
    std::ostringstream ss;
    ss << "usage: " << binary_name << " " << command.cmd << " " << command.usage << (parameters.size() > 0 ? " [options]" : "") << "\n";
    if (extraText != NULL) {
        ss << extraText;
    }
    if (outputFlag == 0xFFFFFFFF) {
        ss << " By " << command.author << "\n";
    }
    ss << "options: ";

    struct {
        const char* title;
        unsigned int category;
    } categories[] = {
            {"prefilter",MMseqsParameter::COMMAND_PREFILTER},
            {"align",    MMseqsParameter::COMMAND_ALIGN},
            {"clust",    MMseqsParameter::COMMAND_CLUST},
            {"kmermatcher", MMseqsParameter::COMMAND_CLUSTLINEAR},
            {"profile",  MMseqsParameter::COMMAND_PROFILE},
            {"misc",     MMseqsParameter::COMMAND_MISC},
            {"common",   MMseqsParameter::COMMAND_COMMON},
            {"expert",   MMseqsParameter::COMMAND_EXPERT}
    };

    const bool printExpert = MMseqsParameter::COMMAND_EXPERT & outputFlag;
    size_t maxParamWidth = 0;
    for (size_t i = 0; i < ARRAY_SIZE(categories); ++i) {
        for (size_t j = 0; j < parameters.size(); j++) {
            const MMseqsParameter *par = parameters[j];
            bool isExpert = (par->category & MMseqsParameter::COMMAND_EXPERT);
            if (par->category & categories[i].category && (printExpert || isExpert == false)) {
                maxParamWidth = std::max(strlen(parameters[j]->name), maxParamWidth);
            }
        }
    }
    // space in front of options
    maxParamWidth += 6;
    std::map<int, bool> alreadyPrintMap;
    // header
    ss << std::setprecision(3) << std::fixed;
    bool printHeader = true;
    std::string paramString;
    paramString.reserve(100);
    bool hasExpert = false;
    for (size_t i = 0; i < ARRAY_SIZE(categories); ++i) {
        bool categoryFound = false;
        for (size_t j = 0; j < parameters.size(); j++) {
            const MMseqsParameter * par = parameters[j];
            bool isExpert = (par->category & MMseqsParameter::COMMAND_EXPERT);
            bool alreadyPrint = alreadyPrintMap[par->uniqid];
            if (par->category & categories[i].category && (printExpert || isExpert == false ) && alreadyPrint == false ) {
//                int others = (par.category ^ categories[i].category);
//                if(others & outputFlag  )
//                    continue;
                categoryFound = true;
                break;
            }
        }
        if (categoryFound) {
            paramString.clear();
            if (outputFlag == 0xFFFFFFFF) {
                paramString.append(categories[i].title).append(":");
            }
            ss << paramString << std::string(maxParamWidth < paramString.size()? 1 : maxParamWidth - paramString.size(), ' ');

            if(printHeader==true){
                //ss << " ";
                //ss << std::left << std::setw(10) << "default" << " ";
                //ss << std::fixed << "description [value range]";
                printHeader = false;
            }
            ss << std::endl;
            // body
            for (size_t j = 0; j < parameters.size(); j++) {
                paramString.clear();
                const MMseqsParameter * par = parameters[j];
                bool isExpert = (par->category & MMseqsParameter::COMMAND_EXPERT);
                hasExpert |= isExpert;
                bool alreadyPrint = alreadyPrintMap[par->uniqid];
                if (par->category & categories[i].category && (printExpert || isExpert == false) && alreadyPrint == false ) {
                    paramString.append(par->name);
                    std::string valueString;
                    if (par->type == typeid(int)) {
                        paramString.append(" INT");
                        valueString = SSTR(*(int *) par->value);
                    } else if (par->type == typeid(size_t)) {
                        paramString.append(" INT");
                        valueString = SSTR(*(size_t *) par->value);
                    } else if (par->type == typeid(float)) {
                        paramString.append(" FLOAT");
                        valueString = SSTR(*(float *) par->value);
                    } else if (par->type == typeid(double)) {
                        paramString.append(" DOUBLE");
                        valueString = SSTR(*(double *) par->value);
                    } else if (par->type == typeid(ByteParser)) {
                        paramString.append(" BYTE");
                        valueString = ByteParser::format(*((size_t *) par->value), 'a', 'h');
                    } else if (par->type == typeid(bool)) {
                        paramString.append(" BOOL");
                        valueString = SSTR(*(bool *)par->value);
                    } else if (par->type == typeid(std::string)) {
                        paramString.append(" STR");
                        valueString = *((std::string *) par->value);
                    } else if (par->type == typeid(MultiParam<NuclAA<std::string>>)) {
                        paramString.append(" TWIN"); //nucl:VAL,aa:VAL"
                        valueString = MultiParam<NuclAA<std::string>>::format(*((MultiParam<NuclAA<std::string>> *) par->value));
                    } else if (par->type == typeid(MultiParam<NuclAA<int>>)) {
                        paramString.append(" TWIN"); //nucl:VAL,aa:VAL"
                        valueString = MultiParam<NuclAA<int>>::format(*((MultiParam<NuclAA<int>> *) par->value));
                    } else if (par->type == typeid(MultiParam<NuclAA<float>>)) {
                        paramString.append(" TWIN"); //nucl:VAL,aa:VAL"
                        valueString = MultiParam<NuclAA<float>>::format(*((MultiParam<NuclAA<float>> *) par->value));
                    } else if (par->type == typeid(MultiParam<SeqProf<int>>)) {
                        paramString.append(" TWIN"); //seq:VAL,prof:VAL"
                        valueString = MultiParam<SeqProf<int>>::format(*((MultiParam<SeqProf<int>> *) par->value));
                    }

                    ss << " " << paramString << std::string(maxParamWidth < paramString.size()? 1 : maxParamWidth - paramString.size(), ' ');

                    ss << " ";
                    const char* description = par->description;
                    while (description != NULL && *description != '\0') {
                        ss.put(*description);
                        if (*description == '\n') {
                            ss << std::string(maxParamWidth + 2, ' ');
                        }
                        description++;
                    }
                    ss << " [" << valueString << "]";
                    ss << std::endl;
                    alreadyPrintMap[par->uniqid] = true;
                }
            }
        }
    }

    if (command.examples) {
        ss << "\nexamples:\n ";
        const char *data = command.examples;
        while (*data != '\0') {
            ss.put(*data);
            if (*data == '\n') {
                ss.put(' ');
            }
            data++;
        }
        if (*(data - 1) != '\n') {
            ss.put('\n');
        }
    }
    if (command.citations > 0) {
        ss << "\nreferences:\n";
        for (unsigned int pos = 0; pos != sizeof(command.citations) * CHAR_BIT; ++pos) {
            unsigned int citation = 1 << pos;
            if (command.citations & citation && citations.find(citation) != citations.end()) {
                ss << " - " << citations.at(citation) << "\n";
            }
        }
    }
    if (printExpert == false && hasExpert) {
        ss << "\n" << "Show an extended list of options by calling '" << binary_name << " " << command.cmd << " -h'.\n";
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

bool parseBool(const std::string &p) {
    if (p == "true" || p == "TRUE" || p == "1") {
        return true;
    } else if (p == "false" || p == "FALSE" || p == "0") {
        return false;
    } else {
        Debug(Debug::ERROR) << "Invalid boolean string " << p << "\n";
        EXIT(EXIT_FAILURE);
    }
}

void Parameters::initMatrices() {
    // set up substituionMatrix
    for(size_t i = 0 ; i < substitutionMatrices.size(); i++) {
        bool isAminoAcid   = scoringMatrixFile.values.aminoacid() == substitutionMatrices[i].name;
        bool isNucleotide  = scoringMatrixFile.values.nucleotide() == substitutionMatrices[i].name;
        bool isSeedAminoAcid   = seedScoringMatrixFile.values.aminoacid() == substitutionMatrices[i].name;
        bool isSeedNucleotide  = seedScoringMatrixFile.values.nucleotide() == substitutionMatrices[i].name;
        if (isAminoAcid || isNucleotide|| isSeedAminoAcid|| isSeedNucleotide) {
            std::string matrixData((const char *)substitutionMatrices[i].subMatData, substitutionMatrices[i].subMatDataLen);
            std::string matrixName = substitutionMatrices[i].name;
            char* serialized = BaseMatrix::serialize(matrixName, matrixData);
            if(isAminoAcid) {
                scoringMatrixFile.values.aminoacid(serialized);
            }
            if(isNucleotide) {
                scoringMatrixFile.values.nucleotide(serialized);
            }
            if(isSeedAminoAcid) {
                seedScoringMatrixFile.values.aminoacid(serialized);
            }
            if(isSeedNucleotide) {
                seedScoringMatrixFile.values.nucleotide(serialized);
            }
            free(serialized);
        }
    }
}

void Parameters::parseParameters(int argc, const char *pargv[], const Command &command, bool printPar, int parseFlags,
                                 int outputFlags) {
    filenames.clear();
    std::vector<MMseqsParameter*> & par = *command.params;

    bool canHandleHelp = false;
    for (size_t parIdx = 0; parIdx < par.size(); parIdx++) {
        if (par[parIdx]->uniqid == PARAM_HELP_ID || par[parIdx]->uniqid == PARAM_HELP_LONG_ID) {
            canHandleHelp = true;
        }
    }

    size_t parametersFound = 0;
    for (int argIdx = 0; argIdx < argc; argIdx++) {
        // it is a parameter if it starts with - or --
        const bool longParameter = (pargv[argIdx][0] == '-' && pargv[argIdx][1] == '-');
        if (longParameter || (pargv[argIdx][0] == '-')) {
            if ((parseFlags & PARSE_REST) && longParameter && pargv[argIdx][2] == '\0') {
                restArgv = pargv + argIdx + 1;
                restArgc = argc - (argIdx + 1);
                break;
            }
            std::string parameter(pargv[argIdx]);
            if (canHandleHelp == false && (parameter.compare("-h") == 0 || parameter.compare("--help") == 0)) {
                printUsageMessage(command, 0xFFFFFFFF);
                EXIT(EXIT_SUCCESS);
            }

            bool hasUnrecognizedParameter = true;
            for (size_t parIdx = 0; parIdx < par.size(); parIdx++) {
                if(parameter.compare(par[parIdx]->name) == 0) {
                    if (typeid(bool) != par[parIdx]->type && argIdx + 1 == argc) {
                        printUsageMessage(command, outputFlags);
                        Debug(Debug::ERROR) << "Missing argument " << par[parIdx]->name << "\n";
                        EXIT(EXIT_FAILURE);
                    }

                    if (par[parIdx]->wasSet) {
                        printUsageMessage(command, outputFlags);
                        Debug(Debug::ERROR) << "Duplicate parameter " << par[parIdx]->name << "\n";
                        EXIT(EXIT_FAILURE);
                    }

                    if (typeid(int) == par[parIdx]->type) {
                        regex_t regex;
                        compileRegex(&regex, par[parIdx]->regex);
                        int nomatch = regexec(&regex, pargv[argIdx+1], 0, NULL, 0);
                        regfree(&regex);
                        // if no match found or two matches found (we want exactly one match)
                        if (nomatch){
                            printUsageMessage(command, 0xFFFFFFFF);
                            Debug(Debug::ERROR) << "Error in argument " << par[parIdx]->name << "\n";
                            EXIT(EXIT_FAILURE);
                        }else{
                            *((int *) par[parIdx]->value) = atoi(pargv[argIdx+1]);
                            par[parIdx]->wasSet = true;
                        }
                        argIdx++;
                    } else if (typeid(size_t) == par[parIdx]->type) {
                        regex_t regex;
                        compileRegex(&regex, par[parIdx]->regex);
                        int nomatch = regexec(&regex, pargv[argIdx+1], 0, NULL, 0);
                        regfree(&regex);
                        // if no match found or two matches found (we want exactly one match)
                        if (nomatch){
                            printUsageMessage(command, 0xFFFFFFFF);
                            Debug(Debug::ERROR) << "Error in argument " << par[parIdx]->name << "\n";
                            EXIT(EXIT_FAILURE);
                        }else{
                            *((size_t *) par[parIdx]->value) = atoi(pargv[argIdx+1]);
                            par[parIdx]->wasSet = true;
                        }
                        argIdx++;
                    } else if (typeid(ByteParser) == par[parIdx]->type) {
                        regex_t regex;
                        compileRegex(&regex, par[parIdx]->regex);
                        int nomatch = regexec(&regex, pargv[argIdx+1], 0, NULL, 0);
                        regfree(&regex);

                        // if no match found or two matches found (we want exactly one match)
                        if (nomatch){
                            printUsageMessage(command, 0xFFFFFFFF);
                            Debug(Debug::ERROR) << "Error in argument regex " << par[parIdx]->name << "\n";
                            EXIT(EXIT_FAILURE);
                        } else {
                            size_t value = ByteParser::parse(pargv[argIdx+1]);
                            if (value == ByteParser::INVALID_SIZE) {
                                printUsageMessage(command, 0xFFFFFFFF);
                                Debug(Debug::ERROR) << "Error in value parsing " << par[parIdx]->name << "\n";
                                EXIT(EXIT_FAILURE);
                            } else {
                                *((size_t *) par[parIdx]->value) = value;
                                par[parIdx]->wasSet = true;
                            }
                        }
                        argIdx++;
                    } else if (typeid(MultiParam<NuclAA<std::string>>) == par[parIdx]->type) {
                        std::string val(pargv[argIdx+1]);
                        if (Util::startWith("b64:", val)) {
                            val = base64_decode(val.c_str() + 4, val.size() - 4);
                        }
                        NuclAA<std::string> value = MultiParam<NuclAA<std::string>>(val.c_str()).values;
                        if (value.first == "INVALID" || value.second == "INVALID") {
                            printUsageMessage(command, 0xFFFFFFFF);
                            Debug(Debug::ERROR) << "Error in value parsing " << par[parIdx]->name << "\n";
                            EXIT(EXIT_FAILURE);
                        } else {
                            *((MultiParam<NuclAA<std::string>> *) par[parIdx]->value) = value;
                            par[parIdx]->wasSet = true;
                        }
                        argIdx++;
                    }else if (typeid(MultiParam<NuclAA<int>>) == par[parIdx]->type) {
                        NuclAA<int> value = MultiParam<NuclAA<int>>(pargv[argIdx+1]).values;
                        if (value.first == INT_MAX || value.second == INT_MAX) {
                            printUsageMessage(command, 0xFFFFFFFF);
                            Debug(Debug::ERROR) << "Error in value parsing " << par[parIdx]->name << "\n";
                            EXIT(EXIT_FAILURE);
                        } else {
                            *((MultiParam<NuclAA<int>> *) par[parIdx]->value) = value;
                            par[parIdx]->wasSet = true;
                        }
                        argIdx++;
                    }else if (typeid(MultiParam<NuclAA<float>>) == par[parIdx]->type) {
                        NuclAA<float> value = MultiParam<NuclAA<float>>(pargv[argIdx + 1]).values;
                        if (value.first == FLT_MAX || value.second == FLT_MAX) {
                            printUsageMessage(command, 0xFFFFFFFF);
                            Debug(Debug::ERROR) << "Error in value parsing " << par[parIdx]->name << "\n";
                            EXIT(EXIT_FAILURE);
                        } else {
                            *((MultiParam<NuclAA<float>> *) par[parIdx]->value) = value;
                            par[parIdx]->wasSet = true;
                        }
                        argIdx++;
                    }else if (typeid(MultiParam<SeqProf<int>>) == par[parIdx]->type) {
                        SeqProf<int> value = MultiParam<SeqProf<int>>(pargv[argIdx+1]).values;
                        *((MultiParam<SeqProf<int>> *) par[parIdx]->value) = value;
                        par[parIdx]->wasSet = true;
                        argIdx++;
                    }else if (typeid(MultiParam<PseudoCounts>) == par[parIdx]->type) {
                        PseudoCounts value = MultiParam<PseudoCounts>(pargv[argIdx + 1]).values;
                        if (value.first == FLT_MAX || value.second == FLT_MAX) {
                            printUsageMessage(command, 0xFFFFFFFF);
                            Debug(Debug::ERROR) << "Error in value parsing " << par[parIdx]->name << "\n";
                            EXIT(EXIT_FAILURE);
                        } else {
                            *((MultiParam<PseudoCounts> *) par[parIdx]->value) = value;
                            par[parIdx]->wasSet = true;
                        }
                        argIdx++;
                    }else if (typeid(float) == par[parIdx]->type) {
                        regex_t regex;
                        compileRegex(&regex, par[parIdx]->regex);
                        int nomatch = regexec(&regex, pargv[argIdx+1], 0, NULL, 0);
                        regfree(&regex);
                        if (nomatch){
                            printUsageMessage(command, 0xFFFFFFFF);
                            Debug(Debug::ERROR) << "Error in argument " << par[parIdx]->name << "\n";
                            EXIT(EXIT_FAILURE);
                        }else{
                            double input = strtod(pargv[argIdx+1], NULL);
                            *((float *) par[parIdx]->value) = static_cast<float>(input);
                            par[parIdx]->wasSet = true;
                        }
                        argIdx++;
                    } else if (typeid(double) == par[parIdx]->type) {
                        regex_t regex;
                        compileRegex(&regex, par[parIdx]->regex);
                        int nomatch = regexec(&regex, pargv[argIdx+1], 0, NULL, 0);
                        regfree(&regex);
                        if (nomatch){
                            printUsageMessage(command, 0xFFFFFFFF);
                            Debug(Debug::ERROR) << "Error in argument " << par[parIdx]->name << "\n";
                            EXIT(EXIT_FAILURE);
                        }else{
                            *((double *) par[parIdx]->value) = strtod(pargv[argIdx+1], NULL);
                            par[parIdx]->wasSet = true;
                        }
                        argIdx++;
                    } else if (typeid(std::string) == par[parIdx]->type) {
                        std::string val(pargv[argIdx+1]);
                        if (Util::startWith("b64:", val)) {
                            val = base64_decode(val.c_str() + 4, val.size() - 4);
                        }
                        std::string* currVal = (std::string*)par[parIdx]->value;
                        currVal->assign(val);
                        par[parIdx]->wasSet = true;
                        argIdx++;
                    } else if (typeid(bool) == par[parIdx]->type) {
                        bool *value = (bool *) par[parIdx]->value;
                        if (argIdx + 1 == argc || pargv[argIdx+1][0] == '-') {
                            *value = !*value;
                        } else {
                            *value = parseBool(pargv[argIdx+1]);
                            argIdx++;
                        }
                        par[parIdx]->wasSet = true;
                    } else {
                        Debug(Debug::ERROR) << "Wrong parameter type in parseParameters. Please inform the developers\n";
                        EXIT(EXIT_FAILURE);
                    }

                    hasUnrecognizedParameter = false;
                    continue;
                }
            }

            if (hasUnrecognizedParameter) {
                printUsageMessage(command, 0xFFFFFFFF);

                // Suggest some parameter that the user might have meant
                std::vector<MMseqsParameter *>::const_iterator index = par.end();
                int maxDistance = 0;
                for (std::vector<MMseqsParameter *>::const_iterator it = par.begin(); it != par.end(); ++it) {
                    int distance = DistanceCalculator::localLevenshteinDistance(parameter, (*it)->name);
                    if (distance > maxDistance) {
                        maxDistance = distance;
                        index = it;
                    }
                }

                Debug(Debug::ERROR) << "Unrecognized parameter \"" << parameter << "\"";
                if (index != par.end()) {
                    Debug(Debug::ERROR) << ". Did you mean \"" << (*index)->name << "\" (" << (*index)->display << ")?\n";
                } else {
                    Debug(Debug::ERROR) << "\n";
                }

                EXIT(EXIT_FAILURE);
            }

            parametersFound++;
        } else {
            // parameter is actually a filename
#ifdef __CYGWIN__
            // normalize windows paths to cygwin unix paths
            const char *path = pargv[argIdx];
            ssize_t size = cygwin_conv_path(CCP_WIN_A_TO_POSIX | CCP_RELATIVE, path, NULL, 0);
            if (size < 0) {
                Debug(Debug::ERROR) << "Could not convert cygwin path!\n";
                EXIT(EXIT_FAILURE);
            } else {
                char *posix = new char[size];
                if (cygwin_conv_path(CCP_WIN_A_TO_POSIX | CCP_RELATIVE, path, posix, size)) {
                    Debug(Debug::ERROR) << "Could not convert cygwin path!\n";
                    EXIT(EXIT_FAILURE);
                }
                filenames.emplace_back(posix);
                delete posix;
            }
#else
            filenames.emplace_back(pargv[argIdx]);
#endif
        }
    }

    if (MMseqsMPI::isMaster()) {
        Debug::setDebugLevel(verbosity);
    }

#ifdef OPENMP
    omp_set_num_threads(threads);
#endif
#ifndef OPENMP
    threads = 1;
#endif


    bool ignorePathCountChecks = command.databases.empty() == false && command.databases[0].specialType & DbType::ZERO_OR_ALL && filenames.size() == 0;
    const size_t MAX_DB_PARAMETER = 6;
    if (ignorePathCountChecks == false && command.databases.size() > MAX_DB_PARAMETER) {
        Debug(Debug::ERROR) << "Use argv if you need more than " << MAX_DB_PARAMETER << " db parameters" << "\n";
        EXIT(EXIT_FAILURE);
    }

    if (ignorePathCountChecks == false && filenames.size() < command.databases.size()){
        printUsageMessage(command, outputFlags);
        Debug(Debug::ERROR) << "Not enough input paths provided. ";
        if (command.databases.size() == 1) {
            Debug(Debug::ERROR) << "1 path is required.\n";
        } else {
            Debug(Debug::ERROR) << command.databases.size() << " paths are required.\n";
        }
        EXIT(EXIT_FAILURE);
    }

    bool isVar = false;
    bool isStartVar = false;
    bool isMiddleVar = false;
    bool isEndVar = false;
    if(command.databases.empty() == false && command.databases[0].validator != NULL) {
        if (command.databases.size() >= 2) {
            for(size_t i = 0; i < command.databases.size();i++){
                if(i == 0){
                    isStartVar |= (command.databases[i].specialType & DbType::VARIADIC);
                } else if(i == command.databases.size() - 1){
                    isEndVar |= (command.databases[i].specialType & DbType::VARIADIC);
                } else {
                    isMiddleVar |= (command.databases[i].specialType & DbType::VARIADIC);
                }

            }
            isVar = isStartVar | isMiddleVar | isEndVar;
        }
        if (ignorePathCountChecks == false && isVar == false && filenames.size() > command.databases.size()) {
            printUsageMessage(command, outputFlags);
            Debug(Debug::ERROR) << "Too many input paths provided. Only " << SSTR(command.databases.size()) << " are allowed\n";
            EXIT(EXIT_FAILURE);
        }
    }
    switch (std::min(filenames.size(), MAX_DB_PARAMETER)) {
        case 6:
            db6 = filenames[5];
            db6Index = db6;
            db6Index.append(".index");
            db6dbtype = db6;
            db6dbtype.append(".dbtype");
            hdr6 = db6;
            hdr6.append("_h");
            hdr6Index = hdr6;
            hdr6Index.append(".index");
            hdr6dbtype = hdr6;
            hdr6dbtype.append(".dbtype");
            // FALLTHROUGH
        case 5:
            db5 = filenames[4];
            db5Index = db5;
            db5Index.append(".index");
            db5dbtype = db5;
            db5dbtype.append(".dbtype");
            hdr5 = db5;
            hdr5.append("_h");
            hdr5Index = hdr5;
            hdr5Index.append(".index");
            hdr5dbtype = hdr5;
            hdr5dbtype.append(".dbtype");
            // FALLTHROUGH
        case 4:
            db4 = filenames[3];
            db4Index = db4;
            db4Index.append(".index");
            db4dbtype = db4;
            db4dbtype.append(".dbtype");
            hdr4 = db4;
            hdr4.append("_h");
            hdr4Index = hdr4;
            hdr4Index.append(".index");
            hdr4dbtype = hdr4;
            hdr4dbtype.append(".dbtype");
            // FALLTHROUGH
        case 3:
            db3 = filenames[2];
            db3Index = db3;
            db3Index.append(".index");
            db3dbtype = db3;
            db3dbtype.append(".dbtype");
            hdr3 = db3;
            hdr3.append("_h");
            hdr3Index = hdr3;
            hdr3Index.append(".index");
            hdr3dbtype = hdr3;
            hdr3dbtype.append(".dbtype");
            // FALLTHROUGH
        case 2:
            db2 = filenames[1];
            db2Index = db2;
            db2Index.append(".index");
            db2dbtype = db2;
            db2dbtype.append(".dbtype");
            hdr2 = db2;
            hdr2.append("_h");
            hdr2Index = hdr2;
            hdr2Index.append(".index");
            hdr2dbtype = hdr2;
            hdr2dbtype.append(".dbtype");
            // FALLTHROUGH
        case 1:
            db1 = filenames[0];
            db1Index = db1;
            db1Index.append(".index");
            db1dbtype = db1;
            db1dbtype.append(".dbtype");
            hdr1 = db1;
            hdr1.append("_h");
            hdr1Index = hdr1;
            hdr1Index.append(".index");
            hdr1dbtype = hdr1;
            hdr1dbtype.append(".dbtype");
            break;
        default:
            // Do not abort execution if we expect a variable amount of parameters
            if (parseFlags & PARSE_VARIADIC)
                break;
            // FALLTHROUGH
        case 0:
            if (parseFlags & PARSE_ALLOW_EMPTY)
                break;
            printUsageMessage(command, outputFlags);
            printParameters(command.cmd, argc, pargv, par);
            Debug(Debug::ERROR) << "Unrecognized parameters!" << "\n";
            EXIT(EXIT_FAILURE);
    }

    initMatrices();

    if (ignorePathCountChecks == false) {
        checkIfDatabaseIsValid(command, argc, pargv, isStartVar, isMiddleVar, isEndVar);
    }

    if (printPar == true) {
        printParameters(command.cmd, argc, pargv, par);
    }
}

std::vector<std::string> Parameters::findMissingTaxDbFiles(const std::string &filename) {
    std::vector<std::string> missingFiles;
    if (FileUtil::fileExists((filename + "_mapping").c_str()) == false) {
        missingFiles.emplace_back(filename + "_mapping");
    } else if (FileUtil::fileExists((filename + "_taxonomy").c_str()) == true) {
        return missingFiles;
    }
    const std::vector<std::string> suffices = {"_nodes.dmp",  "_names.dmp", "_merged.dmp"};
    for (size_t i = 0; i < suffices.size(); ++i) {
        if (FileUtil::fileExists((filename + suffices[i]).c_str()) == false) {
            missingFiles.emplace_back(filename + suffices[i]);
        }
    }
    return missingFiles;
}

void Parameters::printTaxDbError(const std::string &filename, const std::vector<std::string>& missingFiles) {
    Debug(Debug::ERROR) << "Input taxonomy database \"" << filename << "\" is missing files:\n";
    for (size_t i = 0; i < missingFiles.size(); ++i) {
        Debug(Debug::ERROR) << "- " << missingFiles[i] << "\n";
    }
}

void Parameters::checkIfDatabaseIsValid(const Command& command, int argc, const char** argv, bool isStartVar, bool isMiddleVar, bool isEndVar) {
    size_t fileIdx = 0;
    for (size_t dbIdx = 0; dbIdx < command.databases.size(); dbIdx++) {
        const DbType &db = command.databases[dbIdx];

        // special checks
        if (db.accessMode == db.ACCESS_MODE_INPUT) {
            size_t argumentDist = 0;
            if (dbIdx == 0 && isStartVar) {
                argumentDist = (filenames.size() - command.databases.size());
            } else if (dbIdx == command.databases.size() - 1 && isEndVar) {
                argumentDist = (filenames.size() - command.databases.size());
            } else if ((command.databases[dbIdx].specialType & DbType::VARIADIC) && isMiddleVar) {
                argumentDist = (filenames.size() - command.databases.size());
            }

            size_t currFileIdx = fileIdx;
            for (; fileIdx <= currFileIdx + argumentDist; fileIdx++) {
                if (db.validator == NULL) {
                    continue;
                }

                if (filenames[fileIdx] != "stdin" && FileUtil::fileExists((filenames[fileIdx]).c_str()) == false && FileUtil::fileExists((filenames[fileIdx] + ".dbtype").c_str()) == false) {
                    regex_t regex;
                    compileRegex(&regex, "[a-zA-Z][a-zA-Z0-9+-.]*:\\/\\/");
                    int nomatch = regexec(&regex, filenames[fileIdx].c_str(), 0, NULL, 0);
                    regfree(&regex);
                    if (nomatch) {
                        printParameters(command.cmd, argc, argv, *command.params);
                        Debug(Debug::ERROR) << "Input " << filenames[fileIdx] << " does not exist\n";
                        EXIT(EXIT_FAILURE);
                    }
                }
                int dbtype = FileUtil::parseDbType(filenames[fileIdx].c_str());
                if (db.specialType & DbType::NEED_HEADER && FileUtil::fileExists((filenames[fileIdx] + "_h.dbtype").c_str()) == false && Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_INDEX_DB) == false) {
                    printParameters(command.cmd, argc, argv, *command.params);
                    Debug(Debug::ERROR) << "Database " << filenames[fileIdx] << " needs header information\n";
                    EXIT(EXIT_FAILURE);
                }
                if (db.specialType & DbType::NEED_TAXONOMY) {
                    std::vector<std::string> missingFiles = findMissingTaxDbFiles(filenames[fileIdx]);
                    if (missingFiles.empty() == false) {
                        printParameters(command.cmd, argc, argv, *command.params);
                        printTaxDbError(filenames[fileIdx], missingFiles);
                        EXIT(EXIT_FAILURE);
                    }
                }
                if (db.specialType & DbType::NEED_LOOKUP && FileUtil::fileExists((filenames[fileIdx] + ".lookup").c_str()) == false) {
                    printParameters(command.cmd, argc, argv, *command.params);
                    Debug(Debug::ERROR) << "Database " << filenames[fileIdx] << " needs a lookup file\n";
                    EXIT(EXIT_FAILURE);
                }
                bool dbtypeFound = false;
                for (size_t i = 0; i < db.validator->size() && dbtypeFound == false; i++) {
                    int validatorDbtype = db.validator->at(i);
                    if (validatorDbtype == Parameters::DBTYPE_STDIN) {
                        dbtypeFound = (filenames[fileIdx] == "stdin");
                    } else if (validatorDbtype == Parameters::DBTYPE_URI) {
                        regex_t regex;
                        compileRegex(&regex, "[a-zA-Z][a-zA-Z0-9+-.]*:\\/\\/");
                        int nomatch = regexec(&regex, filenames[fileIdx].c_str(), 0, NULL, 0);
                        regfree(&regex);
                        dbtypeFound = nomatch == false;
                    } else if (validatorDbtype == Parameters::DBTYPE_FLATFILE) {
                        dbtypeFound = (FileUtil::fileExists(filenames[fileIdx].c_str()) == true &&
                                       FileUtil::directoryExists(filenames[fileIdx].c_str()) == false);
                    } else if (validatorDbtype == Parameters::DBTYPE_DIRECTORY) {
                        dbtypeFound = FileUtil::directoryExists(filenames[fileIdx].c_str());
                    } else if (Parameters::isEqualDbtype(db.validator->at(i), dbtype)) {
                        dbtypeFound = true;
                    }
                }
                if (dbtypeFound == false) {
                    printParameters(command.cmd, argc, argv, *command.params);
                    Debug(Debug::ERROR) << "Input database \"" << filenames[fileIdx] << "\" has the wrong type ("
                                        << Parameters::getDbTypeName(dbtype) << ")\nAllowed input:\n";
                    for (size_t i = 0; i < db.validator->size(); ++i) {
                        Debug(Debug::ERROR) << "- " << Parameters::getDbTypeName(db.validator->at(i)) << "\n";
                    }
                    EXIT(EXIT_FAILURE);
                }
            }
        } else if (db.accessMode == db.ACCESS_MODE_OUTPUT) {
            if (db.validator == &DbValidator::directory) {
                if (FileUtil::directoryExists(filenames[fileIdx].c_str()) == false) {
                    if (FileUtil::makeDir(filenames[fileIdx].c_str()) == false) {
                        printParameters(command.cmd, argc, argv, *command.params);
                        Debug(Debug::ERROR) << "Cannot create temporary directory " << filenames[fileIdx] << "\n";
                        EXIT(EXIT_FAILURE);
                    } else {
                        Debug(Debug::INFO) << "Create directory " << filenames[fileIdx] << "\n";
                    }
                }
                fileIdx++;
            } else {
                if (FileUtil::fileExists(filenames[fileIdx].c_str()) == true) {
                    Debug(Debug::WARNING) << filenames[fileIdx] << " exists and will be overwritten\n";
                }
                fileIdx++;
            }
        } else {
            fileIdx++;
        }
    }
}


void Parameters::printParameters(const std::string &module, int argc, const char* pargv[],
                                 const std::vector<MMseqsParameter*> &par){
    if (Debug::debugLevel < Debug::INFO) {
        return;
    }

    Debug(Debug::INFO) << module << " ";
    for (int i = 0; i < argc; i++) {
        // don't expose users to the interal b64 masking of whitespace characters
        if (strncmp("b64:", pargv[i], 4) == 0) {
            Debug(Debug::INFO) << "'" << base64_decode(pargv[i] + 4, strlen(pargv[i]) - 4) << "' ";
        } else {
            Debug(Debug::INFO) << pargv[i] << " ";
        }
    }
    Debug(Debug::INFO) << "\n\n";

    if (CommandCaller::getCallDepth() > 0) {
        return;
    }

    size_t maxWidth = 0;
    for(size_t i = 0; i < par.size(); i++) {
        maxWidth = std::max(strlen(par[i]->display), maxWidth);
    }

    std::stringstream ss;
    ss << std::boolalpha;

    ss << std::setw(maxWidth) << std::left  << "MMseqs Version:" << "\t" << version << "\n";


    for (size_t i = 0; i < par.size(); i++) {
        if (par[i]->category & MMseqsParameter::COMMAND_HIDDEN) {
            continue;
        }
        ss << std::setw(maxWidth) << std::left << par[i]->display << "\t";
        if (typeid(int) == par[i]->type ) {
            ss << *((int *)par[i]->value);
        } else if(typeid(size_t) == par[i]->type ){
            ss << *((size_t *)par[i]->value);
        } else if(typeid(ByteParser) == par[i]->type) {
            ss << ByteParser::format(*((size_t *)par[i]->value), 'a', 'h');
        } else if(PARAM_SUB_MAT.uniqid == par[i]->uniqid || PARAM_SEED_SUB_MAT.uniqid == par[i]->uniqid) {
            MultiParam<NuclAA<std::string>>* param = ((MultiParam<NuclAA<std::string>>*) par[i]->value);
            MultiParam<NuclAA<std::string>> tmpPar(NuclAA<std::string>(
                BaseMatrix::unserializeName(param->values.aminoacid().c_str()),
                BaseMatrix::unserializeName(param->values.nucleotide().c_str())
            ));
            ss << MultiParam<NuclAA<std::string>>::format(tmpPar);
        } else if(typeid(MultiParam<NuclAA<std::string>>) == par[i]->type) {
            ss << MultiParam<NuclAA<std::string>>::format(*((MultiParam<NuclAA<std::string>> *)par[i]->value));
        } else if(typeid(MultiParam<NuclAA<int>>) == par[i]->type) {
            ss << MultiParam<NuclAA<int>>::format(*((MultiParam<NuclAA<int>> *)par[i]->value));
        } else if(typeid(MultiParam<NuclAA<float>>) == par[i]->type) {
            ss << MultiParam<NuclAA<float>>::format(*((MultiParam<NuclAA<float>> *)par[i]->value));
        } else if(typeid(MultiParam<SeqProf<int>>) == par[i]->type) {
            ss << MultiParam<SeqProf<int>>::format(*((MultiParam<SeqProf<int>> *)par[i]->value));
        } else if(typeid(MultiParam<PseudoCounts>) == par[i]->type) {
            ss << MultiParam<PseudoCounts>::format(*((MultiParam<PseudoCounts> *)par[i]->value));
        } else if(typeid(float) == par[i]->type) {
            ss << *((float *)par[i]->value);
        } else if(typeid(double) == par[i]->type) {
            ss << *((double *)par[i]->value);
        } else if(typeid(std::string) == par[i]->type) {
            ss << *((std::string *) par[i]->value);
        } else if (typeid(bool) == par[i]->type) {
            ss << *((bool *)par[i]->value);
        }
        ss << "\n";
    }

    Debug(Debug::INFO) << ss.str() << "\n";
}



void Parameters::setDefaults() {
    restArgv = NULL;
    restArgc = 0;

    scoringMatrixFile =  MultiParam<NuclAA<std::string>>(NuclAA<std::string>("blosum62.out", "nucleotide.out"));
    seedScoringMatrixFile = MultiParam<NuclAA<std::string>>(NuclAA<std::string>("VTML80.out", "nucleotide.out"));

    kmerSize =  0;
    targetSearchMode = 0;
    kmerScore.values = INT_MAX;
    alphabetSize = MultiParam<NuclAA<int>>(NuclAA<int>(21,5));
    maxSeqLen = MAX_SEQ_LEN; // 2^16
    maxResListLen = 300;
    sensitivity = 4;
    split = AUTO_SPLIT_DETECTION;
    splitMode = DETECT_BEST_DB_SPLIT;
    splitMemoryLimit = 0;
    diskSpaceLimit = 0;
    splitAA = false;
    spacedKmerPattern = "";
    localTmp = "";

    // search workflow
    numIterations = 1;
    startSens = 4;
    sensSteps = 1;
    exhaustiveSearch = false;
    exhaustiveFilterMsa = 0;
    strand = 1;
    orfFilter = 0;
    orfFilterSens = 2.0;
    orfFilterEval = 100;
    lcaSearch = false;

    greedyBestHits = false;

    threads = 1;
    compressed = WRITER_ASCII_MODE;
#ifdef OPENMP
    char * threadEnv = getenv("MMSEQS_NUM_THREADS");
    if (threadEnv != NULL) {
        threads = (int) Util::fast_atoi<unsigned int>(threadEnv);
    } else {
        #ifdef _SC_NPROCESSORS_ONLN
            threads = sysconf(_SC_NPROCESSORS_ONLN);
        #endif
        if(threads <= 1){
            threads = Util::omp_thread_count();
        }
    }

#endif
    compBiasCorrection = 1;
    compBiasCorrectionScale = 1.0;
    diagonalScoring = true;
    exactKmerMatching = 0;
    maskMode = 1;
    maskProb = 0.9;
    maskLowerCaseMode = 0;
    minDiagScoreThr = 15;
    spacedKmer = true;
    includeIdentity = false;
    alignmentMode = ALIGNMENT_MODE_FAST_AUTO;
    alignmentOutputMode = ALIGNMENT_OUTPUT_ALIGNMENT;
    evalThr = 0.001;
    covThr = 0.0;
    covMode = COV_MODE_BIDIRECTIONAL;
    seqIdMode = SEQ_ID_ALN_LEN;
    maxRejected = INT_MAX;
    maxAccept   = INT_MAX;
    seqIdThr = 0.0;
    alnLenThr = 0;
    altAlignment = 0;
    gapOpen = MultiParam<NuclAA<int>>(NuclAA<int>(11, 5));
    gapExtend = MultiParam<NuclAA<int>>(NuclAA<int>(1, 2));
#ifdef GAP_POS_SCORING
    gapPseudoCount = 10;
#endif
    zdrop = 40;
    addBacktrace = false;
    realign = false;
    clusteringMode = SET_COVER;
    singleStepClustering = false;
    clusterReassignment = 0;
    clusterSteps = 3;
    preloadMode = 0;
    scoreBias = 0.0;
    realignScoreBias = -0.2f;
    realignMaxSeqs = INT_MAX;
    correlationScoreWeight = 0.0;

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
    reuseLatest = false;
    // Clustering workflow
    removeTmpFiles = false;

    // indexdb
    checkCompatible = 0;
    searchType = SEARCH_TYPE_AUTO;
    indexSubset = INDEX_SUBSET_NORMAL;
    indexDbsuffix = "";

    // createdb
    createdbMode = SEQUENCE_SPLIT_MODE_HARD;
    shuffleDatabase = true;
    writeLookup = true;

    // format alignment
    formatAlignmentMode = FORMAT_ALIGNMENT_BLAST_TAB;
    outfmt = "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits";
    dbOut = false;

    // result2msa
    allowDeletion = false;
    summaryPrefix = "cl";
    skipQuery = false;
    msaFormatMode = FORMAT_MSA_FASTADB;
    // convertmsa
    identifierField = 1;

    // msa2profile
    matchMode = 0;
    matchRatio = 0.5;

    // result2profile
    evalProfile = evalThr;
    maskProfile = 1;
    filterMsa = 1;
    filterMaxSeqId = 0.9;
    qid = "0.0";           // default for minimum sequence identity with query
    qsc = -20.0f;        // default for minimum score per column with query
    covMSAThr = 0.0;           // default for minimum coverage threshold
    Ndiff = 1000;        // pick Ndiff most different sequences from alignment
    filterMinEnable = 0;
    wg = false;
    pcmode = PCMODE_SUBSTITUTION_SCORE;
    pca = MultiParam<PseudoCounts>(PseudoCounts(1.1, 1.4));
    pcb = MultiParam<PseudoCounts>(PseudoCounts(4.1, 5.8));

    // sequence2profile
    neff = 1.0;
    tau = 0.9;

    // logging
    verbosity = Debug::INFO;

    //extractorfs
    orfMinLength = 30;
    orfMaxLength = 32734;
    orfMaxGaps = INT_MAX;
    contigStartMode = 2;
    contigEndMode = 2;
    orfStartMode = 1;
    forwardFrames = "1,2,3";
    reverseFrames = "1,2,3";
    useAllTableStarts = false;
    translate = 0;
    createLookup = 0;

    // createdb
    identifierOffset = 0;
    dbType = 0;

    // split sequence
    sequenceOverlap = 0;
    sequenceSplitMode = Parameters::SEQUENCE_SPLIT_MODE_SOFT;
    headerSplitMode = 0;

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
    wrappedScoring = false;
    filterHits = false;
    sortResults = false;

    // filterDb
    filterColumn = 1;
    columnToTake = -1;
    filterColumnRegex = "^.*$";
    positiveFilter = true;
    filteringFile = "";
    filterExpression = "";
    trimToOneColumn = false;
    extractLines = 0;
    compOperator = "";
    compValue = 0;
    sortEntries = 0;
    beatsFirst = false;


    //besthitperset
    simpleBestHit = true;
    alpha = 1;
    shortOutput = false;
    aggregationMode = 0;

    // concatdbs
    preserveKeysB = false;
    takeLargerEntry = false;

    // diff
    useSequenceId = false;

    // prefixid
    prefix = "";
    tsvOut = false;

    // offsetalignment
    chainAlignment = 0;
    mergeQuery = 1;

    // tsv2db
    outputDbType = Parameters::DBTYPE_GENERIC_DB;

    // mergedbs
    mergePrefixes = "";
    mergeStopEmpty = false;

    // summarizetabs
    overlap = 0.0f;
    msaType = 2;

    // createclusterdb
    dbSuffixList = "_h";

    // summarize header
    headerType = Parameters::HEADER_TYPE_UNICLUST;

    // extractalignedregion
    extractMode = EXTRACT_TARGET;

    // convertkb
    kbColumns = "";

    // linearcluster
    kmersPerSequence = 21;
    kmersPerSequenceScale = MultiParam<NuclAA<float>>(NuclAA<float>(0.0, 0.2));
    includeOnlyExtendable = false;
    ignoreMultiKmer = false;
    hashShift = 67;
    pickNbest = 1;
    adjustKmerLength = false;
    resultDirection = Parameters::PARAM_RESULT_DIRECTION_TARGET;
    weightThr = 0.9;
    weightFile = "";

    // result2stats
    stat = "";

    // createtsv
    firstSeqRepr = false;
    fullHeader = false;
    idxSeqSrc = 0;
    targetTsvColumn = 1;

    // createtaxdb
    taxMappingFile = "";
    taxMappingMode = 0;
    taxDbMode = 1;
    ncbiTaxDump = "";

    // filtertaxdb, filtertaxseqdb
    taxonList = "";

    // view
    idList = "";
    idxEntryType = 0;

    // lca
    pickIdFrom = Parameters::EXTRACT_TARGET;

    // createsubdb
    subDbMode = Parameters::SUBDB_MODE_HARD;
    dbIdMode = Parameters::ID_MODE_KEYS;

    // tar2db
    tarInclude = ".*";
    tarExclude = "^$";

    // unpackdb
    unpackSuffix = "";
    unpackNameMode = Parameters::UNPACK_NAME_ACCESSION;

    lcaRanks = "";
    showTaxLineage = 0;
    // bin for all unclassified sequences
    // https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=12908
    // other sequences (plasmids, etc)
    // https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=28384
    blacklist = "12908:unclassified sequences,28384:other sequences";

    // aggregatetax
    majorityThr = 0.5;
    voteMode = AGG_TAX_MINUS_LOG_EVAL;

    // pairaln
    pairdummymode = PAIRALN_DUMMY_MODE_OFF;
    pairmode = PAIRALN_MODE_ALL_PER_SPECIES;

    // taxonomyreport
    reportMode = 0;

    // expandaln
    expansionMode = EXPAND_TRANSFER_EVALUE;
    expandFilterClusters = 0;

    // taxonomy
    taxonomySearchMode = Parameters::TAXONOMY_APPROX_2BLCA;
    taxonomyOutputMode = Parameters::TAXONOMY_OUTPUT_LCA;

    // substituion matrix
    substitutionMatrices = {
            {"nucleotide.out", nucleotide_out, nucleotide_out_len },
            {"blosum62.out", blosum62_out, blosum62_out_len },
            {"VTML80.out", VTML80_out, VTML80_out_len},
            {"VTML40.out", VTML40_out, VTML40_out_len},
            {"PAM30.out", PAM30_out, PAM30_out_len}
    };

    citations = {
            { CITATION_MMSEQS1,  "Hauser M, Steinegger M, Soding J: MMseqs software suite for fast and deep clustering and searching of large protein sequence sets. Bioinformatics, 32(9), 1323-1330 (2016)" },
            { CITATION_MMSEQS2,  "Steinegger M, Soding J: MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, 35(11), 1026-1028 (2017)" },
            { CITATION_UNICLUST, "Mirdita M, von den Driesch L, Galiez C, Martin M, Soding J, Steinegger M: Uniclust databases of clustered and deeply annotated protein sequences and alignments. Nucleic Acids Research 45(D1), D170-D176 (2017)" },
            { CITATION_LINCLUST, "Steinegger M, Soding J: Clustering huge protein sequence sets in linear time. Nature Communications, 9(1), 2542 (2018)" },
            { CITATION_PLASS,    "Steinegger M, Mirdita M, Soding J: Protein-level assembly increases protein sequence recovery from metagenomic samples manyfold. Nature Methods, 16(7), 603-606 (2019)" },
            { CITATION_SERVER,   "Mirdita M, Steinegger M, Soding J: MMseqs2 desktop and local web server app for fast, interactive sequence searches. Bioinformatics, 35(16), 2856–2858 (2019)" },
            { CITATION_TAXONOMY, "Mirdita M, Steinegger M, Breitwieser F, Soding J, Levy Karin E: Fast and sensitive taxonomic assignment to metagenomic contigs. Bioinformatics, btab184 (2021)" },
    };
}

std::vector<MMseqsParameter*> Parameters::combineList(const std::vector<MMseqsParameter*> &par1,
                                                      const std::vector<MMseqsParameter*> &par2) {
    std::vector<MMseqsParameter*> retVec;
    std::vector< std::vector<MMseqsParameter*>> tmp;
    tmp.push_back(par1);
    tmp.push_back(par2);
    for (size_t z = 0; z < tmp.size(); z++) {
        std::vector<MMseqsParameter*> &currPar = tmp[z];
        for (size_t i = 0; i < currPar.size(); i++) {
            bool addPar = true;
            for (size_t j = 0; j < retVec.size(); j++) {
                if (currPar[i]->uniqid == retVec[j]->uniqid){
                    addPar = false;
                }
            }
            if (addPar == true) {
                retVec.emplace_back(currPar[i]);
            }
        }
    }
    return retVec;
}

size_t Parameters::hashParameter(const std::vector<DbType> &dbtypes, const std::vector<std::string> &filenames, const std::vector<MMseqsParameter*> &par){
    std::string hashString;
    hashString.reserve(1024);

    struct stat stat_buf;
    bool stopAfterVariadic = false;
    for (size_t i = 0; i < filenames.size(); ++i){
        hashString.append(filenames[i]);
        if (stopAfterVariadic == false && i < dbtypes.size()) {
            const DbType& type = dbtypes[i];
            if (type.accessMode != DbType::ACCESS_MODE_INPUT) {
                continue;
            }

            if (type.specialType & DbType::VARIADIC) {
                stopAfterVariadic = true;
            }

            if (filenames[i] == "stdin") {
                continue;
            }

            if (::stat(filenames[i].c_str(), &stat_buf) == 0) {
                hashString.append(SSTR(stat_buf.st_size));
#ifdef __APPLE__
                hashString.append(SSTR(stat_buf.st_mtimespec.tv_sec));
#else
                hashString.append(SSTR(stat_buf.st_mtime));
#endif
                continue;
            }

            std::string index(filenames[i]);
            index.append(".index");
            if (::stat(index.c_str(), &stat_buf) == 0) {
                hashString.append(SSTR(stat_buf.st_size));
#ifdef __APPLE__
                hashString.append(SSTR(stat_buf.st_mtimespec.tv_sec));
#else
                hashString.append(SSTR(stat_buf.st_mtime));
#endif
                continue;
            }
        }
    }
    hashString.append(createParameterString(par));
    hashString.append(version);
    for (int i = 0; i < restArgc; ++i) {
        hashString.append(restArgv[i]);
    }
    return Util::hash(hashString.c_str(), hashString.size());
}


std::string Parameters::createParameterString(const std::vector<MMseqsParameter*> &par, bool wasSet) {
    std::ostringstream ss;
    for (size_t i = 0; i < par.size(); ++i) {
        // Never pass the MPI parameters along, they are passed by the environment
        if (par[i]->uniqid == PARAM_RUNNER_ID) {
            continue;
        }
        if(wasSet == true){
            if(par[i]->wasSet==false){
                continue;
            }
        }

        if (typeid(int) == par[i]->type){
            ss << par[i]->name << " ";
            ss << *((int *)par[i]->value) << " ";
        } else if (typeid(size_t) == par[i]->type){
            ss << par[i]->name << " ";
            ss << *((size_t *)par[i]->value) << " ";
        } else if (typeid(ByteParser) == par[i]->type) {
            ss << par[i]->name << " ";
            ss << ByteParser::format(*((size_t *)par[i]->value), 'a', 'h') << " ";
        } else if (typeid(float) == par[i]->type){
            ss << par[i]->name << " ";
            ss << *((float *)par[i]->value) << " ";
        } else if (typeid(double) == par[i]->type){
            ss << par[i]->name << " ";
            ss << *((double *)par[i]->value) << " ";
        } else if (PARAM_SUB_MAT.uniqid == par[i]->uniqid || PARAM_SEED_SUB_MAT.uniqid == par[i]->uniqid){
            MultiParam<NuclAA<std::string>>* param = ((MultiParam<NuclAA<std::string>>*) par[i]->value);
            MultiParam<NuclAA<std::string>> tmpPar(NuclAA<std::string>(
                BaseMatrix::unserializeName(param->values.aminoacid().c_str()),
                BaseMatrix::unserializeName(param->values.nucleotide().c_str())
            ));
            std::string value = MultiParam<NuclAA<std::string>>::format(tmpPar);
            // encode parameters as base64 if it contains whitespaces
            // whitespaces break parameters in the workflow shell scripts
            if (value.find_first_of(" \n\t[]{}^$?|.~!*<>&") != std::string::npos) {
                ss << par[i]->name << " b64:" << base64_encode(value.c_str(), value.size()) << " ";
            } else {
                ss << par[i]->name << " " << value << " ";
            }
        } else if (typeid(std::string) == par[i]->type){
            std::string& value = *((std::string *) par[i]->value);
            if (value != "") {
                // see above
                if (value.find_first_of(" \n\t[]{}^$?|.~!*<>&") != std::string::npos) {
                    ss << par[i]->name << " b64:" << base64_encode(value.c_str(), value.size()) << " ";
                } else {
                    ss << par[i]->name << " " << value << " ";
                }
            }
        } else if (typeid(bool) == par[i]->type){
            bool val = *((bool *)(par[i]->value));
            if (val == true){
                ss << par[i]->name << " 1 ";
            } else {
                ss << par[i]->name << " 0 ";
            }
        } else if (typeid(MultiParam<NuclAA<std::string>>) == par[i]->type) {
            ss << par[i]->name << " ";
            std::string value = MultiParam<NuclAA<std::string>>::format(*((MultiParam<NuclAA<std::string>> *) par[i]->value));
            // see above
            if (value.find_first_of(" \n\t[]{}^$?|.~!*<>&") != std::string::npos) {
                ss << par[i]->name << " b64:" << base64_encode(value.c_str(), value.size()) << " ";
            } else {
                ss << par[i]->name << " " << value << " ";
            }
        } else if (typeid(MultiParam<NuclAA<int>>) == par[i]->type) {
            ss << par[i]->name << " ";
            ss << MultiParam<NuclAA<int>>::format(*((MultiParam<NuclAA<int>> *) par[i]->value)) << " ";
        } else if (typeid(MultiParam<NuclAA<float>>) == par[i]->type) {
            ss << par[i]->name << " ";
            ss << MultiParam<NuclAA<float>>::format(*((MultiParam<NuclAA<float>> *) par[i]->value)) << " ";
        } else if (typeid(MultiParam<SeqProf<int>>) == par[i]->type) {
            ss << par[i]->name << " ";
            ss << MultiParam<SeqProf<int>>::format(*((MultiParam<SeqProf<int>> *) par[i]->value)) << " ";
        } else if (typeid(MultiParam<PseudoCounts>) == par[i]->type) {
            ss << par[i]->name << " ";
            ss << MultiParam<PseudoCounts>::format(*((MultiParam<PseudoCounts> *) par[i]->value)) << " ";
        } else {
            Debug(Debug::ERROR) << "Wrong parameter type. Please inform the developers!\n";
            EXIT(EXIT_FAILURE);
        }
    }

    return ss.str();
}

std::vector<MMseqsParameter*> Parameters::removeParameter(const std::vector<MMseqsParameter*> &par, const MMseqsParameter &x){
    std::vector<MMseqsParameter*> newParamList;
    for (size_t i =0; i < par.size(); i++) {
        MMseqsParameter* p = par[i];
        if (p->name != x.name)
            newParamList.push_back(p);
    }
    return newParamList;
}

void Parameters::overrideParameterDescription(MMseqsParameter& par, const char *description, const char *regex, int category) {
    if (description != NULL) {
        par.description = description;
    }
    if (regex != NULL) {
        par.regex = regex;
    }
    if (category != 0) {
        par.category = category;
    }
}

std::vector<int> Parameters::getOutputFormat(int formatMode, const std::string &outformat, bool &needSequences, bool &needBacktrace, bool &needFullHeaders,
                                             bool &needLookup, bool &needSource, bool &needTaxonomyMapping, bool &needTaxonomy) {
    std::vector<int> formatCodes;
    if (formatMode == Parameters::FORMAT_ALIGNMENT_SAM || formatMode == Parameters::FORMAT_ALIGNMENT_HTML) {
        needSequences = true;
        needBacktrace = true;
        return formatCodes;
    }
    std::vector<std::string> outformatSplit = Util::split(outformat, ",");
    int code = 0;
    for (size_t i = 0; i < outformatSplit.size(); ++i) {
        if(outformatSplit[i].compare("query") == 0){ code = Parameters::OUTFMT_QUERY;}
        else if (outformatSplit[i].compare("target") == 0){ code = Parameters::OUTFMT_TARGET;}
        else if (outformatSplit[i].compare("evalue") == 0){ code = Parameters::OUTFMT_EVALUE;}
        else if (outformatSplit[i].compare("gapopen") == 0){ code = Parameters::OUTFMT_GAPOPEN;}
        else if (outformatSplit[i].compare("pident") == 0){ code = Parameters::OUTFMT_PIDENT;}
        else if (outformatSplit[i].compare("fident") == 0){ code = Parameters::OUTFMT_FIDENT;}
        else if (outformatSplit[i].compare("nident") == 0){ code = Parameters::OUTFMT_NIDENT;}
        else if (outformatSplit[i].compare("qstart") == 0){ code = Parameters::OUTFMT_QSTART;}
        else if (outformatSplit[i].compare("qend") == 0){ code = Parameters::OUTFMT_QEND;}
        else if (outformatSplit[i].compare("qlen") == 0){ code = Parameters::OUTFMT_QLEN;}
        else if (outformatSplit[i].compare("tstart") == 0){ code = Parameters::OUTFMT_TSTART;}
        else if (outformatSplit[i].compare("tend") == 0){ code = Parameters::OUTFMT_TEND;}
        else if (outformatSplit[i].compare("tlen") == 0){ code = Parameters::OUTFMT_TLEN;}
        else if (outformatSplit[i].compare("alnlen") == 0){ code = Parameters::OUTFMT_ALNLEN;}
        else if (outformatSplit[i].compare("raw") == 0){ needSequences = true; code = Parameters::OUTFMT_RAW;}
        else if (outformatSplit[i].compare("bits") == 0){ code = Parameters::OUTFMT_BITS;}
        else if (outformatSplit[i].compare("cigar") == 0){ needBacktrace = true; code = Parameters::OUTFMT_CIGAR;}
        else if (outformatSplit[i].compare("qseq") == 0){ needSequences = true; code = Parameters::OUTFMT_QSEQ;}
        else if (outformatSplit[i].compare("tseq") == 0){ needSequences = true; code = Parameters::OUTFMT_TSEQ;}
        else if (outformatSplit[i].compare("qheader") == 0){ needFullHeaders = true; code = Parameters::OUTFMT_QHEADER;}
        else if (outformatSplit[i].compare("theader") == 0){ needFullHeaders = true; code = Parameters::OUTFMT_THEADER;}
        else if (outformatSplit[i].compare("qaln") == 0){ needBacktrace = true; needSequences = true; code = Parameters::OUTFMT_QALN;}
        else if (outformatSplit[i].compare("taln") == 0){ needBacktrace = true; needSequences = true; code = Parameters::OUTFMT_TALN;}
        else if (outformatSplit[i].compare("qframe") == 0){ code = Parameters::OUTFMT_QFRAME;}
        else if (outformatSplit[i].compare("tframe") == 0){ code = Parameters::OUTFMT_TFRAME;}
        else if (outformatSplit[i].compare("mismatch") == 0){ code = Parameters::OUTFMT_MISMATCH;}
        else if (outformatSplit[i].compare("qcov") == 0){ code = Parameters::OUTFMT_QCOV;}
        else if (outformatSplit[i].compare("tcov") == 0){ code = Parameters::OUTFMT_TCOV;}
        else if (outformatSplit[i].compare("qset") == 0){ needLookup = true; needSource = true; code = Parameters::OUTFMT_QSET;}
        else if (outformatSplit[i].compare("qsetid") == 0){ needLookup = true; needSource = true; code = Parameters::OUTFMT_QSETID;}
        else if (outformatSplit[i].compare("tset") == 0){ needLookup = true; code = Parameters::OUTFMT_TSET;}
        else if (outformatSplit[i].compare("tsetid") == 0){ needLookup = true; needSource = true; code = Parameters::OUTFMT_TSETID;}
        else if (outformatSplit[i].compare("taxid") == 0){ needTaxonomyMapping = true; code = Parameters::OUTFMT_TAXID;}
        else if (outformatSplit[i].compare("taxname") == 0){ needTaxonomyMapping = true; needTaxonomy = true; code = Parameters::OUTFMT_TAXNAME;}
        else if (outformatSplit[i].compare("taxlineage") == 0){ needTaxonomyMapping = true; needTaxonomy = true; code = Parameters::OUTFMT_TAXLIN;}
        else if (outformatSplit[i].compare("qorfstart") == 0){ code = Parameters::OUTFMT_QORFSTART;}
        else if (outformatSplit[i].compare("qorfend") == 0){ code = Parameters::OUTFMT_QORFEND;}
        else if (outformatSplit[i].compare("torfstart") == 0){ code = Parameters::OUTFMT_TORFSTART;}
        else if (outformatSplit[i].compare("torfend") == 0){ code = Parameters::OUTFMT_TORFEND;}
        else if (outformatSplit[i].compare("empty") == 0){ code = Parameters::OUTFMT_EMPTY;}
        else {
            Debug(Debug::ERROR) << "Format code " << outformatSplit[i] << " does not exist.";
            EXIT(EXIT_FAILURE);
        }
        formatCodes.push_back(code);
    }
    return formatCodes;
}
