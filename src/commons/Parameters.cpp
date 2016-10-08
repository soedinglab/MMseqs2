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

Parameters::Parameters():
PARAM_S(PARAM_S_ID,"-s", "Sensitivity","sensitivity: 1.0 faster; 4.0 fast; 6.1 default; 7.5 sensitive [1.0,7.5]", typeid(float), (void *) &sensitivity, "^[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_PREFILTER),
PARAM_K(PARAM_K_ID,"-k", "K-mer size", "k-mer size in the range [6,7] (0: set automatically to optimum)",typeid(int),  (void *) &kmerSize, "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PREFILTER|MMseqsParameter::COMMAND_CLUSTLINEAR),
PARAM_THREADS(PARAM_THREADS_ID,"--threads", "Threads", "number of cores used for the computation (uses all cores by default)",typeid(int), (void *) &threads, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_COMMON),
PARAM_ALPH_SIZE(PARAM_ALPH_SIZE_ID,"--alph-size", "Alphabet size", "alphabet size [2,21]",typeid(int),(void *) &alphabetSize, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PREFILTER|MMseqsParameter::COMMAND_CLUSTLINEAR),
PARAM_MAX_SEQ_LEN(PARAM_MAX_SEQ_LEN_ID,"--max-seq-len","Max. sequence length", "Maximum sequence length [1,32768]",typeid(int), (void *) &maxSeqLen, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_COMMON),
PARAM_PROFILE(PARAM_PROFILE_ID,"--profile", "Profile", "prefilter with query profiles (query DB must be a profile DB)",typeid(bool),(void *) &profile, "", MMseqsParameter::COMMAND_PREFILTER|MMseqsParameter::COMMAND_ALIGN|MMseqsParameter::COMMAND_PROFILE),
//PARAM_NUCL(PARAM_NUCL_ID,"--nucl", "Nucleotide","Nucleotide sequences input",typeid(bool),(void *) &nucl , ""),
PARAM_DIAGONAL_SCORING(PARAM_DIAGONAL_SCORING_ID,"--diag-score", "Diagonal Scoring", "use diagonal score for sorting the prefilter results [0,1]", typeid(int),(void *) &diagonalScoring, "^[0-1]{1}$", MMseqsParameter::COMMAND_PREFILTER),
PARAM_MIN_DIAG_SCORE(PARAM_MIN_DIAG_SCORE_ID,"--min-ungapped-score", "Minimum Diagonal score", "accept only matches with ungapped alignment score above this threshold", typeid(int),(void *) &minDiagScoreThr, "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PREFILTER),
PARAM_K_SCORE(PARAM_K_SCORE_ID,"--k-score", "K-score", "k-mer threshold for generating similar-k-mer lists",typeid(int),(void *) &kmerScore,  "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PREFILTER),
PARAM_MAX_SEQS(PARAM_MAX_SEQS_ID,"--max-seqs", "Max. results per query", "maximum result sequences per query",typeid(int),(void *) &maxResListLen, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_COMMON),
PARAM_SPLIT(PARAM_SPLIT_ID,"--split", "Split DB", "splits target set in n equally distributed chunks. In default the split is automatically set",typeid(int),(void *) &split,  "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PREFILTER),
PARAM_SPLIT_MODE(PARAM_SPLIT_MODE_ID,"--split-mode", "Split mode", "0: split target db; 1: split query db;  2: auto, depending on main memory",typeid(int),(void *) &splitMode,  "^[0-2]{1}$", MMseqsParameter::COMMAND_PREFILTER),
PARAM_SPLIT_AMINOACID(PARAM_SPLIT_AMINOACID_ID,"--split-aa", "Split by amino acid","Try to find the best split for the target database by amino acid count instead",typeid(bool), (void *) &splitAA, "$"),
PARAM_SUB_MAT(PARAM_SUB_MAT_ID,"--sub-mat", "Sub Matrix", "amino acid substitution matrix file",typeid(std::string),(void *) &scoringMatrixFile, "", MMseqsParameter::COMMAND_COMMON),
PARAM_NO_COMP_BIAS_CORR(PARAM_NO_COMP_BIAS_CORR_ID,"--comp-bias-corr", "Compositional bias","correct for locally biased amino acid composition [0,1]",typeid(int), (void *) &compBiasCorrection, "^[0-1]{1}$", MMseqsParameter::COMMAND_PREFILTER|MMseqsParameter::COMMAND_ALIGN|MMseqsParameter::COMMAND_PROFILE),
PARAM_SPACED_KMER_MODE(PARAM_SPACED_KMER_MODE_ID,"--spaced-kmer-mode", "Spaced Kmer", "0: use consecutive positions a k-mers; 1: use spaced k-mers",typeid(int), (void *) &spacedKmer,  "^[0-1]{1}", MMseqsParameter::COMMAND_PREFILTER),
PARAM_REMOVE_TMP_FILES(PARAM_REMOVE_TMP_FILES_ID, "--remove-tmp-files", "Remove Temporary Files" , "Delete temporary files", typeid(bool), (void *) &removeTmpFiles, ""),
PARAM_INCLUDE_IDENTITY(PARAM_INCLUDE_IDENTITY_ID,"--add-self-matches", "Include identical Seq. Id.","artificially add alignments of queries with themselves (for clustering)",typeid(bool), (void *) &includeIdentity, "", MMseqsParameter::COMMAND_PREFILTER|MMseqsParameter::COMMAND_ALIGN),
PARAM_RES_LIST_OFFSET(PARAM_RES_LIST_OFFSET_ID,"--offset-result", "Offset result","Offset result list",typeid(int), (void *) &resListOffset, "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PREFILTER),
// alignment
PARAM_ALIGNMENT_MODE(PARAM_ALIGNMENT_MODE_ID,"--alignment-mode", "Alignment mode", "What to compute: 0: automatic; 1: score+end_pos; 2:+start_pos+cov; 3: +seq.id",typeid(int), (void *) &alignmentMode, "^[0-4]{1}$", MMseqsParameter::COMMAND_ALIGN),
PARAM_E(PARAM_E_ID,"-e", "E-value threshold", "list matches below this E-value [0.0, inf]",typeid(float), (void *) &evalThr, "^([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)|[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_ALIGN),
PARAM_C(PARAM_C_ID,"-c", "Coverage threshold", "list matches above this fraction of aligned (covered) query residues",typeid(float), (void *) &covThr, "^0(\\.[0-9]+)?|1\\.0$", MMseqsParameter::COMMAND_ALIGN| MMseqsParameter::COMMAND_CLUSTLINEAR),
PARAM_MAX_REJECTED(PARAM_MAX_REJECTED_ID,"--max-rejected", "Max Reject", "maximum rejected alignments before alignment calculation for a query is aborted",typeid(int),(void *) &maxRejected, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_ALIGN),
PARAM_ADD_BACKTRACE(PARAM_ADD_BACKTRACE_ID, "-a", "Add backtrace", "add backtrace string (convert to alignments with mmseqs convertalis utility)", typeid(bool), (void *) &addBacktrace, "", MMseqsParameter::COMMAND_ALIGN),
PARAM_REALIGN(PARAM_REALIGN_ID, "--realign", "Realign hit", "compute more conservative, shorter alignments (scores and E-values not changed)", typeid(bool), (void *) &realign, "", MMseqsParameter::COMMAND_ALIGN),
PARAM_MIN_SEQ_ID(PARAM_MIN_SEQ_ID_ID,"--min-seq-id", "Seq. Id Threshold","list matches above this sequence identity (for clustering) [0.0,1.0]",typeid(float), (void *) &seqIdThr, "[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_ALIGN),
PARAM_FRAG_MERGE(PARAM_FRAG_MERGE_ID,"--frag-merge", "Detect fragments", "include alignments with cov > 0.95 and seq. id > 0.90 (for clustering)",typeid(bool), (void *) &fragmentMerge, "", MMseqsParameter::COMMAND_ALIGN),

// clustering
PARAM_CLUSTER_MODE(PARAM_CLUSTER_MODE_ID,"--cluster-mode", "Cluster mode", "0: Setcover, 1: connected component, 2: Greedy clustering by sequence length",typeid(int), (void *) &clusteringMode, "[0-2]{1}$", MMseqsParameter::COMMAND_CLUST),
PARAM_CASCADED(PARAM_CASCADED_ID,"--cascaded", "Cascaded clustering", "start the cascaded instead of simple clustering workflow",typeid(bool), (void *) &cascaded, "", MMseqsParameter::COMMAND_CLUST),
//affinity clustering
PARAM_MAXITERATIONS(PARAM_MAXITERATIONS_ID,"--max-iterations", "Max depth connected component", "maximum depth of breadth first search in connected component",typeid(int), (void *) &maxIteration,  "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_CLUST),
PARAM_SIMILARITYSCORE(PARAM_SIMILARITYSCORE_ID,"--similarity-type", "Similarity type", "type of score used for clustering [1:2]. 1=alignment score. 2=sequence identity ",typeid(int),(void *) &similarityScoreType,  "^[1-2]{1}$", MMseqsParameter::COMMAND_CLUST),
// logging
PARAM_V(PARAM_V_ID,"-v", "Verbosity","verbosity level: 0=nothing, 1: +errors, 2: +warnings, 3: +info",typeid(int), (void *) &verbosity, "^[0-3]{1}$", MMseqsParameter::COMMAND_COMMON),
// create profile (HMM, PSSM)
PARAM_PROFILE_TYPE(PARAM_PROFILE_TYPE_ID,"--profile-type", "Profile type", "0: HMM (HHsuite) 1: PSSM or 3: HMMER3",typeid(int),(void *) &profileMode,  "^[0-2]{1}$"),
// convertalignments
PARAM_FORMAT_MODE(PARAM_FORMAT_MODE_ID,"--format-mode", "Alignment Format", "output format BLAST TAB=0, PAIRWISE=1, or SAM=2 ", typeid(int), (void*) &formatAlignmentMode, "^[0-2]{1}$"),
// result2msa
PARAM_ALLOW_DELETION(PARAM_ALLOW_DELETION_ID,"--allow-deletion", "Allow Deletion", "allow deletions in a MSA", typeid(bool), (void*) &allowDeletion, ""),
PARAM_ADD_INTERNAL_ID(PARAM_ADD_INTERNAL_ID_ID,"--add-iternal-id", "Add internal id", "add internal id as comment to MSA", typeid(bool), (void*) &addInternalId, ""),
PARAM_COMPRESS_MSA(PARAM_COMPRESS_MSA_ID,"--compress", "Compress MSA", "create MSA in ca3m format", typeid(bool), (void*) &compressMSA, ""),
PARAM_SUMMARIZE_HEADER(PARAM_SUMMARIZE_HEADER_ID,"--summarize", "Summarize headers", "summarize cluster headers into a single header description", typeid(bool), (void*) &summarizeHeader, ""),
PARAM_SUMMARY_PREFIX(PARAM_SUMMARY_PREFIX_ID, "--summary-prefix", "Summary prefix","sets the cluster summary prefix",typeid(std::string),(void *) &summaryPrefix, ""),
PARAM_REPSEQ(PARAM_REPSEQ_ID,"--only-rep-seq","Representative sequence", "outputs a ffindex with the representative sequences", typeid(bool), (void*) &onlyRepSeq, ""),
// result2profile
PARAM_E_PROFILE(PARAM_E_PROFILE_ID,"--e-profile", "Profile e-value threshold", "includes sequences matches with < e-value thr. into the profile [>=0.0]", typeid(float), (void *) &evalProfile, "^([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)|([0-9]*(\\.[0-9]+)?)$", MMseqsParameter::COMMAND_PROFILE),
PARAM_FILTER_MAX_SEQ_ID(PARAM_FILTER_MAX_SEQ_ID_ID,"--max-seq-id", "Maximum sequence identity threshold", "reduce redundancy of output MSA using max. pairwise sequence identity [0.0,1.0]", typeid(float), (void*) &filterMaxSeqId, "^[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_PROFILE),
PARAM_FILTER_QSC(PARAM_FILTER_QSC_ID, "--qsc", "Minimum score per column", "reduce diversity of output MSAs using min. score per aligned residue with query sequences [-50.0,100.0]", typeid(float), (void*) &qsc, "^\\-*[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_PROFILE),
PARAM_FILTER_QID(PARAM_FILTER_QID_ID, "--qid", "Minimum seq. id.", "reduce diversity of output MSAs using min.seq. identity with query sequences [0.0,1.0]", typeid(float), (void*) &qid, "^[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_PROFILE),
PARAM_FILTER_COV(PARAM_FILTER_COV_ID, "--cov", "Minimum coverage", "filter output MSAs using min. fraction of query residues covered by matched sequences [0.0,1.0]", typeid(float), (void*) &cov, "^[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_PROFILE),
PARAM_FILTER_NDIFF(PARAM_FILTER_NDIFF_ID, "--diff", "Select n most diverse seqs", "filter MSAs by selecting most diverse set of sequences, keeping at least this many seqs in each MSA block of length 50 (0: filter off)", typeid(int), (void*) &Ndiff, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PROFILE),
PARAM_WG(PARAM_WG_ID, "--wg", "Use global sequence weighting", "use global sequence weighting for profile calculation", typeid(bool), (void*) &wg, "", MMseqsParameter::COMMAND_PROFILE),
PARAM_PCA(PARAM_PCA_ID, "--pca", "Pseudo count a", "pseudo count admixture strength", typeid(float), (void*) &pca, "^[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_PROFILE),
PARAM_PCB(PARAM_PCB_ID, "--pcb", "Pseudo count b", "pseudo counts: Neff at half of maximum admixture (0.0,infinity)", typeid(float), (void*) &pcb, "^[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_PROFILE),
//PARAM_FIRST_SEQ_REP_SEQ(PARAM_FIRST_SEQ_REP_SEQ_ID, "--first-seq-as-repr", "first sequence as respresentative", "Use the first sequence of the clustering result as representative sequence", typeid(bool), (void*) &firstSeqRepr, "", MMseqsParameter::COMMAND_PROFILE),
// result2stats
PARAM_STAT(PARAM_STAT_ID, "--stat", "Statistics to be computed", "can be one of: linecount, mean, doolittle, charges, seqlen.", typeid(std::string), (void*) &stat, ""),
// linearcluster
PARAM_KMER_PER_SEQ(PARAM_KMER_PER_SEQ_ID, "--kmer-per-seq", "Kmer per sequence", "kmer per sequence", typeid(int), (void*) &kmersPerSequence, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_CLUSTLINEAR),
// workflow
PARAM_RUNNER(PARAM_RUNNER_ID, "--mpi-runner", "Sets the MPI runner","use MPI on compute grid with this MPI command (e.g. \"mpirun -np 42\")",typeid(std::string),(void *) &runner, ""),
// search workflow
PARAM_NUM_ITERATIONS(PARAM_NUM_ITERATIONS_ID, "--num-iterations", "Number search iterations","Search iterations",typeid(int),(void *) &numIterations, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_PROFILE),
PARAM_START_SENS(PARAM_START_SENS_ID, "--start-sens", "Start sensitivity","start sensitivity",typeid(int),(void *) &startSens, "^[1-9]{1}$"),
PARAM_SENS_STEP_SIZE(PARAM_SENS_STEP_SIZE_ID, "--sens-step-size", "sensitivity step size","sensitivity step sizes",typeid(int),(void *) &sensStepSize, "^[1-9]{1}$"),
// Orfs
PARAM_ORF_MIN_LENGTH(PARAM_ORF_MIN_LENGTH_ID, "--min-length", "Min codons in orf", "minimum codon number in open reading frames",typeid(int),(void *) &orfMinLength, "^[1-9]{1}[0-9]*$"),
PARAM_ORF_MAX_LENGTH(PARAM_ORF_MAX_LENGTH_ID, "--max-length", "Max codons in length", "maximum codon number in open reading frames",typeid(int),(void *) &orfMaxLength, "^[1-9]{1}[0-9]*$"),
PARAM_ORF_MAX_GAP(PARAM_ORF_MAX_GAP_ID, "--max-gaps", "Max orf gaps", "maximum number of codons with gaps or unknown residues before an open reading frame is rejected",typeid(int),(void *) &orfMaxGaps, "^(0|[1-9]{1}[0-9]*)$"),
PARAM_ORF_SKIP_INCOMPLETE(PARAM_ORF_SKIP_INCOMPLETE_ID,"--skip-incomplete", "Skip incomplete orfs", "Skip orfs that have only an end or only a start codon or neither of those",typeid(bool),(void *) &orfSkipIncomplete, ""),
PARAM_ORF_LONGEST(PARAM_ORF_LONGEST_ID,"--longest-orf", "Find longest orf", "does the first found start codon start an orf (results in the longst possible orf)",typeid(bool),(void *) &orfLongest, ""),
PARAM_ORF_EXTENDMIN(PARAM_ORF_EXTENDMIN_ID,"--extend-min", "Extend short orfs", "if an orf would be rejected because of the min length threshold, allow it to be extended to the next stop codon",typeid(bool),(void *) &orfExtendMin, ""),
PARAM_ORF_FORWARD_FRAMES(PARAM_ORF_FORWARD_FRAMES_ID, "--forward-frames", "Forward Frames", "comma-seperated list of ORF frames on the forward strand to be extracted", typeid(std::string), (void *) &forwardFrames, ""),
PARAM_ORF_REVERSE_FRAMES(PARAM_ORF_REVERSE_FRAMES_ID, "--reverse-frames", "Reverse Frames", "comma-seperated list of ORF frames on the reverse strand to be extracted", typeid(std::string), (void *) &reverseFrames, ""),
// createdb
PARAM_USE_HEADER(PARAM_USE_HEADER_ID,"--use-fasta-header", "Use fasta header", "use the id parsed from the fasta header as the index key instead of using incrementing numeric identifiers",typeid(bool),(void *) &useHeader, ""),
PARAM_ID_OFFSET(PARAM_ID_OFFSET_ID, "--id-offset", "Offset of numeric ids", "numeric ids in index file are offset by this value ",typeid(int),(void *) &identifierOffset, "^(0|[1-9]{1}[0-9]*)$"),
PARAM_DONT_SPLIT_SEQ_BY_LEN(PARAM_DONT_SPLIT_SEQ_BY_LEN_ID,"--dont-split-seq-by-len", "Split Seq. by len", "Dont split sequences by --max-seq-len",typeid(bool),(void *) &splitSeqByLen, ""),
PARAM_USE_HEADER_FILE(PARAM_USE_HEADER_FILE_ID, "--use-header-file", "Use ffindex header", "use the ffindex header file instead of the body to map the entry keys",typeid(bool),(void *) &useHeaderFile, ""),
PARAM_GFF_TYPE(PARAM_GFF_TYPE_ID,"--gff-type", "GFF Type", "type in the GFF file to filter by",typeid(std::string),(void *) &gffType, ""),
PARAM_TRANSLATION_TABLE(PARAM_TRANSLATION_TABLE_ID,"--translation-table", "Translation Table", "1=CANONICAL, 2=VERT_MITOCHONDRIAL, 3=YEAST_MITOCHONDRIAL, 4=MOLD_MITOCHONDRIAL, 5=INVERT_MITOCHONDRIAL, 6=CILIATE, 9=FLATWORM_MITOCHONDRIAL, 10=EUPLOTID, 11=PROKARYOTE, 12=ALT_YEAST, 13=ASCIDIAN_MITOCHONDRIAL, 14=ALT_FLATWORM_MITOCHONDRIAL, 15=BLEPHARISMA, 16=CHLOROPHYCEAN_MITOCHONDRIAL, 21=TREMATODE_MITOCHONDRIAL, 22=SCENEDESMUS_MITOCHONDRIAL, 23=THRAUSTOCHYTRIUM_MITOCHONDRIAL, 24=PTEROBRANCHIA_MITOCHONDRIAL, 25=GRACILIBACTERI (Note gaps between tables)", typeid(int),(void *) &translationTable, "(^[1-6]{1}$|9|10|11|12|13|14|15|16|21|22|23|24|25)"),
PARAM_MIN_SEQUENCES(PARAM_MIN_SEQUENCES_ID,"--min-sequences", "Min Sequences", "minimum number of sequences a cluster may contain", typeid(int),(void *) &minSequences,"^[1-9]{1}[0-9]*$"),
PARAM_MAX_SEQUENCES(PARAM_MAX_SEQUENCES_ID,"--max-sequences", "Max Sequences", "maximum number of sequences a cluster may contain", typeid(int),(void *) &maxSequences,"^[1-9]{1}[0-9]*$"),
PARAM_HH_FORMAT(PARAM_HH_FORMAT_ID,"--hh-format", "HH format", "format entries to use with hhsuite (for singleton clusters)", typeid(bool), (void *) &hhFormat, ""),
// filterdb
PARAM_FILTER_COL(PARAM_FILTER_COL_ID,"--filter-column", "Filter column", "column", typeid(int),(void *) &filterColumn,"^[1-9]{1}[0-9]*$"),
PARAM_FILTER_REGEX(PARAM_FILTER_REGEX_ID,"--filter-regex", "Filter regex", "regex to select column (example float: [0-9]*(.[0-9]+)? int:[1-9]{1}[0-9])", typeid(std::string),(void *) &filterColumnRegex,"^.*$"),
PARAM_FILTER_POS(PARAM_FILTER_POS_ID,"--positive-filter", "Positive filter", "used in conjunction with --filter-file. If true, out  = in \\intersect filter ; if false, out = in - filter", typeid(bool),(void *) &positiveFilter,""),
PARAM_FILTER_FILE(PARAM_FILTER_FILE_ID,"--filter-file", "Filter file", "specify a file that contains the filtering elements", typeid(std::string),(void *) &filteringFile,""),
PARAM_MAPPING_FILE(PARAM_MAPPING_FILE_ID,"--mapping-file", "Mapping file", "specify a file that translates the keys of a result DB to new keys", typeid(std::string),(void *) &mappingFile,""),
PARAM_TRIM_TO_ONE_COL(PARAM_TRIM_TO_ONE_COL_ID,"--trim-to-one-column", "trim the results to one column","Output only the column specified by --filter-column.",typeid(bool), (void *) &trimToOneColumn, ""),
PARAM_EXTRACT_LINES(PARAM_EXTRACT_LINES_ID,"--extract-lines", "Extract n lines", "extract n lines of each entry.",typeid(int), (void *) &extractLines, "^[1-9]{1}[0-9]*$"),
PARAM_COMP_OPERATOR(PARAM_COMP_OPERATOR_ID,"--comparison-operator", "Numerical comparison operator", "compare numerically (le, ge, e) each entry to a comparison value.",typeid(std::string), (void *) &compOperator, ""),
PARAM_COMP_VALUE(PARAM_COMP_VALUE_ID,"--comparison-value", "Numerical comparison value", "compare numerically (le, ge, e) each entry to this comparison value.",typeid(float), (void *) &compValue, ""),
// concatdb
PARAM_PRESERVEKEYS(PARAM_PRESERVEKEYS_ID,"--preserve-keys", "Preserve the keys", "the keys of the two DB should be distinct, and they will be preserved in the concatenation.",typeid(bool), (void *) &preserveKeysB, ""),
// mergedbs
PARAM_MERGE_PREFIXES(PARAM_MERGE_PREFIXES_ID, "--prefixes", "Merge prefixes", "comma separated list of prefixes for each entry", typeid(std::string),(void *) &mergePrefixes,""),
// evaluationscores
PARAM_EVALUATION_ALLVSALL(PARAM_EVALUATION_ALLVSALL_ID, "-a", "All vs all","all cluster members vs all cluster members, otherwise: all against representative",typeid(bool),(void *) &allVsAll, ""),
PARAM_EVALUATION_RANDOMIZEDREPRESENTATIVE(PARAM_EVALUATION_RANDOMIZEDREPRESENTATIVE_ID, "-r", "random representative choice","Instead of first cluster member as representative choose a random one.",typeid(bool),(void *) &randomizedRepresentative, ""),
PARAM_EVALUATION_USE_SEQUENCEHEADER(PARAM_EVALUATION_USE_SEQUENCEHEADER_ID, "-h", "Use sequence db to map numerical ids back to UniProt Id","use sequence db to map numerical ids back to UniProt Id, should always be set except for UniRef",typeid(bool),(void *) &use_sequenceheader, ""),
PARAM_OVERLAP(PARAM_OVERLAP_ID, "--overlap", "Overlap", "maximum overlap", typeid(float), (void*) &overlap, "^[0-9]*(\\.[0-9]+)?$"),
PARAM_MSA_TYPE(PARAM_MSA_TYPE_ID,"--msa-type", "MSA type", "MSA Type: cA3M 0 or A3M 1", typeid(int), (void *) &msaType, "^[0-2]{1}$"),
// extractalignedregion
PARAM_EXTRACT_MODE(PARAM_EXTRACT_MODE_ID,"--extract-mode", "Extract mode", "Query 1, Target 2", typeid(int), (void *) &extractMode, "^[1-2]{1}$"),
// convertkb
PARAM_KB_COLUMNS(PARAM_KB_COLUMNS_ID, "--kb-columns", "UniprotKB Columns", "list of indices of UniprotKB columns to be extracted", typeid(std::string), (void *) &kbColumns, ""),
PARAM_COUNT_CHARACTER(PARAM_COUNT_CHARACTER_ID, "--count-char", "Count Char", "character to count", typeid(std::string), (void *) &countCharacter, "")
{
    
    // alignment
    align.push_back(PARAM_SUB_MAT);
    align.push_back(PARAM_ADD_BACKTRACE);
    align.push_back(PARAM_ALIGNMENT_MODE);
    align.push_back(PARAM_E);
    align.push_back(PARAM_MIN_SEQ_ID);
    align.push_back(PARAM_C);
    align.push_back(PARAM_MAX_SEQ_LEN);
    align.push_back(PARAM_MAX_SEQS);
    align.push_back(PARAM_NO_COMP_BIAS_CORR);
    //    alignment.push_back(PARAM_NUCL);
    align.push_back(PARAM_PROFILE);
    align.push_back(PARAM_REALIGN);
    align.push_back(PARAM_MAX_REJECTED);
    align.push_back(PARAM_FRAG_MERGE);
    align.push_back(PARAM_INCLUDE_IDENTITY);
    align.push_back(PARAM_THREADS);
    align.push_back(PARAM_V);
    
    // prefilter
    prefilter.push_back(PARAM_SUB_MAT);
    prefilter.push_back(PARAM_S);
    prefilter.push_back(PARAM_K);
    prefilter.push_back(PARAM_K_SCORE);
    prefilter.push_back(PARAM_ALPH_SIZE);
    prefilter.push_back(PARAM_MAX_SEQ_LEN);
    prefilter.push_back(PARAM_PROFILE);
    //    prefilter.push_back(PARAM_NUCL);
    prefilter.push_back(PARAM_MAX_SEQS);
    prefilter.push_back(PARAM_RES_LIST_OFFSET);
    prefilter.push_back(PARAM_SPLIT);
    prefilter.push_back(PARAM_SPLIT_MODE);
    prefilter.push_back(PARAM_NO_COMP_BIAS_CORR);
    prefilter.push_back(PARAM_DIAGONAL_SCORING);
    prefilter.push_back(PARAM_MIN_DIAG_SCORE);
    prefilter.push_back(PARAM_INCLUDE_IDENTITY);
    prefilter.push_back(PARAM_SPACED_KMER_MODE);
    prefilter.push_back(PARAM_THREADS);
    prefilter.push_back(PARAM_V);
    
    // clustering
    clust.push_back(PARAM_CLUSTER_MODE);
    clust.push_back(PARAM_MAX_SEQS);
    clust.push_back(PARAM_V);
    clust.push_back(PARAM_MAXITERATIONS);
    clust.push_back(PARAM_SIMILARITYSCORE);
    clust.push_back(PARAM_THREADS);
    
    // find orf
    onlyverbosity.push_back(PARAM_V);
    
    // convertprofiledb
    convertprofiledb.push_back(PARAM_SUB_MAT);
    convertprofiledb.push_back(PARAM_PROFILE_TYPE);
    convertprofiledb.push_back(PARAM_V);
    
    // create fasta
    createFasta.push_back(PARAM_USE_HEADER);
    createFasta.push_back(PARAM_V);
    
    // result2profile
    result2profile.push_back(PARAM_SUB_MAT);
    result2profile.push_back(PARAM_PROFILE);
    result2profile.push_back(PARAM_E_PROFILE);
    result2profile.push_back(PARAM_NO_COMP_BIAS_CORR);
    result2profile.push_back(PARAM_WG);
    result2profile.push_back(PARAM_FILTER_MAX_SEQ_ID);
    result2profile.push_back(PARAM_FILTER_QID);
    result2profile.push_back(PARAM_FILTER_QSC);
    result2profile.push_back(PARAM_FILTER_COV);
    result2profile.push_back(PARAM_FILTER_NDIFF);
    result2profile.push_back(PARAM_PCA);
    result2profile.push_back(PARAM_PCB);
    result2profile.push_back(PARAM_THREADS);
    result2profile.push_back(PARAM_V);
    //result2profile.push_back(PARAM_FIRST_SEQ_REP_SEQ);
    
    //result2stats
    result2stats.push_back(PARAM_STAT);
    
    // format alignment
    convertalignments.push_back(PARAM_FORMAT_MODE);
    //convertalignments.push_back(PARAM_THREADS);
    convertalignments.push_back(PARAM_V);
    // result2msa
    result2msa.push_back(PARAM_SUB_MAT);
    result2msa.push_back(PARAM_PROFILE);
    result2msa.push_back(PARAM_E_PROFILE);
    result2msa.push_back(PARAM_ALLOW_DELETION);
    result2msa.push_back(PARAM_ADD_INTERNAL_ID);
    result2msa.push_back(PARAM_NO_COMP_BIAS_CORR);
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
    result2msa.push_back(PARAM_REPSEQ);
    //result2msa.push_back(PARAM_FIRST_SEQ_REP_SEQ);
    
    // extract orf
    extractorfs.push_back(PARAM_ORF_MIN_LENGTH);
    extractorfs.push_back(PARAM_ORF_MAX_LENGTH);
    extractorfs.push_back(PARAM_ORF_MAX_GAP);
    extractorfs.push_back(PARAM_ORF_SKIP_INCOMPLETE);
    extractorfs.push_back(PARAM_ORF_LONGEST);
    extractorfs.push_back(PARAM_ORF_EXTENDMIN);
    extractorfs.push_back(PARAM_ORF_FORWARD_FRAMES);
    extractorfs.push_back(PARAM_ORF_REVERSE_FRAMES);
    extractorfs.push_back(PARAM_USE_HEADER);
    extractorfs.push_back(PARAM_ID_OFFSET);
    
    // splitdb
    splitdb.push_back(PARAM_SPLIT);
    splitdb.push_back(PARAM_SPLIT_AMINOACID);
    
    // create index
    createindex.push_back(PARAM_SUB_MAT);
    createindex.push_back(PARAM_K);
    createindex.push_back(PARAM_ALPH_SIZE);
    createindex.push_back(PARAM_MAX_SEQ_LEN);
    createindex.push_back(PARAM_SPLIT);
    createindex.push_back(PARAM_SPACED_KMER_MODE);
    createindex.push_back(PARAM_THREADS);
    createindex.push_back(PARAM_V);
    
    // create db
    createdb.push_back(PARAM_MAX_SEQ_LEN);
    createdb.push_back(PARAM_DONT_SPLIT_SEQ_BY_LEN);
    createdb.push_back(PARAM_USE_HEADER);
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
    gff2ffindex.push_back(PARAM_USE_HEADER);
    gff2ffindex.push_back(PARAM_ID_OFFSET);
    gff2ffindex.push_back(PARAM_V);
    
    searchworkflow = combineList(align, prefilter);
    searchworkflow = combineList(searchworkflow, result2profile);
    searchworkflow.push_back(PARAM_NUM_ITERATIONS);
    searchworkflow.push_back(PARAM_START_SENS);
    searchworkflow.push_back(PARAM_SENS_STEP_SIZE);
    searchworkflow.push_back(PARAM_RUNNER);
    
    
    clusteringWorkflow = combineList(prefilter, align);
    clusteringWorkflow = combineList(clusteringWorkflow, clust);
    clusteringWorkflow.push_back(PARAM_CASCADED);
    clusteringWorkflow.push_back(PARAM_REMOVE_TMP_FILES);
    clusteringWorkflow.push_back(PARAM_RUNNER);
    
    clusterUpdateSearch = removeParameter(searchworkflow,PARAM_MAX_SEQS);
    clusterUpdateClust = removeParameter(clusteringWorkflow,PARAM_MAX_SEQS);
    clusterUpdate = combineList(clusterUpdateSearch, clusterUpdateClust);
    
    // translate nucleotide
    translatenucs.push_back(PARAM_TRANSLATION_TABLE);
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
    filterDb.push_back(PARAM_MAPPING_FILE);
    filterDb.push_back(PARAM_THREADS);
    filterDb.push_back(PARAM_V);
    filterDb.push_back(PARAM_TRIM_TO_ONE_COL);
    filterDb.push_back(PARAM_EXTRACT_LINES);
    filterDb.push_back(PARAM_COMP_OPERATOR);
    filterDb.push_back(PARAM_COMP_VALUE);
    
    
    // swapreults
    swapresults.push_back(PARAM_SUB_MAT);
    swapresults.push_back(PARAM_MAX_SEQ_LEN);
    swapresults.push_back(PARAM_THREADS);
    swapresults.push_back(PARAM_V);
    
    // subtractdbs
    subtractdbs.push_back(PARAM_THREADS);
    subtractdbs.push_back(PARAM_E_PROFILE);
    subtractdbs.push_back(PARAM_V);
    
    //evaluationscores
    evaluationscores.push_back(PARAM_EVALUATION_ALLVSALL);
    evaluationscores.push_back(PARAM_EVALUATION_RANDOMIZEDREPRESENTATIVE);
    evaluationscores.push_back(PARAM_EVALUATION_USE_SEQUENCEHEADER);
    
    // clusthash
    clusthash.push_back(PARAM_SUB_MAT);
    clusthash.push_back(PARAM_ALPH_SIZE);
    clusthash.push_back(PARAM_MIN_SEQ_ID);
    clusthash.push_back(PARAM_MAX_SEQ_LEN);
    clusthash.push_back(PARAM_THREADS);
    clusthash.push_back(PARAM_V);
    
    // linearfilter
    linearfilter.push_back(PARAM_SUB_MAT);
    linearfilter.push_back(PARAM_ALPH_SIZE);
    linearfilter.push_back(PARAM_KMER_PER_SEQ);
    linearfilter.push_back(PARAM_K);
    linearfilter.push_back(PARAM_C);
    linearfilter.push_back(PARAM_MAX_SEQ_LEN);
    linearfilter.push_back(PARAM_THREADS);
    linearfilter.push_back(PARAM_V);
    // result2newick
    result2newick.push_back(PARAM_THREADS);
    result2newick.push_back(PARAM_V);
    
    // mergedbs
    mergedbs.push_back(PARAM_MERGE_PREFIXES);
    mergedbs.push_back(PARAM_V);
    
    // summarize
    summarizeheaders.push_back(PARAM_SUMMARY_PREFIX);
    summarizeheaders.push_back(PARAM_V);
    
    // diff
    diff.push_back(PARAM_THREADS);
    diff.push_back(PARAM_V);
    
    // prefixid
    prefixid.push_back(PARAM_MAPPING_FILE);
    prefixid.push_back(PARAM_THREADS);
    prefixid.push_back(PARAM_V);
    
    // annoate
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
    dbconcat.push_back(PARAM_PRESERVEKEYS);
    
    // extractalignedregion
    extractalignedregion.push_back(PARAM_EXTRACT_MODE);
    extractalignedregion.push_back(PARAM_THREADS);
    extractalignedregion.push_back(PARAM_V);
    
    // count
    count.push_back(PARAM_COUNT_CHARACTER);
    count.push_back(PARAM_THREADS);
    count.push_back(PARAM_V);
    
    // convertkb
    convertkb.push_back(PARAM_KB_COLUMNS);
    convertkb.push_back(PARAM_V);
    
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
        if(command.citations & CITATION_MMSEQS2) {
            ss << "Steinegger, M. & Soding, J. Sensitive protein sequence searching for the analysis of massive data sets. bioRxiv 079681 (2016)\n\n";
        }
        
        if(command.citations & CITATION_MMSEQS1) {
            ss << "Hauser, M., Steinegger, M. & Soding, J. MMseqs software suite for fast and deep clustering and searching of large protein sequence sets. Bioinformatics, 32(9), 1323-1330 (2016). \n\n";
        }
        
        if(command.citations & CITATION_UNICLUST) {
            ss << "Mirdita, M., von den Driesch, L., Galiez, C., Martin M., Soding J. & Steinegger M. Uniclust databases of clustered and deeply annotated protein sequences and alignments. In revision. \n\n";
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
        {"clustlinear", MMseqsParameter::COMMAND_CLUSTLINEAR},
        {"profile",  MMseqsParameter::COMMAND_PROFILE},
        {"misc",     MMseqsParameter::COMMAND_MISC},
        {"common",   MMseqsParameter::COMMAND_COMMON},
    };
    
    size_t maxWidth = 0;
    for(size_t i = 0; i < parameters.size(); i++) {
        maxWidth = std::max(strlen(parameters[i].name), maxWidth);
    }
    maxWidth+=2; // space in front of options
    
    // header
    ss << std::setprecision(1) << std::fixed;
    for(size_t i = 0; i < ARRAY_SIZE(categories); ++i) {
        bool categoryFound = false;
        for (size_t j = 0; j < parameters.size(); j++) {
            const MMseqsParameter &par = parameters[j];
            if (par.category & categories[i].category) {
                int others = (par.category ^ categories[i].category);
                if(others & outputFlag)
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
                if(par.category & categories[i].category){
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
                }
            }
            ss << "\n";
        }
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
                                 bool isVariadic,
                                 int outputFlag)
{
    std::vector<std::string> getFilename;
    std::vector<MMseqsParameter>& par = *command.params;
    size_t parametersFound = 0;
    for(int argIdx = 0; argIdx < argc; argIdx++ ){
        // it is a parameter if it starts with - or --
        if ((pargv[argIdx][0] == '-' && pargv[argIdx][1] == '-') || (pargv[argIdx][0] == '-')) {
            std::string parameter(pargv[argIdx]);
            bool hasUnrecognizedParameter = true;
            for(size_t parIdx = 0; parIdx < par.size(); parIdx++){
                if(parameter.compare(par[parIdx].name) == 0) {
                    if (typeid(bool) != par[parIdx].type && argIdx + 1 == argc) {
                        printUsageMessage(command, outputFlag);
                        Debug(Debug::ERROR) << "Missing argument " << par[parIdx].name << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                    
                    if (par[parIdx].wasSet) {
                        printUsageMessage(command, outputFlag);
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
                            printUsageMessage(command, outputFlag);
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
                            printUsageMessage(command, outputFlag);
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
                printUsageMessage(command, outputFlag);
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
            getFilename.push_back(pargv[argIdx]);
        }
    }
    
    Debug::setDebugLevel(verbosity);
#ifdef OPENMP
    omp_set_num_threads(threads);
#endif
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
    
    const size_t MAX_DB_PARAMETER = 5;
    
    if (requiredParameterCount > MAX_DB_PARAMETER) {
        Debug(Debug::ERROR) << "Use argv if you need more than " << MAX_DB_PARAMETER << " db parameters" << "\n";
        EXIT(EXIT_FAILURE);
    }
    
    if (getFilename.size() < requiredParameterCount){
        printUsageMessage(command, outputFlag);
        Debug(Debug::ERROR) << requiredParameterCount << " Database paths are required" << "\n";
        EXIT(EXIT_FAILURE);
    }
    
    switch (std::min(getFilename.size(), MAX_DB_PARAMETER)) {
        case 5:
            db5 = getFilename[4];
            db5Index = db5;
            db5Index.append(".index");
        case 4:
            db4 = getFilename[3];
            db4Index = db4;
            db4Index.append(".index");
        case 3:
            db3 = getFilename[2];
            db3Index = db3;
            db3Index.append(".index");
        case 2:
            db2 = getFilename[1];
            db2Index = db2;
            db2Index.append(".index");
        case 1:
            db1 = getFilename[0];
            db1Index = db1;
            db1Index.append(".index");
            break;
        default:
            // Do not abort execution if we expect a variable amount of parameters
            if(isVariadic)
                break;
        case 0:
            printUsageMessage(command, outputFlag);
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
    
#ifdef GIT_SHA1
#define str2(s) #s
#define str(s) str2(s)
    std::string gitHash(str(GIT_SHA1));
    ss << std::setw(maxWidth) << std::left  << "MMseqs Version:" << "\t" << gitHash << "\n";
#undef str
#undef str2
#endif
    
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

void Parameters::checkSaneEnvironment() {
    bool isInsane = false;
    
    char* mmdirStr = getenv("MMDIR");
    if (mmdirStr == NULL){
        Debug(Debug::ERROR) << "Please set the environment variable $MMDIR to your MMSEQS installation directory.\n";
        isInsane = true;
    }
    
    if(isInsane) {
        EXIT(EXIT_FAILURE);
    }
}

void Parameters::setDefaults() {
    scoringMatrixFile = "blosum62.out";
    
    kmerSize =  0;
    kmerScore = INT_MAX;
    alphabetSize = 21;
    maxSeqLen = MAX_SEQ_LEN; // 2^16
    maxResListLen = 300;
    sensitivity = 4;
    split = AUTO_SPLIT_DETECTION;
    splitMode = DETECT_BEST_DB_SPLIT;
    splitAA = false;
    querySeqType  = Sequence::AMINO_ACIDS;
    targetSeqType = Sequence::AMINO_ACIDS;
    
    // search workflow
    numIterations = 1;
    startSens = 4;
    sensStepSize = 1;
    
    
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
    minDiagScoreThr = 15;
    spacedKmer = true;
    profile = false;
    nucl = false;
    includeIdentity = false;
    alignmentMode = ALIGNMENT_MODE_FAST_AUTO;
    evalThr = 0.001;
    covThr = 0.0;
    fragmentMerge = false;
    maxRejected = INT_MAX;
    seqIdThr = 0.0;
    addBacktrace = false;
    realign = false;
    clusteringMode = SET_COVER;
    validateClustering = 0;
    cascaded = false;
    resListOffset = 0;
    
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
    
    // result2msa
    allowDeletion = false;
    addInternalId = false;
    compressMSA = false;
    summarizeHeader = false;
    summaryPrefix = "cl";
    onlyRepSeq = false;
    compressMSA = false;
    
    // result2profile
    evalProfile = evalThr;
    filterMaxSeqId = 0.9;
    qid = 0.0;           // default for minimum sequence identity with query
    qsc = -20.0f;        // default for minimum score per column with query
    cov = 0.0;           // default for minimum coverage threshold
    Ndiff = 100;         // pick Ndiff most different sequences from alignment
    wg = false;
    pca = 1.0;
    pcb = 1.5;
    useConsensus = true;
    
    // logging
    verbosity = Debug::INFO;
    
    //extractorfs
    orfMinLength = 1;
    orfMaxLength = INT_MAX;
    orfMaxGaps = INT_MAX;
    orfSkipIncomplete = false;
    orfLongest = false;
    orfExtendMin = false;
    forwardFrames = "1,2,3";
    reverseFrames = "1,2,3";
    
    // createdb
    useHeader = false;
    identifierOffset = 0;
    
    // convert2fasta
    useHeaderFile = false;
    
    // translate nucleotide
    translationTable = 1;
    
    // createseqfiledb
    minSequences = 1;
    maxSequences = INT_MAX;
    hhFormat = false;
    
    // filterDb
    filterColumn = 1;
    filterColumnRegex = "^.*$";
    positiveFilter = true;
    filteringFile = "";
    trimToOneColumn = false;
    extractLines = 0;
    compOperator = "";
    compValue = 0;
    
    // concatdbs
    preserveKeysB = false;
    
    // mergedbs
    mergePrefixes = "";
    
    // evaluationscores
    allVsAll = false;
    randomizedRepresentative = false;
    
    // summarizetabs
    overlap = 0.0f;
    msaType = 0;
    
    // extractalignedregion
    extractMode = EXTRACT_TARGET;
    
    // convertkb
    kbColumns = "";
    
    // linearcluster
    kmersPerSequence = 20;
    
    // count
    countCharacter = "\n";
}

std::vector<MMseqsParameter> Parameters::combineList(std::vector<MMseqsParameter> &par1,
                                                     std::vector<MMseqsParameter> &par2) {
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

std::string Parameters::createParameterString(std::vector<MMseqsParameter> &par) {
    std::stringstream ss;
    for (size_t i = 0; i < par.size(); i++) {
        if(typeid(int) == par[i].type ){
            ss << par[i].name << " ";
            ss << *((int *)par[i].value) << " ";
        } else if(typeid(float) == par[i].type ){
            ss << par[i].name << " ";
            ss << *((float *)par[i].value) << " ";
        }else if(typeid(std::string) == par[i].type ){
            if (*((std::string *) par[i].value)  != "")
            {
                ss << par[i].name << " ";
                ss << *((std::string *) par[i].value) << " ";
            }
        }else if (typeid(bool) == par[i].type){
            bool val = *((bool *)(par[i].value));
            if(val == true){
                ss << par[i].name << " ";
            }
        } else {
            Debug(Debug::ERROR) << "Wrong parameter type. Please inform the developers\n";
            EXIT(EXIT_FAILURE);
        }
    }
    return ss.str();
}

std::vector<MMseqsParameter> Parameters::removeParameter(std::vector<MMseqsParameter> par,MMseqsParameter x){
    std::vector<MMseqsParameter> newParamList;
    for (std::vector<MMseqsParameter>::iterator i = par.begin();i!=par.end();i++)
    {
        if (i->name != x.name)
            newParamList.push_back(*i);
    }
    return newParamList;
}
