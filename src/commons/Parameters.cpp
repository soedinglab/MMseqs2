#include "Parameters.h"
#include "Sequence.h"
#include "Debug.h"
#include "Util.h"

#include <iomanip>
#include <regex.h>
#include <string>
#include <sstream>
#include <climits>

#ifdef OPENMP
#include <omp.h>
#endif

Parameters::Parameters():
        PARAM_S(PARAM_S_ID,"-s", "Sensitivity","Sensitivity in the range [1.0:10.0]. From low (1.0) to high (10.0) sensitivity.", typeid(float), (void *) &sensitivity, "^[0-9]*(\\.[0-9]+)?$"),
        PARAM_K(PARAM_K_ID,"-k", "K-mer size", "k-mer size in the range [6,7]",typeid(int),  (void *) &kmerSize, "^[6-7]{1}$"),
        PARAM_THREADS(PARAM_THREADS_ID,"--threads", "Threads", "Number of cores used for the computation (uses all cores by default)",typeid(int), (void *) &threads, "^[1-9]{1}[0-9]*$"),
        PARAM_ALPH_SIZE(PARAM_ALPH_SIZE_ID,"--alph-size", "Alphabet size", "Amino acid alphabet size [2,21]",typeid(int),(void *) &alphabetSize, "^[1-9]{1}[0-9]*$"),
        PARAM_MAX_SEQ_LEN(PARAM_MAX_SEQ_LEN_ID,"--max-seq-len","Max. sequence length", "Maximum sequence length [1,32768]",typeid(int), (void *) &maxSeqLen, "^[1-9]{1}[0-9]*$"),
        PARAM_PROFILE(PARAM_PROFILE_ID,"--profile", "Profile", "HMM Profile input",typeid(bool),(void *) &profile, ""),
        PARAM_NUCL(PARAM_NUCL_ID,"--nucl", "Nucleotid","Nucleotide sequences input",typeid(bool),(void *) &nucl , ""),
        PARAM_DIAGONAL_SCORING(PARAM_DIAGONAL_SCORING_ID,"--diag-score", "Diagonal Scoring", "Use diagonal score for sorting the prefilter results [0,1]", typeid(int),(void *) &diagonalScoring, "^[0-1]{1}$"),
        PARAM_MIN_DIAG_SCORE(PARAM_MIN_DIAG_SCORE_ID,"--min-diag-score", "Minimum Diagonal score", "Accepts only hits with a ungapped diagonal score above the min score threshold", typeid(int),(void *) &minDiagScoreThr, "^[0-9]{1}[0-9]*$"),
        PARAM_K_SCORE(PARAM_K_SCORE_ID,"--k-score", "K-score", "Set the K-mer threshold for the K-mer generation",typeid(int),(void *) &kmerScore,  "^[1-9]{1}[0-9]*$"),
        PARAM_MAX_SEQS(PARAM_MAX_SEQS_ID,"--max-seqs", "Max. results per query", "Maximum result sequences per query",typeid(int),(void *) &maxResListLen, "^[1-9]{1}[0-9]*$"),
        PARAM_SPLIT(PARAM_SPLIT_ID,"--split", "Split DB", "Splits target set in n equally distributed chunks",typeid(int),(void *) &split,  "^[1-9]{1}[0-9]*$"),
        PARAM_SPLIT_MODE(PARAM_SPLIT_MODE_ID,"--split-mode", "Split mode", "MPI Option: Target set: 0 (low memory) or query set: 1 (faster but memory intensive)",typeid(int),(void *) &splitMode,  "^[0-1]{1}$"),
        PARAM_SPLIT_AMINOACID(PARAM_SPLIT_AMINOACID_ID,"--split-aa", "Split by amino acid","Try to find the best split for the target database by amino acid count instead",typeid(bool), (void *) &splitAA, "$"),
        PARAM_SUB_MAT(PARAM_SUB_MAT_ID,"--sub-mat", "Sub Matrix", "Amino acid substitution matrix file",typeid(std::string),(void *) &scoringMatrixFile, ""),
        PARAM_SEARCH_MODE(PARAM_SEARCH_MODE_ID,"--search-mode", "Search mode", "Search mode. Debug: 1 (debug) Normal: 2 (default)",typeid(int), (void *) &searchMode, "^[0-2]{1}$"),
        PARAM_NO_COMP_BIAS_CORR(PARAM_NO_COMP_BIAS_CORR_ID,"--comp-bias-corr", "Compositional bias","Switch off local amino acid composition bias correction[0,1]",typeid(int), (void *) &compBiasCorrection, "^[0-1]{1}$"),
        PARAM_SPACED_KMER_MODE(PARAM_SPACED_KMER_MODE_ID,"--spaced-kmer-mode", "Spaced Kmer", "Spaced kmers mode (use consecutive pattern). Disable: 0, Enable: 1",typeid(int), (void *) &spacedKmer,  "^[0-1]{1}" ),
        PARAM_REMOVE_TMP_FILES(PARAM_REMOVE_TMP_FILES_ID, "--remove-tmp-files", "Remove Temporary Files" , "Delete temporary files", typeid(bool), (void *) &removeTmpFiles, ""),
// alignment
        PARAM_ALIGNMENT_MODE(PARAM_ALIGNMENT_MODE_ID,"--alignment-mode", "Alignment mode", "Alignment mode 0=fastest based on parameters, 1=score; 2=score,cov,start/end pos; 3=score,cov,start/end pos,seq.id",typeid(int), (void *) &alignmentMode, "^[0-4]{1}$"),
        PARAM_E(PARAM_E_ID,"-e", "E-value threshold", "Maximum e-value[0.0,1.0]",typeid(float), (void *) &evalThr, "^[0-9]*(\\.[0-9]+)?$"),
        PARAM_C(PARAM_C_ID,"-c", "Coverage threshold", "Minimum alignment coverage [0.0,1.0]",typeid(float), (void *) &covThr, "^0(\\.[0-9]+)?|1\\.0$"),
        PARAM_FRAG_MERGE(PARAM_FRAG_MERGE_ID,"--frag-merge", "Detect fragments", "Add Hits with cov > 0.95 and seq. id > 0.90",typeid(bool), (void *) &fragmentMerge, ""),
        PARAM_MAX_REJECTED(PARAM_MAX_REJECTED_ID,"--max-rejected", "Max Reject", "Maximum rejected alignments before alignment calculation for a query is aborted",typeid(int),(void *) &maxRejected, "^[1-9]{1}[0-9]*$"),
        PARAM_ADD_BACKTRACE(PARAM_ADD_BACKTRACE_ID, "--add-backtrace", "Add backtrace", "Add backtrace string to results (M=Match, D=deletion, I=insertion)", typeid(bool), (void *) &addBacktrace, ""),
        PARAM_REALIGN(PARAM_REALIGN_ID, "--realign", "Realign hit", "Realign hit with conservative scoring scheme (keeps old evalue and score but overwrites alignment)", typeid(bool), (void *) &realign, ""),
        PARAM_MIN_SEQ_ID(PARAM_MIN_SEQ_ID_ID,"--min-seq-id", "Seq. Id Threshold","Minimum sequence identity of sequences in a cluster [0.0,1.0]",typeid(float), (void *) &seqIdThr, "[0-9]*(\\.[0-9]+)?$"),
// clustering
        PARAM_CLUSTER_MODE(PARAM_CLUSTER_MODE_ID,"--cluster-mode", "Cluster mode", "0 Setcover, 1 connected component, 2 Greedy clustering by sequence length",typeid(int), (void *) &clusteringMode, "[0-2]{1}$"),
        PARAM_CASCADED(PARAM_CASCADED_ID,"--cascaded", "Cascaded clustering", "Start the cascaded instead of simple clustering workflow",typeid(bool), (void *) &cascaded, ""),
//affinity clustering
        PARAM_MAXITERATIONS(PARAM_MAXITERATIONS_ID,"--max-iterations", "Max depth connected component", "Maximum depth of breadth first search in connected component",typeid(int), (void *) &maxIteration,  "^[1-9]{1}[0-9]*$"),
        PARAM_SIMILARITYSCORE(PARAM_SIMILARITYSCORE_ID,"--similarity-type", "Similarity type", "Type of score used for clustering [1:2]. 1=alignment score. 2=sequence identity ",typeid(int),(void *) &similarityScoreType,  "^[1-2]{1}$"),
// logging
        PARAM_V(PARAM_V_ID,"-v", "Verbosity","Verbosity level: 0=NOTHING, 1=ERROR, 2=WARNING, 3=INFO",typeid(int), (void *) &verbosity, "^[0-3]{1}$"),
// create profile (HMM, PSSM)
        PARAM_PROFILE_TYPE(PARAM_PROFILE_TYPE_ID,"--profile-type", "Profile type", "MPI Option: HMM 0 or PSSM",typeid(int),(void *) &profileMode,  "^[0-1]{1}$"),
// formatalignment
        PARAM_FORMAT_MODE(PARAM_FORMAT_MODE_ID,"--format-mode", "Alignment Format", "Output format BLAST TAB=0, PAIRWISE=1, or SAM=2 ", typeid(int), (void*) &formatAlignmentMode, "^[0-2]{1}$"),
// result2msa
        PARAM_ALLOW_DELETION(PARAM_ALLOW_DELETION_ID,"--allow-deletion", "Allow Deletion", "Allow deletions in a MSA", typeid(bool), (void*) &allowDeletion, ""),
        PARAM_ADD_INTERNAL_ID(PARAM_ADD_INTERNAL_ID_ID,"--add-iternal-id", "Add internal id", "Add internal id as comment to MSA", typeid(bool), (void*) &addInternalId, ""),
// result2profile
        PARAM_FILTER_MAX_SEQ_ID(PARAM_FILTER_MAX_SEQ_ID_ID,"--max-seq-id", "Maximum sequence identity threshold", "Maximum sequence identity with all other sequences in alignment [0.0,1.0]", typeid(float), (void*) &filterMaxSeqId, "^[0-9]*(\\.[0-9]+)?$"),
        PARAM_FILTER_QSC(PARAM_FILTER_QSC_ID, "--qsc", "Minimum score per column", "Minimum score per column with master sequence [-50.0,100.0]", typeid(float), (void*) &qsc, "^\\-*[0-9]*(\\.[0-9]+)?$"),
        PARAM_FILTER_QID(PARAM_FILTER_QID_ID, "--qid", "Minimum seq. id.", "Minimum sequence identity with master sequence [0.0,1.0]", typeid(float), (void*) &qid, "^[0-9]*(\\.[0-9]+)?$"),
        PARAM_FILTER_COV(PARAM_FILTER_COV_ID, "--cov", "Minimum coverage", "Minimum coverage with master sequence [0.0,1.0]", typeid(float), (void*) &cov, "^[0-9]*(\\.[0-9]+)?$"),
        PARAM_FILTER_NDIFF(PARAM_FILTER_NDIFF_ID, "--diff", "Select n most diverse seqs", "Filter MSAs by selecting most diverse set of sequences, keeping at least this many seqs in each MSA block of length 50", typeid(int), (void*) &Ndiff, "^[1-9]{1}[0-9]*$"),
        PARAM_WG(PARAM_WG_ID, "--wg", "Use global sequence weighting", "Use global sequence weighting for profile calculation", typeid(bool), (void*) &wg, ""),
        PARAM_PCA(PARAM_PCA_ID, "--pca", "Pseudo count a", "Overall pseudocount admixture", typeid(float), (void*) &pca, "^[0-9]*(\\.[0-9]+)?$"),
        PARAM_PCB(PARAM_PCB_ID, "--pcb", "Pseudo count b", "Admixture paramter b", typeid(float), (void*) &pcb, "^[0-9]*(\\.[0-9]+)?$"),
// workflow
        PARAM_RUNNER(PARAM_RUNNER_ID, "--mpi-runner", "Sets the MPI runner","Sets the MPI runner",typeid(std::string),(void *) &runner, ""),
// search workflow
        PARAM_NUM_ITERATIONS(PARAM_NUM_ITERATIONS_ID, "--num-iterations", "Number search iterations","Search iterations",typeid(int),(void *) &numIterations, "^[1-9]{1}[0-9]*$"),
        PARAM_START_SENS(PARAM_START_SENS_ID, "--start-sens", "Start sensitivity","Start sensitivity",typeid(int),(void *) &startSens, "^[1-9]{1}$"),
        PARAM_SENS_STEP_SIZE(PARAM_SENS_STEP_SIZE_ID, "--sens-step-size", "Sensitivity step size","Sensitivity step sizes",typeid(int),(void *) &sensStepSize, "^[1-9]{1}$"),
        PARAM_USE_INDEX(PARAM_USE_INDEX_ID, "--use-index", "Use index","Use precomputed index for prefilter",typeid(bool),(void *) &useIndex, ""),
// Orfs
        PARAM_ORF_MIN_LENGTH(PARAM_ORF_MIN_LENGTH_ID, "--min-length", "Min codons in orf", "Minimum codon number in open reading frames",typeid(int),(void *) &orfMinLength, "^[1-9]{1}[0-9]*$"),
        PARAM_ORF_MAX_LENGTH(PARAM_ORF_MAX_LENGTH_ID, "--max-length", "Max codons in length", "Maximum codon number in open reading frames",typeid(int),(void *) &orfMaxLength, "^[1-9]{1}[0-9]*$"),
        PARAM_ORF_MAX_GAP(PARAM_ORF_MAX_GAP_ID, "--max-gaps", "Max orf gaps", "Maximum number of codons with gaps or unknown residues before an open reading frame is rejected",typeid(int),(void *) &orfMaxGaps, "^(0|[1-9]{1}[0-9]*)$"),
        PARAM_ORF_SKIP_INCOMPLETE(PARAM_ORF_SKIP_INCOMPLETE_ID,"--skip-incomplete", "Skip incomplete orfs", "Skip orfs that have only an end or only a start codon or neither of those",typeid(bool),(void *) &orfSkipIncomplete, ""),
        PARAM_ORF_LONGEST(PARAM_ORF_LONGEST_ID,"--longest-orf", "Find longest orf", "Does the first found start codon start an orf (results in the longst possible orf)",typeid(bool),(void *) &orfLongest, ""),
        PARAM_ORF_EXTENDMIN(PARAM_ORF_EXTENDMIN_ID,"--extend-min", "Extend short orfs", "If an orf would be rejected because of the min length threshold, allow it to be extended to the next stop codon",typeid(bool),(void *) &orfExtendMin, ""),
        PARAM_ORF_FORWARD_FRAMES(PARAM_ORF_FORWARD_FRAMES_ID, "--forward-frames", "Forward Frames", "Comma-seperated list of ORF frames on the forward strand to be extracted", typeid(std::string), (void *) &forwardFrames, ""),
        PARAM_ORF_REVERSE_FRAMES(PARAM_ORF_REVERSE_FRAMES_ID, "--reverse-frames", "Reverse Frames", "Comma-seperated list of ORF frames on the reverse strand to be extracted", typeid(std::string), (void *) &reverseFrames, ""),
// createdb
        PARAM_USE_HEADER(PARAM_USE_HEADER_ID,"--use-fasta-header", "Use fasta header", "Use the id parsed from the fasta header as the index key instead of using incrementing numeric identifiers",typeid(bool),(void *) &useHeader, ""),
        PARAM_ID_OFFSET(PARAM_ID_OFFSET_ID, "--id-offset", "Offset of numeric ids", "Numeric ids in index file are offset by this value ",typeid(int),(void *) &identifierOffset, "^(0|[1-9]{1}[0-9]*)$"),
        PARAM_DONT_SPLIT_SEQ_BY_LEN(PARAM_DONT_SPLIT_SEQ_BY_LEN_ID,"--dont-split-seq-by-len", "Split Seq. by len", "Dont split sequences by --max-seq-len",typeid(bool),(void *) &splitSeqByLen, ""),
        PARAM_USE_HEADER_FILE(PARAM_USE_HEADER_FILE_ID, "--use-header-file", "Use ffindex header", "Use the ffindex header file instead of the body to map the entry keys",typeid(bool),(void *) &useHeaderFile, ""),
        PARAM_GFF_TYPE(PARAM_GFF_TYPE_ID,"--gff-type", "GFF Type", "Type in the GFF file to filter by",typeid(std::string),(void *) &gffType, ""),
        PARAM_TRANSLATION_TABLE(PARAM_TRANSLATION_TABLE_ID,"--translation-table", "Translation Table", "1=CANONICAL, 2=VERT_MITOCHONDRIAL, 3=YEAST_MITOCHONDRIAL, 4=MOLD_MITOCHONDRIAL, 5=INVERT_MITOCHONDRIAL, 6=CILIATE, 9=FLATWORM_MITOCHONDRIAL, 10=EUPLOTID, 11=PROKARYOTE, 12=ALT_YEAST, 13=ASCIDIAN_MITOCHONDRIAL, 14=ALT_FLATWORM_MITOCHONDRIAL, 15=BLEPHARISMA, 16=CHLOROPHYCEAN_MITOCHONDRIAL, 21=TREMATODE_MITOCHONDRIAL, 22=SCENEDESMUS_MITOCHONDRIAL, 23=THRAUSTOCHYTRIUM_MITOCHONDRIAL, 24=PTEROBRANCHIA_MITOCHONDRIAL, 25=GRACILIBACTERI (Note gaps between tables)", typeid(int),(void *) &translationTable, "(^[1-6]{1}$|9|10|11|12|13|14|15|16|21|22|23|24|25)"),
        PARAM_MIN_SEQUENCES(PARAM_MIN_SEQUENCES_ID,"--min-sequences", "Min Sequences", "Minimum number of sequences a cluster may contain", typeid(int),(void *) &minSequences,"^[1-9]{1}[0-9]*$"),
        PARAM_FILTER_COL(PARAM_FILTER_COL_ID,"--filter-column", "Filter column", "Column", typeid(int),(void *) &filterColumn,"^[1-9]{1}[0-9]*$"),
        PARAM_FILTER_REGEX(PARAM_FILTER_REGEX_ID,"--filter-regex", "Filter regex", "Regex to select column (example float: [0-9]*(.[0-9]+)? int:[1-9]{1}[0-9])", typeid(std::string),(void *) &filterColumnRegex,"^.*$"),
// evaluationscores
                PARAM_EVALUATION_ALLVSALL(PARAM_EVALUATION_ALLVSALL_ID, "-a", "All vs all","All cluster members vs all cluster members, otherwise: all against representative",typeid(bool),(void *) &allVsAll, ""),
        PARAM_EVALUATION_RANDOMIZEDREPRESENTATIVE(PARAM_EVALUATION_RANDOMIZEDREPRESENTATIVE_ID, "-r", "Random representative choice","Instead of first cluster member as representative choose a random one.",typeid(bool),(void *) &randomizedRepresentative, "")





{
    // alignment
    alignment.push_back(PARAM_SUB_MAT);
    alignment.push_back(PARAM_ALIGNMENT_MODE);
    alignment.push_back(PARAM_E);
    alignment.push_back(PARAM_C);
    alignment.push_back(PARAM_FRAG_MERGE);
    alignment.push_back(PARAM_NO_COMP_BIAS_CORR);
    alignment.push_back(PARAM_MIN_SEQ_ID);
    alignment.push_back(PARAM_MAX_SEQ_LEN);
    alignment.push_back(PARAM_MAX_SEQS);
    alignment.push_back(PARAM_MAX_REJECTED);
    alignment.push_back(PARAM_NUCL);
    alignment.push_back(PARAM_PROFILE);
    alignment.push_back(PARAM_ADD_BACKTRACE);
    alignment.push_back(PARAM_REALIGN);
    alignment.push_back(PARAM_THREADS);
    alignment.push_back(PARAM_V);

    // prefilter
    prefilter.push_back(PARAM_SUB_MAT);
    prefilter.push_back(PARAM_S);
    prefilter.push_back(PARAM_K);
    prefilter.push_back(PARAM_K_SCORE);
    prefilter.push_back(PARAM_ALPH_SIZE);
    prefilter.push_back(PARAM_MAX_SEQ_LEN);
    prefilter.push_back(PARAM_PROFILE);
    prefilter.push_back(PARAM_NUCL);
    prefilter.push_back(PARAM_MAX_SEQS);
    prefilter.push_back(PARAM_SPLIT);
    prefilter.push_back(PARAM_SPLIT_MODE);
    prefilter.push_back(PARAM_SEARCH_MODE);
    prefilter.push_back(PARAM_NO_COMP_BIAS_CORR);
    prefilter.push_back(PARAM_DIAGONAL_SCORING);
    prefilter.push_back(PARAM_MIN_DIAG_SCORE);
    prefilter.push_back(PARAM_SPACED_KMER_MODE);
    prefilter.push_back(PARAM_THREADS);
    prefilter.push_back(PARAM_V);

    // clustering
    clustering.push_back(PARAM_CLUSTER_MODE);
    clustering.push_back(PARAM_MAX_SEQS);
    clustering.push_back(PARAM_V);
    clustering.push_back(PARAM_MAXITERATIONS);
    clustering.push_back(PARAM_SIMILARITYSCORE);
    clustering.push_back(PARAM_THREADS);

    // find orf
    onlyverbosity.push_back(PARAM_V);

    // createprofiledb
    createprofiledb.push_back(PARAM_SUB_MAT);
    createprofiledb.push_back(PARAM_PROFILE_TYPE);
    createprofiledb.push_back(PARAM_V);

    // result2profile
    result2profile.push_back(PARAM_SUB_MAT);
    result2profile.push_back(PARAM_PROFILE);
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

    // format alignment
    formatalignment.push_back(PARAM_FORMAT_MODE);
    //formatalignment.push_back(PARAM_THREADS);
    formatalignment.push_back(PARAM_V);
    // result2msa
    result2msa.push_back(PARAM_SUB_MAT);
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

    // extract orf
    extractorf.push_back(PARAM_ORF_MIN_LENGTH);
    extractorf.push_back(PARAM_ORF_MAX_LENGTH);
    extractorf.push_back(PARAM_ORF_MAX_GAP);
    extractorf.push_back(PARAM_ORF_SKIP_INCOMPLETE);
    extractorf.push_back(PARAM_ORF_LONGEST);
    extractorf.push_back(PARAM_ORF_EXTENDMIN);
    extractorf.push_back(PARAM_ORF_FORWARD_FRAMES);
    extractorf.push_back(PARAM_ORF_REVERSE_FRAMES);
    extractorf.push_back(PARAM_USE_HEADER);
    extractorf.push_back(PARAM_ID_OFFSET);

    // splitffindex
    splitffindex.push_back(PARAM_SPLIT);
    splitffindex.push_back(PARAM_SPLIT_AMINOACID);

    // create index
    createindex.push_back(PARAM_SUB_MAT);
    createindex.push_back(PARAM_K);
    createindex.push_back(PARAM_ALPH_SIZE);
    createindex.push_back(PARAM_MAX_SEQ_LEN);
    createindex.push_back(PARAM_SPLIT);
    createindex.push_back(PARAM_SEARCH_MODE);
    createindex.push_back(PARAM_SPACED_KMER_MODE);
    createindex.push_back(PARAM_V);

    // create db
    createdb.push_back(PARAM_MAX_SEQ_LEN);
    createdb.push_back(PARAM_DONT_SPLIT_SEQ_BY_LEN);
    createdb.push_back(PARAM_USE_HEADER);
    createdb.push_back(PARAM_ID_OFFSET);
    createdb.push_back(PARAM_V);

    // rebuildfasta
    rebuildfasta.push_back(PARAM_USE_HEADER_FILE);
    rebuildfasta.push_back(PARAM_V);

    // gff2ffindex
    gff2ffindex.push_back(PARAM_GFF_TYPE);
    gff2ffindex.push_back(PARAM_USE_HEADER);
    gff2ffindex.push_back(PARAM_ID_OFFSET);
    gff2ffindex.push_back(PARAM_V);

    searchworkflow = combineList(alignment, prefilter);
    searchworkflow = combineList(searchworkflow, result2profile);
    searchworkflow.push_back(PARAM_NUM_ITERATIONS);
    searchworkflow.push_back(PARAM_START_SENS);
    searchworkflow.push_back(PARAM_SENS_STEP_SIZE);
    searchworkflow.push_back(PARAM_USE_INDEX);
    searchworkflow.push_back(PARAM_RUNNER);

    clusteringWorkflow = combineList(prefilter, alignment);
    clusteringWorkflow = combineList(clusteringWorkflow, clustering);
    clusteringWorkflow.push_back(PARAM_CASCADED);
    clusteringWorkflow.push_back(PARAM_REMOVE_TMP_FILES);
    clusteringWorkflow.push_back(PARAM_RUNNER);

    clusterUpdate = combineList(alignment, prefilter);
    clusterUpdate = combineList(clusterUpdate, clustering);

    // translate nucleotide
    translateNucleotide.push_back(PARAM_TRANSLATION_TABLE);
    translateNucleotide.push_back(PARAM_V);

    // addSequences
    addSequences.push_back(PARAM_MIN_SEQUENCES);
    addSequences.push_back(PARAM_V);

    // filterDb
    filterDb.push_back(PARAM_FILTER_COL);
    filterDb.push_back(PARAM_FILTER_REGEX);
    filterDb.push_back(PARAM_THREADS);
    filterDb.push_back(PARAM_V);

    // swapreults
    swapresults.push_back(PARAM_SPLIT);
    swapresults.push_back(PARAM_V);

    // substractresult
    substractresult.push_back(PARAM_THREADS);
    substractresult.push_back(PARAM_V);

    //evaluationscores
    evaluationscores.push_back(PARAM_EVALUATION_ALLVSALL);
    evaluationscores.push_back(PARAM_EVALUATION_RANDOMIZEDREPRESENTATIVE);
    // detectredundancy
    detectredundancy.push_back(PARAM_SUB_MAT);
    detectredundancy.push_back(PARAM_ALPH_SIZE);
    detectredundancy.push_back(PARAM_MIN_SEQ_ID);
    detectredundancy.push_back(PARAM_MAX_SEQ_LEN);
    detectredundancy.push_back(PARAM_THREADS);
    detectredundancy.push_back(PARAM_V);

    checkSaneEnvironment();
    setDefaults();
}

void Parameters::printUsageMessage(const std::string &programUsageHeader,
                                   const std::vector<MMseqsParameter> &parameters){
    std::ostringstream ss;
    ss << programUsageHeader << std::endl;

    size_t maxWidth = 0;
    for(size_t i = 0; i < parameters.size(); i++) {
        maxWidth = std::max(strlen(parameters[i].name), maxWidth);
    }

    // header
    ss << std::left << std::setw(maxWidth) << "Parameter Name" << "\t";
    ss << std::left << std::setw(16) << "Type & Value" << "\t";
    ss << "Description" << std::endl;

    // body
    for(size_t i = 0; i < parameters.size(); i++) {
        const MMseqsParameter& par = parameters[i];
        ss << std::left << std::setw(maxWidth) << par.name << "\t";
        ss << std::boolalpha << std::left << std::setw(6);
        if (par.type == typeid(int)) {
            ss << "[int:" << std::right << std::setw(10) << *((int *) par.value) << "]";
        }else if(par.type == typeid(float)){
            ss << "[real:" << std::right << std::setw(10) << *((float *) par.value) << "]";
        }else if(par.type == typeid(bool)) {
            ss << "[bool:" << std::right << std::setw(10) << *((bool *) par.value) << "]";
        }else if (par.type == typeid(std::string)) {
            std::string& out = *((std::string *) par.value);
            size_t index = out.find(mmdir);
            if(index != std::string::npos) {
                out.replace(index, mmdir.length(), "$MMDIR");
            }
            ss << "[text:" << std::right << std::setw(10) << out << "]";
        }
        ss << "\t";
        ss << std::left << std::setw(60) << par.description << std::endl;
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
                                 const std::string &programUsageHeader,
                                 std::vector<MMseqsParameter> &par,
                                 size_t requiredParameterCount,
                                 bool printPar,
                                 bool isVariadic)
{
    std::vector<std::string> getFilename;
    size_t parametersFound = 0;
    for(int argIdx = 0; argIdx < argc; argIdx++ ){
        // it is a parameter if it starts with - or --
        if ((pargv[argIdx][0] == '-' && pargv[argIdx][1] == '-') || (pargv[argIdx][0] == '-')) {
            std::string parameter(pargv[argIdx]);
            bool hasUnrecognizedParameter = true;
            for(size_t parIdx = 0; parIdx < par.size(); parIdx++){
                if(parameter.compare(par[parIdx].name) == 0) {
                    if (typeid(bool) != par[parIdx].type && argIdx + 1 == argc) {
                        printUsageMessage(programUsageHeader, par);
                        Debug(Debug::ERROR) << "Missing argument " << par[parIdx].name << "\n";
                        EXIT(EXIT_FAILURE);
                    }

                    if (par[parIdx].wasSet) {
                        printUsageMessage(programUsageHeader, par);
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
                            printUsageMessage(programUsageHeader, par);
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
                            printUsageMessage(programUsageHeader, par);
                            Debug(Debug::ERROR) << "Error in argument " << par[parIdx].name << "\n";
                            EXIT(EXIT_FAILURE);
                        }else{
                            *((float *) par[parIdx].value) = atof(pargv[argIdx+1]);
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
                printUsageMessage(programUsageHeader, par);
                Debug(Debug::ERROR) << "Unrecognized parameter " << parameter << "\n";
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
        printUsageMessage(programUsageHeader, par);
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
            // Do not abort execution if we exect a variable amount of parameters
            if(isVariadic)
                break;
        case 0:
            printUsageMessage(programUsageHeader, par);
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
    mmdir = getenv("MMDIR");
    scoringMatrixFile = mmdir;
    scoringMatrixFile.append("/data/blosum62.out");

    kmerSize =  7;
    kmerScore = INT_MAX;
    alphabetSize = 21;
    maxSeqLen = 32000; // 2^15
    maxResListLen = 300;
    sensitivity = 4;
    split = 1;
    splitMode = TARGET_DB_SPLIT;
    splitAA = false;
    querySeqType  = Sequence::AMINO_ACIDS;
    targetSeqType = Sequence::AMINO_ACIDS;

    // search workflow
    numIterations = 1;
    startSens = 4;
    sensStepSize = 1;
    useIndex = false;


    threads = 1;
#ifdef OPENMP
    threads = omp_get_max_threads();
#endif
    compBiasCorrection = 1;
    diagonalScoring = 1;
    minDiagScoreThr = 30;
    spacedKmer = true;
    searchMode = SEARCH_LOCAL_FAST;
    profile = false;
    nucl = false;

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

    // affinity clustering
    maxIteration=1000;
    similarityScoreType=APC_SEQID;

    // workflow
    runner = "";

    // Clustering workflow
    removeTmpFiles = false;

    // createprofiledb
    profileMode = PROFILE_MODE_HMM;

    // createdb
    splitSeqByLen = true;

    // format alignment
    formatAlignmentMode = FORMAT_ALIGNMENT_BLAST_TAB;

    // result2msa
    allowDeletion = false;
    addInternalId = false;
    // result2profile
    filterMaxSeqId = 0.9;
    qid = 0.0;           // default for minimum sequence identity with query
    qsc = -20.0f;        // default for minimum score per column with query
    cov = 0.0;           // default for minimum coverage threshold
    Ndiff = 100;         // pick Ndiff most different sequences from alignment
    wg = false;
    pca = 1.4;
    pcb = 2.5;
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

    // rebuildfasta
    useHeaderFile = false;

    // translate nucleotide
    translationTable = 1;

    // addSequences
    minSequences = 1;

    // filterDb
    filterColumn = 1;
    filterColumnRegex = "^.*$";

    // evaluationscores
    allVsAll = false;
    randomizedRepresentative = false;

}

std::vector<MMseqsParameter> Parameters::combineList(std::vector<MMseqsParameter> &par1,
                                                     std::vector<MMseqsParameter> &par2) {
    std::vector<MMseqsParameter> retVec;
    std::vector< std::vector<MMseqsParameter>> tmp;
    tmp.push_back(par1);
    tmp.push_back(par2);
    for(size_t i = 0; i < tmp.size(); i++) {
        std::vector<MMseqsParameter> currPar = tmp[i];
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
            ss << par[i].name << " ";
            ss << *((std::string *) par[i].value) << " ";
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
