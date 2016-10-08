#include "Debug.h"
#include "Command.h"
#include "CommandDeclarations.h"
#include "Util.h"
#include "Parameters.h"
#include "DistanceCalculator.h"

#include <iomanip>
#include <CpuInfo.h>

Parameters& par = Parameters::getInstance();

static struct Command commands[] = {
// Main tools  (for non-experts)
        {"createdb",             createdb,             &par.createdb,             COMMAND_MAIN,
            "Convert protein sequence set in a FASTA file to MMseqs’ sequence DB format",
            "converts a protein sequence set in a FASTA formatted file to MMseqs’ sequence DB format. This format is needed as input to mmseqs search and many other tools.",
            "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
            "<i:fastaFile>  <o:sequenceDB> [mappingFasta]",
            CITATION_MMSEQS2},
        {"search",               search,               &par.searchworkflow,       COMMAND_MAIN,
            "Search with query sequence or profile DB (iteratively) through target sequence DB",
            "Searches with the sequences or profiles query DB through the target sequence DB by running the prefilter tool and the align tool for Smith-Waterman alignment. For each query a results file with sequence matches is written as entry into a database of search results (“alignmentDB”).\nIn iterative profile search mode, the detected sequences satisfying user-specified criteria are aligned to the query MSA, and the resulting query profile is used for the next search iteration. Iterative profile searches are usually much more sensitive than (and at least as sensitive as) searches with single query sequences.",
            "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
            "<i:queryDB> <i:targetDB> <o:alignmentDB> <tmpDir>",
            CITATION_MMSEQS2},
        {"cluster",              clusteringworkflow,   &par.clusteringWorkflow,   COMMAND_MAIN,
            "Compute clustering of a sequence DB (quadratic time)",
            "Clusters sequences by similarity. It compares all sequences in the sequence DB with each other using mmseqs search, filters alignments according to user-specified criteria (max. E-value, min. coverage,...),   and runs mmseqs clust to group similar sequences together into clusters.",
            "Martin Steinegger <martin.steinegger@mpibpc.mpg.de> & Lars von den Driesch",
            "<i:sequenceDB> <o:clusterDB> <tmpDir>",
            CITATION_MMSEQS2|CITATION_MMSEQS1},
        {"createindex",          createindex,          &par.createindex,          COMMAND_MAIN,
            "Precompute index table of sequence DB for faster searches",
            "Precomputes an index table for the sequence DB. Handing over the precomputed index table as input to mmseqs search or mmseqs prefilter eliminates the computational overhead of building the index table on the fly.",
            "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
            "<i:sequenceDB> <o:indexDB> <tmpDir>",
            CITATION_MMSEQS2},
// Utility tools for format conversions
        {"createtsv",            createtsv,            &par.onlyverbosity,        COMMAND_FORMAT_CONVERSION,
            "Create tab-separated flat file from prefilter DB, alignment DB, or cluster DB",
            NULL,
            "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>r",
            "<i:queryDB> <i:targetDB> <i:resultDB> <o:tsvFile>",
            CITATION_MMSEQS2},
        {"convertalis",          convertalignments,    &par.convertalignments,    COMMAND_FORMAT_CONVERSION,
            "Convert alignment DB to BLAST-tab format, SAM flat file, or to raw pairwise alignments",
            NULL,
            "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
            "<i:queryDb> <i:targetDb> <i:alignmentDB> <o:alignmentFile>",
            CITATION_MMSEQS2},
        {"convertprofiledb",     convertprofiledb,     &par.convertprofiledb,     COMMAND_FORMAT_CONVERSION,
            "Convert ffindex DB of HMM/HMMER3/PSSM files to MMseqs profile DB",
            NULL,
            "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
            "<i:ffindexProfileDB> <o:profileDB>",
            CITATION_MMSEQS2},
        {"convert2fasta",        convert2fasta,        &par.convert2fasta,        COMMAND_FORMAT_CONVERSION,
            "Convert sequence DB to FASTA format",
            NULL,
            "Milot Mirdita <milot@mirdita.de>",
            "<i:sequenceDB> <o:fastaFile>",
            CITATION_MMSEQS2},
        {"result2flat",          result2flat,          &par.result2flat,          COMMAND_FORMAT_CONVERSION,
            "Create a FASTA-like flat file from prefilter DB, alignment DB, or cluster DB",
            NULL,
            "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
            "<i:queryDB> <i:targetDB> <i:resultDB> <o:fastaDB>",
            CITATION_MMSEQS2},
// Utility tools for clustering
        {"clusterupdate",        clusterupdate,        &par.clusterUpdate,        COMMAND_CLUSTER,
            "Update clustering of old sequence DB to clustering of new sequence DB",
            NULL,
            "Clovis Galiez & Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
            "<i:oldSequenceDB> <i:newSequenceDB> <i:oldClustResultDB> <o:newClustResultDB> <tmpDir>",
            CITATION_MMSEQS2|CITATION_MMSEQS1},
        {"createseqfiledb",      createseqfiledb,      &par.createseqfiledb,      COMMAND_CLUSTER,
            "Create DB of unaligned FASTA files (1 per cluster) from sequence DB and cluster DB",
            NULL,
            "Milot Mirdita <milot@mirdita.de>",
            "<i:sequenceDB> <i:clusterDB> <o:fastaDB>",
            CITATION_MMSEQS2},
        {"mergeclusters",        mergeclusters,        &par.onlyverbosity,        COMMAND_CLUSTER,
            "Merge multiple cluster DBs into single cluster DB",
            NULL,
            "Maria Hauser & Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
            "<i:sequenceDB> <o:clusterDB> <i:clusterDB1> ... <i:clusterDBn>",
            CITATION_MMSEQS2},
// Expert tools (for advanced users)
        {"prefilter",            prefilter,            &par.prefilter,            COMMAND_EXPERT,
            "Search with query sequence / profile DB through target DB (k-mer matching + ungapped alignment)",
            "Searches with the sequences or profiles in query DB through the target sequence DB in two consecutive stages: a very fast k-mer matching stage (double matches on same diagonal) and a subsequent ungapped alignment stage. For each query a results file with sequence matches is written as entry into the prefilter DB.",
            "Martin Steinegger <martin.steinegger@mpibpc.mpg.de> & Maria Hauser",
            "<i:queryDB> <i:targetDB> <o:prefilterDB>",
            CITATION_MMSEQS2},
        {"align",                align,                &par.align,                COMMAND_EXPERT,
            "Compute Smith-Waterman alignments for previous results (e.g. prefilter DB, cluster DB)",
            "Calculates Smith-Waterman alignment scores between all sequences in the query database and the sequences of the target database which passed the prefiltering.",
            "Martin Steinegger <martin.steinegger@mpibpc.mpg.de> & Maria Hauser",
            "<i:queryDB> <i:targetDB> <i:prefilterDB> <o:alignmentDB>",
            CITATION_MMSEQS2},
        {"clust",                clust,                &par.clust,                COMMAND_EXPERT,
            "Cluster sequence DB from alignment DB (e.g. created by searching DB against itself)",
            "Computes a clustering of a sequence DB based on the alignment DB containing for each query sequence or profile the Smith Waterman alignments generated by mmseqs align. (When given a prefilter DB as input the tool will use the ungapped alignment scores scores for the clustering.) The tool reads the search results DB,  constructs a similarity graph based on the matched sequences in alignment DB, and applies one of several clustering algorithms. The first, representative sequence of each cluster is connected by an edge to each cluster member. Its names are used as ID in the resulting cluster DB, the entries contain the names of all member sequences.",
            "Martin Steinegger <martin.steinegger@mpibpc.mpg.de> & Lars von den Driesch & Maria Hauser",
            "<i:sequenceDB> <i:alignmentDB> <o:clusterDB>",
            CITATION_MMSEQS2|CITATION_MMSEQS1},
        {"clustlinear",          clustlinear,          &par.prefilter,            COMMAND_EXPERT,
            "Cluster sequences of >70% sequence identity *in linear time*",
            "Detects redundant sequences based on reduced alphabet and k-mer sorting.",
            "Martin Steinegger <martin.steinegger@mpibpc.mpg.de> ",
            "<i:sequenceDB> <o:alignmentDB>",
            CITATION_MMSEQS2},
        {"clusthash",            clusthash,            &par.clusthash,            COMMAND_EXPERT,
            "Cluster sequences of same length and >90% sequence identity *in linear time*",
            "Detects redundant sequences based on reduced alphabet hashing and hamming distance.",
            "Martin Steinegger <martin.steinegger@mpibpc.mpg.de> ",
            "<i:sequenceDB> <o:alignmentDB>",
            CITATION_MMSEQS2},
// Utility tools to manipulate DBs
        {"extractorfs",          extractorfs,          &par.extractorfs,          COMMAND_DB,
            "Extract open reading frames from all six frames from nucleotide sequence DB",
            NULL,
            "Milot Mirdita <milot@mirdita.de>",
            "<i:sequenceDB> <o:sequenceDB>",
            CITATION_MMSEQS2},
        {"translatenucs",        translatenucs,        &par.translatenucs,        COMMAND_DB,
            "Translate nucleotide sequence DB into protein sequence DB",
            NULL,
            "Milot Mirdita <milot@mirdita.de>",
            "<i:sequenceDB> <o:sequenceDB>",
            CITATION_MMSEQS2},
        {"swapresults",          swapresults,          &par.empty,                COMMAND_DB,
            "Reformat prefilter/alignment/cluster DB as if target DB had been searched through query DB",
            NULL,
            "Martin Steinegger <martin.steinegger@mpibpc.mpg.de> & Clovis Galiez",
            "<i:queryDB> <i:targetDB> <i:resultDB> <o:resultDB>",
            CITATION_MMSEQS2},
        {"mergedbs",             mergedbs,             &par.mergedbs,             COMMAND_DB,
            "Merge multiple DBs into a single DB, based on IDs (names) of entries",
            NULL,
            "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
            "<i:sequenceDB> <o:resultDB> <i:resultDB1> ... <i:resultDBn>",
            CITATION_MMSEQS2},
        {"splitdb",              splitdb,              &par.splitdb,              COMMAND_DB,
            "Split a mmseqs DB into multiple DBs",
            NULL,
            "Milot Mirdita <milot@mirdita.de>",
            "<i:sequenceDB> <o:sequenceDB_1..N>",
            CITATION_MMSEQS2},
        {"subtractdbs",          subtractdbs,          &par.subtractdbs,          COMMAND_DB,
            "Generate a DB with entries of first DB not occurring in second DB",
            NULL,
            "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
            "<i:resultDBLeft> <i:resultDBRight> <o:resultDB>",
            CITATION_MMSEQS2},
        {"filterdb",             filterdb,             &par.filterDb,             COMMAND_DB,
            "Filter a DB by conditioning (regex, numerical, ...) on one of its whitespace-separated columns",
            NULL,
            "Clovis Galiez & Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
            "<i:resultDB> <o:resultDB>",
            CITATION_MMSEQS2},
        {"createsubdb",          createsubdb,          &par.onlyverbosity,        COMMAND_DB,
            "Create a subset of a DB from a file of IDs of entries",
            NULL,
            "Milot Mirdita <milot@mirdita.de>",
            "<i:subsetFile> <i:resultDB> <o:resultDB>",
            CITATION_MMSEQS2},
        {"result2profile",       result2profile,       &par.result2profile,       COMMAND_DB,
            "Compute profile and consensus DB from a prefilter, alignment or cluster DB",
            NULL,
            "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
            "<i:queryDB> <targetDB> <i:resultDB> <o:profileDB>",
            CITATION_MMSEQS2},
        {"result2msa",           result2msa,           &par.result2msa,           COMMAND_DB,
            "Generate MSAs for queries by locally aligning their matched targets in prefilter/alignment/cluster DB",
            NULL,
            "Martin Steinegger (martin.steinegger@mpibpc.mpg.de) & Milot Mirdita <milot@mirdita.de> & Clovis Galiez",
            "<i:queryDB> <i:targetDB> <i:resultDB> <o:msaDB>",
            CITATION_MMSEQS2},
        {"result2stats",         result2stats,         &par.result2stats,         COMMAND_DB,
            "Compute statistics for each entry in a sequence, prefilter, alignment or cluster DB",
            NULL,
            "Clovis Galiez & Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
            "<i:queryDB> <i:targetDB> <i:resultDB> <o:statsDB>",
            CITATION_MMSEQS2},
// Special-purpose utilities
        {"diffseqdbs",           diffseqdbs,           &par.onlyverbosity,        COMMAND_SPECIAL,
            "Find IDs of sequences kept, added and removed between two versions of sequence DB",
            "It creates 3 filtering files, that can be used in cunjunction with \"createsubdb\" tool.\nThe first file contains the keys that has been removed from DBold to DBnew.\nThe second file maps the keys of the kept sequences from DBold to DBnew.\nThe third file contains the keys of the sequences that have been added in DBnew.",
            "Clovis Galiez & Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
            "<i:oldSequenceDB> <i:newSequenceDB> <o:rmSeqKeysFile> <o:keptSeqKeysFile> <o:newSeqKeysFile>",
            CITATION_MMSEQS2},
        {"concatdbs",            concatdbs,            &par.onlyverbosity,        COMMAND_SPECIAL,
            "Concatenate two DBs, giving new IDs to entries from second input DB",
            NULL,
            "Clovis Galiez & Martin Steinegger (martin.steinegger@mpibpc.mpg.de)",
            "<i:resultDB> <i:resultDB> <o:resultDB>",
            CITATION_MMSEQS2},
        {"summarizetabs",        summarizetabs,        &par.summarizetabs,        COMMAND_SPECIAL,
            "Extract annotations from HHblits BAST-tab-formatted results",
            NULL,
            "Milot Mirdita <milot@mirdita.de> & Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
            "<blastTabDB> <lengthFile> <outDB>",
            CITATION_MMSEQS2|CITATION_UNICLUST},
        {"gff2db",               gff2db,               &par.gff2ffindex,          COMMAND_SPECIAL,
            "Turn a gff3 (generic feature format) file into a gff3 DB",
            NULL,
            "Milot Mirdita <milot@mirdita.de>",
            "<i:gff3File> <i:sequenceDB> <o:sequenceDB>",
            CITATION_MMSEQS2},
        {"maskbygff",            maskbygff,            &par.gff2ffindex,          COMMAND_SPECIAL,
            "X out sequence regions in a sequence DB by features in a gff3 file",
            NULL,
            "Milot Mirdita <milot@mirdita.de>",
            "<i:gff3File> <i:sequenceDB> <o:sequenceDB>",
            CITATION_MMSEQS2},
        {"prefixid",             prefixid,             &par.prefixid,             COMMAND_SPECIAL,
            "For each entry in a DB prepend the entry ID to the entry itself",
            NULL,
            "Milot Mirdita <milot@mirdita.de>",
            "<i:resultDB> <o:resultDB>",
            CITATION_MMSEQS2},
        {"convertkb",            convertkb,            &par.convertkb,            COMMAND_SPECIAL,
            "Convert UniProt knowledge flat file into knowledge DB for the selected column types",
            NULL,
            "Milot Mirdita <milot@mirdita.de>",
            "<uniprotkbFile> <uniprotkbDB>",
            CITATION_MMSEQS2},
        {"summarizeheaders",     summarizeheaders,     &par.summarizeheaders,     COMMAND_SPECIAL,
            "Return a new summarized header DB from the UniProt headers of a cluster DB",
            NULL,
            "Milot Mirdita <milot@mirdita.de>",
            "<i:queryHeaderDB> <i:targetHeaderDB> <i:clusterDB> <o:headerDB>",
            CITATION_MMSEQS2|CITATION_UNICLUST},
        {"extractalignedregion", extractalignedregion, &par.extractalignedregion, COMMAND_SPECIAL,
            "Extract aligned sequence region",
            NULL,
            "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
            "<i:queryDB> <i:targetDB> <i:resultDB> <o:domainDB>",
            CITATION_MMSEQS2},
        {"extractdomains",       extractdomains,       &par.extractdomains,       COMMAND_SPECIAL,
            "Extract highest scoring alignment region for each sequence from BLAST-tab file",
            NULL,
            "Milot Mirdita <milot@mirdita.de> & Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
            "<i:domainDB> <i:msaDB> <o:domainDB>",
            CITATION_MMSEQS2|CITATION_UNICLUST},
        {"shellcompletion",      shellcompletion,      &par.empty,                COMMAND_HIDDEN,
            "",
            NULL,
            "",
            "",
            CITATION_MMSEQS2},
        {"computeGOscore",       computeGOscore,       &par.evaluationscores,     COMMAND_HIDDEN,
            "Compute GO scores for a result of clustering",
            NULL,
            "Lars von den Driesch",
            "<gofolder> <clustering_file> <prefix> <outputfolder>",
            CITATION_MMSEQS2},
};

void checkCpu();

void printUsage() {
    std::stringstream usage;
    usage << "MMseqs2 (Many against Many sequence searching) is an open-source software suite for very fast, \n"
            "parallelizable protein sequence searches and clustering of huge protein sequence data sets.\n\n";
    usage << "Please cite: M. Steinegger and J. Soding. Sensitive protein sequence searching for the analysis of massive data sets. bioRxiv 079681 (2016).\n\n";
#ifdef GIT_SHA1
#define str2(s) #s
#define str(s) str2(s)
	std::string gitHash(str(GIT_SHA1));
	usage << "MMseqs Version: " << gitHash << "\n";
#undef str
#undef str2
#endif
    usage << "© Martin Steinegger (martin.steinegger@mpibpc.mpg.de)\n";

    struct {
        const char* title;
        CommandMode mode;
    } categories[] = {
            {"Main tools  (for non-experts)",  COMMAND_MAIN},
            {"Utility tools for format conversions",   COMMAND_FORMAT_CONVERSION},
            {"Utility tools for clustering",     COMMAND_CLUSTER},
            {"Core tools (for advanced users)",     COMMAND_EXPERT},
            {"Utility tools to manipulate DBs",     COMMAND_DB},
            {"Special-purpose utilities",     COMMAND_SPECIAL},
    };

    for(size_t i = 0; i < ARRAY_SIZE(categories); ++i) {
        usage << "\n" << std::setw(20) << categories[i].title << "\n";
        for (size_t j = 0; j < ARRAY_SIZE(commands); j++) {
            struct Command *p = commands + j;
            if (p->mode == categories[i].mode) {
                usage << std::left << std::setw(20) << "  " + std::string(p->cmd) << "\t" << p->shortDescription << "\n";
            }
        }
    }

    usage << "\nBash completion for tools and parameters can be installed by adding \"source path/to/mmseqs/util/bash-completion.sh\" to your \"$HOME/.bash_profile\".\n"
            "Include the location of the MMseqs binaries is in your \"$PATH\" environment variable.";

    Debug(Debug::INFO) << usage.str() << "\n";
}


int isCommand(const char *s) {
    for (size_t i = 0; i < ARRAY_SIZE(commands); i++) {
        struct Command *p = commands + i;
        if (!strcmp(s, p->cmd))
            return 1;
    }
    return 0;
}

int runCommand(const Command &p, int argc, const char **argv) {
    int status = p.commandFunction(argc, argv, p);
    if (status)
        return status;
    return 0;
}

int shellcompletion(int argc, const char** argv, const Command& command) {
    // mmseqs programs
    if(argc == 0) {
        for (size_t i = 0; i < ARRAY_SIZE(commands); i++) {
            struct Command *p = commands + i;
            if(p->mode == COMMAND_HIDDEN)
                continue;
            Debug(Debug::INFO) << p->cmd << " ";
        }
        Debug(Debug::INFO) << "\n";
    }

    // mmseqs parameters for given program
    if(argc == 1) {
        for (size_t i = 0; i < ARRAY_SIZE(commands); i++) {
            struct Command *p = commands + i;
            if(strcmp(p->cmd, argv[0]) != 0)
                continue;
            if(!p->params)
                continue;
            for(std::vector<MMseqsParameter>::const_iterator it = p->params->begin(); it != p->params->end(); ++it) {
                Debug(Debug::INFO) << it->name << " ";
            }
            Debug(Debug::INFO) << "\n";
            break;
        }
        Debug(Debug::INFO) << "\n";
    }

    return EXIT_SUCCESS;
}

int main(int argc, const char **argv) {
    checkCpu();
    if (argc < 2) {
        printUsage();
        EXIT(EXIT_FAILURE);
    }
    setenv("MMSEQS", argv[0], true);
    if (isCommand(argv[1])) {
        for (size_t i = 0; i < ARRAY_SIZE(commands); i++) {
            const struct Command *p = commands + i;
            if (strcmp(p->cmd, argv[1]))
                continue;

            EXIT(runCommand(*p, argc - 2, argv + 2));
        }
    } else {
        printUsage();
        Debug(Debug::ERROR) << "\nInvalid Command: " << argv[1] << "\n";

        // Suggest some command that the user might have meant
        size_t index = SIZE_MAX;
        size_t minDistance = SIZE_MAX;
        for (size_t i = 0; i < ARRAY_SIZE(commands); ++i) {
            struct Command *p = commands + i;
            if(p->mode == COMMAND_HIDDEN)
                continue;

            size_t distance = DistanceCalculator::levenshteinDistance(argv[1], p->cmd);
            if(distance < minDistance) {
                minDistance = distance;
                index = i;
            }
        }

        if(index != SIZE_MAX) {
            Debug(Debug::ERROR) << "Did you mean \"mmseqs " << (commands + index)->cmd << "\"?\n";
        }

        EXIT(EXIT_FAILURE);
    }

    return 0;
}

void checkCpu() {
    CpuInfo info;
    if(info.HW_x64 == false) {
        Debug(Debug::ERROR) << "64 bit system is required to run MMseqs.\n";
        EXIT(EXIT_FAILURE);
    }
#ifdef SEE
    if(info.HW_SSE41 == false) {
        Debug(Debug::ERROR) << "SSE4.1 is required to run MMseqs.\n";
        EXIT(EXIT_FAILURE);
    }
#endif
#ifdef AVX2
    if(info.HW_AVX2 == false){
        Debug(Debug::ERROR) << "Your machine does not support AVX2.\n";
        if(info.HW_SSE41 == true) {
            Debug(Debug::ERROR) << "Please compile with SSE4.1 cmake -DHAVE_SSE4_1=1 \n";
        }else{
            Debug(Debug::ERROR) << "SSE 4.1 is the minimum requirement to run MMseqs.\n";
        }
        EXIT(EXIT_FAILURE);
    }
#endif
}
