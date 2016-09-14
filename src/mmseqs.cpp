#include "Debug.h"
#include "CommandDeclarations.h"
#include "Util.h"
#include "Parameters.h"

#include <iomanip>
#include <CpuInfo.h>

Parameters par;

enum CommandMode {
    COMMAND_MAIN = 0,
    COMMAND_FORMAT_CONVERSION,
    COMMAND_CLUSTER,
    COMMAND_DB,
    COMMAND_EXPERT,
    COMMAND_SPECIAL,
    COMMAND_HIDDEN
};

struct Command {
    const char *cmd;
    int (*commandFunction)(int, const char **);
    std::vector<MMseqsParameter>* params;
    CommandMode mode;
    const char *description;
};

static struct Command commands[] = {
// Main tools  (for non-experts)
        {"createdb",            createdb,               &par.createdb,              COMMAND_MAIN,
                "Convert protein sequence set in a FASTA file to MMseqs’ sequence DB format"},
        {"search",              search,                 &par.searchworkflow,        COMMAND_MAIN,
                "Search with query sequence or profile DB (iteratively) through target sequence DB"},
        {"cluster",  clusteringworkflow,     &par.clusteringWorkflow,    COMMAND_MAIN,
                "Compute clustering of a sequence DB (quadratic time)"},
        {"createindex",         createindex,            &par.createindex,           COMMAND_MAIN,
                "Precompute index table of sequence DB for faster searches"},
// Utility tools for format conversions
        {"createtsv",           createtsv,              &par.onlyverbosity,         COMMAND_FORMAT_CONVERSION,
                "Create tab-separated flat file from prefilter DB, alignment DB, or cluster DB"},
        {"convertalignments",     convertalignments,        &par.convertalignments,       COMMAND_FORMAT_CONVERSION,
                "Convert alignment DB to BLAST-tab format, SAM flat file, or to raw pairwise alignments"},
        {"convertprofiledb",     convertprofiledb,        &par.convertprofiledb,       COMMAND_FORMAT_CONVERSION,
                "Convert ffindex DB of HMM/HMMER3/PSSM files to MMseqs profile DB"},
        {"convert2fasta",        convert2fasta,           &par.convert2fasta,          COMMAND_FORMAT_CONVERSION,
                "Convert sequence DB to FASTA format"},
        {"result2flat",         result2flat,            &par.onlyverbosity,         COMMAND_FORMAT_CONVERSION,
                "Create a FASTA-like flat file from prefilter DB, alignment DB, or cluster DB"},
// Utility tools for clustering
        {"clusterupdate",       clusterupdate,          &par.clusterUpdate,         COMMAND_CLUSTER,
                "Update clustering of old sequence DB to clustering of new sequence DB"},
        {"createseqfiledb",        createseqfiledb,           &par.createseqfiledb,          COMMAND_CLUSTER,
                "Create DB of unaligned FASTA files, one per cluster in cluster DB"},
        {"mergeclusters",        mergeclusters,           &par.onlyverbosity,         COMMAND_CLUSTER,
                "Merge multiple cluster DBs into single cluster DB"},
// Expert tools (for advanced users)
        {"prefilter",           prefilter,              &par.clustlinear,             COMMAND_EXPERT,
                "Search with query sequence / profile DB through target DB (k-mer matching + ungapped alignment)"},
        {"align",           align,              &par.align,             COMMAND_EXPERT,
                "Compute Smith-Waterman alignments for previous results (e.g. prefilter DB, cluster DB)"},
        {"clust",             clust,                &par.clust,            COMMAND_EXPERT,
                "Cluster sequence DB from alignment DB (e.g. created by searching DB against itself)"},
        {"clustlinear",           clustlinear,       &par.clustlinear,             COMMAND_EXPERT,
                "Cluster sequences of >70% sequence identity *in linear time*"},
        {"clusthash",   clusthash,       &par.clusthash,             COMMAND_EXPERT,
                "Cluster sequences of same length and >90% sequence identity *in linear time*"},
// Utility tools to manipulate DBs
        {"extractorfs",          extractorfs,             &par.extractorfs,            COMMAND_DB,
                "Extract open reading frames from all six frames from nucleotide sequence DB"},
        {"translatenucs", translatenucs,    &par.translatenucs,   COMMAND_DB,
                "Translate nucleotide sequence DB into protein sequence DB"},
        {"swapresults",         swapresults,            &par.empty,                 COMMAND_DB,
                "Reformat prefilter/alignment/cluster DB as if target DB had been searched through query DB"},
        {"mergedbs",        mergedbs,           &par.empty,                 COMMAND_DB,
                "Merge multiple DBs into a single DB, based on IDs (names) of entries"},
        {"splitdb",        splitdb,           &par.splitdb,          COMMAND_DB,
                "Split a mmseqs DB into multiple DBs"},
        {"subtractdbs",     subtractdbs,        &par.subtractdbs,       COMMAND_DB,
                "Generate a DB with entries of first DB not occurring in second DB"},
        {"filterdb",            filterdb,               &par.filterDb,              COMMAND_DB,
                "Filter a DB by conditioning (regex, numerical, …) on one of its whitespace-separated columns"},
        {"createsubdb",         createsubdb,                  &par.onlyverbosity,         COMMAND_DB,
                "Create a subset of a DB from a file of IDs of entries"},
        {"result2profile",      result2profile,         &par.result2profile,        COMMAND_DB,
                "Compute profile and consensus DB from a prefilter, alignment or cluster DB"},
        {"result2msa",          result2msa,             &par.result2msa,            COMMAND_DB,
                "Generate MSAs for queries by locally aligning their matched targets in prefilter/alignment/cluster DB"},
        {"result2stats",          result2stats,             &par.result2stats,            COMMAND_DB,
                "Compute statistics for each entry in a sequence, prefilter, alignment or cluster DB"},
// Special-purpose utilities
        {"diffseqdbs",           diffseqdbs,                   &par.onlyverbosity,         COMMAND_SPECIAL,
                "Find IDs of sequences kept, added and removed between two versions of sequence DB"},
        {"concatdbs",            concatdbs,               &par.onlyverbosity,         COMMAND_SPECIAL,
                "Concatenate two DBs, giving new IDs to entries from second input DB"},
        {"summarizetabs",       summarizetabs,          &par.summarizetabs,         COMMAND_SPECIAL,
                "Extract annotations from HHblits BAST-tab-formatted results"},
        {"gff2db",         gff2db ,           &par.gff2ffindex,           COMMAND_SPECIAL,
                "Turn a gff3 (generic feature format) file into a gff3 DB"},
        {"maskbygff",           maskbygff,              &par.gff2ffindex,           COMMAND_SPECIAL,
                "X out sequence regions in a sequence DB by features in a gff3 file"},
        {"prefixid",            prefixid,               &par.prefixid,              COMMAND_SPECIAL,
                "For each entry in a DB prepend the entry ID to the entry itself"},
        {"convertkb",           convertkb,              &par.convertkb,             COMMAND_SPECIAL,
            "Convert UniProt knowledge flat file into knowledge DB for the selected column types"},
        {"summarizeheaders",    summarizeheaders,       &par.summarizeheaders,      COMMAND_SPECIAL,
                "Return a new summarized header DB from the UniProt headers of a cluster DB"},
        {"extractalignedregion",extractalignedregion,   &par.extractalignedregion,  COMMAND_SPECIAL,
                "Extract aligned sequence region"},
        {"extractdomains",      extractdomains,         &par.extractdomains,        COMMAND_SPECIAL,
                "Extract highest scoring alignment region for each sequence from BLAST-tab file"},
        {"shellcompletion",     shellcompletion,        &par.empty,                 COMMAND_HIDDEN, ""},
        {"computeGOscore",      computeGOscore,         &par.evaluationscores,      COMMAND_HIDDEN,
                "Compute GO scores for a result of clustering"},
};


void checkCpu();

void printUsage() {
    std::stringstream usage;
    usage << "MMseqs2 (Many against Many sequence searching) is an open-source software suite for very fast, \n"
            "parallelizable protein sequence searches and clustering of huge protein sequence data sets.\n\n";
    usage << "Please cite: M. Steinegger and J. Soding. Sensitive protein sequence searching for the analysis of massive data sets. bioRxiv XXXX (2016).\n\n";
#ifdef GIT_SHA1
#define str2(s) #s
#define str(s) str2(s)
	std::string gitHash(str(GIT_SHA1));
	usage << "MMseqs Version: " << gitHash << "\n";
#undef str
#undef str2
#endif
    usage << "(C) Martin Steinegger (martin.steinegger@mpibpc.mpg.de) & Maria Hauser)\n";

    struct {
        const char* title;
        CommandMode mode;
    } categories[] = {
            {"Main tools  (for non-experts)",  COMMAND_MAIN},
            {"Utility tools for format conversions",   COMMAND_FORMAT_CONVERSION},
            {"Utility tools for clustering",     COMMAND_CLUSTER},
            {"Expert tools (for advanced users)",     COMMAND_EXPERT},
            {"Utility tools to manipulate DBs",     COMMAND_DB},
            {"Special-purpose utilities",     COMMAND_SPECIAL},
    };

    for(size_t i = 0; i < ARRAY_SIZE(categories); ++i) {
        usage << "\n" << std::setw(20) << categories[i].title << "\n";
        for (size_t j = 0; j < ARRAY_SIZE(commands); j++) {
            struct Command *p = commands + j;
            if (p->mode == categories[i].mode)
                usage << std::setw(20) << p->cmd << "\t" << p->description << "\n";
        }
    }

    usage << "\nBash completion for tools and parameters can be installed by adding a line \"util/bash-completion.sh\" in your \"$HOME/.bash_profile\".\n"
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

int runCommand(Command *p, int argc, const char **argv) {
    int status = p->commandFunction(argc, argv);
    if (status)
        return status;
    return 0;
}

int shellcompletion(int argc, const char** argv) {
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
            struct Command *p = commands + i;
            if (strcmp(p->cmd, argv[1]))
                continue;
            EXIT(runCommand(p, argc - 2, argv + 2));
        }
    } else {
        printUsage();
        Debug(Debug::ERROR) << "Invalid Command: " << argv[1] << "\n";
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
