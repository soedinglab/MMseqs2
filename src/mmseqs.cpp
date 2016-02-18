#include "Debug.h"
#include "CommandDeclarations.h"
#include "Util.h"
#include "Parameters.h"

#include <iomanip>

Parameters par;

enum CommandMode {
    COMMAND_MAIN = 0,
    COMMAND_WORKFLOW,
    COMMAND_HELPER,
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
        {"prefilter",           prefilter,              &par.prefilter,             COMMAND_MAIN,
                "Calculates similarity scores between all sequences in the query db and all sequences in the target db"},
        {"alignment",           alignment,              &par.alignment,             COMMAND_MAIN,
                "Calculates Smith-Waterman alignment scores from prefilter output"},
        {"cluster",             cluster,                &par.clustering,            COMMAND_MAIN,
                "Calculates clustering of a sequence database based on alignment output with set cover algorithm"},
        {"search",              search,                 &par.searchworkflow,        COMMAND_WORKFLOW,
                "Searches protein sequences in a protein database"},
        {"clusteringworkflow",  clusteringworkflow,     &par.clusteringWorkflow,    COMMAND_WORKFLOW,
                "Calculates cascaded clustering of a ffindex sequence database. (Prefiltering -> Alignment -> Cluster)*n"},
        {"clusterupdate",       clusterupdate,          &par.clusterUpdate,         COMMAND_WORKFLOW,
                "Updates the existing clustering of the previous database version with new sequences from the current version"},
        {"createdb",            createdb,               &par.createdb,              COMMAND_HELPER,
                "Convert fasta to ffindex (all programs need ffindex as input)"},
        {"createindex",         createindex,            &par.createindex,           COMMAND_HELPER,
                "Convert ffindex to fast index for prefiltering"},
        {"createfasta",         createfasta,            &par.onlyverbosity,         COMMAND_HELPER,
                "Convert ffindex to fasta"},
        {"createtsv",           createtsv,              &par.onlyverbosity,         COMMAND_HELPER,
                "Convert ffindex to tsv"},
        {"createprofiledb",     createprofiledb,        &par.createprofiledb,       COMMAND_HELPER,
                "Convert ffindex profile databse (HMM/PSSM) to MMseqs ffindex profile database"},
        {"filterdb",            filterdb,               &par.filterDb,              COMMAND_HELPER,
                "Filter a database by column regex"},
        {"formatalignment",     formatalignment,        &par.formatalignment,       COMMAND_HELPER,
                "Convert a ffindex alignment database to BLAST tab or SAM flat file"},
        {"swapresults",         swapresults,            &par.empty,                 COMMAND_HELPER,
                "Swaps results from the mapping A -> (A,B,C) to A -> A, B -> A, C -> A"},
        {"addsequences",        addsequences,           &par.addSequences,          COMMAND_HELPER,
                "Adds sequences in fasta format to a mmseqs clustering"},
        {"mergeffindex",        mergeffindex,           &par.empty,                 COMMAND_HELPER,
                "Merge multiple ffindex files based on ids into one file"},
        {"splitffindex",        splitffindex,           &par.splitffindex,          COMMAND_HELPER,
                "Splits a ffindex database into multiple ffindex databases"},
        {"mergecluster",        mergecluster,           &par.onlyverbosity,         COMMAND_HELPER,
                "Merge multiple cluster result files into one"},
        {"substractresult",     substractresult,        &par.substractresult,       COMMAND_HELPER,
                "Removes all entries with same ID (out = left - right)"},
        {"result2profile",      result2profile,         &par.result2profile,        COMMAND_HELPER,
                "Calculates profiles from a clustering"},
        {"result2msa",          result2msa,             &par.result2msa,            COMMAND_HELPER,
                "Calculates MSAs from a clustering"},
        {"rebuildfasta",        rebuildfasta,           &par.rebuildfasta,          COMMAND_HELPER,
                "Rebuild a fasta file from a ffindex database"},
        {"extractorf",          extractorf,             &par.extractorf,            COMMAND_HELPER,
                "Extract all open reading frames from a nucleotide ffindex into a second ffindex database"},
        {"translatenucleotide", translatenucleotide,    &par.translateNucleotide,   COMMAND_HELPER,
                "Translate nucleotide sequences into aminoacid sequences in a ffindex database"},
        {"maskbygff",           maskbygff,              &par.gff2ffindex,           COMMAND_HELPER,
                "Masks the sequences in an ffindex database by the selected rows in a gff file"},
        {"gff2ffindex",         gff2ffindex ,           &par.gff2ffindex,           COMMAND_HELPER,
                "Turn a GFF3 file into a ffindex database"},
        {"convertkb",           convertkb,              &par.onlyverbosity,         COMMAND_HELPER,
            "Converts Uniprot flat file into ffindex database"},
        {"kbtotsv",             kbtotsv,                &par.onlyverbosity,         COMMAND_HELPER,
                "Turns an UniprotKB file into separate TSV tables"},
        {"result2newick",       result2newick,          &par.onlyverbosity,         COMMAND_HELPER,
            "Extracts clustering relationship from clustering steps into Newick trees"},
        {"timetest",            timetest,               &par.empty,                 COMMAND_HIDDEN, ""},
        {"shellcompletion",     shellcompletion,        &par.empty,                 COMMAND_HIDDEN, ""},
};


void printUsage() {
    std::stringstream usage;
    usage << "All available MMseqs commands\n";
    usage << "Written by Martin Steinegger (martin.steinegger@mpibpc.mpg.de) & Maria Hauser (mhauser@genzentrum.lmu.de)\n";

    struct {
        const char* title;
        CommandMode mode;
    } categories[] = {
            {"Main Tools",  COMMAND_MAIN},
            {"Workflows",   COMMAND_WORKFLOW},
            {"Helpers",     COMMAND_HELPER},
    };

    for(size_t i = 0; i < ARRAY_SIZE(categories); ++i) {
        usage << "\n" << std::setw(20) << categories[i].title << "\n";
        for (size_t j = 0; j < ARRAY_SIZE(commands); j++) {
            struct Command *p = commands + j;
            if (p->mode == categories[i].mode)
                usage << std::setw(20) << p->cmd << "\t" << p->description << "\n";
        }
    }

    usage << "\nBash completion for subcommands and parameters can be installed by adding the following lines to your ~/.bash_profile:\n";
    usage << "if [ -f $MMDIR/util/bash-completion.sh ]; then\n";
    usage << "\t. $MMDIR/util/bash-completion.sh\n";
    usage << "fi\n";

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
    if (argc < 2) {
        printUsage();
        EXIT(EXIT_FAILURE);
    }
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
