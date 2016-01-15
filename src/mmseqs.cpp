#include "Debug.h"
#include "CommandDeclarations.h"
#include "Util.h"
#include "Parameters.h"

Parameters par;

struct Command {
    const char *cmd;

    int (*commandFunction)(int, const char **);

    std::vector<MMseqsParameter>* params;
};

static struct Command commands[] = {
    {"alignment",           alignment,              &par.alignment},
    {"cluster",             cluster,                &par.clustering},
    {"search",              search,                 &par.searchworkflow},
    {"clusteringworkflow",  clusteringworkflow,     &par.clusteringWorkflow},
    {"clusterupdate",       clusterupdate,          &par.clusterUpdate},
    {"prefilter",           prefilter,              &par.prefilter},
    {"createdb",            createdb,               &par.createdb},
    {"createfasta",         createfasta,            &par.onlyverbosity},
    {"createtsv",           createtsv,              &par.onlyverbosity},
    {"formatalignment",     formatalignment,        &par.formatalignment},
    {"createindex",         createindex,            &par.formatalignment},
    {"mergeffindex",        mergeffindex,           &par.empty},
    {"mergecluster",        mergecluster,           &par.onlyverbosity},
    {"clustertofastadb",    clusteringtofastadb,    &par.empty},
    {"swapresults",         swapresults,            &par.empty},
    {"extractorf",          extractorf,             &par.extractorf},
    {"createprofiledb",     createprofiledb,        &par.createprofiledb},
    {"translatenucleotide", translatenucleotide,    &par.onlyverbosity},
    {"timetest",            timetest,               &par.empty},
    {"legacycs219",         legacycs219,            &par.onlyverbosity},
    {"findsorf",            findsorf,               &par.onlyverbosity},
    {"resulttoprofiledb",   result2profile,         &par.createprofiledb},
    {"rebuildfasta",        rebuildfasta,           &par.rebuildfasta},
    {"splitffindex",        splitffindex,           &par.splitffindex},
    {"gff2ffindex",         gff2ffindex ,           &par.gff2ffindex},
    {"shellcompletion",     shellcompletion,        &par.empty}
};


void printUsage() {
    std::string usage("\nAll possible mmseqs commands\n");
    usage.append("Written by Martin Steinegger (martin.steinegger@mpibpc.mpg.de) & Maria Hauser (mhauser@genzentrum.lmu.de)\n\n");
    usage.append("Main tools: \n"
            "prefilter          \tCalculates similarity scores between all sequences in the query db and all sequences in the target db\n"
            "alignment          \tCalculates Smith-Waterman alignment scores from prefilter output\n"
            "cluster            \tCalculates clustering of a sequence database based on alignment output with set cover algorithm\n"
            "clusteringworkflow \tCalculates cascaded clustering of a ffindex sequence database. (Prefiltering -> Alignment -> cluster)*n \n"
            "clusterupdate      \tUpdates the existing clustering of the previous database version with new sequences from the current version\n"
            "\nHelper: \n"
            "createdb           \tConvert fasta to ffindex (all programs need ffindex as input)\n"
            "createindex        \tConvert ffindex to fast index for prefiltering\n"
            "createfasta        \tConvert ffindex to fasta\n"
            "createtsv          \tConvert ffindex to tsv\n"
            "formatalignment    \tConvert a ffindex alignment database to BLAST tab or SAM flat file.\n"
            "createprofiledb    \tConvert ffindex profile databse (HMM/PSSM) to MMseqs ffindex profile database.\n"
            "swapresults        \tSwaps results from the mapping A->A,B,C to A -> A, B -> A, C -> A\n"
            "clustertofastadb   \tConvert Convert mmseqs clustering to ffindex indexed fasta format\n"
            "clustertoprofiledb \tCalculates profile from clustering\n"
            "mergeffindex       \tMerge multiple ffindex files based on similar id into one file\n"
            "splitffindex       \tSplits a ffindex database into multiple ffindex databases.\n"
            "extractorf         \tExtract all open reading frames from a nucleotide fasta file into a ffindex database\n"
            "translatenucleotide\tTranslate nucleotide sequences into aminoacid sequences in a ffindex database\n"
            "legacycs219        \tTranslates a cs219 ffindex database into its legacy format. This tool is part of the mmseqs-based HH-suite database pipeline\n"
            "rebuildfasta       \tRebuild a fasta file from a ffindex database\n"
            "gff2ffindex        \tTurn a GFF3 file into a ffindex database\n"
    );
    Debug(Debug::INFO) << usage << "\n";
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
            if(!strcmp(p->cmd, "shellcompletion"))
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
