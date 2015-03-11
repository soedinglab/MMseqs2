#include "Debug.h"
#include "CommandDeclarations.h"
#include "Util.h"


struct Command {
    const char *cmd;

    int (*commandFunction)(int, const char **);
};


static struct Command commands[] = {
        {"alignment", alignment},
        {"cluster", cluster},
        {"search", search},
        {"clusteringworkflow", clusteringworkflow},
        {"clusterupdate", clusterupdate},
        {"prefilter", prefilter},
        {"createdb", createdb},
        {"createfasta", createfasta},
        {"createindex", createindex},
        {"mergeffindex", mergeffindex},
        {"clusteringtofastadb", clusteringtofastadb},
        {"swapresults", swapresults},
        {"extractorf", extractorf},
        {"createprofiledb", createprofiledb},
        {"translatenucleotide", translatenucleotide},
        {"timetest", timetest}
};


void printUsage() {
    std::string usage("\nAll possible mmseqs command\n");
    usage.append("Written by Martin Steinegger (Martin.Steinegger@campus.lmu.de) & Maria Hauser (mhauser@genzentrum.lmu.de)\n\n");
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
            "createprofiledb    \tConvert ffindex profile databse (HMM/PSSM) to MMseqs ffindex profile database.\n"
            "swapresults        \tSwaps results from the mapping A->A,B,C to A -> A, B -> A, C -> A\n"
            "clusteringtofastadb\tCConvert Convert mmseqs clustering to ffindex indexed fasta format\n"
            "mergeffindex       \tMerge multiple ffindex files based on similar id into one file\n"
            "extractorf         \tExtract all open reading frames from a nucleotide fasta file into a ffindex database\n"
            "translatenucleotide\tTranslate nucleotide sequences into aminoacid sequences in a ffindex database\n"
    );
    Debug(Debug::INFO) << usage;
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
            EXIT(runCommand(p, argc - 1, argv + 1));
        }
    } else {
        printUsage();
        Debug(Debug::ERROR) << "Invalid Command: " << argv[1] << "\n";
        EXIT(EXIT_FAILURE);
    }
    return 0;
}
