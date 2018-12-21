#ifndef MMSEQS_COMMAND_H
#define MMSEQS_COMMAND_H

#include <vector>

const int CITATION_MMSEQS2  = 1 << 0;
const int CITATION_MMSEQS1  = 1 << 1;
const int CITATION_UNICLUST = 1 << 2;
const int CITATION_LINCLUST = 1 << 3;
const int CITATION_PLASS    = 1 << 4;
const int CITATION_SERVER   = 1 << 5;

struct MMseqsParameter;

enum CommandMode {
    COMMAND_MAIN = 0,
    COMMAND_FORMAT_CONVERSION,
    COMMAND_CLUSTER,
    COMMAND_TAXONOMY,
    COMMAND_MULTIHIT,
    COMMAND_DB,
    COMMAND_EXPERT,
    COMMAND_SPECIAL,
    COMMAND_HIDDEN,
    COMMAND_EASY
};

struct Command {
    const char *cmd;
    int (*commandFunction)(int, const char **, const Command&);
    std::vector<MMseqsParameter*>* params;
    CommandMode mode;
    const char *shortDescription;
    const char *longDescription;
    const char *author;
    const char *usage;
    int citations;
};

struct Categories {
    const char* title;
    CommandMode mode;
};

#endif
