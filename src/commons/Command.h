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



struct DbValidator {
    static std::vector<int> sequenceDb;
    static std::vector<int> nuclDb;
    static std::vector<int> aaDb;
    static std::vector<int> prefAlnResDb;
    static std::vector<int> taxSequenceDb;
    static std::vector<int> nuclAaDb;
    static std::vector<int> alignmentDb;
    static std::vector<int> prefilterDb;
    static std::vector<int> clusterDb;
    static std::vector<int> resultDb;
    static std::vector<int> ca3mDb;
    static std::vector<int> msaDb;
    static std::vector<int> genericDb;
    static std::vector<int> profileDb;
    static std::vector<int> csDb;
    static std::vector<int> indexDb;
    static std::vector<int> allDb;
    static std::vector<int> allDbAndFlat;
    static std::vector<int> taxResult;
    static std::vector<int> directory;
    static std::vector<int> flatfile;
};


struct DbType{
    static const int ACCESS_MODE_INPUT = 1;
    static const int ACCESS_MODE_OUTPUT = 2;
    static const int NEED_DATA = 0;
    static const int NEED_HEADER = 1;
    static const int NEED_LOOKUP = 2;
    static const int NEED_TAXONOMY = 4;
    static const int VARIADIC = 8;

    const char *usageText;
    int accessMode;
    int specialType;
    std::vector<int> * validator;
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
    std::vector<DbType> databases;
};

struct Categories {
    const char* title;
    CommandMode mode;
};

#endif
