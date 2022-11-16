#ifndef MMSEQS_COMMAND_H
#define MMSEQS_COMMAND_H

#include <vector>

const unsigned int CITATION_MMSEQS2  = 1U << 0;
const unsigned int CITATION_MMSEQS1  = 1U << 1;
const unsigned int CITATION_UNICLUST = 1U << 2;
const unsigned int CITATION_LINCLUST = 1U << 3;
const unsigned int CITATION_PLASS    = 1U << 4;
const unsigned int CITATION_SERVER   = 1U << 5;
const unsigned int CITATION_TAXONOMY = 1U << 6;

// Make sure this is always the last bit
// citations from inheriting modules will start from here
const unsigned int CITATION_END      = CITATION_TAXONOMY << 1;

struct MMseqsParameter;

typedef const unsigned int CommandMode;

CommandMode COMMAND_MAIN              = 1U << 1;
CommandMode COMMAND_FORMAT_CONVERSION = 1U << 2;
CommandMode COMMAND_TAXONOMY          = 1U << 3;
CommandMode COMMAND_MULTIHIT          = 1U << 4;
CommandMode COMMAND_DB                = 1U << 5;
CommandMode COMMAND_SPECIAL           = 1U << 6;
CommandMode COMMAND_HIDDEN            = 1U << 7;
CommandMode COMMAND_EASY              = 1U << 8;
CommandMode COMMAND_DATABASE_CREATION = 1U << 9;
CommandMode COMMAND_STORAGE           = 1U << 10;
CommandMode COMMAND_SET               = 1U << 11;
CommandMode COMMAND_SEQUENCE          = 1U << 12;
CommandMode COMMAND_RESULT            = 1U << 13;
CommandMode COMMAND_PREFILTER         = 1U << 14;
CommandMode COMMAND_ALIGNMENT         = 1U << 15;
CommandMode COMMAND_CLUSTER           = 1U << 16;
CommandMode COMMAND_PROFILE           = 1U << 17;
CommandMode COMMAND_PROFILE_PROFILE   = 1U << 18;

CommandMode COMMAND_EXPERT            = 1U << 31;



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
    static std::vector<int> ppResultDb;
    static std::vector<int> ca3mDb;
    static std::vector<int> msaDb;
    static std::vector<int> genericDb;
    static std::vector<int> profileDb;
    static std::vector<int> csDb;
    static std::vector<int> indexDb;
    static std::vector<int> allDb;
    static std::vector<int> allDbAndFlat;
    static std::vector<int> taxResult;
    static std::vector<int> taxonomyReportInput;
    static std::vector<int> directory;
    static std::vector<int> flatfile;
    static std::vector<int> flatfileAndStdin;
    static std::vector<int> flatfileStdinAndGeneric;
    static std::vector<int> flatfileStdinGenericUri;
    static std::vector<int> empty;
};


struct DbType{
    static const int ACCESS_MODE_INPUT = 1;
    static const int ACCESS_MODE_OUTPUT = 2;
    static const int NEED_DATA = 0;
    static const int NEED_HEADER = 1;
    static const int NEED_LOOKUP = 2;
    static const int NEED_TAXONOMY = 4;
    static const int VARIADIC = 8;
    static const int ZERO_OR_ALL = 16;

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
    const char *description;
    const char *examples;
    const char *author;
    const char *usage;
    unsigned int citations;
    std::vector<DbType> databases;
};

struct Categories {
    const char* title;
    CommandMode mode;
};

#endif
