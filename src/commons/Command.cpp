#include "Command.h"

std::vector<Categories> categories = {
        {"Easy workflows (for non-experts)",     COMMAND_EASY},
        {"Main tools  (for non-experts)",        COMMAND_MAIN},
        {"Utility tools for format conversions", COMMAND_FORMAT_CONVERSION},
        {"Taxonomy tools",                       COMMAND_TAXONOMY},
        {"Multi-hit search tools",               COMMAND_MULTIHIT},
        {"Utility tools for clustering",         COMMAND_CLUSTER},
        {"Core tools (for advanced users)",      COMMAND_EXPERT},
        {"Utility tools to manipulate DBs",      COMMAND_DB},
        {"Special-purpose utilities",            COMMAND_SPECIAL},
};
