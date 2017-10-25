#include "Command.h"

std::vector<Categories> categories = {
        {"Main tools  (for non-experts)",        COMMAND_MAIN},
        {"Utility tools for format conversions", COMMAND_FORMAT_CONVERSION},
        {"Taxonomy tools",                       COMMAND_TAXONOMY},
        {"Utility tools for clustering",         COMMAND_CLUSTER},
        {"Core tools (for advanced users)",      COMMAND_EXPERT},
        {"Utility tools to manipulate DBs",      COMMAND_DB},
        {"Special-purpose utilities",            COMMAND_SPECIAL},
};
