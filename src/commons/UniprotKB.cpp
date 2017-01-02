#include <string>
#include <cstring>
#include <vector>
#include <sstream>
#include <algorithm>

#include "UniprotKB.h"
#include "Debug.h"
#include "Util.h"

const std::string UniprotKB::columnNames[] = {
        "ID", "AC",
        "DT", "DE", "GN",
        "OS", "OG", "OC",
        "OX", "OH", "REF",
        "CC", "DR", "PE",
        "KW", "FT", "SEQ"
};

std::string removeAfterFirstColon(std::string in) {
    in.erase(in.find_first_of(":"));
    return in;
}

enum uniprot_prefix_type_t {
    UNIPROT_ID = 0,
    UNIPROT_AC,
    UNIPROT_DT,
    UNIPROT_DE,
    UNIPROT_GN,
    UNIPROT_OS,
    UNIPROT_OG,
    UNIPROT_OC,
    UNIPROT_OX,
    UNIPROT_OH,
    UNIPROT_RN,
    UNIPROT_RP,
    UNIPROT_RC,
    UNIPROT_RX,
    UNIPROT_RG,
    UNIPROT_RA,
    UNIPROT_RT,
    UNIPROT_RL,
    UNIPROT_CC,
    UNIPROT_DR,
    UNIPROT_PE,
    UNIPROT_KW,
    UNIPROT_FT,
    UNIPROT_SQ,
    UNIPROT_SEQ,
    UNIPROT_END
};

enum occurence_t {
    OCCURENCE_REQUIRED = 0,
    OCCURENCE_OPTIONAL
};

enum lines_t {
    LINES_SINGLE = 0,
    LINES_MULTIPLE,
    LINES_MULTIPLE_FOLD,
    LINES_MULTIPLE_CONCAT
};

struct uniprotkb {
    uniprot_prefix_type_t type;
    char prefix[3];
    std::string description;
    occurence_t occurence;
    lines_t lines;
    bool includeInDB;
    int dbColumn;
    std::string (*transform)(std::string);
};

const uniprotkb uniprotkb_prefix[] = {
        {UNIPROT_ID,  "ID", "Identification",               OCCURENCE_REQUIRED, LINES_SINGLE,          true,  0,  Util::removeAfterFirstSpace},
        {UNIPROT_AC,  "AC", "Accession number(s)",          OCCURENCE_REQUIRED, LINES_MULTIPLE_CONCAT, false, 1,  Util::removeWhiteSpace},
        {UNIPROT_DT,  "DT", "Date",                         OCCURENCE_REQUIRED, LINES_MULTIPLE,        true,  2,  NULL},
        {UNIPROT_DE,  "DE", "Description",                  OCCURENCE_REQUIRED, LINES_MULTIPLE,        true,  3,  NULL},
        {UNIPROT_GN,  "GN", "Gene name(s)",                 OCCURENCE_OPTIONAL, LINES_MULTIPLE_FOLD,   true,  4,  NULL},
        {UNIPROT_OS,  "OS", "Organism species",             OCCURENCE_REQUIRED, LINES_MULTIPLE,        true,  5,  NULL},
        {UNIPROT_OG,  "OG", "Organelle",                    OCCURENCE_OPTIONAL, LINES_MULTIPLE,        true,  6,  NULL},
        {UNIPROT_OC,  "OC", "Organism classification",      OCCURENCE_REQUIRED, LINES_MULTIPLE_FOLD,   true,  7,  NULL},
        {UNIPROT_OX,  "OX", "Taxonomy cross-reference",     OCCURENCE_REQUIRED, LINES_SINGLE,          true,  8,  NULL},
        {UNIPROT_OH,  "OH", "Organism host",                OCCURENCE_OPTIONAL, LINES_MULTIPLE,        true,  9,  NULL},
        {UNIPROT_RN,  "RN", "Reference number",             OCCURENCE_REQUIRED, LINES_MULTIPLE,        true,  10, NULL},
        {UNIPROT_RP,  "RP", "Reference position",           OCCURENCE_REQUIRED, LINES_MULTIPLE,        false, 10, NULL},
        {UNIPROT_RC,  "RC", "Reference comment(s)",         OCCURENCE_OPTIONAL, LINES_MULTIPLE,        false, 10, NULL},
        {UNIPROT_RX,  "RX", "Reference cross-reference(s)", OCCURENCE_OPTIONAL, LINES_MULTIPLE,        false, 10, NULL},
        {UNIPROT_RG,  "RG", "Reference group",              OCCURENCE_OPTIONAL, LINES_MULTIPLE,        false, 10, NULL},
        {UNIPROT_RA,  "RA", "Reference authors",            OCCURENCE_OPTIONAL, LINES_MULTIPLE,        false, 10, NULL},
        {UNIPROT_RT,  "RT", "Reference title",              OCCURENCE_OPTIONAL, LINES_MULTIPLE,        false, 10, NULL},
        {UNIPROT_RL,  "RL", "Reference location",           OCCURENCE_REQUIRED, LINES_MULTIPLE,        false, 10, NULL},
        {UNIPROT_CC,  "CC", "Comments or notes",            OCCURENCE_OPTIONAL, LINES_MULTIPLE,        true,  11, NULL},
        {UNIPROT_DR,  "DR", "Database cross-references",    OCCURENCE_OPTIONAL, LINES_MULTIPLE,        true,  12, NULL},
        {UNIPROT_PE,  "PE", "Protein existence",            OCCURENCE_REQUIRED, LINES_SINGLE,          true,  13, removeAfterFirstColon},
        {UNIPROT_KW,  "KW", "Keywords",                     OCCURENCE_OPTIONAL, LINES_MULTIPLE_FOLD,   true,  14, NULL},
        {UNIPROT_FT,  "FT", "Feature table data",           OCCURENCE_OPTIONAL, LINES_MULTIPLE,        true,  15, NULL},
        {UNIPROT_SQ,  "SQ", "Sequence header",              OCCURENCE_REQUIRED, LINES_SINGLE,          false, -1, NULL},
        {UNIPROT_SEQ, "  ", "Sequence data",                OCCURENCE_REQUIRED, LINES_MULTIPLE_CONCAT, true,  16, Util::removeWhiteSpace},
        {UNIPROT_END, "//", "Termination line",             OCCURENCE_REQUIRED, LINES_SINGLE,          false, -1, NULL}
};

bool UniprotKB::readLine(const char* line) {
    hasEntry = false;
    if (strncmp(uniprotkb_prefix[UNIPROT_ID].prefix, line, 2) == 0) {
        for (size_t i = 0; i < dbColumns; ++i) {
            streams[i].str("");
            streams[i].clear();
        }
        isInEntry = true;
    }

    if (isInEntry) {
        for (size_t i = 0; i < UNIPROT_END; ++i) {
            if (strncmp(uniprotkb_prefix[i].prefix, line, 2) == 0) {
                if(uniprotkb_prefix[i].dbColumn == -1)
                    continue;

                const char *start = line + 5;
                std::ostringstream& stream = streams[uniprotkb_prefix[i].dbColumn];
                if (uniprotkb_prefix[i].transform != NULL) {
                    stream << uniprotkb_prefix[i].transform(start);
                } else {
                    stream << start;
                }

                switch (uniprotkb_prefix[i].lines) {
                    case LINES_MULTIPLE:
                        stream << "\n";
                        break;
                    case LINES_MULTIPLE_FOLD:
                        streams[uniprotkb_prefix[i].dbColumn] << " ";
                        break;
                    default:
                        break;
                }
            }
        }
    }

    if (strncmp(uniprotkb_prefix[UNIPROT_END].prefix, line, 2) == 0) {
        isInEntry = false;
        hasEntry = true;
        return true;
    }

    return false;
}

std::string UniprotKB::getColumn(size_t column) {
    if (hasEntry == false) {
        Debug(Debug::WARNING) << "No UniprotKB entry ready to read!\n";
        return "";
    }

    if (column >= dbColumns) {
        Debug(Debug::WARNING) << "Invalid column selected\n";
        return std::string();
    }

    return streams[column].str();
}
