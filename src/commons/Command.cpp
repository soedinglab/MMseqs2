#include "Command.h"
#include "Parameters.h"

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


std::vector<int> DbValidator::sequenceDb = {Parameters::DBTYPE_INDEX_DB, Parameters::DBTYPE_NUCLEOTIDES,
                                            Parameters::DBTYPE_HMM_PROFILE, Parameters::DBTYPE_AMINO_ACIDS};
std::vector<int> DbValidator::nuclDb = {Parameters::DBTYPE_NUCLEOTIDES};
std::vector<int> DbValidator::aaDb = {Parameters::DBTYPE_AMINO_ACIDS};
std::vector<int> DbValidator::prefAlnResDb =  {Parameters::DBTYPE_ALIGNMENT_RES, Parameters::DBTYPE_PREFILTER_RES};
std::vector<int> DbValidator::taxSequenceDb = {Parameters::DBTYPE_SEQTAXDB, Parameters::DBTYPE_INDEX_DB, Parameters::DBTYPE_NUCLEOTIDES,
                                               Parameters::DBTYPE_HMM_PROFILE, Parameters::DBTYPE_AMINO_ACIDS};
std::vector<int> DbValidator::allDb = {Parameters::DBTYPE_SEQTAXDB, Parameters::DBTYPE_INDEX_DB, Parameters::DBTYPE_NUCLEOTIDES,
                                      Parameters::DBTYPE_HMM_PROFILE, Parameters::DBTYPE_AMINO_ACIDS, Parameters::DBTYPE_ALIGNMENT_RES,
                                      Parameters::DBTYPE_PREFILTER_RES, Parameters::DBTYPE_PREFILTER_REV_RES, Parameters::DBTYPE_CLUSTER_RES,
                                      Parameters::DBTYPE_OFFSETDB, Parameters::DBTYPE_GENERIC_DB, Parameters::DBTYPE_TAXONOMICAL_RESULT};
std::vector<int> DbValidator::allDbAndFlat = {Parameters::DBTYPE_SEQTAXDB, Parameters::DBTYPE_INDEX_DB, Parameters::DBTYPE_NUCLEOTIDES,
                                              Parameters::DBTYPE_HMM_PROFILE, Parameters::DBTYPE_AMINO_ACIDS, Parameters::DBTYPE_ALIGNMENT_RES,
                                              Parameters::DBTYPE_PREFILTER_RES, Parameters::DBTYPE_CLUSTER_RES, Parameters::DBTYPE_GENERIC_DB,
                                              Parameters::DBTYPE_FLATFILE};
std::vector<int> DbValidator::csDb = {Parameters::DBTYPE_PROFILE_STATE_SEQ};
std::vector<int> DbValidator::ca3mDb = {Parameters::DBTYPE_CA3M_DB};
std::vector<int> DbValidator::msaDb = {Parameters::DBTYPE_MSA_DB};
std::vector<int> DbValidator::genericDb = {Parameters::DBTYPE_GENERIC_DB};
std::vector<int> DbValidator::profileDb = {Parameters::DBTYPE_HMM_PROFILE};
std::vector<int> DbValidator::prefilterDb = {Parameters::DBTYPE_PREFILTER_RES};
std::vector<int> DbValidator::clusterDb = {Parameters::DBTYPE_CLUSTER_RES};
std::vector<int> DbValidator::indexDb = {Parameters::DBTYPE_INDEX_DB};
std::vector<int> DbValidator::taxResult = {Parameters::DBTYPE_TAXONOMICAL_RESULT};
std::vector<int> DbValidator::nuclAaDb = {Parameters::DBTYPE_NUCLEOTIDES, Parameters::DBTYPE_AMINO_ACIDS};
std::vector<int> DbValidator::alignmentDb = {Parameters::DBTYPE_ALIGNMENT_RES};
std::vector<int> DbValidator::directory = {Parameters::DBTYPE_DIRECTORY};
std::vector<int> DbValidator::flatfile = {Parameters::DBTYPE_FLATFILE};
std::vector<int> DbValidator::resultDb =  {Parameters::DBTYPE_ALIGNMENT_RES, Parameters::DBTYPE_PREFILTER_RES, Parameters::DBTYPE_PREFILTER_REV_RES, Parameters::DBTYPE_CLUSTER_RES};
