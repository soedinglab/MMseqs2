#include "Command.h"
#include "Parameters.h"

std::vector<Categories> categories = {
        {"Easy workflows for plain text input/output",   COMMAND_EASY},
        {"Main workflows for database input/output",     COMMAND_MAIN},
        {"Input database creation",                      COMMAND_DATABASE_CREATION},
        {"Handle databases on storage and memory",       COMMAND_STORAGE},
        {"Unite and intersect databases",                COMMAND_SET},
        {"Format conversion for downstream processing",  COMMAND_FORMAT_CONVERSION},
        {"Sequence manipulation/transformation",         COMMAND_SEQUENCE},
        {"Result manipulation",                          COMMAND_RESULT},
        {"Taxonomy assignment",                          COMMAND_TAXONOMY},
        {"Multi-hit search",                             COMMAND_MULTIHIT},
        {"Prefiltering",                                 COMMAND_PREFILTER},
        {"Alignment",                                    COMMAND_ALIGNMENT},
        {"Clustering",                                   COMMAND_CLUSTER},
        {"Profile databases",                            COMMAND_PROFILE},
        {"Profile-profile databases",                    COMMAND_PROFILE_PROFILE},
        {"Utility modules to manipulate DBs",            COMMAND_DB},
        {"Special-purpose utilities",                    COMMAND_SPECIAL},
};


std::vector<int> DbValidator::sequenceDb = {Parameters::DBTYPE_INDEX_DB, Parameters::DBTYPE_NUCLEOTIDES,
                                            Parameters::DBTYPE_HMM_PROFILE, Parameters::DBTYPE_AMINO_ACIDS};
std::vector<int> DbValidator::nuclDb = {Parameters::DBTYPE_NUCLEOTIDES};
std::vector<int> DbValidator::aaDb = {Parameters::DBTYPE_AMINO_ACIDS};
std::vector<int> DbValidator::prefAlnResDb =  {Parameters::DBTYPE_ALIGNMENT_RES, Parameters::DBTYPE_PREFILTER_RES};
std::vector<int> DbValidator::taxSequenceDb = {Parameters::DBTYPE_SEQTAXDB, Parameters::DBTYPE_INDEX_DB, Parameters::DBTYPE_NUCLEOTIDES,
                                               Parameters::DBTYPE_HMM_PROFILE, Parameters::DBTYPE_AMINO_ACIDS};
std::vector<int> DbValidator::allDb = {Parameters::DBTYPE_SEQTAXDB, Parameters::DBTYPE_INDEX_DB, Parameters::DBTYPE_NUCLEOTIDES, Parameters::DBTYPE_MSA_DB,
                                      Parameters::DBTYPE_HMM_PROFILE, Parameters::DBTYPE_AMINO_ACIDS, Parameters::DBTYPE_ALIGNMENT_RES,
                                      Parameters::DBTYPE_PREFILTER_RES, Parameters::DBTYPE_PREFILTER_REV_RES, Parameters::DBTYPE_CLUSTER_RES,
                                      Parameters::DBTYPE_OFFSETDB, Parameters::DBTYPE_GENERIC_DB, Parameters::DBTYPE_TAXONOMICAL_RESULT};
std::vector<int> DbValidator::allDbAndFlat = {Parameters::DBTYPE_SEQTAXDB, Parameters::DBTYPE_INDEX_DB, Parameters::DBTYPE_NUCLEOTIDES, Parameters::DBTYPE_MSA_DB,
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
std::vector<int> DbValidator::flatfileAndStdin = {Parameters::DBTYPE_FLATFILE, Parameters::DBTYPE_STDIN};
std::vector<int> DbValidator::resultDb =  {Parameters::DBTYPE_ALIGNMENT_RES, Parameters::DBTYPE_PREFILTER_RES, Parameters::DBTYPE_PREFILTER_REV_RES, Parameters::DBTYPE_CLUSTER_RES};
std::vector<int> DbValidator::empty = {};
