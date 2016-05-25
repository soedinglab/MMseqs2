#ifndef MMSEQS_COMPAREGOTERMS_H
#define MMSEQS_COMPAREGOTERMS_H

#include "DBReader.h"

#include <string>
#include <set>
#include <list>
#include <map>

class CompareGOTerms {
public:
    void init();

    CompareGOTerms(const std::string &go_ffindex, const std::string &go_ffindex_indexfile,
                   const std::string &protid_go_ffindex, const std::string &protid_go_ffindex_indexfile,
                   const std::string &evaluationfolder, const std::string &sequencedb,
                   bool usesequencedb);

//runmodes
    void all_against_all_comparison();

    void all_against_all_comparison_proteinset();

    double compare_protein_ids(const char *prot1, const char *prot2);

    void run_evaluation_mmseqsclustering(const std::string &cluster_ffindex,
                                         const std::string &cluster_ffindex_indexfile,
                                         const std::string &fileprefix,
                                         const std::string &filesuffix,
                                         bool allagainstall, bool randomized);

    ~CompareGOTerms();

private:

    DBReader<std::string> *go_ffindex_reader;
    DBReader<std::string> *protid_go_ffindex_reader;
    DBReader<std::string> *targetdb_header;

    bool usesequencedb;

    std::string evaluationfolder;

    size_t total_go_number;

    /* int * * go_to_index;
     int * go_to_indexData;

     int * index_to_go;
     size_t number_of_goterms;
 */
    //go graph
    size_t is_a_relation_number;
    int **is_a_relation;
    int *is_a_relation_data;
    int *is_a_relation_size;
    //counts of goterms
    int *count_goterm;

    int *count_accumulated_goterm;
    double count_goterm_total_sum;

    std::set<int> *parentsets;
    //mapping
    //gomapping
    std::map<int, int> goterm_to_index;
    std::map<int, int> index_togoterm;
    //proteinmapping

    //utils
    int convert_GOterm_to_index(int Goterm);

    int convert_index_toGOterm(int index);

    void compute_parentnodes(std::set<int> &result, int id);

    //scorecomputation
    double similarity(int id1, int id2);

    int most_specific_parent(int id1, int id2);

    double similarity_of_list(const std::list<int> &ids1, const std::list<int> &ids2);

    std::list<int> getGOListforProtein(const char *protid);

    std::string getProteinNameForID(const char *dbKey);
};

#endif //MMSEQS_COMPAREGOTERMS_H
