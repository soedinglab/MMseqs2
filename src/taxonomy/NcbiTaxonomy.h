// Ported from blast2lca
// https://github.com/emepyc/Blast2lca
// Originally licensed under GPLv2 or later

#ifndef MMSEQS_NCBITAXONOMY_H
#define MMSEQS_NCBITAXONOMY_H

#include <map>
#include <vector>
#include <string>

struct TaxonNode {
    int id;
    int taxon;
    int parentTaxon;
    std::string rank;
    std::string name;
    std::vector<int> children;

    TaxonNode(int taxon, int parentTaxon, const std::string& rank)
            : id(taxon), taxon(taxon), parentTaxon(parentTaxon), rank(rank), name("") {};
};

class NcbiTaxonomy {
public:
    NcbiTaxonomy(const std::string &namesFile,  const std::string &nodesFile,
                 const std::string &mergedFile, const std::string &delnodesFile);
    ~NcbiTaxonomy();

    TaxonNode* LCA(const std::vector<int>& taxa);
    std::vector<std::string> AtRanks(TaxonNode *node, const std::vector<std::string> &levels);
    std::map<std::string, std::string> AllRanks(TaxonNode *node);
    bool IsAncestor(int ancestor, int child);

    TaxonNode* findNode(int taxonId);


private:
    void InitLevels();
    void loadNodes(const std::string &nodesFile);
    void updateIndices(int node, int *index);
    void restoreRelations();
    void loadNames(const std::string &namesFile);
    void elh(int node, int level);
    void InitRangeMinimumQuery();
    void loadMerged(const std::string &mergedFile);
    void loadDelnodes(const std::string &delnodesFile);

    int RangeMinimumQuery(int i, int j);
    int lcaHelper(int i, int j);
    TaxonNode* Parent(int parentTaxon);


    std::map<int, TaxonNode> taxonTree;
    std::map<int, int> D;
    std::vector<int> E;
    std::vector<int> L;
    int* H;
    int **M;
    size_t maxNodes;

    std::map<std::string, int> sortedLevels;
};

#endif
