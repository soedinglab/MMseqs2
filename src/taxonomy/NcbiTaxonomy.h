#ifndef MMSEQS_NCBITAXONOMY_H
#define MMSEQS_NCBITAXONOMY_H

#include "StringBlock.h"

#include <map>
#include <unordered_map>
#include <vector>
#include <string>

typedef int TaxID;

struct TaxonNode {
public:
    int id;
    TaxID taxId;
    TaxID parentTaxId;
    size_t rankIdx;
    size_t nameIdx;

    TaxonNode() {};

    TaxonNode(int id, TaxID taxId, TaxID parentTaxId, size_t rankIdx, size_t nameIdx)
            : id(id), taxId(taxId), parentTaxId(parentTaxId), rankIdx(rankIdx), nameIdx(nameIdx) {};
};

const double MAX_TAX_WEIGHT = 1000;
struct WeightedTaxHit {
    WeightedTaxHit(const TaxID taxon, const float evalue, const int weightVoteMode);

    TaxID taxon;
    double weight;
};

struct WeightedTaxResult {
    WeightedTaxResult(TaxID taxon, size_t assignedSeqs, size_t unassignedSeqs, size_t seqsAgreeWithSelectedTaxon, double selectedPercent)
            : taxon(taxon), assignedSeqs(assignedSeqs), unassignedSeqs(unassignedSeqs), seqsAgreeWithSelectedTaxon(seqsAgreeWithSelectedTaxon), selectedPercent(selectedPercent) {};

    TaxID  taxon;
    size_t assignedSeqs;
    size_t unassignedSeqs;
    size_t seqsAgreeWithSelectedTaxon;
    double selectedPercent;
};

struct TaxonCounts {
    unsigned int taxCount;       // number of reads/sequences matching to taxa
    unsigned int cladeCount;     // number of reads/sequences matching to taxa or its children
    std::vector<TaxID> children; // list of children
};

static const std::map<std::string, int> NcbiRanks = {{ "forma", 1 },
                                                     { "varietas", 2 },
                                                     { "subspecies", 3 },
                                                     { "species", 4 },
                                                     { "species subgroup", 5 },
                                                     { "species group", 6 },
                                                     { "subgenus", 7 },
                                                     { "genus", 8 },
                                                     { "subtribe", 9 },
                                                     { "tribe", 10 },
                                                     { "subfamily", 11 },
                                                     { "family", 12 },
                                                     { "superfamily", 13 },
                                                     { "parvorder", 14 },
                                                     { "infraorder", 15 },
                                                     { "suborder", 16 },
                                                     { "order", 17 },
                                                     { "superorder", 18 },
                                                     { "infraclass", 19 },
                                                     { "subclass", 20 },
                                                     { "class", 21 },
                                                     { "superclass", 22 },
                                                     { "subphylum", 23 },
                                                     { "phylum", 24 },
                                                     { "superphylum", 25 },
                                                     { "subkingdom", 26 },
                                                     { "kingdom", 27 },
                                                     { "superkingdom", 28 }};

static const std::map<std::string, char> NcbiShortRanks = {{ "species", 's' },
                                                           { "genus", 'g' },
                                                           { "family", 'f' },
                                                           { "order", 'o' },
                                                           { "class", 'c' },
                                                           { "phylum", 'p' },
                                                           { "kingdom", 'k' },
                                                           { "superkingdom", 'd' }};

class NcbiTaxonomy {
public:
    static NcbiTaxonomy* openTaxonomy(const std::string &database);
    NcbiTaxonomy(const std::string &namesFile,  const std::string &nodesFile, const std::string &mergedFile);
    ~NcbiTaxonomy();

    TaxonNode const * LCA(const std::vector<TaxID>& taxa) const;
    TaxID LCA(TaxID taxonA, TaxID taxonB) const;
    std::vector<std::string> AtRanks(TaxonNode const * node, const std::vector<std::string> &levels) const;
    std::map<std::string, std::string> AllRanks(TaxonNode const *node) const;
    std::string taxLineage(TaxonNode const *node, bool infoAsName = true);

    static std::vector<std::string> parseRanks(const std::string& ranks);
    static int findRankIndex(const std::string& rank);
    static char findShortRank(const std::string& rank);

    bool IsAncestor(TaxID ancestor, TaxID child);
    TaxonNode const* taxonNode(TaxID taxonId, bool fail = true) const;
    bool nodeExists(TaxID taxId) const;

    std::unordered_map<TaxID, std::vector<TaxID>> getParentToChildren() const;
    std::unordered_map<TaxID, TaxonCounts> getCladeCounts(const std::unordered_map<TaxID, unsigned int>& taxonCounts, const std::unordered_map<TaxID, std::vector<TaxID>>& parentToChildren) const;

    WeightedTaxResult weightedMajorityLCA(const std::vector<WeightedTaxHit> &setTaxa, const float majorityCutoff);

    const char* getString(size_t blockIdx) const;

    static std::pair<char*, size_t> serialize(const NcbiTaxonomy& taxonomy);
    static NcbiTaxonomy* unserialize(char* data);

    TaxonNode* taxonNodes;
    size_t maxNodes;
    int maxTaxID;
private:
    size_t loadNodes(std::vector<TaxonNode> &tmpNodes, const std::string &nodesFile);
    size_t loadMerged(const std::string &mergedFile);
    void loadNames(std::vector<TaxonNode> &tmpNodes, const std::string &namesFile);
    void elh(std::vector<std::vector<TaxID>> const & children, int node, int level, std::vector<int> &tmpE, std::vector<int> &tmpL);
    void computeSparseTable();
    int nodeId(TaxID taxId) const;

    int RangeMinimumQuery(int i, int j) const;
    int lcaHelper(int i, int j) const;

    NcbiTaxonomy(TaxonNode* taxonNodes, size_t maxNodes, int maxTaxID, int *D, int *E, int *L, int *H, int **M, StringBlock<unsigned int> *block)
        : taxonNodes(taxonNodes), maxNodes(maxNodes), maxTaxID(maxTaxID), D(D), E(E), L(L), H(H), M(M), block(block), externalData(true), mmapData(NULL), mmapSize(0) {};
    int *D; // maps from taxID to node ID in taxonNodes
    int *E; // for Euler tour sequence (size 2N-1)
    int *L; // Level of nodes in tour sequence (size 2N-1)
    int *H;
    int **M;
    StringBlock<unsigned int>* block;

    bool externalData;
    char* mmapData;
    size_t mmapSize;

    static const int SERIALIZATION_VERSION;
};

#endif
