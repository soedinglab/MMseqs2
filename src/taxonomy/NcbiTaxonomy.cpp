// Ported from blast2lca
// https://github.com/emepyc/Blast2lca
// Originally licensed under GPLv2 or later

#include "NcbiTaxonomy.h"
#include "FileUtil.h"
#include "MathUtil.h"
#include "Debug.h"
#include "Util.h"

#include <fstream>
#include <algorithm>
#include <cassert>

int **makeMatrix(size_t maxNodes) {
    size_t dimension = maxNodes * 2;
    int **M = new int*[dimension];
    int k = (int)(MathUtil::flog2(dimension)) + 1;
    for (size_t i = 0; i < dimension; ++i) {
        M[i] = new int[k]();
    }
    return M;
}

void deleteMatrix(int** M, size_t maxNodes) {
    for (size_t i = 0; i < (maxNodes * 2); ++i) {
        delete[] M[i];
    }
    delete[] M;
}

NcbiTaxonomy::NcbiTaxonomy(const std::string &namesFile,  const std::string &nodesFile,
                           const std::string &mergedFile, const std::string &delnodesFile) {
    maxNodes = FileUtil::countLines(nodesFile.c_str());

    InitLevels();

    loadNodes(nodesFile);
    loadNames(namesFile);

    E.reserve(maxNodes * 2);
    L.reserve(maxNodes * 2);
    H = new int[maxNodes];
    std::fill(H, H + (maxNodes), 0);
    elh(1, 0);
    E.resize(maxNodes * 2, 0);
    L.resize(maxNodes * 2, 0);

    M = makeMatrix(maxNodes);
    InitRangeMinimumQuery();

    loadMerged(mergedFile);
    loadDelnodes(delnodesFile);
}

NcbiTaxonomy::~NcbiTaxonomy() {
    delete[] H;
    deleteMatrix(M, maxNodes);
}

void NcbiTaxonomy::InitLevels() {
    sortedLevels["forma"] = 1;
    sortedLevels["varietas"] = 2;
    sortedLevels["subspecies"] = 3;
    sortedLevels["species"] = 4;
    sortedLevels["species subgroup"] = 5;
    sortedLevels["species group"] = 6;
    sortedLevels["subgenus"] = 7;
    sortedLevels["genus"] = 8;
    sortedLevels["subtribe"] = 9;
    sortedLevels["tribe"] = 10;
    sortedLevels["subfamily"] = 11;
    sortedLevels["family"] = 12;
    sortedLevels["superfamily"] = 13;
    sortedLevels["parvorder"] = 14;
    sortedLevels["infraorder"] = 15;
    sortedLevels["suborder"] = 16;
    sortedLevels["order"] = 17;
    sortedLevels["superorder"] = 18;
    sortedLevels["infraclass"] = 19;
    sortedLevels["subclass"] = 20;
    sortedLevels["class"] = 21;
    sortedLevels["superclass"] = 22;
    sortedLevels["subphylum"] = 23;
    sortedLevels["phylum"] = 24;
    sortedLevels["superphylum"] = 25;
    sortedLevels["subkingdom"] = 26;
    sortedLevels["kingdom"] = 27;
    sortedLevels["superkingdom"] = 28;
}

std::vector<std::string> splitByDelimiter(const std::string &s, const std::string &delimiter, int maxCol) {
    std::vector<std::string> result;
    size_t prev = 0, pos = 0;
    int i = 0;
    do {
        pos = s.find(delimiter, prev);
        if (pos == std::string::npos) pos = s.length();
        result.emplace_back(s.substr(prev, pos - prev));
        prev = pos + delimiter.length();
        i++;
    } while (pos < s.length() && prev < s.length() && i < maxCol);

    return result;
}

void NcbiTaxonomy::loadNodes(const std::string &nodesFile) {
    std::ifstream ss(nodesFile);
    if (ss.fail()) {
        Debug(Debug::ERROR) << "File " << nodesFile << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(ss, line)) {
        std::vector<std::string> result = splitByDelimiter(line, "\t|\t", 3);
        int currentId = (int) strtol(result[0].c_str(), NULL, 10);
        int parentId = (int) strtol(result[1].c_str(), NULL, 10);

        // no circular references
        if (currentId == parentId) {
            continue;
        }

        std::map<int, TaxonNode>::iterator it = taxonTree.find(currentId);
        if (it != taxonTree.end()) {
            (*it).second.parentTaxon = parentId;
            (*it).second.rank = result[2];
        } else {
            taxonTree.emplace(currentId, TaxonNode(currentId, parentId, result[2]));
        }

        it = taxonTree.find(parentId);
        if (it != taxonTree.end()) {
            (*it).second.children.emplace_back(currentId);
        } else {
            TaxonNode node(parentId, 0, "");
            node.children.emplace_back(currentId);
            taxonTree.emplace(parentId, node);
        }
    }

//    for (std::map<int, TaxonNode>::iterator it = taxonTree.begin(); it != taxonTree.end(); ++it) {
//        TaxonNode& node = (*it).second;
//        assert(node.taxon > 0);
//        assert(node.parentTaxon > 0 || node.taxon == 1);
//    }

    int index = 1;
    updateIndices(1,  &index);
    restoreRelations();
}

void NcbiTaxonomy::updateIndices(int node, int *index) {
    std::map<int, TaxonNode>::iterator it = taxonTree.find(node);
    if (it == taxonTree.end()) {
        Debug(Debug::ERROR) << "Missing node!\n";
        EXIT(EXIT_FAILURE);
    }

    int newIndex = (*index)++;
    it->second.id = newIndex;

    std::map<int, int>::iterator it2 = D.find(node);
    if (it2 == D.end()) {
        D.emplace(node, newIndex);
    }

    for (std::vector<int>::iterator jt = it->second.children.begin();
         jt != it->second.children.end(); ++jt) {
        updateIndices(*jt, index);
    }
}

void NcbiTaxonomy::restoreRelations() {
    std::map<int, TaxonNode> helperTree;

    for (std::map<int, TaxonNode>::iterator it = taxonTree.begin(); it != taxonTree.end(); ++it) {
        TaxonNode node = it->second;
        node.parentTaxon = D[node.parentTaxon];
        for (size_t i = 0; i < node.children.size(); ++i) {
            node.children[i] = D[node.children[i]];
        }

        helperTree.emplace(node.id, node);
    }

    taxonTree.clear();
    taxonTree = helperTree;
}

std::pair<int, std::string> parseName(const std::string &line) {
    std::vector<std::string> result = splitByDelimiter(line, "\t|\t", 2);
    if (result.size() != 2) {
        Debug(Debug::ERROR) << "Invalid name entry!\n";
        EXIT(EXIT_FAILURE);
    }
    return std::make_pair((int)strtol(result[0].c_str(), NULL, 10), result[1]);
}

void NcbiTaxonomy::loadNames(const std::string &namesFile) {
    std::ifstream ss(namesFile);
    if (ss.fail()) {
        Debug(Debug::ERROR) << "File " << namesFile << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(ss, line)) {
        if (line.find("scientific name") == std::string::npos) {
            continue;
        }

        std::pair<int, std::string> entry = parseName(line);

        std::map<int, int>::iterator it1 = D.find(entry.first);
        if (it1 == D.end()) {
            Debug(Debug::ERROR) << "Invalid node!\n";
            EXIT(EXIT_FAILURE);
        }

        std::map<int, TaxonNode>::iterator it2 = taxonTree.find(it1->second);
        if (it2 == taxonTree.end()) {
            Debug(Debug::ERROR) << "Invalid node!\n";
            EXIT(EXIT_FAILURE);
        }

        assert(it1->first == entry.first);

        it2->second.name = entry.second;
    }
}

void NcbiTaxonomy::elh(int node, int level) {
    std::map<int, TaxonNode>::iterator it = taxonTree.find(node);
    if (it == taxonTree.end()) {
        Debug(Debug::ERROR) << "Invalid node!\n";
        EXIT(EXIT_FAILURE);
    }

    int id = it->first;

    assert (id > 0);

    if (H[id - 1] == 0) {
        H[id - 1] = E.size();
    }

    E.emplace_back(id);
    L.emplace_back(level);


    for (std::vector<int>::iterator jt = it->second.children.begin(); jt != it->second.children.end(); ++jt) {
        elh((*jt), level + 1);
    }
    E.emplace_back(it->second.parentTaxon);
    L.emplace_back(level - 1);
}

void NcbiTaxonomy::InitRangeMinimumQuery() {
    for (unsigned int i = 0; i < (maxNodes * 2); ++i) {
        M[i][0] = i;
    }

    for (unsigned int j = 1; (1ul << j) <= (maxNodes * 2); ++j) {
        for (unsigned int i = 0; (i + (1ul << j) - 1) < (maxNodes * 2); ++i) {
            int A = M[i][j - 1];
            int B = M[i + (1ul << (j - 1))][j - 1];
            if (L[A] < L[B]) {
                M[i][j] = A;
            } else {
                M[i][j] = B;
            }
        }
    }
}

int NcbiTaxonomy::RangeMinimumQuery(int i, int j) {
    assert(j >= i);
    int k = (int)MathUtil::flog2(j - i + 1);
    int A = M[i][k];
    int B = M[j - MathUtil::ipow<int>(2, k) + 1][k];
    if (L[A] <= L[B]) {
        return A;
    }
    return B;
}

int NcbiTaxonomy::lcaHelper(int i, int j) {
    assert(i > 0);
    assert(j > 0);
    if (i == j) {
        return i;
    }
    int v1 = H[i - 1];
    int v2 = H[j - 1];
    if (v1 > v2) {
        int tmp = v1;
        v1 = v2;
        v2 = tmp;
    }
    int rmq = RangeMinimumQuery(v1, v2);
    return E[rmq];
}

bool NcbiTaxonomy::IsAncestor(int ancestor, int child) {
    std::map<int, int>::iterator jt = D.find(ancestor);
    if (jt == D.end()) {
        Debug(Debug::ERROR) << "Invalid taxon tree: Could not find node " << ancestor << "!\n";
        EXIT(EXIT_FAILURE);
    } else {
        int value = jt->second;
        // -1 nodes was deleted (in delnodes)
        if (value == -1) {
            return false;
        }
        ancestor = value;
    }

    jt = D.find(child);
    if (jt == D.end()) {
        Debug(Debug::ERROR) << "Invalid taxon tree: Could not find node " << ancestor << "!\n";
        EXIT(EXIT_FAILURE);
    } else {
        int value = jt->second;
        if (value == -1) {
            return false;
        }
        child = value;
    }

    return lcaHelper(child, ancestor) == ancestor;
}

TaxonNode* NcbiTaxonomy::LCA(const std::vector<int>& taxa) {
    std::vector<int> indices;
    for (std::vector<int>::const_iterator it = taxa.begin(); it != taxa.end(); ++it) {
        std::map<int, int>::iterator jt = D.find(*it);
        if (jt == D.end()) {
            Debug(Debug::ERROR) << "Invalid taxon tree: Could not find node " << (*it) << "!\n";
            EXIT(EXIT_FAILURE);
        } else {
            int value = jt->second;
            // -1 nodes was deleted (in delnodes)
            if (value == -1) {
                continue;
            }
            indices.emplace_back(value);
        }
    }

    if (indices.empty()) {
        return NULL;
    }

    std::vector<int>::iterator it = indices.begin();
    int red = (*it);
    std::advance(it, 1);
    for (; it != indices.end(); ++it) {
        red = lcaHelper(red, *it);
    }

    std::map<int, TaxonNode>::iterator it2 = taxonTree.find(red);
    if (it2 == taxonTree.end()) {
        Debug(Debug::ERROR) << "Invalid taxon tree: Could not find mapping "  << red << "!\n";
        EXIT(EXIT_FAILURE);
    }

    return &(it2->second);
}


// AtRanks returns a slice of slices having the taxons at the specified taxonomic levels
std::vector<std::string> NcbiTaxonomy::AtRanks(TaxonNode *node, const std::vector<std::string> &levels) {
    std::vector<std::string> result;
    std::map<std::string, std::string> allRanks = AllRanks(node);
    int baseRankIndex = sortedLevels[node->rank];
    std::string baseRank = "uc_" + node->name;
    for (std::vector<std::string>::const_iterator it = levels.begin(); it != levels.end(); ++it) {
        std::map<std::string, std::string>::iterator jt = allRanks.find(*it);
        if (jt != allRanks.end()) {
            result.emplace_back(jt->second);
            continue;
        }

        // If not ... 2 possible causes: i) too low level ("uc_")
        if (sortedLevels[*it] < baseRankIndex) {
            result.emplace_back(baseRank);
            continue;
        }

        // ii) No taxon for the LCA at the required level -- give the first known upstream
        result.emplace_back("unknown");
    }
    return result;
}

TaxonNode* NcbiTaxonomy::Parent(int parentTaxon) {
    std::map<int, TaxonNode>::iterator it = taxonTree.find(parentTaxon);
    if (it == taxonTree.end()) {
        Debug(Debug::ERROR) << "Invalid Node!\n";
        EXIT(EXIT_FAILURE);
    }

    return &(it->second);
}


TaxonNode* NcbiTaxonomy::findNode(int taxonId) {
    std::map<int, int>::iterator it1 = D.find(taxonId);
    if (it1 == D.end()) {
        Debug(Debug::ERROR) << "Invalid node!\n";
        EXIT(EXIT_FAILURE);
    }

    if (it1->second == -1) {
        return NULL;
    }
    std::map<int, TaxonNode>::iterator it2 = taxonTree.find(it1->second);

    if (it2 == taxonTree.end()) {
        Debug(Debug::ERROR) << "Invalid Node " << taxonId << "!\n";
        EXIT(EXIT_FAILURE);
    }

    return &(it2->second);
}

std::map<std::string, std::string> NcbiTaxonomy::AllRanks(TaxonNode *node) {
    std::map<std::string, std::string> result;
    while (true) {
        if (node->taxon == 1) {
            result.emplace(node->rank, node->name);
            return result;
        }

        if (node->rank != "no_rank") {
            result.emplace(node->rank, node->name);
        }

        node = Parent(node->parentTaxon);
    }
}

void NcbiTaxonomy::loadMerged(const std::string &mergedFile) {
    std::ifstream ss(mergedFile);
    if (ss.fail()) {
        Debug(Debug::ERROR) << "File " << mergedFile << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(ss, line)) {
        std::vector<std::string> result = splitByDelimiter(line, "\t|\t", 2);
        if (result.size() != 2) {
            Debug(Debug::ERROR) << "Invalid name entry!\n";
            EXIT(EXIT_FAILURE);
        }

        unsigned int oldId = (unsigned int)strtoul(result[0].c_str(), NULL, 10);
        unsigned int mergedId = (unsigned int)strtoul(result[1].c_str(), NULL, 10);
        std::map<int, int>::iterator it = D.find(mergedId);
        if (it == D.end()) {
            Debug(Debug::ERROR) << "Invalid taxon tree: Could not map node " << mergedId << "!\n";
            EXIT(EXIT_FAILURE);
        }

        D.emplace(oldId, D[mergedId]);
    }
}

void NcbiTaxonomy::loadDelnodes(const std::string &delnodesFile) {
    std::ifstream ss(delnodesFile);
    if (ss.fail()) {
        Debug(Debug::ERROR) << "File " << delnodesFile << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(ss, line)) {
        unsigned int oldId = (unsigned int)strtoul(line.c_str(), NULL, 10);
        D.emplace(oldId, -1);
    }
}
