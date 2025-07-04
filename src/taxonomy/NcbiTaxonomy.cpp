#include "NcbiTaxonomy.h"
#include "FileUtil.h"
#include "MathUtil.h"
#include "Debug.h"
#include "Util.h"
#include "sys/mman.h"

#include <fstream>
#include <algorithm>
#include <cassert>

const int NcbiTaxonomy::SERIALIZATION_VERSION = 2;

int **makeMatrix(size_t maxNodes) {
    size_t dimension = maxNodes * 2;
    int **M = new int*[dimension];
    int k = (int)(MathUtil::flog2(dimension)) + 1;
    M[0] = new int[dimension * k]();
    for(size_t i = 1; i < dimension; i++) {
        M[i] = M[i-1] + k;
    }

    return M;
}

void deleteMatrix(int** M) {
    delete[] M[0];
    delete[] M;
}

NcbiTaxonomy::NcbiTaxonomy(const std::string &namesFile, const std::string &nodesFile, const std::string &mergedFile) : externalData(false) {
    block = new StringBlock<unsigned int>();
    std::vector<TaxonNode> tmpNodes;
    loadNodes(tmpNodes, nodesFile);
    loadMerged(mergedFile);
    loadNames(tmpNodes, namesFile);

    maxNodes = tmpNodes.size();
    taxonNodes = new TaxonNode[maxNodes];
    std::copy(tmpNodes.begin(), tmpNodes.end(), taxonNodes);

    std::vector<int> tmpE;
    tmpE.reserve(maxNodes * 2);

    std::vector<int> tmpL;
    tmpL.reserve(maxNodes * 2);

    H = new int[maxNodes];
    std::fill(H, H + maxNodes, 0);

    std::vector<std::vector<TaxID>> children(tmpNodes.size());
    for (std::vector<TaxonNode>::const_iterator it = tmpNodes.begin(); it != tmpNodes.end(); ++it) {
        if (it->parentTaxId != it->taxId) {
            children[nodeId(it->parentTaxId)].push_back(it->taxId);
        }
    }

    elh(children, 1, 0, tmpE, tmpL);
    tmpE.resize(maxNodes * 2, 0);
    tmpL.resize(maxNodes * 2, 0);

    E = new int[maxNodes * 2];
    std::copy(tmpE.begin(), tmpE.end(), E);
    L = new int[maxNodes * 2];
    std::copy(tmpL.begin(), tmpL.end(), L);

    M = makeMatrix(maxNodes);
    computeSparseTable();

    mmapData = NULL;
    mmapSize = 0;
}

NcbiTaxonomy::~NcbiTaxonomy() {
    if (externalData) {
        delete[] M;
    } else {
        delete[] taxonNodes;
        delete[] H;
        delete[] D;
        delete[] E;
        delete[] L;
        deleteMatrix(M);
    }
    delete block;
    if (mmapData != NULL) {
        munmap(mmapData, mmapSize);
    }
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

size_t NcbiTaxonomy::loadNodes(std::vector<TaxonNode> &tmpNodes, const std::string &nodesFile) {
    Debug(Debug::INFO) << "Loading nodes file ...";
    std::ifstream ss(nodesFile);
    if (ss.fail()) {
        Debug(Debug::ERROR) << "File " << nodesFile << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    std::map<TaxID, int> Dm; // temporary map TaxID -> internal ID;
    maxTaxID = 0;
    int currentId = 0;
    std::string line;
    while (std::getline(ss, line)) {
        std::vector<std::string> result = splitByDelimiter(line, "\t|\t", 3);
        TaxID taxId = (TaxID) strtol(result[0].c_str(), NULL, 10);
        TaxID parentTaxId = (TaxID) strtol(result[1].c_str(), NULL, 10);
        if (taxId > maxTaxID) {
            maxTaxID = taxId;
        }
        size_t rankIdx = block->append(result[2].c_str(), result[2].size());
        tmpNodes.emplace_back(currentId, taxId, parentTaxId, rankIdx, (size_t)-1);
        Dm.emplace(taxId, currentId);
        ++currentId;
    }

    D = new int[maxTaxID + 1];
    std::fill_n(D, maxTaxID + 1, -1);
    for (std::map<TaxID, int>::iterator it = Dm.begin(); it != Dm.end(); ++it) {
        assert(it->first <= maxTaxID);
        D[it->first] = it->second;
    }

    // Loop over taxonNodes and check all parents exist
    for (std::vector<TaxonNode>::iterator it = tmpNodes.begin(); it != tmpNodes.end(); ++it) {
        if (!nodeExists(it->parentTaxId)) {
            Debug(Debug::ERROR) << "Inconsistent nodes.dmp taxonomy file! Cannot find parent taxon with ID " << it->parentTaxId << "!\n";
            EXIT(EXIT_FAILURE);
        }
    }

    Debug(Debug::INFO) << " Done, got " << tmpNodes.size() << " nodes\n";
    return tmpNodes.size();
}

std::pair<int, std::string> parseName(const std::string &line) {
    std::vector<std::string> result = splitByDelimiter(line, "\t|\t", 2);
    if (result.size() != 2) {
        Debug(Debug::ERROR) << "Invalid name entry!\n";
        EXIT(EXIT_FAILURE);
    }
    return std::make_pair((int)strtol(result[0].c_str(), NULL, 10), result[1]);
}

void NcbiTaxonomy::loadNames(std::vector<TaxonNode> &tmpNodes, const std::string &namesFile) {
    Debug(Debug::INFO) << "Loading names file ...";
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
        if (!nodeExists(entry.first)) {
            Debug(Debug::ERROR) << "loadNames: Taxon " << entry.first << " not present in nodes file!\n";
            EXIT(EXIT_FAILURE);
        }
        tmpNodes[nodeId(entry.first)].nameIdx = block->append(entry.second.c_str(), entry.second.size());
    }
    Debug(Debug::INFO) << " Done\n";
}

// Euler traversal of tree
void NcbiTaxonomy::elh(std::vector<std::vector<TaxID>> const & children, TaxID taxId, int level, std::vector<int> &tmpE, std::vector<int> &tmpL) {
    assert (taxId > 0);
    int id = nodeId(taxId);

    if (H[id] == 0) {
        H[id] = tmpE.size();
    }

    tmpE.emplace_back(id);
    tmpL.emplace_back(level);

    for (std::vector<TaxID>::const_iterator child_it = children[id].begin(); child_it != children[id].end(); ++child_it) {
        elh(children, *child_it, level + 1, tmpE, tmpL);
    }
    tmpE.emplace_back(nodeId(taxonNodes[id].parentTaxId));
    tmpL.emplace_back(level - 1);
}

void NcbiTaxonomy::computeSparseTable() {
    Debug(Debug::INFO) << "Init computeSparseTable ...";
    // sparse table M has N rows and log(N) columns.
    // M[i][j] refers to the subarray L[i..2^j]
    // M[i][j] holds the index of the minimal value in the subarray
    size_t N = maxNodes * 2; // TO DO - I think this can actually be changed to maxNodes!!!
    // Debug(Debug::INFO) << "N: " << N << "\n";

    // size_t helperCount = 0;

    // initialize all rows for column 0
    for (size_t row_ind = 0; row_ind < N; row_ind++) {
        M[row_ind][0] = row_ind;
        // helperCount++;
    }

    // fill in column after column
    size_t col_ind = 1;
    size_t exp_prev_col_ind = 1ul; // 2 ^ 0
    size_t exp_col_ind = (1ul << col_ind);

    while (exp_col_ind <= N) {   
        size_t row_ind = 0;
        while (row_ind + exp_col_ind - 1 < N) {
            int min_ind_first_half = M[row_ind][col_ind - 1];
            int min_ind_second_half = M[row_ind + exp_prev_col_ind][col_ind - 1];
            if (L[min_ind_first_half] < L[min_ind_second_half]) {
                M[row_ind][col_ind] = min_ind_first_half;
                // helperCount++;
            } else {
                M[row_ind][col_ind] = min_ind_second_half;
                // helperCount++;
            }
            // increase row_ind
            row_ind = row_ind + 1;
        }
        // increase col_ind
        col_ind = col_ind + 1;
        exp_prev_col_ind = exp_col_ind;
        exp_col_ind = (1ul << col_ind);
    }
    // Debug(Debug::INFO) << "updated cells of M: " << helperCount << "\n";
    // Debug(Debug::INFO) << "last used exponent: " << exp_prev_col_ind << "\n";
    // Debug(Debug::INFO) << "last unused exponent: " << exp_col_ind << "\n";
    // Debug(Debug::INFO) << "col_ind: " << col_ind << "\n";
    Debug(Debug::INFO) << "Done\n";
}

int NcbiTaxonomy::RangeMinimumQuery(int i, int j) const {
    assert(j >= i);
    int k = (int)MathUtil::flog2(j - i + 1);
    int A = M[i][k];
    int B = M[j - MathUtil::ipow<int>(2, k) + 1][k];
    if (L[A] <= L[B]) {
        return A;
    }
    return B;
}

int NcbiTaxonomy::lcaHelper(int i, int j) const {
    if (i == 0 || j == 0) {
        return 0;
    }
    assert(i > 0);
    assert(j > 0);
    if (i == j) {
        return i;
    }
    int v1 = H[i];
    int v2 = H[j];
    if (v1 > v2) {
        int tmp = v1;
        v1 = v2;
        v2 = tmp;
    }
    int rmq = RangeMinimumQuery(v1, v2);
    assert(E[rmq] >= 0);
    return E[rmq];
}

bool NcbiTaxonomy::IsAncestor(TaxID ancestor, TaxID child) {
    if (ancestor == child) {
        return true;
    }

    if (ancestor == 0 || child == 0) {
        return false;
    }

    if (!nodeExists(child)) {
        return false;
    }

    if (!nodeExists(ancestor)) {
        return false;
    }

    return lcaHelper(nodeId(child), nodeId(ancestor)) == nodeId(ancestor);
}


TaxID NcbiTaxonomy::LCA(TaxID taxonA, TaxID taxonB) const {
    if (!nodeExists(taxonA)) {
        return taxonB;
    } else if (!nodeExists(taxonB)) {
        return taxonA;
    }
    return taxonNodes[lcaHelper(nodeId(taxonA), nodeId(taxonB))].taxId;
}


TaxonNode const * NcbiTaxonomy::LCA(const std::vector<TaxID>& taxa) const {
    std::vector<int>::const_iterator it = taxa.begin();
    while (it != taxa.end() && !nodeExists(*it)) {
        Debug(Debug::WARNING) << "No node for taxID " << *it << ", ignoring it.\n";
        ++it;
    }
    if (it == taxa.end()) { return NULL; }
    int red = nodeId(*it++);
    for (; it != taxa.end(); ++it) {
        if (nodeExists(*it)) {
            red = lcaHelper(red, nodeId(*it));
        } else {
            Debug(Debug::WARNING) << "No node for taxID " << *it << ", ignoring it.\n";
        }
    }

    assert(red >= 0 && static_cast<unsigned int>(red) < maxNodes);

    return &(taxonNodes[red]);
}


// AtRanks returns a slice of slices having the taxons at the specified taxonomic levels
std::vector<std::string> NcbiTaxonomy::AtRanks(TaxonNode const *node, const std::vector<std::string> &levels) const {
    std::vector<std::string> result;
    std::map<std::string, std::string> allRanks = AllRanks(node);
    // map does not include "no rank" nor "no_rank"
    const char* rank = getString(node->rankIdx);
    int baseRankIndex = findRankIndex(rank);
    std::string baseRank = "uc_";
    baseRank.append(getString(node->nameIdx));
    for (std::vector<std::string>::const_iterator it = levels.begin(); it != levels.end(); ++it) {
        std::map<std::string, std::string>::iterator jt = allRanks.find(*it);
        if (jt != allRanks.end()) {
            result.emplace_back(jt->second);
            continue;
        }

        // If not ... 2 possible causes: i) too low level ("uc_")
        if (NcbiRanks.at(*it) < baseRankIndex) {
            result.emplace_back(baseRank);
            continue;
        }

        // ii) No taxon for the LCA at the required level -- give the first known upstream
        result.emplace_back("unknown");
    }
    return result;
}

std::vector<std::string> NcbiTaxonomy::parseRanks(const std::string& ranks) {
    std::vector<std::string> temp = Util::split(ranks, ",");
    for (size_t i = 0; i < temp.size(); ++i) {
        if (findRankIndex(temp[i]) == -1) {
            Debug(Debug::ERROR) << "Invalid taxonomic rank " << temp[i] << "given\n";
            EXIT(EXIT_FAILURE);
        }
    }
    return temp;
}

int NcbiTaxonomy::findRankIndex(const std::string& rank) {
    std::map<std::string, int>::const_iterator it;
    if ((it = NcbiRanks.find(rank)) != NcbiRanks.end()) {
        return it->second;
    }
    return -1;
}

char NcbiTaxonomy::findShortRank(const std::string& rank) {
    std::map<std::string, char>::const_iterator it;
    if ((it = NcbiShortRanks.find(rank)) != NcbiShortRanks.end()) {
        return it->second;
    }
    return '-';
}

std::string NcbiTaxonomy::taxLineage(TaxonNode const *node, bool infoAsName) {
    std::vector<TaxonNode const *> taxLineageVec;
    std::string taxLineage;
    taxLineage.reserve(4096);
    do {
        taxLineageVec.push_back(node);
        node = taxonNode(node->parentTaxId);
    } while (node->parentTaxId != node->taxId);

    for (int i = taxLineageVec.size() - 1; i >= 0; --i) {
        if (infoAsName) {
            taxLineage += findShortRank(getString(taxLineageVec[i]->rankIdx));
            taxLineage += '_';
            taxLineage += getString(taxLineageVec[i]->nameIdx);
        } else {
            taxLineage += SSTR(taxLineageVec[i]->taxId);
        }

        if (i > 0) {
            taxLineage += ";";
        }
    }
    return taxLineage;
}

int NcbiTaxonomy::nodeId(TaxID taxonId) const {
    if (taxonId < 0 || !nodeExists(taxonId)) {
        Debug(Debug::ERROR) << "Invalid node " << taxonId << "!\n";
        EXIT(EXIT_FAILURE);
    }
    return D[taxonId];
}

bool NcbiTaxonomy::nodeExists(TaxID taxonId) const {
    return taxonId <= maxTaxID && D[taxonId] != -1;
}

TaxonNode const * NcbiTaxonomy::taxonNode(TaxID taxonId, bool fail) const {
    if (taxonId == 0 || (!fail && !nodeExists(taxonId))) {
        return NULL;
    }
    return &(taxonNodes[nodeId(taxonId)]);
}

std::map<std::string, std::string> NcbiTaxonomy::AllRanks(TaxonNode const *node) const {
    std::map<std::string, std::string> result;
    while (true) {
        std::string rank = getString(node->rankIdx);
        std::string name = getString(node->nameIdx);
        if (node->taxId == 1) {
            result.emplace(rank, name);
            return result;
        }

        if ((rank != "no_rank") && (rank != "no rank")) {
            result.emplace(rank, name);
        }

        node = taxonNode(node->parentTaxId);
    }
}

size_t NcbiTaxonomy::loadMerged(const std::string &mergedFile) {
    Debug(Debug::INFO) << "Loading merged file ...";
    std::ifstream ss(mergedFile);
    if (ss.fail()) {
        Debug(Debug::ERROR) << "File " << mergedFile << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    std::unordered_map<TaxID, TaxID> mergedMap;
    TaxID localMaxTaxID = maxTaxID;
    std::string line;
    while (std::getline(ss, line)) {
        std::vector<std::string> result = splitByDelimiter(line, "\t|\t", 2);
        if (result.size() != 2) {
            Debug(Debug::ERROR) << "Invalid name entry!\n";
            EXIT(EXIT_FAILURE);
        }

        TaxID oldId = (TaxID) strtoul(result[0].c_str(), NULL, 10);
        TaxID mergedId = (TaxID) strtoul(result[1].c_str(), NULL, 10);

        // Only update if the oldId doesn't exist yet AND the mergedId does exist
        if (!nodeExists(oldId) && nodeExists(mergedId)) {
            if (oldId > localMaxTaxID) {
                localMaxTaxID = oldId;
            }
            if (mergedId > localMaxTaxID) {
                localMaxTaxID = mergedId;
            }
            mergedMap[oldId] = mergedId;
        }
    }

    // realloc D if we find a higher maxTaxID
    if (localMaxTaxID > maxTaxID) {
        int* newD = new int[localMaxTaxID + 1];
        std::copy(D, D + maxTaxID + 1, newD);
        std::fill_n(newD + maxTaxID + 1, localMaxTaxID - maxTaxID, -1);
        delete[] D;
        D = newD;
        maxTaxID = localMaxTaxID;
    }

    size_t count = 0;
    for (std::unordered_map<TaxID, TaxID>::iterator it = mergedMap.begin(); it != mergedMap.end(); ++it) {
        D[it->first] = D[it->second];
        ++count;
    }
    Debug(Debug::INFO) << " Done, added " << count << " merged nodes.\n";
    return count;
}

std::unordered_map<TaxID, std::vector<TaxID>> NcbiTaxonomy::getParentToChildren() const {
    std::unordered_map<TaxID, std::vector<TaxID>> result;
    result.reserve(maxNodes);
    
    // Build the adjacency (parent -> children)
    for (size_t i = 0; i < maxNodes; ++i) {
        const TaxonNode& tn = taxonNodes[i];
        if (tn.parentTaxId == tn.taxId) {
            continue;
        }
        result[tn.parentTaxId].push_back(tn.taxId);
    }

    return result;
}

std::unordered_map<TaxID, TaxonCounts> NcbiTaxonomy::getCladeCounts(const std::unordered_map<TaxID, unsigned int>& taxonCounts, const std::unordered_map<TaxID, std::vector<TaxID>>& parentToChildren) const {
    std::unordered_map<TaxID, TaxonCounts> cladeCounts;

    for (std::unordered_map<TaxID, unsigned int>::const_iterator it = taxonCounts.begin(); it != taxonCounts.end(); ++it) {
        cladeCounts[it->first].taxCount = it->second;
        cladeCounts[it->first].cladeCount += it->second;
        if (nodeExists(it->first)) {
            TaxonNode const* taxon = taxonNode(it->first);
            while (taxon->parentTaxId != taxon->taxId && nodeExists(taxon->parentTaxId)) {
                taxon = taxonNode(taxon->parentTaxId);
                cladeCounts[taxon->taxId].cladeCount += it->second;
            }
        }
    }

   for (std::unordered_map<TaxID, TaxonCounts>::iterator it = cladeCounts.begin(); it != cladeCounts.end(); ++it) {
        TaxID parentTaxId = it->first;
        TaxonCounts& taxCounts = it->second;
        std::unordered_map<TaxID, std::vector<TaxID>>::const_iterator ptcIt = parentToChildren.find(parentTaxId);
        if (ptcIt != parentToChildren.end()) {
            taxCounts.children = ptcIt->second;
        }
    }

    return cladeCounts;
}

NcbiTaxonomy * NcbiTaxonomy::openTaxonomy(const std::string &database){
    std::string binFile = database + "_taxonomy";
    if (FileUtil::fileExists(binFile.c_str())) {
        FILE* handle = fopen(binFile.c_str(), "r");
        struct stat sb;
        if (fstat(fileno(handle), &sb) < 0) {
            Debug(Debug::ERROR) << "Failed to fstat file " << binFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        char* data = (char*)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fileno(handle), 0);
        if (data == MAP_FAILED){
            Debug(Debug::ERROR) << "Failed to mmap file " << binFile << " with error " << errno << "\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(handle);
        NcbiTaxonomy* t = NcbiTaxonomy::unserialize(data);
        if (t != NULL) {
            t->mmapData = data;
            t->mmapSize = sb.st_size;
            return t;
        } else {
            Debug(Debug::WARNING) << "Outdated taxonomy information, please recreate with createtaxdb.\n";
        }
    }
    Debug(Debug::INFO) << "Loading NCBI taxonomy\n";
    std::string nodesFile = database + "_nodes.dmp";
    std::string namesFile = database + "_names.dmp";
    std::string mergedFile = database + "_merged.dmp";
    if (FileUtil::fileExists(nodesFile.c_str())
        && FileUtil::fileExists(namesFile.c_str())
        && FileUtil::fileExists(mergedFile.c_str())) {
    } else if (FileUtil::fileExists("nodes.dmp")
               && FileUtil::fileExists("names.dmp")
               && FileUtil::fileExists("merged.dmp")) {
        nodesFile = "nodes.dmp";
        namesFile = "names.dmp";
        mergedFile = "merged.dmp";
    } else {
        Debug(Debug::ERROR) << "names.dmp, nodes.dmp, merged.dmp from NCBI taxdump could not be found!\n";
        EXIT(EXIT_FAILURE);
    }
    return new NcbiTaxonomy(namesFile, nodesFile, mergedFile);
}

const TaxID ROOT_TAXID = 1;
const int ROOT_RANK = INT_MAX;

struct TaxNode {
    TaxNode(const double weight, const bool isCandidate, const TaxID childTaxon)
            : weight(weight), isCandidate(isCandidate), childTaxon(childTaxon) {}

    void update(const double weightToAdd, const TaxID & childTaxonInput) {
        if (childTaxon != childTaxonInput) {
            isCandidate = true;
            childTaxon = childTaxonInput;
        }
        weight += weightToAdd;
    }

    double weight;
    bool isCandidate;
    TaxID childTaxon;
};

const char* NcbiTaxonomy::getString(size_t blockIdx) const {
    return block->getString(blockIdx);
}

WeightedTaxHit::WeightedTaxHit(const TaxID taxon, const float evalue, const int weightVoteMode) : taxon(taxon) {
    switch (weightVoteMode) {
        case Parameters::AGG_TAX_UNIFORM:
            weight = 1.0;
            break;
        case Parameters::AGG_TAX_MINUS_LOG_EVAL:
            weight = evalue;
            if (evalue != FLT_MAX) {
                if (evalue > 0) {
                    weight = -log(evalue);
                } else {
                    weight = MAX_TAX_WEIGHT;
                }
            }
            break;
        case Parameters::AGG_TAX_SCORE:
            weight = evalue;
            break;
        default:
            Debug(Debug::ERROR) << "Invalid weight vote mode\n";
            EXIT(EXIT_FAILURE);
    }
}

WeightedTaxResult NcbiTaxonomy::weightedMajorityLCA(const std::vector<WeightedTaxHit> &setTaxa, const float majorityCutoff) {
    // count num occurences of each ancestor, possibly weighted
    std::map<TaxID, TaxNode> ancTaxIdsCounts;

    // initialize counters and weights
    size_t assignedSeqs = 0;
    size_t unassignedSeqs = 0;
    double totalAssignedSeqsWeights = 0.0;

    for (size_t i = 0; i < setTaxa.size(); ++i) {
        TaxID currTaxId = setTaxa[i].taxon;
        double currWeight = setTaxa[i].weight;
        // ignore unassigned sequences
        if (currTaxId == 0) {
            unassignedSeqs++;
            continue;
        }
        TaxonNode const *node = taxonNode(currTaxId, false);
        if (node == NULL) {
            unassignedSeqs++;
            continue;
        }
        totalAssignedSeqsWeights += currWeight;
        assignedSeqs++;

        // each start of a path due to an orf is a candidate
        std::map<TaxID, TaxNode>::iterator it;
        if ((it = ancTaxIdsCounts.find(currTaxId)) != ancTaxIdsCounts.end()) {
            it->second.update(currWeight, 0);
        } else {
            TaxNode current(currWeight, true, 0);
            ancTaxIdsCounts.emplace(currTaxId, current);
        }

        // iterate all ancestors up to root (including). add currWeight and candidate status to each
        TaxID currParentTaxId = node->parentTaxId;
        while (currParentTaxId != currTaxId) {
            if ((it = ancTaxIdsCounts.find(currParentTaxId)) != ancTaxIdsCounts.end()) {
                it->second.update(currWeight, currTaxId);
            } else {
                TaxNode parent(currWeight, false, currTaxId);
                ancTaxIdsCounts.emplace(currParentTaxId, parent);
            }
            // move up
            currTaxId = currParentTaxId;
            node = taxonNode(currParentTaxId, false);
            currParentTaxId = node->parentTaxId;
        }
    }

    TaxID selctedTaxon = 0;
    if (totalAssignedSeqsWeights == 0) {
        return WeightedTaxResult(selctedTaxon, assignedSeqs, unassignedSeqs, 0, 0.0);
    }

    // select the lowest ancestor that meets the cutoff
    int minRank = INT_MAX;
    double selectedPercent = 0;
    for (std::map<TaxID, TaxNode>::iterator it = ancTaxIdsCounts.begin(); it != ancTaxIdsCounts.end(); it++) {
        // consider only candidates
        if (it->second.isCandidate == false) {
            continue;
        }

        double currPercent = it->second.weight / totalAssignedSeqsWeights;
        if (currPercent >= majorityCutoff) {
            // iterate all ancestors to find lineage min rank (the candidate is a descendant of a node with this rank)
            TaxID currTaxId = it->first;
            TaxonNode const *node = taxonNode(currTaxId, false);
            int currMinRank = ROOT_RANK;
            TaxID currParentTaxId = node->parentTaxId;
            while (currParentTaxId != currTaxId) {
                int currRankInd = NcbiTaxonomy::findRankIndex(getString(node->rankIdx));
                if ((currRankInd > 0) && (currRankInd < currMinRank)) {
                    currMinRank = currRankInd;
                    // the rank can only go up on the way to the root, so we can break
                    break;
                }
                // move up:
                currTaxId = currParentTaxId;
                node = taxonNode(currParentTaxId, false);
                currParentTaxId = node->parentTaxId;
            }

            if ((currMinRank < minRank) || ((currMinRank == minRank) && (currPercent > selectedPercent))) {
                selctedTaxon = it->first;
                minRank = currMinRank;
                selectedPercent = currPercent;
            }
        }
    }

    // count the number of seqs who have selectedTaxon in their ancestors (agree with selection):
    if (selctedTaxon == ROOT_TAXID) {
        // all agree with "root"
        return WeightedTaxResult(selctedTaxon, assignedSeqs, unassignedSeqs, assignedSeqs, selectedPercent);
    }
    if (selctedTaxon == 0) {
        // nothing informative
        return WeightedTaxResult(selctedTaxon, assignedSeqs, unassignedSeqs, 0, selectedPercent);
    }
    size_t seqsAgreeWithSelectedTaxon = 0;
    // otherwise, iterate over all seqs
    for (size_t i = 0; i < setTaxa.size(); ++i) {
        TaxID currTaxId = setTaxa[i].taxon;
        // ignore unassigned sequences
        if (currTaxId == 0) {
            continue;
        }
        TaxonNode const *node = taxonNode(currTaxId, false);
        if (node == NULL) {
            continue;
        }

        // iterate all ancestors up to the root
        TaxID currParentTaxId = node->parentTaxId;
        while (currParentTaxId != currTaxId) {
            if (currTaxId == selctedTaxon) {
                seqsAgreeWithSelectedTaxon++;
                break;
            }
            currTaxId = currParentTaxId;
            node = taxonNode(currParentTaxId, false);
            currParentTaxId = node->parentTaxId;
        }
    }

    return WeightedTaxResult(selctedTaxon, assignedSeqs, unassignedSeqs, seqsAgreeWithSelectedTaxon, selectedPercent);
}

std::pair<char*, size_t> NcbiTaxonomy::serialize(const NcbiTaxonomy& t) {
    t.block->compact();
    size_t matrixDim = (t.maxNodes * 2);
    size_t matrixK = (int)(MathUtil::flog2(matrixDim)) + 1;
    size_t matrixSize = matrixDim * matrixK * sizeof(int);
    size_t blockSize = StringBlock<unsigned int>::memorySize(*t.block);
    size_t memSize = sizeof(int) // SERIALIZATION_VERSION
        + sizeof(size_t) // maxNodes
        + sizeof(int) // maxTaxID
        + t.maxNodes * sizeof(TaxonNode) // taxonNodes
        + (t.maxTaxID + 1) * sizeof(int) // D
        + 2 * (t.maxNodes * 2) * sizeof(int) // E,L
        + t.maxNodes * sizeof(int) // H
        + matrixSize // M
        + blockSize; // block

    char* mem = (char*) malloc(memSize);
    char* p = mem;
    memcpy(p, &t.SERIALIZATION_VERSION, sizeof(int));
    p += sizeof(int);
    memcpy(p, &t.maxNodes, sizeof(size_t));
    p += sizeof(size_t);
    memcpy(p, &t.maxTaxID, sizeof(int));
    p += sizeof(int);
    memcpy(p, t.taxonNodes, t.maxNodes * sizeof(TaxonNode));
    p += t.maxNodes * sizeof(TaxonNode);
    memcpy(p, t.D, (t.maxTaxID + 1) * sizeof(int));
    p += (t.maxTaxID + 1) * sizeof(int);
    memcpy(p, t.E, (t.maxNodes * 2) * sizeof(int));
    p += (t.maxNodes * 2) * sizeof(int);
    memcpy(p, t.L, (t.maxNodes * 2) * sizeof(int));
    p += (t.maxNodes * 2) * sizeof(int);
    memcpy(p, t.H, t.maxNodes * sizeof(int));
    p += t.maxNodes * sizeof(int);
    memcpy(p, t.M[0], matrixSize);
    p += matrixSize;
    char* blockData = StringBlock<unsigned int>::serialize(*t.block);
    memcpy(p, blockData, blockSize);
    p += blockSize;
    free(blockData);
    return std::make_pair(mem, memSize);
}

NcbiTaxonomy* NcbiTaxonomy::unserialize(char* mem) {
    const char* p = mem;
    int version = *((int*)p);
    p += sizeof(int);
    if (version != NcbiTaxonomy::SERIALIZATION_VERSION) {
        return NULL;
    }
    size_t maxNodes = *((size_t*)p);
    p += sizeof(size_t);
    int maxTaxID = *((int*)p);
    p += sizeof(int);
    TaxonNode* taxonNodes = (TaxonNode*)p;
    p += maxNodes * sizeof(TaxonNode);
    int* D = (int*)p;
    p += (maxTaxID + 1) * sizeof(int);
    int* E = (int*)p;
    p += (maxNodes * 2) * sizeof(int);
    int* L = (int*)p;
    p += (maxNodes * 2) * sizeof(int);
    int* H = (int*)p;
    p += maxNodes * sizeof(int);
    size_t matrixDim = (maxNodes * 2);
    size_t matrixK = (int)(MathUtil::flog2(matrixDim)) + 1;
    size_t matrixSize = matrixDim * matrixK * sizeof(int);
    int** M = new int*[matrixDim];
    M[0] = (int*)p;
    for(size_t i = 1; i < matrixDim; i++) {
        M[i] = M[i-1] + matrixK;
    }
    p += matrixSize;
    StringBlock<unsigned int>* block = StringBlock<unsigned int>::unserialize(p);
    return new NcbiTaxonomy(taxonNodes, maxNodes, maxTaxID, D, E, L, H, M, block);
}
