#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "Debug.h"
#include "NcbiTaxonomy.h"
#include "FastSort.h"
#include "MemoryMapped.h"
#include "GzReader.h"

#ifdef OPENMP
#include <omp.h>
#endif

static bool compareToFirstString(const std::pair<std::string, TaxID>& lhs, const std::pair<std::string, TaxID>& rhs){
    return (lhs.first <= rhs.first);
}

static bool sortMappingByDbKey(const std::pair<unsigned int, TaxID>& lhs, const std::pair<unsigned int, TaxID>& rhs){
    return (lhs.first <= rhs.first);
}

static bool sortByFirstString(const std::pair<std::string, TaxID>& lhs, const std::pair<std::string, TaxID>& rhs){
    return (lhs.first < rhs.first);
}

struct SortByName {
    SortByName(NcbiTaxonomy* taxonomy) : taxonomy(taxonomy) {}
    bool operator() (const TaxonNode& lhs, const TaxonNode& rhs) const {
        return strcmp(taxonomy->getString(lhs.nameIdx), taxonomy->getString(rhs.nameIdx)) <  0;
    }
    const NcbiTaxonomy* taxonomy;
};

TaxID lookupTaxID(const std::vector<std::pair<std::string, TaxID>>& mapping, const std::string& value) {
    std::pair<std::string, TaxID> val;
    val.first = value;
    std::vector<std::pair<std::string, TaxID>>::const_iterator mappingIt
        = std::upper_bound(mapping.begin(), mapping.end(), val, compareToFirstString);
    if (mappingIt == mapping.end() || mappingIt->first != val.first) {
        return 0;
    } else {
        return mappingIt->second;
    }
}

int nrtotaxmapping(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string resultDbData = par.filenames.back();
    std::string resultDbIndex = resultDbData + ".index";
    par.filenames.pop_back();

    std::string seqDbData = par.filenames.back();
    std::string seqHdrData = seqDbData + "_h";
    std::string seqHdrIndex = seqHdrData + ".index";
    par.filenames.pop_back();

    Debug::Progress progress;

    std::vector<std::pair<std::string, TaxID>> accessionMapping;
    for (size_t i = 0; i < par.filenames.size(); i++) {
        GzReader kbIn(par.filenames[i]);
        if (kbIn.fail()) {
            Debug(Debug::ERROR) << "File " << par.filenames[i] << " not found\n";
            EXIT(EXIT_FAILURE);
        }

        std::string line;
        const char *entry[255];
        while (kbIn.getline(line)) {
            progress.updateProgress();
            const size_t columns = Util::getWordsOfLine(line.c_str(), entry, 255);
            if (columns < 4) {
                Debug(Debug::ERROR) << "Invalid accession2taxid file " << par.filenames[i] << "\n";
                EXIT(EXIT_FAILURE);
            }
            std::string accession(entry[0], entry[1] - entry[0] - 1);
            unsigned int taxid = Util::fast_atoi<unsigned int>(entry[2]);
            accessionMapping.emplace_back(accession, taxid);
        }
    }
    SORT_PARALLEL(accessionMapping.begin(), accessionMapping.end(), sortByFirstString);

    NcbiTaxonomy* taxonomy = NcbiTaxonomy::openTaxonomy(seqDbData);

    // make sure to create a copy since taxonNodes is still used later
    std::vector<TaxonNode> nodesCopy(taxonomy->taxonNodes, taxonomy->taxonNodes + taxonomy->maxNodes);
    SORT_PARALLEL(nodesCopy.begin(), nodesCopy.end(), SortByName(taxonomy));

    // get a sorted list of taxa that uniquely point to a taxid
    std::vector<std::pair<std::string, TaxID>> uniqueNames;
    size_t nodesSize = nodesCopy.size();
    if (nodesSize >= 2 && taxonomy->getString(nodesCopy[0].nameIdx) != taxonomy->getString(nodesCopy[1].nameIdx)) {
        uniqueNames.emplace_back(taxonomy->getString(nodesCopy[0].nameIdx), nodesCopy[0].taxId);
    }
    for (size_t i = 1; i < (nodesSize - 1); ++i) {
        if ((taxonomy->getString(nodesCopy[i - 1].nameIdx) != taxonomy->getString(nodesCopy[i].nameIdx)) && (taxonomy->getString(nodesCopy[i].nameIdx) != taxonomy->getString(nodesCopy[i + 1].nameIdx))) {
            uniqueNames.emplace_back(taxonomy->getString(nodesCopy[i].nameIdx), nodesCopy[i].taxId);
        }
    }
    if (nodesSize > 2 && (taxonomy->getString(nodesCopy[nodesSize - 1].nameIdx) != taxonomy->getString(nodesCopy[nodesSize - 2].nameIdx))) {
        uniqueNames.emplace_back(taxonomy->getString(nodesCopy[nodesSize - 1].nameIdx), nodesCopy[nodesSize - 1].taxId);
    }
    nodesCopy.clear();

    DBReader<unsigned int> reader(seqHdrData.c_str(), seqHdrIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(resultDbData.c_str(), resultDbIndex.c_str(), par.threads, false, Parameters::DBTYPE_OMIT_FILE);
    writer.open();

    size_t processed = 0;
    size_t entries = reader.getSize();
    progress.reset(entries);
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        std::vector<TaxID> taxa;
        taxa.reserve(64);

        std::string result;

#pragma omp for schedule(dynamic, 100) reduction(+:processed)
        for (size_t i = 0; i < entries; ++i) {
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            char* data = reader.getData(i, thread_idx);

            char* start = data;
            bool isInAccession = true;
            char* startName = NULL;
            char* endName = NULL;
            bool isInSpeciesName = false;
            bool needSpeciesName = false;
            bool done = false;
            while (done == false) {
                switch (*data) {
                    case '\n':
                        // FALLTHROUGH
                    case '\0':
                        done = true;
                        // FALLTHROUGH
                    case '\1':
                        if (needSpeciesName == true && isInSpeciesName == true) {
                            std::string species(startName, endName - startName);
                            TaxID taxID = lookupTaxID(uniqueNames, species);
                            if (taxID != 0) {
                                taxa.emplace_back(taxID);
                            }
                        }
                        start = ++data;
                        isInAccession = true;
                        needSpeciesName = false;
                        isInSpeciesName = false;
                        break;
                    case '[':
                        // take last bracket with space before instead of first
                        // takes care of protein names with brackets
                        if (*(data - 1) == ' ') {
                            startName = ++data;
                            endName = data;
                            isInSpeciesName = true;
                        }
                        break;
                    case ']':
                        endName = data;
                        break;
                    // strip NR accession version
                    case '.':
                        // FALLTHROUGH
                    case ' ':
                        if (isInAccession) {
                            std::string accession(start, data - start);
                            TaxID taxID = lookupTaxID(accessionMapping, accession);
                            if (taxID != 0) {
                                taxa.emplace_back(taxID);
                            } else {
                                needSpeciesName = true;
                            }
                            isInAccession = false;
                        }
                        break;
                }
                ++data;
            }

            const TaxonNode* node = taxonomy->LCA(taxa);
            if (node != NULL) {
                result.append(SSTR(key));
                result.append(1, '\t');
                result.append(SSTR(node->taxId));
                result.append(1, '\n');
                writer.writeData(result.c_str(), result.length(), key, thread_idx, false);
                result.clear();
                processed++;
            }
            taxa.clear();
        }
    }
    writer.close(true);
    FileUtil::remove(resultDbIndex.c_str());
    reader.close();
    uniqueNames.clear();
    accessionMapping.clear();
    delete taxonomy;

    // rewrite mapping to be sorted to avoid future on-the-fly sorting
    MemoryMapped mappingUnsorted(resultDbData, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    if (!mappingUnsorted.isValid()){
        Debug(Debug::ERROR) << "Could not open mapping file " << resultDbData << "\n";
        EXIT(EXIT_FAILURE);
    }
    char* data = (char *) mappingUnsorted.getData();
    std::vector<std::pair<unsigned int, TaxID>> mapping;
    mapping.reserve(processed);
    const char *entry[255];
    progress.reset(processed);
    while (*data != '\0') {
        progress.updateProgress();
        const size_t columns = Util::getWordsOfLine(data, entry, 255);
        data = Util::skipLine(data);
        if (columns < 2) {
            Debug(Debug::ERROR) << "Invalid mapping file " << resultDbData << "\n";
            EXIT(EXIT_FAILURE);
        }
        unsigned int dbKey = Util::fast_atoi<unsigned int>(entry[0]);
        unsigned int taxId = Util::fast_atoi<unsigned int>(entry[1]);
        mapping.emplace_back(dbKey, taxId);
    }
    mappingUnsorted.close();
    SORT_PARALLEL(mapping.begin(), mapping.end(), sortMappingByDbKey);
    FILE* handle = FileUtil::openFileOrDie(resultDbData.c_str(), "w", true);
    if (handle == NULL) {
        Debug(Debug::ERROR) << "Could not write to mapping file " << resultDbData << "\n";
        EXIT(EXIT_FAILURE);
    }

    std::string result;
    result.reserve(128);
    progress.reset(processed);
    for (size_t i = 0; i < mapping.size(); ++i) {
        progress.updateProgress();
        result.append(SSTR(mapping[i].first));
        result.append(1, '\t');
        result.append(SSTR(mapping[i].second));
        result.append(1, '\n');
        int written = fwrite(result.c_str(), sizeof(char), result.size(), handle);
        if (written != (int) result.size()) {
            Debug(Debug::ERROR) << "Could not write to mapping file " << resultDbData << "\n";
            EXIT(EXIT_FAILURE);
        }
        result.clear();
    }

    if (fclose(handle) != 0) {
        Debug(Debug::ERROR) << "Could not close mapping file " << resultDbData << "\n";
        EXIT(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
