#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "IndexReader.h"
#include "MemoryMapped.h"
#include "NcbiTaxonomy.h"
#include "MappingReader.h"
#include "algorithm"
#include <iterator>
#include <map>

#ifdef OPENMP
#include <omp.h>
#endif

#ifndef SIZE_T_MAX
#define SIZE_T_MAX ((size_t) -1)
#endif


#define ABS_DIFF(a, b)  ( ((a) > (b)) ? ((a) - (b)) : ((b) - (a)) )

class UniProtConverter {
public:
    size_t toStructuredNumber(std::string &uniprotId) const {
        if (uniprotId.compare(0, 6,  "UniRef") == 0) {
            std::vector<std::string> parts = Util::split(uniprotId, "_");
            if (parts.size() > 1) {
                uniprotId = parts[1];
            }
        }
        if (uniprotId.find('-') != std::string::npos) {
            uniprotId = uniprotId.substr(0, uniprotId.find('-'));
        }

        if (uniprotId.empty()) return 0;

        const size_t len = uniprotId.length();
        const char firstChar = static_cast<char>(std::toupper(uniprotId[0]));

        if (len == 6 && (firstChar == 'O' || firstChar == 'P' || firstChar == 'Q')) {
            return convertOpqPattern(uniprotId);
        }

        if ((len == 6 || len == 10)) {
            return convertAnrzPattern(uniprotId);
        }

        if (uniprotId[0] == 'U' && uniprotId[1] == 'P' && uniprotId[2] == 'I') {
            static constexpr size_t UPI_OFFSET = 1000000000000000ULL;
            std::string hex_part(uniprotId.substr(3));
            return UPI_OFFSET + std::stoll(hex_part, nullptr, 16);
        }
        //TODO need to do this better
        return 0;
    }

private:
    size_t convertOpqPattern(const std::string & id) const {
        size_t number = 0;
        size_t multiplier = 1;

        for (int i = 5; i >= 0; --i) {
            const char currentChar = static_cast<char>(std::toupper(id[i]));
            int charVal = -1;
            int radix = 0;
            switch (i) {
                case 0: charVal = getOpqValue(currentChar); radix = 3; break;
                case 1: case 5: charVal = getDigitValue(currentChar); radix = 10; break;
                case 2: case 3: case 4: charVal = getAlphanumValue(currentChar); radix = 36; break;
            }
            if (charVal == -1) {
                return 0;
            }
            number += static_cast<size_t>(charVal) * multiplier;
            multiplier *= radix;
        }
        return number;
    }

    size_t convertAnrzPattern(const std::string & id) const {
        size_t number = 0;
        size_t multiplier = 1;
        size_t len = id.length();

        for (int i = len - 1; i >= 0; --i) {
            const char currentChar = static_cast<char>(std::toupper(id[i]));
            int charVal = -1;
            int radix = 0;

            // --- Logic using a switch statement ---
            switch (i) {
                case 0:
                    charVal = getAnrzValue(currentChar); radix = 23;
                    break;
                case 1: case 5: case 9:
                    charVal = getDigitValue(currentChar); radix = 10;
                    break;
                case 2: case 6:
                    charVal = getAlphaValue(currentChar); radix = 26;
                    break;
                case 3: case 4: case 7: case 8:
                    charVal = getAlphanumValue(currentChar); radix = 36;
                    break;
                default:
                    return 0;
            }
            if (charVal == -1) return 0;

            number += static_cast<size_t>(charVal) * multiplier;
            multiplier *= radix;
        }
        return number;
    }

    // --- Static Character-to-Value Helpers ---
    static int getAlphanumValue(char c) { if (c >= '0' && c <= '9') return c - '0'; if (c >= 'A' && c <= 'Z') return c - 'A' + 10; return -1; }
    static int getAlphaValue(char c) { if (c >= 'A' && c <= 'Z') return c - 'A'; return -1; }
    static int getDigitValue(char c) { if (c >= '0' && c <= '9') return c - '0'; return -1; }
    static int getOpqValue(char c) { if (c == 'O') return 0; if (c == 'P') return 1; if (c == 'Q') return 2; return -1; }
    static int getAnrzValue(char c) { if (c >= 'A' && c <= 'N') return c - 'A'; if (c >= 'R' && c <= 'Z') return c - 'A' - 3; return -1; }
};

struct CompareUniProt {
    bool operator()(const Matcher::result_t& res, uint64_t num) const {
        return getUniProtNumber(res) < num;
    }
    static size_t getUniProtNumber(const Matcher::result_t& res){
        return (static_cast<uint64_t>(res.queryOrfStartPos) << 32) | static_cast<uint32_t>(res.queryOrfEndPos);
    };
};

struct CompareByTaxon {
    bool operator()(const Matcher::result_t& res, int tax) const { return res.dbOrfStartPos < tax; }
    bool operator()(int tax, const Matcher::result_t& res) const { return tax < res.dbOrfStartPos; }
};

size_t findNearestPartner(const Matcher::result_t& query, const std::vector<Matcher::result_t>& results2)
{
    typedef std::vector<Matcher::result_t>::const_iterator ResultConstIterator;
//    std::pair<ResultConstIterator, ResultConstIterator> range;
//    range = std::equal_range(results2.begin(), results2.end(), taxon_id, CompareByTaxon());
//
//    if (range.first == range.second) {
//        return SIZE_T_MAX;
//    }

    size_t bestIdx = SIZE_T_MAX;
    size_t min_dist = SIZE_T_MAX;

    size_t query_uniprot_num = CompareUniProt::getUniProtNumber(query);
    ResultConstIterator it2 = std::lower_bound(results2.begin(), results2.end(), query_uniprot_num, CompareUniProt());
    if (it2 != results2.end()) {
        size_t target_uniprot_num = CompareUniProt::getUniProtNumber(*it2);
        size_t dist = ABS_DIFF(target_uniprot_num, query_uniprot_num);

        if (dist < min_dist) {
            min_dist = dist;
            bestIdx = std::distance(results2.begin(), it2);
        }
    }

    // Check the element just before the found position.
    if (it2 != results2.begin()) {
        ResultConstIterator prev_it2 = it2;
        --prev_it2; // Use pre-decrement to get the previous iterator.

        size_t target_uniprot_num = CompareUniProt::getUniProtNumber(*prev_it2);
        size_t dist = ABS_DIFF(query_uniprot_num, target_uniprot_num);

        if (dist < min_dist) {
            bestIdx = std::distance(results2.begin(), prev_it2);
        }
    }
    return bestIdx;
}

static bool compareByTaxId(const Matcher::result_t &first, const Matcher::result_t &second) {
    return (first.dbOrfStartPos < second.dbOrfStartPos);
}

static bool compareByTaxIdAndUniProtNum(const Matcher::result_t &first, const Matcher::result_t &second) {
//    if (first.dbOrfStartPos != second.dbOrfStartPos) {
//        return first.dbOrfStartPos < second.dbOrfStartPos;
//    }

    if (first.queryOrfStartPos != second.queryOrfStartPos) {
        return first.queryOrfStartPos < second.queryOrfStartPos;
    }
    return static_cast<unsigned int>(first.queryOrfEndPos) < static_cast<unsigned int>(second.queryOrfEndPos);
}

int pairaln(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> qdbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_LOOKUP_REV);
    qdbr.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int>::LookupEntry* lookup = qdbr.getLookup();
    unsigned int maxFileNumber = 0;
    for (unsigned int i = 0; i < qdbr.getLookupSize(); i++) {
        maxFileNumber = std::max(maxFileNumber, lookup[i].fileNumber);
    }
    //build a mapping from file number to ids from lookup
    std::vector<std::vector<unsigned int>> fileToIds(maxFileNumber + 1, std::vector<unsigned int>());
    for (size_t i = 0; i < qdbr.getLookupSize(); ++i) {
        fileToIds[lookup[i].fileNumber].push_back(lookup[i].id);
    }
    IndexReader *targetHeaderReaderIdx = NULL;
    if(par.pairfilter == Parameters::PAIRALN_FILTER_PROXIMITY) {
        uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(FileUtil::parseDbType(par.db3.c_str()));
        bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
        targetHeaderReaderIdx = new IndexReader(par.db2, par.threads,
                                                extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC
                                                ? IndexReader::SRC_HEADERS : IndexReader::HEADERS,
                                                (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    }

    std::string db2NoIndexName = PrefilteringIndexReader::dbPathWithoutIndex(par.db2);
    MappingReader* mapping = new MappingReader(db2NoIndexName);

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::NOSORT);

    size_t localThreads = 1;
#ifdef OPENMP
    localThreads = std::max(std::min((size_t)par.threads, fileToIds.size()), (size_t)1);
#endif

    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), localThreads, par.compressed, alnDbr.getDbtype());
    resultWriter.open();

    Debug::Progress progress(fileToIds.size());
#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        std::vector<Matcher::result_t> result;
        std::vector<Matcher::result_t> compatible;
//        std::vector<Matcher::result_t> bestCompatible;

        result.reserve(100000);
        std::unordered_map<unsigned int, size_t> findPair;
        std::vector<unsigned int> taxonToPair;
        std::string output;
        output.reserve(100000);
        bool hasBacktrace = false;
        UniProtConverter converter;
        unsigned int minResultDbKey = UINT_MAX;
        Matcher::result_t emptyResult(UINT_MAX, 0, 0, 0, 0, 0,
                                      0, UINT_MAX, 0, 0, UINT_MAX, 0, 0, "1M");
#pragma omp for schedule(dynamic, 1)
        for (size_t fileNumber = 0; fileNumber < fileToIds.size(); fileNumber++) {
            char buffer[1024 + 32768 * 4];
            findPair.clear();
            taxonToPair.clear();
            progress.updateProgress();

            // find intersection between all proteins
            for (size_t i = 0; i < fileToIds[fileNumber].size(); i++) {
                result.clear();
                size_t id = fileToIds[fileNumber][i];
                Matcher::readAlignmentResults(result, alnDbr.getData(id, thread_idx), true);
                for (size_t resIdx = 0; resIdx < result.size(); ++resIdx) {
                    hasBacktrace = result[resIdx].backtrace.size() > 0;
                    unsigned int taxon = mapping->lookup(result[resIdx].dbKey);
                    // we don't want to introduce a new field, reuse existing unused field here
                    result[resIdx].dbOrfStartPos = taxon;

                    minResultDbKey = std::min(result[resIdx].dbKey, minResultDbKey);
                }
                std::stable_sort(result.begin(), result.end(), compareByTaxId);
                unsigned int prevTaxon = UINT_MAX;
                // find pairs
                for (size_t resIdx = 0; resIdx < result.size(); ++resIdx) {
                    //  dbOrfStartPos is the taxon here, see above.
                    unsigned int taxon = result[resIdx].dbOrfStartPos;
                    if (taxon == prevTaxon) {
                        continue;
                    }
                    if (findPair.find(taxon) != findPair.end()) {
                        findPair[taxon]++;
                    } else {
                        findPair.emplace(taxon, 1);
                    }
                    prevTaxon = taxon;
                }
            }
            emptyResult.dbKey = minResultDbKey;
            // fill taxonToPair vector
            std::unordered_map<unsigned int, size_t>::iterator it;
            for (it = findPair.begin(); it != findPair.end(); ++it) {
                size_t thresholdToPair = (par.pairmode == Parameters::PAIRALN_MODE_ALL_PER_SPECIES) ? 1 : fileToIds[fileNumber].size() - 1;
                if (it->second > thresholdToPair) {
                    taxonToPair.emplace_back(it->first);
                }
            }
            std::sort(taxonToPair.begin(), taxonToPair.end());
            if(par.pairfilter == Parameters::PAIRALN_FILTER_PROXIMITY){
                std::vector<std::vector<Matcher::result_t>> resultPerId(fileToIds[fileNumber].size());
                std::vector<std::string> outputs(fileToIds[fileNumber].size());
                std::vector<size_t> foundIds;

                for (size_t i = 0; i < fileToIds[fileNumber].size(); i++) {
                    resultPerId[i].clear();
                    size_t id = fileToIds[fileNumber][i];
                    Matcher::readAlignmentResults(resultPerId[i], alnDbr.getData(id, thread_idx), true);
                    // find pairs
                    for (size_t resIdx = 0; resIdx < resultPerId[i].size(); ++resIdx) {
                        unsigned int taxon = mapping->lookup(resultPerId[i][resIdx].dbKey);
                        // we don't want to introduce a new field, reuse existing unused field here
                        resultPerId[i][resIdx].dbOrfStartPos = taxon;
                        size_t headerId = targetHeaderReaderIdx->sequenceReader->getId(resultPerId[i][resIdx].dbKey);
                        char *headerData = targetHeaderReaderIdx->sequenceReader->getData(headerId, thread_idx);
                        std::string targetAccession = Util::parseFastaHeader(headerData);
                        size_t uniProtNumber = converter.toStructuredNumber(targetAccession);
                        resultPerId[i][resIdx].queryOrfStartPos = static_cast<int>(uniProtNumber >> 32);
                        resultPerId[i][resIdx].queryOrfEndPos = static_cast<int>(uniProtNumber & 0xFFFFFFFF);
                    }
                    std::stable_sort(resultPerId[i].begin(), resultPerId[i].end(), compareByTaxIdAndUniProtNum);
                }
                outputs.resize(resultPerId.size());
                for(size_t i = 0; i < resultPerId.size(); i++) {
                    outputs[i].clear();
                }

                // iterate over taxonToPair
//                for (size_t taxonInList: taxonToPair) {
                    //bestCompatible.clear();
                    //bestCompatible.resize(resultPerId.size(), Matcher::result_t());
//                    int bestCompatibleSize = 0;
//                    typedef std::vector<Matcher::result_t>::const_iterator ResultIterator;
//                    std::pair<ResultIterator, ResultIterator> range;
//                    range = std::equal_range(resultPerId[0].cbegin(), resultPerId[0].cend(), taxonInList, CompareByTaxon());
//                    size_t startIndex = std::distance(resultPerId[0].cbegin(), range.first);
//                    size_t endIndex = std::distance(resultPerId[0].cbegin(), range.second);
                    for(size_t resIdx = 0; resIdx < resultPerId[0].size(); ++resIdx) {
                        // check if pairable
                        compatible.clear();
                        compatible.resize(resultPerId.size(), Matcher::result_t());
                        compatible[0] = resultPerId[0][resIdx];
                        int compatibleSize = 1;
                        for (size_t i = 1; i < resultPerId.size(); ++i) {
                            size_t partnerIdx = findNearestPartner(compatible[0],
                                                                   resultPerId[i]);
                            if (partnerIdx == SIZE_T_MAX) { // no partner found
                                if (par.pairdummymode == Parameters::PAIRALN_DUMMY_MODE_ON) {
                                    compatible[i] = emptyResult;
                                }
                                continue;
                            }

                            size_t currNum = CompareUniProt::getUniProtNumber(resultPerId[i][partnerIdx]);
                            bool isCompatible = false;

                            for (size_t j = 0; j < compatible.size(); ++j) {
                                if (compatible[j].dbKey == UINT_MAX) continue;   // not set yet
                                size_t prevNum = CompareUniProt::getUniProtNumber(compatible[j]);
                                size_t diff = ABS_DIFF(currNum, prevNum);
                                if (diff <= static_cast<size_t>(par.pairProximityDistance)) {
                                    isCompatible = true;
                                    break;
                                }
                            }

                            if (isCompatible) {
                                compatible[i] = resultPerId[i][partnerIdx];
                                compatibleSize++;
                            } else {
                                compatible[i] = emptyResult;
                            }
                        }

                        if((par.pairmode == Parameters::PAIRALN_MODE_COVER_ALL_CHAINS && compatibleSize != static_cast<int>(resultPerId.size()))
                            || compatibleSize == 1) {
                            continue;
                        }

                        for (size_t i = 0; i < compatible.size(); i++) {
                            if (compatible[i].dbKey == UINT_MAX &&
                                par.pairdummymode != Parameters::PAIRALN_DUMMY_MODE_ON) {
                                continue;
                            }
                            size_t len = Matcher::resultToBuffer(buffer, compatible[i], hasBacktrace, false,
                                                                 false);
                            outputs[i].append(buffer, len);
                        }
                    }

//                }
                for(size_t i = 0; i < resultPerId.size(); i++) {
                    resultWriter.writeData(outputs[i].c_str(), outputs[i].length(),
                                           alnDbr.getDbKey(fileToIds[fileNumber][i]), thread_idx);
                }
            } else {
                for (size_t i = 0; i < fileToIds[fileNumber].size(); i++) {
                    result.clear();
                    output.clear();
                    size_t id = fileToIds[fileNumber][i];
                    Matcher::readAlignmentResults(result, alnDbr.getData(id, thread_idx), true);
                    // find pairs
                    for (size_t resIdx = 0; resIdx < result.size(); ++resIdx) {
                        unsigned int taxon = mapping->lookup(result[resIdx].dbKey);
                        // we don't want to introduce a new field, reuse existing unused field here
                        result[resIdx].dbOrfStartPos = taxon;
                    }

                    // stable sort is required to assure that best hit is first per taxon
                    std::stable_sort(result.begin(), result.end(), compareByTaxId);
                    // find hit in proximity based on uniProtNumber within a taxon and sort by compareByTaxId + closest pairs first.


                    unsigned int prevTaxon = UINT_MAX;
                    size_t resIdxStart = 0;
                    for (size_t taxonInList: taxonToPair) {
                        bool taxonFound = false;
                        for (size_t resIdx = resIdxStart; resIdx < result.size(); ++resIdx) {
                            unsigned int taxon = result[resIdx].dbOrfStartPos;
                            // check if this taxon has enough information to pair
                            if (taxonInList != taxon) {
                                continue;
                            }
                            bool bestTaxonHit = (taxon != prevTaxon);
                            taxonFound = true;
                            if (bestTaxonHit) {
                                size_t len = Matcher::resultToBuffer(buffer, result[resIdx], hasBacktrace, false,
                                                                     false);
                                output.append(buffer, len);
                                resIdxStart = resIdx + 1;
                                break;
                            }
                            prevTaxon = taxon;
                        }
                        if (taxonFound == false &&
                            par.pairdummymode == Parameters::PAIRALN_DUMMY_MODE_ON) { // par.addDummyPairedAlignment
                            // write an Matcher::result_t with UINT_MAX as dbKey
                            size_t len = Matcher::resultToBuffer(buffer, emptyResult, hasBacktrace, false, false);
                            output.append(buffer, len);
                        }
                    }

                    resultWriter.writeData(output.c_str(), output.length(), alnDbr.getDbKey(id), thread_idx);
                }
            }
        }
    }
    resultWriter.close();
    qdbr.close();
    // clean up
    delete mapping;
    if(par.pairfilter == Parameters::PAIRALN_FILTER_PROXIMITY) {
        delete targetHeaderReaderIdx;
    }
    alnDbr.close();
    return EXIT_SUCCESS;
}

#undef ABS_DIFF
#undef SIZE_T_MAX