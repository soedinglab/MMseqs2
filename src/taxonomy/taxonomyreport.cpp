#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "FastSort.h"
#include "MappingReader.h"

#include <unordered_map>

#include "krona_prelude.html.h"

#ifdef OPENMP
#include <omp.h>
#endif

unsigned int cladeCountVal(const std::unordered_map<TaxID, TaxonCounts> &map, TaxID key) {
    std::unordered_map<TaxID, TaxonCounts>::const_iterator it = map.find(key);
    if (it == map.end()) {
        return 0;
    } else {
        return it->second.cladeCount;
    }
}

void taxReport(
    DBWriter& writer,
    unsigned int thread_idx,
    const NcbiTaxonomy &taxDB,
    const std::unordered_map<TaxID, TaxonCounts> &cladeCounts,
    unsigned long totalReads,
    TaxID taxID = 0,
    int depth = 0
) {
    std::unordered_map<TaxID, TaxonCounts>::const_iterator it = cladeCounts.find(taxID);
    unsigned int cladeCount = it == cladeCounts.end() ? 0 : it->second.cladeCount;
    unsigned int taxCount = it == cladeCounts.end() ? 0 : it->second.taxCount;
    char buffer[4096];
    int len = 0;
    if (taxID == 0) {
        if (cladeCount > 0) {
            len = snprintf(buffer, sizeof(buffer), "%.4f\t%i\t%i\tno rank\t0\tunclassified\n",
                           100 * cladeCount / double(totalReads),
                           cladeCount, taxCount);
            writer.writeAdd(buffer, static_cast<size_t>(len), thread_idx);
        }
        taxReport(writer, thread_idx, taxDB, cladeCounts, totalReads, 1);
    } else {
        if (cladeCount == 0) {
            return;
        }
        const TaxonNode *taxon = taxDB.taxonNode(taxID);
        std::string indent = std::string(2 * depth, ' ');
        len = snprintf(buffer, sizeof(buffer), "%.4f\t%i\t%i\t%s\t%i\t%s%s\n",
                       100 * cladeCount / double(totalReads),
                       cladeCount, taxCount,
                       taxDB.getString(taxon->rankIdx),
                       taxID,
                       indent.c_str(),
                       taxDB.getString(taxon->nameIdx));
        writer.writeAdd(buffer, static_cast<size_t>(len), thread_idx);
        std::vector<TaxID> children = it->second.children;
        SORT_SERIAL(children.begin(), children.end(), [&](int a, int b) { return cladeCountVal(cladeCounts, a) > cladeCountVal(cladeCounts, b); });
        for (size_t i = 0; i < children.size(); ++i) {
            TaxID childTaxId = children[i];
            if (cladeCounts.count(childTaxId)) {
                taxReport(writer, thread_idx, taxDB, cladeCounts, totalReads, childTaxId, depth + 1);
            } else {
                break;
            }
        }
    }
}

std::string escapeAttribute(const std::string &data) {
    std::string buffer;
    buffer.reserve(data.size() * 1.1);
    for (size_t i = 0; i < data.size(); ++i) {
        switch (data[i]) {
            case '&':
                buffer.append("&amp;");
                break;
            case '\"':
                buffer.append("&quot;");
                break;
            case '\'':
                buffer.append("&apos;");
                break;
            case '<':
                buffer.append("&lt;");
                break;
            case '>':
                buffer.append("&gt;");
                break;
            default:
                buffer.append(1, data[i]);
                break;
        }
    }
    return buffer;
}

void kronaReport(
    DBWriter& writer,
    unsigned int thread_idx,
    const NcbiTaxonomy &taxDB,
    const std::unordered_map<TaxID, TaxonCounts> &cladeCounts,
    unsigned long totalReads,
    TaxID taxID = 0,
    int depth = 0
) {
    std::unordered_map<TaxID, TaxonCounts>::const_iterator it = cladeCounts.find(taxID);
    unsigned int cladeCount = it == cladeCounts.end() ? 0 : it->second.cladeCount;
    char buffer[1024];
    int len = 0;
    if (taxID == 0) {
        if (cladeCount > 0) {
            len = snprintf(buffer, sizeof(buffer), "<node name=\"unclassified\"><magnitude><val>%d</val></magnitude></node>", cladeCount);
            writer.writeAdd(buffer, static_cast<size_t>(len), thread_idx);
        }
        kronaReport(writer, thread_idx, taxDB, cladeCounts, totalReads, 1);
    } else {
        if (cladeCount == 0) {
            return;
        }
        const TaxonNode *taxon = taxDB.taxonNode(taxID);
        std::string escapedName = escapeAttribute(taxDB.getString(taxon->nameIdx));
        len = snprintf(buffer, sizeof(buffer), "<node name=\"%s\"><magnitude><val>%d</val></magnitude>", escapedName.c_str(), cladeCount);
        writer.writeAdd(buffer, static_cast<size_t>(len), thread_idx);
        std::vector<TaxID> children = it->second.children;
        SORT_SERIAL(children.begin(), children.end(), [&](int a, int b) { return cladeCountVal(cladeCounts, a) > cladeCountVal(cladeCounts, b); });
        for (size_t i = 0; i < children.size(); ++i) {
            TaxID childTaxId = children[i];
            if (cladeCounts.count(childTaxId)) {
                kronaReport(writer, thread_idx, taxDB, cladeCounts, totalReads, childTaxId, depth + 1);
            } else {
                break;
            }
        }
        len = snprintf(buffer, sizeof(buffer), "</node>");
        writer.writeAdd(buffer, static_cast<size_t>(len), thread_idx);
    }
}

void mergeMaps(std::unordered_map<TaxID, unsigned int>& target, const std::unordered_map<TaxID, unsigned int>& source) {
    for (std::unordered_map<TaxID, unsigned int>::const_iterator it = source.cbegin(); it != source.cend(); ++it) {
        std::unordered_map<TaxID, unsigned int>::iterator found = target.find(it->first);
        if (found != target.end()) {
            found->second += it->second;
        } else {
            target.emplace(it->first, it->second);
        }
    }
}

int taxonomyreport(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    if (par.reportMode == Parameters::REPORT_MODE_SKIP) {
        Debug(Debug::ERROR) << "Report mode " << Parameters::REPORT_MODE_SKIP << " can only be use in workflows\n";
        EXIT(EXIT_FAILURE);
    }

    // allow reading any kind of sequence database
    const int readerDbType = FileUtil::parseDbType(par.db2.c_str());
    const bool isSequenceDB = Parameters::isEqualDbtype(readerDbType, Parameters::DBTYPE_HMM_PROFILE)
                             || Parameters::isEqualDbtype(readerDbType, Parameters::DBTYPE_AMINO_ACIDS)
                             || Parameters::isEqualDbtype(readerDbType, Parameters::DBTYPE_NUCLEOTIDES);
    if (par.reportMode == Parameters::REPORT_MODE_KRAKENDB && isSequenceDB == true) {
        Debug(Debug::ERROR) << "Cannot use Kraken DB report mode with sequence db input\n";
        EXIT(EXIT_FAILURE);
    }
    int dataMode = DBReader<unsigned int>::USE_INDEX;
    if (isSequenceDB == false) {
        dataMode |= DBReader<unsigned int>::USE_DATA;
    }
    DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str(), par.threads, dataMode);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    // support reading both LCA databases and result databases (e.g. alignment)
    const bool isTaxonomyInput = Parameters::isEqualDbtype(reader.getDbtype(), Parameters::DBTYPE_TAXONOMICAL_RESULT);
    NcbiTaxonomy *taxDB = NcbiTaxonomy::openTaxonomy(par.db1);
    MappingReader* mapping = NULL;
    if (isTaxonomyInput == false) {
        mapping = new MappingReader(par.db1);
    }

    int mode = Parameters::DBTYPE_OMIT_FILE;
    unsigned int localThreads = 1;
    if (par.reportMode == Parameters::REPORT_MODE_KRAKENDB) {
        mode = Parameters::DBTYPE_GENERIC_DB;
        localThreads = 1;
#ifdef OPENMP
        localThreads = std::max(std::min((size_t)par.threads, reader.getSize()), (size_t)1);
#endif
    }
    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), localThreads, false, mode);
    writer.open();

    std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = taxDB->getParentToChildren();

    std::unordered_map<TaxID, unsigned int> taxCounts;
    Debug::Progress progress(reader.getSize());
#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

        std::unordered_map<TaxID, unsigned int> localTaxCounts;
#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();

            if (isSequenceDB == true) {
                unsigned int taxon = mapping->lookup(reader.getDbKey(i));
                if (taxon != 0) {
                    ++localTaxCounts[taxon];
                }
                continue;
            }

            char *data = reader.getData(i, thread_idx);
            size_t entryCount = 0;
            while (*data != '\0') {
                if (isTaxonomyInput) {
                    TaxID taxon = Util::fast_atoi<int>(data);
                    ++localTaxCounts[taxon];
                } else {
                    // match dbKey to its taxon based on mapping
                    unsigned int taxon = mapping->lookup(Util::fast_atoi<unsigned int>(data));
                    if (taxon != 0) {
                        ++localTaxCounts[taxon];
                    }
                }
                entryCount++;
                data = Util::skipLine(data);
            }
            if (par.reportMode == Parameters::REPORT_MODE_KRAKENDB) {
                std::unordered_map<TaxID, TaxonCounts> cladeCounts = taxDB->getCladeCounts(localTaxCounts, parentToChildren);
                writer.writeStart(thread_idx);
                taxReport(writer, thread_idx, *taxDB, cladeCounts, entryCount);
                writer.writeEnd(reader.getDbKey(i), thread_idx);
                localTaxCounts.clear();
            }
        }
        if (par.reportMode != Parameters::REPORT_MODE_KRAKENDB) {
#pragma omp critical
{
            mergeMaps(taxCounts, localTaxCounts);
}
        }
    }

    int status = EXIT_SUCCESS;
    if (par.reportMode == Parameters::REPORT_MODE_KRAKENDB) {
        reader.close();
        writer.close(true);
    } else {
        Debug(Debug::INFO) << "Found " << taxCounts.size() << " different taxa for " << reader.getSize() << " different reads\n";
        unsigned int unknownCnt = (taxCounts.find(0) != taxCounts.end()) ? taxCounts.at(0) : 0;
        Debug(Debug::INFO) << unknownCnt << " reads are unclassified\n";
        const size_t entryCount = reader.getSize();
        reader.close();

        Debug(Debug::INFO) << "Calculating clade counts ... ";
        std::unordered_map<TaxID, TaxonCounts> cladeCounts = taxDB->getCladeCounts(taxCounts, parentToChildren);
        Debug(Debug::INFO) << " Done\n";
        if (par.reportMode == Parameters::REPORT_MODE_KRAKEN) {
            writer.writeStart(0);
            taxReport(writer, 0, *taxDB, cladeCounts, entryCount);
            writer.writeEnd(0, 0, false, false);
        } else if (par.reportMode == Parameters::REPORT_MODE_KRONA) {
            writer.writeStart(0);
            writer.writeAdd(reinterpret_cast<const char*>(krona_prelude_html), krona_prelude_html_len, 0);
            char buffer[1024];
            int len = snprintf(buffer, sizeof(buffer), "<node name=\"all\"><magnitude><val>%zu</val></magnitude>", entryCount);
            writer.writeAdd(buffer, static_cast<size_t>(len), 0);
            kronaReport(writer, 0, *taxDB, cladeCounts, entryCount);
            len = snprintf(buffer, sizeof(buffer), "</node></krona></div></body></html>");
            writer.writeAdd(buffer, static_cast<size_t>(len), 0);
            writer.writeEnd(0, 0, false, false);
        } else {
            Debug(Debug::ERROR) << "Invalid report mode " << par.reportMode << "\n";
            status = EXIT_FAILURE;
        }
        writer.close(true);
        FileUtil::remove(writer.getIndexFileName());
    }
    delete taxDB;
    return status;
}

