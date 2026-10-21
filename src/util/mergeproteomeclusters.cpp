#include "Util.h"
#include "Parameters.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "FastSort.h"
#include "itoa.h"

#include <algorithm>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

struct ProteomeCluster {
    DBKeyType referenceKey;
    std::vector<DBKeyType> proteomeKeys;
    std::string entry;
};

struct PassedScore {
    size_t clusterIdx;
    size_t proteomeIdx;
    float uniScore;
    float biScore;

    static bool compareByClusterIdx(const PassedScore &a, const PassedScore &b) {
        return (a.clusterIdx != b.clusterIdx) ? a.clusterIdx < b.clusterIdx : a.proteomeIdx < b.proteomeIdx;
    }
};

int mergeproteomeclusters(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    DBReader<DBKeyType> seqDbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<DBKeyType>::USE_INDEX|DBReader<DBKeyType>::USE_LOOKUP|DBReader<DBKeyType>::USE_SOURCE);
    seqDbr.open(DBReader<DBKeyType>::NOSORT);
    std::vector<unsigned int> proteinCount(seqDbr.getSourceSize(), 0);
    for (size_t i = 0; i < seqDbr.getLookupSize(); i++) {
        proteinCount[static_cast<size_t>(seqDbr.getLookupFileNumber(i))]++;
    }

    std::vector<ProteomeCluster> clusters;
    std::vector<size_t> clusterOf(proteinCount.size(), SIZE_MAX);
    char keyBuffer[255];
    for (size_t i = 3; i < par.filenames.size(); i++) {
        DBReader<DBKeyType> cluDbr(par.filenames[i].c_str(), (par.filenames[i] + ".index").c_str(), par.threads, DBReader<DBKeyType>::USE_INDEX|DBReader<DBKeyType>::USE_DATA);
        cluDbr.open(DBReader<DBKeyType>::LINEAR_ACCCESS);
        for (size_t id = 0; id < cluDbr.getSize(); id++) {
            ProteomeCluster cluster;
            cluster.referenceKey = cluDbr.getDbKey(id);
            char *data = cluDbr.getData(id, 0);
            cluster.entry = data;
            while (*data != '\0') {
                Util::parseKey(data, keyBuffer);
                cluster.proteomeKeys.push_back(Util::fast_atoi<DBKeyType>(keyBuffer));
                data = Util::skipLine(data);
            }
            clusterOf[static_cast<size_t>(cluster.referenceKey)] = clusters.size();
            clusters.push_back(cluster);
        }
        cluDbr.close();
    }

    DBReader<DBKeyType> scoreDbr(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<DBKeyType>::USE_INDEX|DBReader<DBKeyType>::USE_DATA);
    scoreDbr.open(DBReader<DBKeyType>::NOSORT);

    std::vector<std::vector<PassedScore>> passedScores(clusters.size());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        const char *entry[255];
        std::vector<PassedScore> passed;
#pragma omp for schedule(dynamic, 1)
        for (size_t clusterIdx = 0; clusterIdx < clusters.size(); clusterIdx++) {
            const std::vector<DBKeyType> &proteomeKeys = clusters[clusterIdx].proteomeKeys;
            passed.clear();
            for (size_t proteomeIdx = 0; proteomeIdx < proteomeKeys.size(); proteomeIdx++) {
                const float proteomeSize = proteinCount[static_cast<size_t>(proteomeKeys[proteomeIdx])];
                char *data = scoreDbr.getDataByDBKey(proteomeKeys[proteomeIdx], thread_idx);
                while (data != NULL && *data != '\0') {
                    Util::getWordsOfLine(data, entry, 255);
                    const DBKeyType referenceKey = Util::fast_atoi<DBKeyType>(entry[0]);
                    const size_t referenceClusterIdx = clusterOf[static_cast<size_t>(referenceKey)];
                    const float uniScore = strtof(entry[1], NULL);
                    const float referenceSize = proteinCount[static_cast<size_t>(referenceKey)];
                    const float biScore = 2.0f * std::min(uniScore * proteomeSize, referenceSize) / (referenceSize + proteomeSize);
                    if (referenceClusterIdx < clusterIdx && uniScore >= par.proteomeSimThr && biScore >= par.proteomeRelativeSimThr) {
                        passed.push_back({referenceClusterIdx, proteomeIdx, uniScore, biScore});
                    }
                    data = Util::skipLine(data);
                }
            }
            SORT_SERIAL(passed.begin(), passed.end(), PassedScore::compareByClusterIdx);
            for (size_t i = 0, j = 0; i < passed.size(); i = j) {
                while (j < passed.size() && passed[j].clusterIdx == passed[i].clusterIdx) {
                    j++;
                }
                if (j - i == proteomeKeys.size()) {
                    passedScores[clusterIdx].insert(passedScores[clusterIdx].end(), passed.begin() + i, passed.begin() + j);
                }
            }
        }
    }

    std::vector<size_t> mergedInto(clusters.size(), SIZE_MAX);
    size_t mergedCount = 0;
    char buffer[1024];
    for (size_t clusterIdx = 0; clusterIdx < clusters.size(); clusterIdx++) {
        const std::vector<PassedScore> &passed = passedScores[clusterIdx];
        const size_t proteomeCount = clusters[clusterIdx].proteomeKeys.size();
        for (size_t i = 0; i < passed.size(); i += proteomeCount) {
            const size_t referenceClusterIdx = passed[i].clusterIdx;
            if (mergedInto[referenceClusterIdx] != SIZE_MAX) {
                continue;
            }
            for (size_t k = i; k < i + proteomeCount; k++) {
                char *tmpBuffer = Itoa::u64toa_sse2(static_cast<uint64_t>(clusters[clusterIdx].proteomeKeys[passed[k].proteomeIdx]), buffer);
                *(tmpBuffer - 1) = '\t';
                tmpBuffer = Util::fastSeqIdToBuffer(passed[k].uniScore, tmpBuffer);
                *(tmpBuffer - 1) = '\t';
                tmpBuffer = Util::fastSeqIdToBuffer(passed[k].biScore, tmpBuffer);
                *(tmpBuffer - 1) = '\n';
                clusters[referenceClusterIdx].entry.append(buffer, tmpBuffer - buffer);
            }
            mergedInto[clusterIdx] = referenceClusterIdx;
            mergedCount++;
            break;
        }
    }

    int proteomeDBType = DBReader<DBKeyType>::setExtendedDbtype(Parameters::DBTYPE_GENERIC_DB, Parameters::DBTYPE_EXTENDED_SET);
    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), 1, par.compressed, proteomeDBType);
    writer.open();
    for (size_t clusterIdx = 0; clusterIdx < clusters.size(); clusterIdx++) {
        if (mergedInto[clusterIdx] == SIZE_MAX) {
            writer.writeData(clusters[clusterIdx].entry.c_str(), clusters[clusterIdx].entry.size(), clusters[clusterIdx].referenceKey, 0);
        }
    }
    writer.close();
    Debug(Debug::INFO) << "Merged " << mergedCount << " of " << clusters.size() << " proteome clusters into earlier clusters\n";

    scoreDbr.close();
    seqDbr.close();
    return EXIT_SUCCESS;
}
