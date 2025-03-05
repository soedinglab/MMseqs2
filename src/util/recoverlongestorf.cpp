#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "Orf.h"

#include <unordered_map>
#include <unordered_set>

#ifdef OPENMP
#include <omp.h>
#endif

int recoverlongestorf(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> headerReader(par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    headerReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBReader<unsigned int> resultReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), 1, false, Parameters::DBTYPE_OMIT_FILE);
    writer.open();

    Debug::Progress progress(resultReader.getSize());
    std::unordered_map<unsigned int, std::pair<unsigned int, size_t>> contigToLongest;
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());

#endif 
        std::unordered_map<unsigned int, std::pair<unsigned int, size_t>> localContigToLongest;

#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < headerReader.getSize(); ++id) {
            progress.updateProgress();
            unsigned int orfKey = headerReader.getDbKey(id);
            char *orfHeader = headerReader.getData(id, thread_idx);
            Orf::SequenceLocation orf = Orf::parseOrfHeader(orfHeader);
            unsigned int contigKey = orf.id; 
            size_t orfLen = std::max(orf.from, orf.to) - std::min(orf.from, orf.to) + 1;
            std::unordered_map<unsigned int, std::pair<unsigned int, size_t>>::iterator it = localContigToLongest.find(contigKey);
            if (it != localContigToLongest.end()) {
                std::pair<unsigned int, size_t> orfKeyToLength = it->second;
                if (orfLen > orfKeyToLength.second) {
                    it->second = std::make_pair(orfKey, orfLen);
                }
            } else {
                localContigToLongest.emplace(contigKey, std::make_pair(orfKey, orfLen));
            }
        }
        
#pragma omp critical
        {
            for (auto &entry : localContigToLongest) {
                auto &contigKey = entry.first;
                auto &localData = entry.second;
                auto it = contigToLongest.find(contigKey);
                if (it != contigToLongest.end()) {
                    if (localData.second > it->second.second) {
                        it->second = localData;
                    }
                } else {
                    contigToLongest[contigKey] = localData;
                }
            }
        }
    }

    progress.reset(resultReader.getSize());
    std::unordered_set<unsigned int> acceptedContigs;    
    std::unordered_set<unsigned int> eliminatedContigs;
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        std::string resultBuffer;
        resultBuffer.reserve(1024 * 1024);

        std::unordered_set<unsigned int> localAcceptedContigs;
        std::unordered_set<unsigned int> localEliminatedContigs;
#pragma omp for schedule(dynamic, 5)
        for (size_t i = 0; i < resultReader.getSize(); ++i) {
            progress.updateProgress();

            unsigned int key = resultReader.getDbKey(i);
            size_t entryLength = resultReader.getEntryLen(i);
            if (entryLength > 1) {
                size_t id = headerReader.getId(key);
                char *orfHeader = headerReader.getData(id, thread_idx);
                Orf::SequenceLocation orf = Orf::parseOrfHeader(orfHeader);
                unsigned int contigKey = orf.id; 
                localAcceptedContigs.emplace(contigKey);
            }

            size_t id = headerReader.getId(key);
            char *orfHeader = headerReader.getData(id, thread_idx);
            Orf::SequenceLocation orf = Orf::parseOrfHeader(orfHeader);
            unsigned int contigKey = orf.id; 
            localEliminatedContigs.emplace(contigKey);
        }

#pragma omp critical
        {
            acceptedContigs.insert(localAcceptedContigs.begin(), localAcceptedContigs.end());
            eliminatedContigs.insert(localEliminatedContigs.begin(), localEliminatedContigs.end());
        }
    }

    for (auto it = eliminatedContigs.begin(); it != eliminatedContigs.end(); ) {
        if (acceptedContigs.find(*it) != acceptedContigs.end()) {
            it = eliminatedContigs.erase(it);
        } else {
            ++it;
        }
    }

    std::string resultBuffer;
    resultBuffer.reserve(1024 * 1024);
    for (auto contigIt = eliminatedContigs.begin(); contigIt != eliminatedContigs.end(); ++contigIt) {
        unsigned int contigKey = *contigIt;
        std::unordered_map<unsigned int, std::pair<unsigned int, size_t>>::iterator it = contigToLongest.find(contigKey);
        if (it != contigToLongest.end()) {
            unsigned int orfKey = it->second.first;
            resultBuffer.append(SSTR(orfKey));
            resultBuffer.append(1, '\n');
            writer.writeData(resultBuffer.c_str(), resultBuffer.length(), orfKey, 0, false, false);
            resultBuffer.clear();
        } else {
            Debug(Debug::ERROR) << "Missing contig " << contigKey << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

    writer.close(true);
    FileUtil::remove(writer.getIndexFileName());
    headerReader.close();
    resultReader.close();

    return EXIT_SUCCESS;
}
