#include "DBWriter.h"
#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "ReducedMatrix.h"
#include "DistanceCalculator.h"
#include "Orf.h"
#include "omptl/omptl_algorithm"

#ifdef OPENMP
#include <omp.h>
#endif

int clusthash(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.alphabetSize = MultiParam<int>(Parameters::CLUST_HASH_DEFAULT_ALPH_SIZE,5);
    par.seqIdThr = (float)Parameters::CLUST_HASH_DEFAULT_MIN_SEQ_ID/100.0f;
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        reader.readMmapedDataInMemory();
    }

    const bool isNuclInput = Parameters::isEqualDbtype(reader.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES);
    BaseMatrix *subMat = NULL;
    if (isNuclInput == false) {
        SubstitutionMatrix sMat(par.scoringMatrixFile.aminoacids, 2.0, -0.2);
        subMat = new ReducedMatrix(sMat.probMatrix, sMat.subMatrixPseudoCounts, sMat.aa2num, sMat.num2aa, sMat.alphabetSize, par.alphabetSize.aminoacids, 2.0);
    }

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    writer.open();
    Debug(Debug::INFO) << "Hashing sequences...\n";
    std::pair<size_t, unsigned int> *hashSeqPair = new std::pair<size_t, unsigned int>[reader.getSize() + 1];
    // needed later to check if one of array
    hashSeqPair[reader.getSize()] = std::make_pair(UINT_MAX, 0);
    Debug::Progress progress(reader.getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        if (isNuclInput) {
#pragma omp for schedule(dynamic, 100)
            for (size_t id = 0; id < reader.getSize(); ++id) {
                progress.updateProgress();
                char *data = reader.getData(id, thread_idx);
                size_t length = reader.getSeqLen(id);
                const size_t INITIAL_VALUE = 0;
                const size_t A = 31;
                size_t h1 = INITIAL_VALUE;
                size_t h2 = INITIAL_VALUE;
                for (size_t i = 0; i < length; ++i){
                    h1 = ((h1*A) + data[i]);
                    h2 = ((h2*A) + Orf::complement(data[length - i - 1]));
                }
                hashSeqPair[id] = std::make_pair(std::min(h1, h2), id);
            }
        } else {
            Sequence seq(par.maxSeqLen, reader.getDbtype(), subMat, 0, false, false);
#pragma omp for schedule(dynamic, 100)
            for (size_t id = 0; id < reader.getSize(); ++id) {
                progress.updateProgress();
                char *data = reader.getData(id, thread_idx);
                size_t length = reader.getSeqLen(id);
                size_t seqHash;
                seq.mapSequence(id, 0, data, length);
                seqHash = Util::hash(seq.numSequence, seq.L);
                hashSeqPair[id] = std::make_pair(seqHash, id);
            }
        }
    }

    // sort by hash and set up the pointer for parallel processing
    omptl::sort(hashSeqPair, hashSeqPair + reader.getSize());

    size_t uniqHashes = 1;
    size_t prevHash = hashSeqPair[0].first;
    for (size_t id = 0; id < reader.getSize(); id++) {
        if (prevHash != hashSeqPair[id].first) {
            uniqHashes++;
        }
        prevHash = hashSeqPair[id].first;
    }
    std::pair<size_t, unsigned int> **hashLookup = new std::pair<size_t, unsigned int>*[uniqHashes];
    hashLookup[0] = hashSeqPair;
    size_t currKey = 1;
    prevHash = hashSeqPair[0].first;
    for (size_t id = 0; id < reader.getSize(); ++id) {
        if (prevHash != hashSeqPair[id].first) {
            hashLookup[currKey] = (hashSeqPair + id);
            currKey++;
        }
        prevHash = hashSeqPair[id].first;
    }

    Debug(Debug::INFO) << "Found " << uniqHashes << " unique hashes\n";
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        std::vector<unsigned int> setIds;
        setIds.reserve(300);
        std::vector<bool> found;
        found.reserve(300);
        std::string result;
        result.reserve(1024);
        char buffer[64];

#pragma omp for schedule(dynamic, 10)
        for (size_t hashId = 0; hashId < uniqHashes; ++hashId) {
            progress.updateProgress();

            size_t initHash = hashLookup[hashId]->first;
            size_t pos = 0;
            while (hashLookup[hashId][pos].first == initHash) {
                setIds.push_back(hashLookup[hashId][pos].second);
                found.push_back(false);
                pos++;
            }
            for (size_t i = 0; i < setIds.size(); i++) {
                unsigned int queryKey = reader.getDbKey(setIds[i]);
                unsigned int queryLength = reader.getSeqLen(setIds[i]);
                const char *querySeq = reader.getData(setIds[i], thread_idx);
                result.append(SSTR(queryKey));
                result.append("\t255\t1.00\t0\t0\t");
                result.append(SSTR(queryLength - 1));
                result.append(1, '\t');
                result.append(SSTR(queryLength));
                result.append("\t0\t");
                result.append(SSTR(queryLength - 1));
                result.append(1, '\t');
                result.append(SSTR(queryLength));
                result.append(1, '\n');
                if (found[i] == true) {
                    goto outer;
                }

                for (size_t j = 0; j < setIds.size(); j++) {
                    if (found[j] == true) {
                        continue;
                    }
                    unsigned int targetLength = reader.getSeqLen(setIds[j]);
                    if (i != j && queryLength == targetLength) {
                        const char *targetSeq = reader.getData(setIds[j], thread_idx);
                        unsigned int distance = DistanceCalculator::computeInverseHammingDistance(querySeq, targetSeq, queryLength);
                        const float seqId = (static_cast<float>(distance)) / static_cast<float>(queryLength);
                        if (seqId >= par.seqIdThr) {
                            result.append(SSTR(reader.getDbKey(setIds[j])));
                            result.append("\t255\t");
                            Util::fastSeqIdToBuffer(seqId, buffer);
                            result.append(buffer);
                            result.append("\t0\t0\t");
                            result.append(SSTR(queryLength - 1));
                            result.append(1, '\t');
                            result.append(SSTR(queryLength));
                            result.append("\t0\t");
                            result.append(SSTR(queryLength - 1));
                            result.append(1, '\t');
                            result.append(SSTR(queryLength));
                            result.append(1, '\n');
                            found[j] = true;
                        }
                    }
                }
                outer:
                writer.writeData(result.c_str(), result.length(), queryKey, thread_idx);
                result.clear();
            }
            setIds.clear();
            found.clear();
        }
    }
    writer.close();
    reader.close();
    delete[] hashLookup;
    delete[] hashSeqPair;

    if (subMat != NULL) {
        delete subMat;
    }

    return EXIT_SUCCESS;
}
