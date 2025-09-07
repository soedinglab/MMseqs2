#include "DBConcat.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "itoa.h"
#include "Util.h"
#include "Debug.h"
#include "FileUtil.h"
#include "Parameters.h"

#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif

DBConcat::DBConcat(const std::string &dataFileNameA, const std::string &indexFileNameA,
                   const std::string &dataFileNameB, const std::string &indexFileNameB,
                   const std::string &dataFileNameC, const std::string &indexFileNameC,
                   unsigned int threads, bool write, bool preserveKeysA, bool preserveKeysB, bool takeLargerEntry, size_t trimRight) {
    sameDatabase = dataFileNameA == dataFileNameB;

    bool shouldConcatMapping = false;
    bool shouldConcatLookup = false;
    bool shouldConcatSource = false;
    if (write == true) {
        if (FileUtil::fileExists((dataFileNameA + "_mapping").c_str()) && FileUtil::fileExists((dataFileNameB + "_mapping").c_str())) {
            shouldConcatMapping = true;
        }
        if (FileUtil::fileExists((dataFileNameA + ".lookup").c_str()) && FileUtil::fileExists((dataFileNameB + ".lookup").c_str())) {
            shouldConcatLookup = true;
        }
        if (FileUtil::fileExists((dataFileNameA + ".source").c_str()) && FileUtil::fileExists((dataFileNameB + ".source").c_str())) {
            shouldConcatSource = true;
        }
    }

    int mode = DBReader<KeyType>::USE_INDEX;
    if (write == true) {
        mode |= DBReader<KeyType>::USE_DATA;
    }
    if (shouldConcatLookup) {
        mode |= DBReader<KeyType>::USE_LOOKUP;
    }
    DBReader<KeyType> dbA(dataFileNameA.c_str(), indexFileNameA.c_str(), threads, mode);
    DBReader<KeyType> dbB(dataFileNameB.c_str(), indexFileNameB.c_str(), threads, mode);
    dbA.open(DBReader<KeyType>::NOSORT);
    dbB.open(DBReader<KeyType>::NOSORT);
    indexSizeA = dbA.getSize();
    indexSizeB = dbB.getSize();

    // keys paris are like : (key,i) where key is the ith key in the database
    keysA = new std::pair<KeyType, KeyType>[indexSizeA];
    keysB = new std::pair<KeyType, KeyType>[indexSizeB];

    DBWriter* concatWriter = NULL;
    if (write == true) {
        concatWriter = new DBWriter(dataFileNameC.c_str(), indexFileNameC.c_str(), threads, Parameters::WRITER_ASCII_MODE, dbA.getDbtype());
        concatWriter->open();
    }

    Debug::Progress progress(indexSizeA);
    // where the new key numbering of B should start
    const bool writeNull = trimRight > 0;
    KeyType maxKeyA = 0;
#pragma omp parallel num_threads(threads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for schedule(dynamic, 10) reduction(max:maxKeyA)
        for (size_t id = 0; id < indexSizeA; id++) {
            progress.updateProgress();

            KeyType newKey;
            if (preserveKeysA) {
                newKey = dbA.getDbKey(id);
            } else {
                newKey = static_cast<KeyType>(id);
            }

            if (write) {
                char *data = dbA.getData(id, thread_idx);
                size_t dataSizeA = std::max(dbA.getEntryLen(id), trimRight) - trimRight;
                if (takeLargerEntry == true) {
                    KeyType idB = dbB.getId(newKey);
                    size_t dataSizeB = std::max(dbB.getEntryLen(idB), trimRight) - trimRight;
                    if (dataSizeA >= dataSizeB) {
                        concatWriter->writeData(data, dataSizeA, newKey, thread_idx, writeNull);
                    }
                } else {
                    concatWriter->writeData(data, dataSizeA, newKey, thread_idx, writeNull);
                }
            }

            // need to store the index, because it'll be sorted out by keys later
            keysA[id] = std::make_pair(dbA.getDbKey(id), newKey);
            maxKeyA = std::max(maxKeyA, newKey);
        }
    }
    maxKeyA++;

    progress.reset(indexSizeB);
#pragma omp parallel num_threads(threads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for schedule(dynamic, 10)
        for (size_t id = 0; id < indexSizeB; id++) {
            progress.updateProgress();

            KeyType newKey;
            if (preserveKeysB) {
                newKey = dbB.getDbKey(id);
            } else {
                newKey = static_cast<KeyType>(id) + maxKeyA;
            }

            if (write) {
                char *data = dbB.getData(id, thread_idx);
                size_t dataSizeB = std::max(dbB.getEntryLen(id), trimRight) - trimRight;
                if (takeLargerEntry) {
                    KeyType idB = dbA.getId(newKey);
                    size_t dataSizeA = std::max(dbA.getEntryLen(idB), trimRight) - trimRight;
                    if (dataSizeB > dataSizeA) {
                        concatWriter->writeData(data, dataSizeB, newKey, thread_idx, writeNull);
                    }
                } else {
                    concatWriter->writeData(data, dataSizeB, newKey, thread_idx, writeNull);
                }
            }

            // need to store the index, because it'll be sorted out by keys later
            keysB[id] = std::make_pair(dbB.getDbKey(id), id + maxKeyA);
        }
    }

    //sort by key
    std::stable_sort(keysA, keysA + indexSizeA, compareFirstEntry());
    std::stable_sort(keysB, keysB + indexSizeB, compareFirstEntry());

    if (write) {
        concatWriter->close(true);
        delete concatWriter;
    }
    dbA.close();
    dbB.close();

    // handle mapping
    if (shouldConcatMapping) {
        char buffer[1024];
        std::vector<std::pair<KeyType, unsigned int>> mappingA;
        Util::readMapping((dataFileNameA + "_mapping"), mappingA);
        std::vector<std::pair<KeyType, unsigned int>> mappingB;
        Util::readMapping((dataFileNameB + "_mapping"), mappingB);

        FILE* mappingFilePtr = fopen((dataFileNameC + "_mapping").c_str(), "w");

        for(size_t i = 0; i < mappingA.size(); ++i) {
            KeyType prevKeyA = mappingA[i].first;
            unsigned int taxidA = mappingA[i].second;
            KeyType newKeyA = dbAKeyMap(prevKeyA);

            char * basePos = buffer;
            constexpr bool keyIsU32 = std::is_same<KeyType, unsigned int>::value;
            char * tmpBuff = keyIsU32
                             ? Itoa::u32toa_sse2(newKeyA, buffer)
                             : Itoa::u64toa_sse2(newKeyA, buffer);

            *(tmpBuff-1) = '\t';
            tmpBuff = Itoa::u32toa_sse2(taxidA, tmpBuff);;
            *(tmpBuff-1) = '\n';
            size_t length = tmpBuff - basePos;

            size_t written = fwrite(buffer, sizeof(char), length, mappingFilePtr);
            if (written != length) {
                Debug(Debug::ERROR) << "Cannot write to data file " << dataFileNameC << "_mapping\n";
                EXIT(EXIT_FAILURE);
            }
        }

        for(size_t i = 0; i < mappingB.size(); ++i) {
            KeyType prevKeyB = mappingB[i].first;
            unsigned int taxidB = mappingB[i].second;
            KeyType newKeyB = dbBKeyMap(prevKeyB);

            char * basePos = buffer;
            constexpr bool keyIsU32 = std::is_same<KeyType, unsigned int>::value;
            char * tmpBuff = keyIsU32
                             ? Itoa::u32toa_sse2(static_cast<uint64_t>(newKeyB), buffer)
                             : Itoa::u64toa_sse2(static_cast<uint64_t>(newKeyB), buffer);
            *(tmpBuff-1) = '\t';
            tmpBuff = Itoa::u32toa_sse2(taxidB, tmpBuff);
            *(tmpBuff-1) = '\n';
            size_t length = tmpBuff - basePos;

            size_t written = fwrite(buffer, sizeof(char), length, mappingFilePtr);
            if (written != length) {
                Debug(Debug::ERROR) << "Cannot write to data file " << dataFileNameC << "_mapping\n";
                EXIT(EXIT_FAILURE);
            }
        }
        if (fclose(mappingFilePtr) != 0) {
            Debug(Debug::ERROR) << "Cannot close data file " << dataFileNameC << "_mapping\n";
            EXIT(EXIT_FAILURE);
        }
    }

    KeyType maxSetIdA = 0;
    // handle lookup
    if (shouldConcatLookup) {
        DBReader<KeyType> lookupReaderA(dataFileNameA.c_str(), indexFileNameA.c_str(), 1, DBReader<KeyType>::USE_LOOKUP);
        lookupReaderA.open(DBReader<KeyType>::NOSORT);
        DBReader<KeyType>::LookupEntry* lookupA = lookupReaderA.getLookup();

        FILE* lookupFilePtr = fopen((dataFileNameC + ".lookup").c_str(), "w");

        char buffer[1024];
        std::string line;

        for (size_t i = 0; i < lookupReaderA.getLookupSize(); ++i) {
            KeyType prevKeyA = lookupA[i].id;
            std::string accA = lookupA[i].entryName;
            KeyType setIdA = lookupA[i].fileNumber;
            if (setIdA > maxSetIdA) {
                maxSetIdA = setIdA;
            }

            KeyType newKeyA = dbAKeyMap(prevKeyA);
            constexpr bool keyIsU32 = std::is_same<KeyType, unsigned int>::value;
            char * tmpBuff = keyIsU32
                             ? Itoa::u32toa_sse2(newKeyA, buffer)
                             : Itoa::u64toa_sse2(newKeyA, buffer);

            line.append(buffer, tmpBuff - buffer - 1);
            line.append(1, '\t');
            line.append(accA);
            line.append(1, '\t');
            tmpBuff = keyIsU32
                      ? Itoa::u32toa_sse2(setIdA, buffer)
                      : Itoa::u64toa_sse2(setIdA, buffer);
            line.append(buffer, tmpBuff - buffer - 1);
            line.append(1, '\n');

            size_t written = fwrite(line.c_str(), sizeof(char), line.size(), lookupFilePtr);
            if (written != line.size()) {
                Debug(Debug::ERROR) << "Cannot write to data file " << dataFileNameC << ".lookup\n";
                EXIT(EXIT_FAILURE);
            }
            line.clear();
        }
        lookupReaderA.close();

        // for B we compute: newSetIdB = maxSetIdA + 1 + setIdB
        DBReader<KeyType> lookupReaderB(dataFileNameB.c_str(), indexFileNameB.c_str(), 1, DBReader<KeyType>::USE_LOOKUP);
        lookupReaderB.open(DBReader<KeyType>::NOSORT);
        DBReader<KeyType>::LookupEntry* lookupB = lookupReaderB.getLookup();
        for (size_t i = 0; i < lookupReaderB.getLookupSize(); ++i) {
            KeyType prevKeyB = lookupB[i].id;
            std::string accB = lookupB[i].entryName;
            KeyType setIdB = lookupB[i].fileNumber;

            KeyType newKeyB = dbBKeyMap(prevKeyB);
            KeyType newSetIdB = maxSetIdA + 1 + setIdB;

            constexpr bool keyIsU32 = std::is_same<KeyType, unsigned int>::value;
            char * tmpBuff = keyIsU32
                             ? Itoa::u32toa_sse2(newKeyB, buffer)
                             : Itoa::u64toa_sse2(newKeyB, buffer);

            line.append(buffer, tmpBuff - buffer - 1);
            line.append(1, '\t');
            line.append(accB);
            line.append(1, '\t');
            tmpBuff =  keyIsU32
                       ? Itoa::u32toa_sse2(newSetIdB, buffer)
                       : Itoa::u64toa_sse2(newSetIdB, buffer);
            line.append(buffer, tmpBuff - buffer - 1);
            line.append(1, '\n');

            size_t written = fwrite(line.c_str(), sizeof(char), line.size(), lookupFilePtr);
            if (written != line.size()) {
                Debug(Debug::ERROR) << "Cannot write to data file " << dataFileNameC << ".lookup\n";
                EXIT(EXIT_FAILURE);
            }

            line.clear();
        }
        if (fclose(lookupFilePtr) != 0) {
            Debug(Debug::ERROR) << "Cannot close file " << dataFileNameC << ".lookup\n";
            EXIT(EXIT_FAILURE);
        }
        lookupReaderB.close();
    }

    // handle source
    if (shouldConcatSource) {
        KeyType sourceMaxSetIdA = 0;
        std::map<KeyType, std::string> sourceMapA = Util::readLookup((dataFileNameA + ".source"), false);
        std::map<KeyType, std::string>::iterator itA;

        char buffer[1024];
        std::string line;

        FILE* sourceFilePtr = fopen((dataFileNameC + ".source").c_str(), "w");
        for (itA = sourceMapA.begin(); itA != sourceMapA.end(); itA++) {
            KeyType setIdA = itA->first;
            std::string fileNameA = itA->second;
            if (setIdA > sourceMaxSetIdA) {
                sourceMaxSetIdA = setIdA;
            }

            constexpr bool keyIsU32 = std::is_same<KeyType, unsigned int>::value;
            char * tmpBuff = keyIsU32
                             ? Itoa::u32toa_sse2(setIdA, buffer)
                             : Itoa::u64toa_sse2(setIdA, buffer);
            line.append(buffer, tmpBuff - buffer - 1);
            line.append(1, '\t');
            line.append(fileNameA);
            line.append(1, '\n');

            size_t written = fwrite(line.c_str(), sizeof(char), line.size(), sourceFilePtr);
            if (written != line.size()) {
                Debug(Debug::ERROR) << "Cannot write to data file " << dataFileNameC << ".source\n";
                EXIT(EXIT_FAILURE);
            }
            line.clear();
        }

        // if lookup was concatenated - make sure maxSetId there is consistent with sourceMaxSetIdA
        if (shouldConcatLookup && (sourceMaxSetIdA != maxSetIdA)) {
            Debug(Debug::ERROR) << "The maxSetId in " << dataFileNameA << ".lookup is " << maxSetIdA << " and in " << dataFileNameA << ".source is " << sourceMaxSetIdA << "\n";
            EXIT(EXIT_FAILURE);
        }

        std::map<KeyType, std::string> sourceMapB = Util::readLookup((dataFileNameB + ".source"), false);
        std::map<KeyType, std::string>::iterator itB;

        for (itB = sourceMapB.begin(); itB != sourceMapB.end(); itB++) {
            KeyType setIdB = itB->first;
            std::string fileNameB = itB->second;

            KeyType newSetIdB = sourceMaxSetIdA + 1 + setIdB;

            constexpr bool keyIsU32 = std::is_same<KeyType, unsigned int>::value;
            char * tmpBuff = keyIsU32
                             ? Itoa::u32toa_sse2(newSetIdB, buffer)
                             : Itoa::u64toa_sse2(newSetIdB, buffer);

            line.append(buffer, tmpBuff - buffer - 1);
            line.append(1, '\t');
            line.append(fileNameB);
            line.append(1, '\n');

            size_t written = fwrite(line.c_str(), sizeof(char), line.size(), sourceFilePtr);
            if (written != line.size()) {
                Debug(Debug::ERROR) << "Cannot write to data file " << dataFileNameC << ".source\n";
                EXIT(EXIT_FAILURE);
            }
            line.clear();
        }
        if (fclose(sourceFilePtr) != 0) {
            Debug(Debug::ERROR) << "Cannot close file " << dataFileNameC << ".source\n";
            EXIT(EXIT_FAILURE);
        }
    }
}

KeyType DBConcat::dbAKeyMap(KeyType key) {
    if (sameDatabase)
        return key;

    std::pair<KeyType, KeyType> *originalMap = std::upper_bound(keysA, keysA + indexSizeA, key, compareKeyToFirstEntry());
    return originalMap->second;
}

KeyType DBConcat::dbBKeyMap(KeyType key) {
    if (sameDatabase)
        return key;

    std::pair<KeyType, KeyType> *originalMap = std::upper_bound(keysB, keysB + indexSizeB, key, compareKeyToFirstEntry());
    return originalMap->second;
}

DBConcat::~DBConcat() {
    if (sameDatabase) {
        return;
    }
    delete[] keysA;
    delete[] keysB;
}

void setDbConcatDefault(Parameters *par) {
    par->threads = 1;
}

int concatdbs(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    setDbConcatDefault(&par);
    par.parseParameters(argc, argv, command, true, 0, 0);

    // TODO check equal db type
    DBConcat outDB(par.db1.c_str(), par.db1Index.c_str(),
                   par.db2.c_str(), par.db2Index.c_str(),
                   par.db3.c_str(), par.db3Index.c_str(),
                   static_cast<unsigned int>(par.threads), true, true, par.preserveKeysB, par.takeLargerEntry);

    return EXIT_SUCCESS;
}
