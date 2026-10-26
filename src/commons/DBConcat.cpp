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

    int mode = DBReader<DBKeyType>::USE_INDEX;
    if (write == true) {
        mode |= DBReader<DBKeyType>::USE_DATA;
    }
    if (shouldConcatLookup) {
        mode |= DBReader<DBKeyType>::USE_LOOKUP;
    }
    DBReader<DBKeyType> dbA(dataFileNameA.c_str(), indexFileNameA.c_str(), threads, mode);
    DBReader<DBKeyType> dbB(dataFileNameB.c_str(), indexFileNameB.c_str(), threads, mode);
    dbA.open(DBReader<DBKeyType>::NOSORT);
    dbB.open(DBReader<DBKeyType>::NOSORT);
    indexSizeA = dbA.getSize();
    indexSizeB = dbB.getSize();

    // keys paris are like : (key,i) where key is the ith key in the database
    keysA = new std::pair<DBKeyType, DBKeyType>[indexSizeA];
    keysB = new std::pair<DBKeyType, DBKeyType>[indexSizeB];

    DBWriter* concatWriter = NULL;
    if (write == true) {
        concatWriter = new DBWriter(dataFileNameC.c_str(), indexFileNameC.c_str(), threads, Parameters::WRITER_ASCII_MODE, dbA.getDbtype());
        concatWriter->open();
    }

    Debug::Progress progress(indexSizeA);
    // where the new key numbering of B should start
    const bool writeNull = trimRight > 0;
    DBKeyType maxKeyA = 0;
#pragma omp parallel num_threads(threads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for schedule(dynamic, 10) reduction(max:maxKeyA)
        for (size_t id = 0; id < indexSizeA; id++) {
            progress.updateProgress();

            DBKeyType newKey;
            if (preserveKeysA) {
                newKey = dbA.getDbKey(id);
            } else {
                newKey = static_cast<DBKeyType>(id);
            }

            if (write) {
                char *data = dbA.getData(id, thread_idx);
                size_t dataSizeA = std::max(dbA.getEntryLen(id), trimRight) - trimRight;
                if (takeLargerEntry == true) {
                    size_t idB = dbB.getId(newKey);
                    size_t dataSizeB = (idB == DB_ENTRY_NOT_FOUND) ? 0 : std::max(dbB.getEntryLen(idB), trimRight) - trimRight;
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

            DBKeyType newKey;
            if (preserveKeysB) {
                newKey = dbB.getDbKey(id);
            } else {
                newKey = static_cast<DBKeyType>(id) + maxKeyA;
            }

            if (write) {
                char *data = dbB.getData(id, thread_idx);
                size_t dataSizeB = std::max(dbB.getEntryLen(id), trimRight) - trimRight;
                if (takeLargerEntry) {
                    size_t idB = dbA.getId(newKey);
                    size_t dataSizeA = (idB == DB_ENTRY_NOT_FOUND) ? 0 : std::max(dbA.getEntryLen(idB), trimRight) - trimRight;
                    if (dataSizeB > dataSizeA) {
                        concatWriter->writeData(data, dataSizeB, newKey, thread_idx, writeNull);
                    }
                } else {
                    concatWriter->writeData(data, dataSizeB, newKey, thread_idx, writeNull);
                }
            }

            // need to store the index, because it'll be sorted out by keys later
            keysB[id] = std::make_pair(dbB.getDbKey(id), newKey);
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
        std::vector<std::pair<DBKeyType, DBKeyType>> mappingA;
        Util::readMappingDBKey((dataFileNameA + "_mapping"), mappingA);
        std::vector<std::pair<DBKeyType, DBKeyType>> mappingB;
        Util::readMappingDBKey((dataFileNameB + "_mapping"), mappingB);

        FILE* mappingFilePtr = fopen((dataFileNameC + "_mapping").c_str(), "w");

        for(size_t i = 0; i < mappingA.size(); ++i) {
            DBKeyType prevKeyA = mappingA[i].first;
            DBKeyType taxidA = mappingA[i].second;
            DBKeyType newKeyA = dbAKeyMap(prevKeyA);

            char * basePos = buffer;
            char * tmpBuff = Itoa::u64toa_sse2(static_cast<uint64_t>(newKeyA), buffer);
            *(tmpBuff-1) = '\t';
            tmpBuff = Itoa::u64toa_sse2(static_cast<uint64_t>(taxidA), tmpBuff);;
            *(tmpBuff-1) = '\n';
            size_t length = tmpBuff - basePos;

            int written = fwrite(buffer, sizeof(char), length, mappingFilePtr);
            if (written != (int) length) {
                Debug(Debug::ERROR) << "Cannot write to data file " << dataFileNameC << "_mapping\n";
                EXIT(EXIT_FAILURE);
            }
        }

        for(size_t i = 0; i < mappingB.size(); ++i) {
            DBKeyType prevKeyB = mappingB[i].first;
            DBKeyType taxidB = mappingB[i].second;
            DBKeyType newKeyB = dbBKeyMap(prevKeyB);

            char * basePos = buffer;
            char * tmpBuff = Itoa::u64toa_sse2(static_cast<uint64_t>(newKeyB), buffer);
            *(tmpBuff-1) = '\t';
            tmpBuff = Itoa::u64toa_sse2(static_cast<uint64_t>(taxidB), tmpBuff);;
            *(tmpBuff-1) = '\n';
            size_t length = tmpBuff - basePos;

            int written = fwrite(buffer, sizeof(char), length, mappingFilePtr);
            if (written != (int) length) {
                Debug(Debug::ERROR) << "Cannot write to data file " << dataFileNameC << "_mapping\n";
                EXIT(EXIT_FAILURE);
            }
        }
        if (fclose(mappingFilePtr) != 0) {
            Debug(Debug::ERROR) << "Cannot close data file " << dataFileNameC << "_mapping\n";
            EXIT(EXIT_FAILURE);
        }
    }

    DBKeyType maxSetIdA = 0;
    // handle lookup
    if (shouldConcatLookup) {
        DBReader<DBKeyType> lookupReaderA(dataFileNameA.c_str(), indexFileNameA.c_str(), 1, DBReader<DBKeyType>::USE_LOOKUP);
        lookupReaderA.open(DBReader<DBKeyType>::NOSORT);
        DBReader<DBKeyType>::LookupEntry* lookupA = lookupReaderA.getLookup();

        FILE* lookupFilePtr = fopen((dataFileNameC + ".lookup").c_str(), "w");

        char buffer[1024];
        std::string line;

        for (size_t i = 0; i < lookupReaderA.getLookupSize(); ++i) {
            DBKeyType prevKeyA = lookupA[i].id;
            std::string accA = lookupA[i].entryName;
            DBKeyType setIdA = lookupA[i].fileNumber;
            if (setIdA > maxSetIdA) {
                maxSetIdA = setIdA;
            }

            DBKeyType newKeyA = dbAKeyMap(prevKeyA);

            char *tmpBuff = Itoa::u64toa_sse2(static_cast<uint64_t>(newKeyA), buffer);
            line.append(buffer, tmpBuff - buffer - 1);
            line.append(1, '\t');
            line.append(accA);
            line.append(1, '\t');
            tmpBuff = Itoa::u64toa_sse2(static_cast<uint64_t>(setIdA), buffer);
            line.append(buffer, tmpBuff - buffer - 1);
            line.append(1, '\n');
            
            int written = fwrite(line.c_str(), sizeof(char), line.size(), lookupFilePtr);
            if (written != (int) line.size()) {
                Debug(Debug::ERROR) << "Cannot write to data file " << dataFileNameC << ".lookup\n";
                EXIT(EXIT_FAILURE);
            }
            line.clear();
        }
        lookupReaderA.close();

        // for B we compute: newSetIdB = maxSetIdA + 1 + setIdB
        DBReader<DBKeyType> lookupReaderB(dataFileNameB.c_str(), indexFileNameB.c_str(), 1, DBReader<DBKeyType>::USE_LOOKUP);
        lookupReaderB.open(DBReader<DBKeyType>::NOSORT);
        DBReader<DBKeyType>::LookupEntry* lookupB = lookupReaderB.getLookup();
        for (size_t i = 0; i < lookupReaderB.getLookupSize(); ++i) {
            DBKeyType prevKeyB = lookupB[i].id;
            std::string accB = lookupB[i].entryName;
            DBKeyType setIdB = lookupB[i].fileNumber;
            
            DBKeyType newKeyB = dbBKeyMap(prevKeyB);
            DBKeyType newSetIdB = maxSetIdA + 1 + setIdB;

            char *tmpBuff = Itoa::u64toa_sse2(static_cast<uint64_t>(newKeyB), buffer);
            line.append(buffer, tmpBuff - buffer - 1);
            line.append(1, '\t');
            line.append(accB);
            line.append(1, '\t');
            tmpBuff = Itoa::u64toa_sse2(static_cast<uint64_t>(newSetIdB), buffer);
            line.append(buffer, tmpBuff - buffer - 1);
            line.append(1, '\n');
            
            int written = fwrite(line.c_str(), sizeof(char), line.size(), lookupFilePtr);
            if (written != (int) line.size()) {
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
        DBKeyType sourceMaxSetIdA = 0;
        std::map<DBKeyType, std::string> sourceMapA = Util::readLookup((dataFileNameA + ".source"), false);
        std::map<DBKeyType, std::string>::iterator itA;
        
        char buffer[1024];
        std::string line;

        FILE* sourceFilePtr = fopen((dataFileNameC + ".source").c_str(), "w");
        for (itA = sourceMapA.begin(); itA != sourceMapA.end(); itA++) {
            DBKeyType setIdA = itA->first;
            std::string fileNameA = itA->second;
            if (setIdA > sourceMaxSetIdA) {
                sourceMaxSetIdA = setIdA;
            }

            char *tmpBuff = Itoa::u64toa_sse2(static_cast<uint64_t>(setIdA), buffer);
            line.append(buffer, tmpBuff - buffer - 1);
            line.append(1, '\t');
            line.append(fileNameA);
            line.append(1, '\n');

            int written = fwrite(line.c_str(), sizeof(char), line.size(), sourceFilePtr);
            if (written != (int) line.size()) {
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

        std::map<DBKeyType, std::string> sourceMapB = Util::readLookup((dataFileNameB + ".source"), false);
        std::map<DBKeyType, std::string>::iterator itB;

        for (itB = sourceMapB.begin(); itB != sourceMapB.end(); itB++) {
            DBKeyType setIdB = itB->first;
            std::string fileNameB = itB->second;

            DBKeyType newSetIdB = sourceMaxSetIdA + 1 + setIdB;

            char *tmpBuff = Itoa::u64toa_sse2(static_cast<uint64_t>(newSetIdB), buffer);
            line.append(buffer, tmpBuff - buffer - 1);
            line.append(1, '\t');
            line.append(fileNameB);
            line.append(1, '\n');

            int written = fwrite(line.c_str(), sizeof(char), line.size(), sourceFilePtr);
            if (written != (int) line.size()) {
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

DBKeyType DBConcat::dbAKeyMap(DBKeyType key) {
    if (sameDatabase)
        return key;

    std::pair<DBKeyType, DBKeyType> *originalMap = std::upper_bound(keysA, keysA + indexSizeA, key, compareKeyToFirstEntry());
    return originalMap->second;
}

DBKeyType DBConcat::dbBKeyMap(DBKeyType key) {
    if (sameDatabase)
        return key;

    std::pair<DBKeyType, DBKeyType> *originalMap = std::upper_bound(keysB, keysB + indexSizeB, key, compareKeyToFirstEntry());
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
