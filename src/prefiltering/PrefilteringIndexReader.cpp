#include "PrefilteringIndexReader.h"
#include "DBWriter.h"
#include "Prefiltering.h"
#include "IndexTableGlobal.h"
#include "IndexTableLocal.h"

const char *PrefilteringIndexReader::CURRENT_VERSION="2.0.0";
const char *PrefilteringIndexReader::VERSION = "0";
const char *PrefilteringIndexReader::ENTRIES = "1";
const char *PrefilteringIndexReader::ENTRIESIZES = "2";
const char *PrefilteringIndexReader::ENTRIESNUM = "3";
const char *PrefilteringIndexReader::TABLESIZE = "4";
const char *PrefilteringIndexReader::META = "5";

bool PrefilteringIndexReader::checkIfIndexFile(DBReader *reader) {
    char * version = reader->getDataByDBKey((char *) VERSION);
    if(version == NULL){
        return false;
    }
    return (strncmp(version, CURRENT_VERSION, strlen(CURRENT_VERSION)) == 0 ) ? true : false;
}

void PrefilteringIndexReader::createIndexFile(std::string outDB, std::string outDBIndex, DBReader *dbr, Sequence *seq,
                                              int split, int alphabetSize, int kmerSize, int skip, bool hasSpacedKmer,
                                              int searchMode) {
    DBWriter writer(outDB.c_str(), outDBIndex.c_str(), DBWriter::BINARY_MODE);
    writer.open();
    int stepCnt = split;

    for (int step = 0; step < stepCnt; step++) {
        size_t splitStart = 0;
        size_t splitSize  = 0;
        Util::decomposeDomainByAminoaAcid(dbr->getAminoAcidDBSize(), dbr->getSeqLens(), dbr->getSize(),
                step, stepCnt, &splitStart, &splitSize);
        IndexTable *indexTable;
        if (searchMode == Parameters::SEARCH_LOCAL || searchMode == Parameters::SEARCH_LOCAL_FAST) {
            indexTable = new IndexTableLocal(alphabetSize, kmerSize, skip);
        }else{
            indexTable = new IndexTableGlobal(alphabetSize, kmerSize, skip);
        }
        Prefiltering::countKmersForIndexTable(dbr, seq, indexTable, splitStart, splitStart + splitSize);
        Debug(Debug::INFO) << "\n";
        Prefiltering::maskLowComplexKmer(indexTable, kmerSize, alphabetSize, seq->int2aa);
        Debug(Debug::INFO) << "\n";
        Prefiltering::fillDatabase(dbr, seq, indexTable, splitStart, splitStart + splitSize);

        // save the entries
        std::string entries_key = std::string(ENTRIES) + SSTR(step);
        Debug(Debug::WARNING) << "Write " << entries_key << "\n";
        char *entries = (char *) indexTable->getEntries();
        writer.write(entries,
                indexTable->getTableEntriesNum() * indexTable->getSizeOfEntry(),
                (char *) entries_key.c_str(), 0);
        indexTable->deleteEntries();

        // save the size
        std::string entriesizes_key = std::string(ENTRIESIZES) + SSTR(step);
        Debug(Debug::WARNING) << "Write " << entriesizes_key << "\n";
        char **sizes = indexTable->getTable();
        size_t *size = new size_t[indexTable->getTableSize()];
        for (size_t i = 0; i < indexTable->getTableSize(); i++) {
            const ptrdiff_t diff =  (sizes[i + 1] - sizes[i]) / indexTable->getSizeOfEntry();
            size[i] = (size_t) diff;
        }
        writer.write((char *) size, indexTable->getTableSize() * sizeof(size_t), (char *) entriesizes_key.c_str(), 0);
        delete[] size;
        // meta data
        // ENTRIESNUM
        std::string entriesnum_key = std::string(ENTRIESNUM) + SSTR(step);
        Debug(Debug::WARNING) << "Write " << entriesnum_key << "\n";
        int64_t entriesNum = indexTable->getTableEntriesNum();
        char *entriesNumPtr = (char *) &entriesNum;
        writer.write(entriesNumPtr, 1 * sizeof(int64_t), (char *) entriesnum_key.c_str(), 0);
        // TABLESIZE
        std::string tablesize_key = std::string(TABLESIZE) + SSTR(step);
        Debug(Debug::WARNING) << "Write " << tablesize_key << "\n";
        size_t tablesize = {indexTable->getSize()};
        char *tablesizePtr = (char *) &tablesize;
        writer.write(tablesizePtr, 1 * sizeof(size_t), (char *) tablesize_key.c_str(), 0);

        delete indexTable;
    }
    Debug(Debug::WARNING) << "Write " << META << "\n";
    int local = (searchMode) ? 1 : 0;
    int spacedKmer = (hasSpacedKmer) ? 1 : 0;
    int metadata[] = {kmerSize, alphabetSize, skip, split, local, spacedKmer};
    char *metadataptr = (char *) &metadata;
    writer.write(metadataptr, 6 * sizeof(int), (char *) META, 0);

    Debug(Debug::WARNING) << "Write " << VERSION << "\n";
    writer.write((char*)CURRENT_VERSION, strlen(CURRENT_VERSION) * sizeof(char), (char *) VERSION, 0);

    Debug(Debug::WARNING) << "Write MMSEQSFFINDEX \n";
    std::ifstream src(dbr->getIndexFileName(), std::ios::binary);
    std::ofstream dst(std::string(outDB + ".mmseqsindex").c_str(), std::ios::binary);
    dst << src.rdbuf();
    src.close();
    dst.close();
    writer.close();
    Debug(Debug::WARNING) << "Done. \n";
}

DBReader *PrefilteringIndexReader::openNewReader(DBReader *dbr) {
    std::string filePath(dbr->getDataFileName());
    std::string fullPath = filePath + std::string(".mmseqsindex");
    DBReader *reader = new DBReader("", fullPath.c_str(), DBReader::INDEXONLY);
    reader->open(DBReader::SORT);
    return reader;
}

IndexTable *PrefilteringIndexReader::generateIndexTable(DBReader *dbr, int split) {
    std::string entriesizes_key = std::string(ENTRIESIZES) + SSTR(split);
    size_t *entrieSizes         = (size_t *) dbr->getDataByDBKey((char *) entriesizes_key.c_str());
    std::string entries_key     = std::string(ENTRIES) + SSTR(split);
    char *entries               = (char *) dbr->getDataByDBKey((char *) entries_key.c_str());
    std::string entriesnum_key  = std::string(ENTRIESNUM) + SSTR(split);
    int64_t *entriesNum         = (int64_t *) dbr->getDataByDBKey((char *) entriesnum_key.c_str());
    std::string tablesize_key   = std::string(TABLESIZE) + SSTR(split);
    size_t *tableSize           = (size_t *) dbr->getDataByDBKey((char *) tablesize_key.c_str());

    PrefilteringIndexData data = getMetadata(dbr);
    IndexTable *retTable;
    if (data.local)
        retTable = new IndexTableLocal(data.alphabetSize, data.kmerSize, data.skip);
    else
        retTable = new IndexTableGlobal(data.alphabetSize, data.kmerSize, data.skip);

    retTable->initTableByExternalData(entriesNum[0], entrieSizes, entries, tableSize[0]);
    return retTable;
}


PrefilteringIndexData PrefilteringIndexReader::getMetadata(DBReader *dbr) {
    PrefilteringIndexData prefData;
    int *version_tmp = (int *) dbr->getDataByDBKey((char *) VERSION);
    Debug(Debug::WARNING) << "Index version: " << version_tmp[0] << "\n";
    int *metadata_tmp = (int *) dbr->getDataByDBKey((char *) META);

    Debug(Debug::WARNING) << "KmerSize:     " << metadata_tmp[0] << "\n";
    Debug(Debug::WARNING) << "AlphabetSize: " << metadata_tmp[1] << "\n";
    Debug(Debug::WARNING) << "Skip:         " << metadata_tmp[2] << "\n";
    Debug(Debug::WARNING) << "Split:        " << metadata_tmp[3] << "\n";
    Debug(Debug::WARNING) << "Type:         " << metadata_tmp[4] << "\n";
    Debug(Debug::WARNING) << "Spaced:       " << metadata_tmp[5] << "\n";

    prefData.kmerSize = metadata_tmp[0];
    prefData.alphabetSize = metadata_tmp[1];
    prefData.skip = metadata_tmp[2];
    prefData.split = metadata_tmp[3];
    prefData.local = metadata_tmp[4];
    prefData.spacedKmer = metadata_tmp[5];

    return prefData;
}
