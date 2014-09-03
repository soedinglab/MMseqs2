#include "PrefilteringIndexReader.h"
#include "Debug.h"
#include "DBWriter.h"
#include "Prefiltering.h"
#include "IndexTableGlobal.h"

const char * PrefilteringIndexReader::VERSION       = "MMSEQSVERSION";
const char * PrefilteringIndexReader::ENTRIES       = "MMSEQSENTRIES";
const char * PrefilteringIndexReader::ENTRIESIZES   = "MMSEQSENTRIESIZES";
const char * PrefilteringIndexReader::ENTRIESNUM    = "MMSEQSENTRIESNUM";
const char * PrefilteringIndexReader::TABLESIZE     = "MMSEQSTABLESIZE";
const char * PrefilteringIndexReader::META          = "MMSEQSMETA";
    
bool PrefilteringIndexReader::checkIfIndexFile(DBReader * reader){
    return (reader->getDataByDBKey((char*)VERSION) == NULL) ? false : true;
}
    
void PrefilteringIndexReader::createIndexFile(std::string outDB, DBReader * dbr, Sequence * seq,
                                int split, int alphabetSize, int kmerSize, int skip )
{
    DBWriter writer(outDB.c_str(), std::string( outDB +".index").c_str(), DBWriter::BINARY_MODE);
    writer.open();
    int stepCnt = split;
    for(int step = 0; step < stepCnt; step++){
        int splitStart, splitSize;
        Util::decomposeDomainByAminoaAcid(dbr->getAminoAcidDBSize(), dbr->getSeqLens(), dbr->getSize(),
                                      step, stepCnt, &splitStart, &splitSize);
        IndexTable * indexTable = new IndexTableGlobal(alphabetSize, kmerSize, skip);

        Prefiltering::countKmersForIndexTable(dbr, seq, indexTable, splitStart, splitStart + splitSize);
        
        std::string entriesizes_key = std::string(ENTRIESIZES) + SSTR(step);
        Debug(Debug::WARNING) << "Write " << entriesizes_key << "\n";
        //char * sizes = (char *) indexTable->getSizes();
        //writer.write(sizes, indexTable->getTableSize() * sizeof(short),(char*) entriesizes_key.c_str(), 0);
        
        Prefiltering::fillDatabase(dbr, seq, indexTable, splitStart, splitStart + splitSize );
        
        std::string entries_key = std::string(ENTRIES) + SSTR(step);
        Debug(Debug::WARNING) << "Write " << entries_key << "\n";
        char * entries = (char *) indexTable->getEntries();
        writer.write(entries, (indexTable->getTableEntriesNum() + indexTable->getTableSize()) * sizeof(int),(char*) entries_key.c_str(), 0);
        
        std::string entriesnum_key = std::string(ENTRIESNUM) + SSTR(step);
        Debug(Debug::WARNING) << "Write " << entriesnum_key << "\n";
        int64_t entriesNum = { indexTable->getTableEntriesNum() };
        char * entriesNumPtr = (char *) &entriesNum;
        writer.write(entriesNumPtr, 1 * sizeof(int64_t),(char *) entriesnum_key.c_str(), 0);
        
        std::string tablesize_key = std::string(TABLESIZE) + SSTR(step);
        Debug(Debug::WARNING) << "Write " << tablesize_key << "\n";
        unsigned int tablesize = { indexTable->getSize() };
        char * tablesizePtr = (char *) &tablesize;
        writer.write(tablesizePtr, 1 * sizeof(int),(char *) tablesize_key.c_str(), 0);
        
        delete indexTable;
    }
    Debug(Debug::WARNING) << "Write " << META << "\n";
    int metadata[] = {kmerSize, alphabetSize, skip, split};
    char * metadataptr = (char *) &metadata;
    writer.write(metadataptr, 4 * sizeof(int),(char *) META, 0);
    
    Debug(Debug::WARNING) << "Write " <<  VERSION << "\n";
    int version = 1;
    char * versionptr = (char *) &version;
    writer.write(versionptr, sizeof(int),(char *) VERSION, 0);
    
    
    Debug(Debug::WARNING) << "Write MMSEQSFFINDEX \n";
    std::ifstream  src(dbr->getIndexFileName(),         std::ios::binary);
    std::ofstream  dst(std::string(outDB+".mmseqsindex").c_str(),   std::ios::binary);
    dst << src.rdbuf();
    src.close();
    dst.close();
    writer.close();
    Debug(Debug::WARNING) << "Done. \n";
}

DBReader * PrefilteringIndexReader::openNewReader(DBReader * dbr){
    std::string filePath(dbr->getDataFileName());
    std::string fullPath = filePath + std::string(".mmseqsindex");
    DBReader * reader = new DBReader("", fullPath.c_str(), DBReader::INDEXONLY );
    reader->open(DBReader::SORT);
    return reader;
}

IndexTable * PrefilteringIndexReader::generateIndexTable(DBReader * dbr, int split){
    std::string entriesizes_key = std::string(ENTRIESIZES) + SSTR(split);
    unsigned short * entrie_sizes   = (unsigned short  *) dbr->getDataByDBKey((char*)entriesizes_key.c_str());
    std::string entries_key = std::string(ENTRIES) + SSTR(split);
    unsigned int * entries        = (unsigned int *) dbr->getDataByDBKey((char*)entries_key.c_str());
    std::string entriesnum_key = std::string(ENTRIESNUM) + SSTR(split);
    int64_t * entriesNum   = (int64_t *) dbr->getDataByDBKey((char*)entriesnum_key.c_str());
    std::string tablesize_key = std::string(TABLESIZE) + SSTR(split);
    unsigned int * tableSize   = (unsigned int *) dbr->getDataByDBKey((char*)tablesize_key.c_str());

    PrefilteringIndexData data = getMetadata(dbr);
    IndexTable * retTable = new IndexTableGlobal(data.alphabetSize, data.kmerSize, data.skip);
    retTable->initTableByExternalData(entriesNum[0], entrie_sizes, entries, tableSize[0]);
    return retTable;
}


PrefilteringIndexData PrefilteringIndexReader::getMetadata(DBReader * dbr){
    PrefilteringIndexData prefData;
    int * version_tmp = (int *) dbr->getDataByDBKey((char*)VERSION);
    Debug(Debug::WARNING) << "Index version: " << version_tmp[0] << "\n";
    int * metadata_tmp = (int *) dbr->getDataByDBKey((char *)META);
    Debug(Debug::WARNING) << "KmerSize:     " << metadata_tmp[0] << "\n";
    Debug(Debug::WARNING) << "AlphabetSize: " << metadata_tmp[1] << "\n";
    Debug(Debug::WARNING) << "Skip:         " << metadata_tmp[2] << "\n";
    Debug(Debug::WARNING) << "Split:        " << metadata_tmp[3] << "\n";

    prefData.kmerSize         = metadata_tmp[0];
    prefData.alphabetSize     = metadata_tmp[1];
    prefData.skip             = metadata_tmp[2];
    prefData.split            = metadata_tmp[3];
    return prefData;
}
