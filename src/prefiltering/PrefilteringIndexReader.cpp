#include "PrefilteringIndexReader.h"
#include "../commons/Debug.h"
#include "../commons/DBWriter.h"
#include "Prefiltering.h"

const char * PrefilteringIndexReader::VERSION       = "MMSEQSVERSION";
const char * PrefilteringIndexReader::ENTRIES       = "MMSEQSENTRIES";
const char * PrefilteringIndexReader::ENTRIESIZES   = "MMSEQSENTRIESIZES";
const char * PrefilteringIndexReader::META          = "MMSEQSMETA";
    
bool PrefilteringIndexReader::checkIfIndexFile(DBReader * reader){
    return (reader->getDataByDBKey((char*)VERSION) == NULL) ? false : true;
}
    
void PrefilteringIndexReader::createIndexFile(std::string outDB, DBReader * dbr, Sequence * seq,
                                int split, int alphabetSize, int kmerSize, int skip )
{
    int dbFrom = 0;
    size_t dbSize = dbr->getSize();
        
    IndexTable * indexTable = Prefiltering::generateIndexTable(dbr, seq, alphabetSize, kmerSize, dbFrom, dbFrom + dbSize , skip);
    DBWriter writer(outDB.c_str(), std::string( outDB +".index").c_str(), DBWriter::BINARY_MODE);
    writer.open();
        
    Debug(Debug::WARNING) << "Write  " << ENTRIES << "\n";
    char * entries = (char *) indexTable->getEntries();
    writer.write(entries, indexTable->tableEntriesNum * sizeof(int),(char*) ENTRIES, 0);
        
        
    Debug(Debug::WARNING) << "Write " << ENTRIESIZES << "\n";
    char * sizes = (char *) indexTable->getSizes();
    writer.write(sizes, indexTable->tableSize * sizeof(int),(char*) ENTRIESIZES, 0);
    
    Debug(Debug::WARNING) << "Write " << META << "\n";
    int64_t metadata[] = {kmerSize, alphabetSize, skip, indexTable->tableEntriesNum};
    char * metadataptr = (char *) &metadata;
    writer.write(metadataptr, 4 * sizeof(int64_t),(char *) META, 0);
    
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
    delete indexTable;
    Debug(Debug::WARNING) << "Done. \n";
}

DBReader * PrefilteringIndexReader::openNewReader(DBReader * dbr){
    std::string filePath(dbr->getDataFileName());
    std::string fullPath = filePath + std::string(".mmseqsindex");
    DBReader * reader = new DBReader("", fullPath.c_str(), DBReader::INDEXONLY );
    reader->open(DBReader::SORT);
    return reader;
}

IndexTable * PrefilteringIndexReader::generateIndexTable(DBReader * dbr){
    int * entrie_sizes   = (int *) dbr->getDataByDBKey((char*)ENTRIESIZES);
    int * entries        = (int *) dbr->getDataByDBKey((char*)ENTRIES);
    PrefilteringIndexData data = getMetadata(dbr);
    IndexTable * retTable = new IndexTable(data.alphabetSize, data.kmerSize, data.skip);
    retTable->initTableByExternalData(data.entrieSize, entrie_sizes, entries);
    return retTable;
}


PrefilteringIndexData PrefilteringIndexReader::getMetadata(DBReader * dbr){
    PrefilteringIndexData prefData;
    int * version_tmp = (int *) dbr->getDataByDBKey((char*)VERSION);
    Debug(Debug::WARNING) << "Index version: " << version_tmp[0] << "\n";
    int64_t * metadata_tmp = (int64_t *) dbr->getDataByDBKey((char *)META);
    Debug(Debug::WARNING) << "KmerSize:     " << metadata_tmp[0] << "\n";
    Debug(Debug::WARNING) << "AlphabetSize: " << metadata_tmp[1] << "\n";
    Debug(Debug::WARNING) << "Skip:         " << metadata_tmp[2] << "\n";
    prefData.kmerSize         = (int) metadata_tmp[0];
    prefData.alphabetSize     = (int) metadata_tmp[1];
    prefData.skip             = (int) metadata_tmp[2];
    prefData.entrieSize       = metadata_tmp[3];
    return prefData;
}
