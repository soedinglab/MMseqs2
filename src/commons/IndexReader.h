#ifndef MMSEQS_INDEXREADER_H
#define MMSEQS_INDEXREADER_H

#include "DBReader.h"
#include "Debug.h"
#include "PrefilteringIndexReader.h"
#include "FileUtil.h"

class IndexReader {
public:
    const static unsigned int PRELOAD_NO = 0;
    const static unsigned int PRELOAD_DATA = 1;
    const static unsigned int PRELOAD_INDEX = 2;

    IndexReader(
            const std::string &dataName,
            int threads,
            unsigned int databaseType = SEQUENCES | HEADERS,
            unsigned int preloadMode = false,
            int dataMode = DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
            std::string failSuffix = ""
    ) : sequenceReader(NULL), index(NULL) {
        int targetDbtype = FileUtil::parseDbType(dataName.c_str());
        if (Parameters::isEqualDbtype(targetDbtype, Parameters::DBTYPE_INDEX_DB)) {
            index = new DBReader<unsigned int>(dataName.c_str(), (dataName + ".index").c_str(), 1, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
            index->open(DBReader<unsigned int>::NOSORT);
            if (PrefilteringIndexReader::checkIfIndexFile(index)) {
                PrefilteringIndexReader::printSummary(index);
                PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(index);
                seqType = data.seqType;
                bool touchIndex = preloadMode & PRELOAD_INDEX;
                bool touchData = preloadMode & PRELOAD_DATA;
                if (databaseType & USER_SELECT) {
                    sequenceReader = PrefilteringIndexReader::openNewReader(
                            index,
                            (databaseType & ~USER_SELECT) + 1,
                            databaseType & ~USER_SELECT,
                            dataMode & DBReader<unsigned int>::USE_DATA, threads, touchIndex, touchData
                    );
                } else if (databaseType & SRC_SEQUENCES) {
                    sequenceReader = PrefilteringIndexReader::openNewReader(index,
                                                                            PrefilteringIndexReader::DBR2DATA, PrefilteringIndexReader::DBR2INDEX, dataMode & DBReader<unsigned int>::USE_DATA, threads, touchIndex, touchData);
                } else if (databaseType & SEQUENCES) {
                    sequenceReader = PrefilteringIndexReader::openNewReader(index,
                                                                            PrefilteringIndexReader::DBR1DATA, PrefilteringIndexReader::DBR1INDEX, dataMode & DBReader<unsigned int>::USE_DATA, threads, touchIndex, touchData);
                } else if (databaseType & SRC_HEADERS) {

                    sequenceReader = PrefilteringIndexReader::openNewHeaderReader(index,
                                                                                  PrefilteringIndexReader::HDR2DATA,
                                                                                  PrefilteringIndexReader::HDR2INDEX,
                                                                                  threads, touchIndex, touchData);
                } else if (databaseType & HEADERS) {
                    sequenceReader = PrefilteringIndexReader::openNewHeaderReader(index,
                                                                                  PrefilteringIndexReader::HDR1DATA,
                                                                                  PrefilteringIndexReader::HDR1INDEX,
                                                                                  threads, touchIndex, touchData);
                } else if (databaseType & ALIGNMENTS) {
                    sequenceReader = PrefilteringIndexReader::openNewReader(index,
                                                                            PrefilteringIndexReader::ALNDATA,
                                                                            PrefilteringIndexReader::ALNINDEX,
                                                                            dataMode & DBReader<unsigned int>::USE_DATA,
                                                                            threads, touchIndex, touchData);
                }
                if (sequenceReader == NULL) {
                    Debug(Debug::INFO) << "Index does not contain plain sequences. Using normal database instead.\n";
                }
                seqType = Parameters::DBTYPE_INDEX_DB;
            } else {
                Debug(Debug::WARNING) << "Outdated index version. Please recompute with 'createindex'!\n";
                index->close();
                delete index;
                index = NULL;
            }
        }

        if (sequenceReader == NULL) {
            if((databaseType & USER_SELECT) == false && failSuffix == "" ){
                if (databaseType & HEADERS) {
                    failSuffix = "_h";
                } else if (databaseType & SRC_HEADERS) {
                    failSuffix = "_seq_h";
                    if(FileUtil::fileExists((dataName + "_seq_h.dbtype").c_str())==false){
                        failSuffix = "_h";
                    }
                } else if (databaseType & SRC_SEQUENCES) {
                    failSuffix = "_seq";
                    if(FileUtil::fileExists((dataName + "_seq.dbtype").c_str())==false){
                        failSuffix = "";
                    }
                } else if (databaseType & SEQUENCES) {
                    failSuffix = "";
                } else if (databaseType & ALIGNMENTS) {
                    if(FileUtil::fileExists((dataName + "_clu.dbtype").c_str())){
                        failSuffix = "_clu";
                    }else{
                        failSuffix = "_aln";
                    }
                }
            }
            sequenceReader = new DBReader<unsigned int>(
                    (dataName + failSuffix).c_str(), (dataName + failSuffix + ".index").c_str(),
                    threads, dataMode
            );
            sequenceReader->open(DBReader<unsigned int>::NOSORT);
            bool touchData = preloadMode & PRELOAD_DATA;
            if (touchData) {
                sequenceReader->readMmapedDataInMemory();
            }
            seqType = sequenceReader->getDbtype();
        }
    }

    static const unsigned int SEQUENCES = 1;
    static const unsigned int HEADERS   = 2;
    static const unsigned int SRC_HEADERS = 4;
    static const unsigned int SRC_SEQUENCES =  8;
    static const unsigned int ALIGNMENTS = 16;
    static const unsigned int USER_SELECT = 1 << 31;

    static unsigned int makeUserDatabaseType(unsigned int baseKey) {
        return baseKey | USER_SELECT;
    }

    int getDbtype() const {
        return seqType;
    }

    ~IndexReader() {
        if (sequenceReader != NULL) {
            sequenceReader->close();
            delete sequenceReader;
        }

        if (index != NULL) {
            index->close();
            delete index;
        }
    }

    DBReader<unsigned int> *sequenceReader;
    DBReader<unsigned int> *index;

private:
    int seqType;
};

#endif
