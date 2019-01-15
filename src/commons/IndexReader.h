#ifndef MMSEQS_INDEXREADER_H
#define MMSEQS_INDEXREADER_H

#include "DBReader.h"
#include "Debug.h"
#include "PrefilteringIndexReader.h"

class IndexReader {
public:
    IndexReader(const std::string &dataName, int threads, int mode = SEQUENCES | HEADERS, bool preload = false, int dataMode=(DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA))
            : sequenceReader(NULL), index(NULL) {
        int targetDbtype = DBReader<unsigned int>::parseDbType(dataName.c_str());
        if (Parameters::isEqualDbtype(targetDbtype, Parameters::DBTYPE_INDEX_DB)) {
            Debug(Debug::INFO) << "Use index " << dataName << "\n";
            index = new DBReader<unsigned int>(dataName.c_str(), (dataName + ".index").c_str(), 1, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
            index->open(DBReader<unsigned int>::NOSORT);
            if (PrefilteringIndexReader::checkIfIndexFile(index)) {
                PrefilteringIndexReader::printSummary(index);
                PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(index);
                seqType = data.seqType;

                if (mode & SRC_SEQUENCES) {
                    sequenceReader = PrefilteringIndexReader::openNewReader(index,
                            PrefilteringIndexReader::DBR2DATA, PrefilteringIndexReader::DBR2INDEX, mode & SEQUENCES, threads, preload);

                } else if (mode & SEQUENCES) {
                    sequenceReader = PrefilteringIndexReader::openNewReader(index,
                            PrefilteringIndexReader::DBR1DATA, PrefilteringIndexReader::DBR1INDEX, mode & SEQUENCES, threads, preload);
                } else if (mode & SRC_HEADERS) {

                    sequenceReader = PrefilteringIndexReader::openNewHeaderReader(index,
                                                                                  PrefilteringIndexReader::HDR2DATA,
                                                                                  PrefilteringIndexReader::HDR2INDEX,
                                                                                  threads, preload);
                }else if(mode & HEADERS) {
                    sequenceReader = PrefilteringIndexReader::openNewHeaderReader(index,
                                                                                PrefilteringIndexReader::HDR1DATA, PrefilteringIndexReader::HDR1INDEX, threads, preload);
                }

                if (sequenceReader == NULL) {
                    Debug(Debug::INFO) << "Index does not contain plain sequences. Using normal database instead.\n";
                }
            } else {
                Debug(Debug::WARNING) << "Outdated index version. Please recompute with 'createindex'!\n";
                index->close();
                delete index;
                index = NULL;
            }
        }

        if (sequenceReader == NULL) {
            if(mode & (HEADERS | SRC_HEADERS)){
                sequenceReader = new DBReader<unsigned int>((dataName+"_h").c_str(), ((dataName+"_h") + ".index").c_str(), threads, dataMode);
            }else{
                sequenceReader = new DBReader<unsigned int>(dataName.c_str(), (dataName + ".index").c_str(), threads, dataMode);
            }
            sequenceReader->open(DBReader<unsigned int>::NOSORT);
            if (preload) {
                sequenceReader->readMmapedDataInMemory();
            }
            seqType = sequenceReader->getDbtype();
        }
    }

    static const int SEQUENCES = 1;
    static const int HEADERS   = 2;
    static const int SRC_HEADERS = 4;
    static const int SRC_SEQUENCES =  8;

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

private:
    DBReader<unsigned int> *index;
    int seqType;
};

#endif
