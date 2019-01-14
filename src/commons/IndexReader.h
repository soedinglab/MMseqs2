#ifndef MMSEQS_INDEXREADER_H
#define MMSEQS_INDEXREADER_H

#include "DBReader.h"
#include "Debug.h"
#include "PrefilteringIndexReader.h"

class IndexReader {
public:
    IndexReader(const std::string &dataName, int threads, int mode = NEED_SEQUENCES | NEED_HEADERS, bool preload = false)
            : sequenceReader(NULL), headerReader(NULL), index(NULL) {
        int targetDbtype = DBReader<unsigned int>::parseDbType(dataName.c_str());
        if (Parameters::isEqualDbtype(targetDbtype, Parameters::DBTYPE_INDEX_DB)) {
            Debug(Debug::INFO) << "Use index " << dataName << "\n";
            index = new DBReader<unsigned int>(dataName.c_str(), (dataName + ".index").c_str(), 1, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
            index->open(DBReader<unsigned int>::NOSORT);
            if (PrefilteringIndexReader::checkIfIndexFile(index)) {
                PrefilteringIndexReader::printSummary(index);
                PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(index);
                seqType = data.seqType;

                if (mode & NEED_SRC_SEQUENCES) {
                    sequenceReader = PrefilteringIndexReader::openNewReader(index,
                            PrefilteringIndexReader::DBR2DATA, PrefilteringIndexReader::DBR2INDEX, mode & NEED_SEQUENCES, threads, preload);

                } else {
                    sequenceReader = PrefilteringIndexReader::openNewReader(index,
                            PrefilteringIndexReader::DBR1DATA, PrefilteringIndexReader::DBR1INDEX, mode & NEED_SEQUENCES, threads, preload);
                }
                if (sequenceReader == NULL) {
                    Debug(Debug::INFO) << "Index does not contain plain sequences. Using normal database instead.\n";
                }

                if (mode & NEED_HEADERS) {
                    if (data.headers2 == 1 && mode & NEED_ALT_HEADERS) {
                        headerReader = PrefilteringIndexReader::openNewHeaderReader(index,
                               PrefilteringIndexReader::HDR1DATA, PrefilteringIndexReader::HDR1INDEX, threads, preload);
                    } else if (data.headers1 == 1) {
                        headerReader = PrefilteringIndexReader::openNewHeaderReader(index,
                               PrefilteringIndexReader::HDR2DATA, PrefilteringIndexReader::HDR2INDEX, threads, preload);
                    } else {
                        Debug(Debug::INFO) << "Index does not contain headers. Using normal database instead.\n";
                    }
                }
            } else {
                Debug(Debug::WARNING) << "Outdated index version. Please recompute with 'createindex'!\n";
                index->close();
                delete index;
                index = NULL;
            }
        }

        if (sequenceReader == NULL) {
            const int dataMode = DBReader<unsigned int>::USE_INDEX | ((mode & NEED_SEQUENCES) ? DBReader<unsigned int>::USE_DATA : 0);
            sequenceReader = new DBReader<unsigned int>(dataName.c_str(), (dataName + ".index").c_str(), threads, dataMode);
            sequenceReader->open(DBReader<unsigned int>::NOSORT);
            if (preload) {
                sequenceReader->readMmapedDataInMemory();
            }
            seqType = sequenceReader->getDbtype();
        }

        if ((mode & NEED_HEADERS) && headerReader == NULL) {
            headerReader = new DBReader<unsigned int>((dataName + "_h").c_str(), (dataName + "_h.index").c_str(), threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
            headerReader->open(DBReader<unsigned int>::NOSORT);
            if (preload) {
                headerReader->readMmapedDataInMemory();
            }
        }
    }

    static const int NEED_SEQ_INDEX = 0;
    static const int NEED_SEQUENCES = 1;
    static const int NEED_HEADERS   = 2;
    static const int NEED_ALT_HEADERS = 2 | 4;
    static const int NEED_SRC_SEQUENCES = 1 | 8;

    int getDbtype() const {
        return seqType;
    }

    ~IndexReader() {
        if (sequenceReader != NULL) {
            sequenceReader->close();
            delete sequenceReader;
        }

        if (headerReader != NULL) {
            headerReader->close();
            delete headerReader;
        }

        if (index != NULL) {
            index->close();
            delete index;
        }
    }

    DBReader<unsigned int> *sequenceReader;
    DBReader<unsigned int> *headerReader;

private:
    DBReader<unsigned int> *index;
    int seqType;
};

#endif
