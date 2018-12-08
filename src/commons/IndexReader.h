#ifndef MMSEQS_INDEXREADER_H
#define MMSEQS_INDEXREADER_H

#include "DBReader.h"
#include "Debug.h"
#include "PrefilteringIndexReader.h"

class IndexReader {
public:
    IndexReader(const std::string &dataName, int threads, int mode = NEED_SEQUENCES | NEED_HEADERS,
                    bool preload = false)
            : sequenceReader(NULL), headerReader(NULL), index(NULL) {
        std::string indexDB = PrefilteringIndexReader::searchForIndex(dataName.c_str());
        if (indexDB != "") {
            Debug(Debug::INFO) << "Use index " << indexDB << "\n";
            index = new DBReader<unsigned int>(indexDB.c_str(), (indexDB + ".index").c_str(), 1, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
            index->open(DBReader<unsigned int>::NOSORT);
            if (PrefilteringIndexReader::checkIfIndexFile(index)) {
                PrefilteringIndexReader::printSummary(index);
                PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(index);
                seqType = data.seqType;

                sequenceReader = PrefilteringIndexReader::openNewReader(index, mode & NEED_SEQUENCES, preload);
                if (sequenceReader == NULL) {
                    Debug(Debug::INFO) << "Index does not contain plain sequences. Using normal database instead.\n";
                }

                if (mode & NEED_HEADERS) {
                    if (data.headers2 == 1 && mode & NEED_ALT_HEADERS) {
                        headerReader = PrefilteringIndexReader::openNewHeaderReader(index,
                                PrefilteringIndexReader::HDR1INDEX, PrefilteringIndexReader::HDR1DATA, preload);
                    } else if (data.headers1 == 1) {
                        headerReader = PrefilteringIndexReader::openNewHeaderReader(index,
                                PrefilteringIndexReader::HDR2INDEX, PrefilteringIndexReader::HDR2DATA, preload);
                    } else {
                        Debug(Debug::INFO) << "Index does not contain headers. Using normal database instead.\n";
                    }
                }
            } else {
                Debug(Debug::WARNING) << "Outdated index version. Please recompute it with 'createindex'!\n";
                index->close();
                delete index;
                index = NULL;
            }
        }

        if (sequenceReader == NULL) {
            const int dataMode = DBReader<unsigned int>::USE_INDEX | ((mode & NEED_SEQUENCES) ? DBReader<unsigned int>::USE_DATA : 0);
            //TODO threads
            sequenceReader = new DBReader<unsigned int>(dataName.c_str(), (dataName + ".index").c_str(), threads, dataMode);
            sequenceReader->open(DBReader<unsigned int>::NOSORT);
            if (preload) {
                sequenceReader->readMmapedDataInMemory();
            }
            seqType = sequenceReader->getDbtype();
        }

        if ((mode & NEED_HEADERS) && headerReader == NULL) {
            //TODO threads
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
