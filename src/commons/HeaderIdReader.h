//
// Created by Martin Steinegger on 19.10.18.
//

#ifndef MMSEQS_HEADERIDREADER_H
#define MMSEQS_HEADERIDREADER_H


#include <string>
#include "DBReader.h"
#include "Debug.h"
#include "PrefilteringIndexReader.h"


class HeaderIdReader {
public:
    HeaderIdReader(const std::string &dataName, const unsigned int headerIdx, const unsigned int dataIdx, bool preload)
            : reader(NULL), index(NULL) {
        std::string indexDB = PrefilteringIndexReader::searchForIndex(dataName.c_str());
        if (indexDB != "") {
            Debug(Debug::INFO) << "Use index  " << indexDB << "\n";
            int dataMode = DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA;
            index = new DBReader<unsigned int>(indexDB.c_str(), (indexDB + ".index").c_str(), dataMode);
            index->open(DBReader<unsigned int>::NOSORT);
            bool templateDBIsIndex = PrefilteringIndexReader::checkIfIndexFile(index);
            if (templateDBIsIndex == true) {
                PrefilteringIndexReader::printSummary(index);
                PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(index);

                if (data.headers == 1) {
                    reader = PrefilteringIndexReader::openNewHeaderReader(index, headerIdx, dataIdx, preload);
                } else {
                    Debug(Debug::INFO) << "Index does not contain headers. Using normal database instead.\n";
                }
            } else {
                Debug(Debug::WARNING) << "Outdated index version. Please recompute it with 'createindex'!\n";
                index->close();
                delete index;
                index = NULL;
            }
        }

        if (reader == NULL) {
            reader = new DBReader<unsigned int>((dataName + "_h").c_str(), (dataName + "_h.index").c_str());
            reader->open(DBReader<unsigned int>::NOSORT);


            if (preload) {
                reader->readMmapedDataInMemory();
            }
        }
    }

    std::string getId(unsigned int key) {
        size_t id = reader->getId(key);
        const char *data = reader->getData(id);
        return Util::parseFastaHeader(data);
    }

    ~HeaderIdReader() {
        reader->close();
        delete reader;

        if (index != NULL) {
            index->close();
            delete index;
        }
    }
    DBReader<unsigned int> * getReader(){
        return reader;
    }
private:
    DBReader<unsigned int> *reader;
    DBReader<unsigned int> *index;
};

#endif //MMSEQS_HEADERIDREADER_H
