#include <cstdio>
#include <Parameters.h>

#include "DBReader.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "MMseqsNode.h"

#ifdef OPENMP
#include <omp.h>
#endif


int result2flat(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    const bool isLinfinity = DBReader<DBKeyType>::getExtendedDbtype(FileUtil::parseDbType(par.db1.c_str())) & Parameters::DBTYPE_EXTENDED_LINFINITY;
    // --createdb-mode 3: every node writes the representatives of its cluster database <clu>.<node> with their sequences to <out>.<node>
    std::string nodeSuffix;
    if (isLinfinity) {
        MMseqsNode::init();
        nodeSuffix = "." + SSTR(MMseqsNode::rank);
    }

    DBReader<DBKeyType> querydb_header(par.hdr1.c_str(), isLinfinity ? (par.db1 + ".nodeinfo").c_str() : par.hdr1Index.c_str(), par.threads, DBReader<DBKeyType>::USE_INDEX|DBReader<DBKeyType>::USE_DATA|DBReader<DBKeyType>::USE_LINFINITY);
    querydb_header.open(DBReader<DBKeyType>::NOSORT);
    querydb_header.readMmapedDataInMemory();

    DBReader<DBKeyType> targetdb_header(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, DBReader<DBKeyType>::USE_INDEX|DBReader<DBKeyType>::USE_DATA);
    DBReader<DBKeyType> seqDbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<DBKeyType>::USE_INDEX|DBReader<DBKeyType>::USE_DATA|DBReader<DBKeyType>::USE_LINFINITY);
    if (isLinfinity) {
        seqDbr.open(DBReader<DBKeyType>::NOSORT);
    } else {
        targetdb_header.open(DBReader<DBKeyType>::NOSORT);
        targetdb_header.readMmapedDataInMemory();
    }

    DBReader<DBKeyType> dbr_data((par.db3 + nodeSuffix).c_str(), (par.db3Index + nodeSuffix).c_str(), par.threads, DBReader<DBKeyType>::USE_INDEX|DBReader<DBKeyType>::USE_DATA|DBReader<DBKeyType>::USE_LINFINITY);
    dbr_data.open(DBReader<DBKeyType>::LINEAR_ACCCESS);

    FILE *fastaFP = fopen((par.db4 + nodeSuffix).c_str(), "w");

    // one split of all entries; --createdb-mode 3: splits of the keys of this node whose sequences, names and entries fit the memory, loaded and prefetched
    std::vector<size_t> splitEntryOffsets(1, 0);
    std::vector<size_t> splitKeyOffsets;
    if (isLinfinity) {
        const std::vector<size_t> nodeKeyOffsets = seqDbr.getFileKeyOffsets(MMseqsNode::count);
        const size_t nodeFirstKey = nodeKeyOffsets[MMseqsNode::rank];
        const size_t nodeEndKey = nodeKeyOffsets[MMseqsNode::rank + 1];
        // a key also holds its name with its offset and its share of the entries with their cells and keys
        const size_t keyBytes = (querydb_header.getNamesSize(MMseqsNode::rank) + dbr_data.getDataSize() + 4 * sizeof(uint64_t) * dbr_data.getSize()) / std::max<size_t>(1, nodeEndKey - nodeFirstKey) + sizeof(uint64_t);
        const size_t splitBytes = std::max<size_t>(1, Util::computeMemory(par.splitMemoryLimit) / 8);
        splitKeyOffsets = seqDbr.getSplitKeyOffsets(nodeFirstKey, nodeEndKey, splitBytes, keyBytes);
        for (size_t split = 1; split + 1 < splitKeyOffsets.size(); split++) {
            splitEntryOffsets.push_back(dbr_data.findEntryByKey(splitKeyOffsets[split]));
        }
    }
    splitEntryOffsets.push_back(dbr_data.getSize());

    const size_t chunkEntries = 1024;
    for (size_t split = 0; split + 1 < splitEntryOffsets.size(); split++) {
        if (isLinfinity) {
            querydb_header.loadData(splitKeyOffsets[split], splitKeyOffsets[split + 1]);
            seqDbr.loadData(splitKeyOffsets[split], splitKeyOffsets[split + 1]);
            dbr_data.loadData(splitEntryOffsets[split], splitEntryOffsets[split + 1]);
            if (split + 2 < splitEntryOffsets.size()) {
                querydb_header.prefetchData(splitKeyOffsets[split + 1], splitKeyOffsets[split + 2]);
                seqDbr.prefetchData(splitKeyOffsets[split + 1], splitKeyOffsets[split + 2]);
                dbr_data.prefetchData(splitEntryOffsets[split + 1], splitEntryOffsets[split + 2]);
            }
        }
        const size_t splitBegin = splitEntryOffsets[split];
        const size_t splitEnd = splitEntryOffsets[split + 1];
        // the threads format chunks of entries and write them in entry order
#pragma omp parallel
        {
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            std::string buffer;
#pragma omp for ordered schedule(dynamic, 1)
            for (size_t chunkBegin = splitBegin; chunkBegin < splitEnd; chunkBegin += chunkEntries) {
                buffer.clear();
                for (size_t i = chunkBegin; i < std::min(splitEnd, chunkBegin + chunkEntries); i++) {

                    // Write the header, taken from the original queryDB
                    buffer.push_back('>');
                    DBKeyType key = dbr_data.getDbKey(i);
                    char *header_data = isLinfinity ? querydb_header.getData(key, thread_idx) : querydb_header.getDataByDBKey(key, thread_idx);

                    std::string headerStr;
                    if (par.useHeader == true){
                        size_t lineLen = Util::skipLine(header_data) - header_data;
                        headerStr.assign(header_data, lineLen);

                        if (headerStr.length() > 0) {
                            if (headerStr[headerStr.length() - 1] == '\n') {
                                headerStr[headerStr.length() - 1] = ' ';
                            }
                        }
                    }else{
                        headerStr=Util::parseFastaHeader(header_data);
                    }
                    buffer.append(headerStr);
                    buffer.push_back('\n');

                    // --createdb-mode 3: the sequence of the representative with its newline, without its null byte
                    if (isLinfinity) {
                        buffer.append(seqDbr.getData(key, thread_idx), seqDbr.getSeqLen(key) + 1);
                        continue;
                    }

                    // write data
                    char *data = dbr_data.getData(i, thread_idx);
                    std::string dbKeyBuffer;
                    const char * words[2];
                    while (*data != '\0') {
                        // dbKeyBuffer can contain sequence
                        Util::getWordsOfLine(data, words, 2);
                        //Util::parseKey(data, dbKeyBuffer);
                        char *target_header_data = NULL;
                        size_t keyLen = 0;
                        for(size_t  i = 0; i < DbValidator::resultDb.size(); i++){
                            if(Parameters::isEqualDbtype(dbr_data.getDbtype(), DbValidator::resultDb[i])  ) {
                                keyLen = (words[1] - words[0]);
                                dbKeyBuffer.append(words[0], keyLen);
                                const DBKeyType dbKey = Util::fast_atoi<DBKeyType>(dbKeyBuffer.c_str());
                                target_header_data = targetdb_header.getDataByDBKey(dbKey, thread_idx);
                            }
                        }
                        std::string dataStr;
                        if (par.useHeader == true && target_header_data != NULL)
                        {
                            dataStr = Util::parseFastaHeader(target_header_data);
                            char *endLenData = Util::skipLine(data);
                            char *dataWithoutKey = data + keyLen;
                            size_t dataToCopySize = endLenData - dataWithoutKey;
                            dataStr.append(dataWithoutKey, dataToCopySize);
                        } else {
                            char *startLine = data;
                            char *endLine = Util::skipLine(data);
                            size_t n = endLine - startLine;
                            dataStr = std::string(startLine, n);
                        }

                        // newline at the end
                        if (dataStr.length() > 0) {
                            if (dataStr[dataStr.length() - 1] != '\n') {
                                dataStr.push_back('\n');
                            }
                        }

//                        std::cout << dataStr << std::endl;
                        buffer.append(dataStr);
                        dataStr.clear();
                        data = Util::skipLine(data);
                    }
                }
#pragma omp ordered
                fwrite(buffer.c_str(), sizeof(char), buffer.length(), fastaFP);
            }
        }
    }

    if (fclose(fastaFP) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << par.db4 << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (isLinfinity) {
        seqDbr.close();
    } else {
        targetdb_header.close();
    }
    querydb_header.close();
    dbr_data.close();

    return EXIT_SUCCESS;
}
