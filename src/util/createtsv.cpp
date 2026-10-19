#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "IndexReader.h"
#include "FileUtil.h"
#include "MMseqsNode.h"
#include "FastSort.h"

#ifdef OPENMP
#include <omp.h>

#endif

#ifndef SIZE_T_MAX
#define SIZE_T_MAX ((size_t) -1)
#endif

struct MemberAccession {
    DBKeyType repKey;
    DBKeyType memberKey;
    size_t accessionOffset;
    unsigned int accessionLength;

    static bool compareByRep(const MemberAccession &first, const MemberAccession &second) {
        const ClusterPair firstPair = { first.repKey, first.memberKey };
        const ClusterPair secondPair = { second.repKey, second.memberKey };
        return ClusterPair::compareByRep(firstPair, secondPair);
    }
};

static const size_t MEMBER_ACCESSION_SIZE = 2 * sizeof(DBKeyType) + sizeof(unsigned int);

struct QueuedMemberAccession {
    DBKeyType repKey;
    DBKeyType memberKey;
    const char *accession;
    unsigned int accessionLength;
    size_t block;
};

struct CompareQueuedMemberAccession {
    bool operator()(const QueuedMemberAccession &first, const QueuedMemberAccession &second) const {
        const ClusterPair firstPair = { first.repKey, first.memberKey };
        const ClusterPair secondPair = { second.repKey, second.memberKey };
        return ClusterPair::compareByRep(secondPair, firstPair);
    }
};
typedef std::priority_queue<QueuedMemberAccession, std::vector<QueuedMemberAccession>, CompareQueuedMemberAccession> MemberAccessionQueue;

static void writeMemberAccessions(std::vector<MemberAccession> &memberAccessions, const std::string &accessions, const std::vector<size_t> &nodeKeyOffsets,
                                  FILE *file, const std::string &fileName, std::string &indexData, size_t &writtenBytes) {
    SORT_PARALLEL(memberAccessions.begin(), memberAccessions.end(), MemberAccession::compareByRep);
    std::vector<size_t> nodeOffsets(MMseqsNode::count + 1, 0);
    std::string writeBuffer;
    size_t node = 0;
    for (size_t i = 0; i < memberAccessions.size(); i++) {
        const MemberAccession &memberAccession = memberAccessions[i];
        for (; node < MMseqsNode::count && nodeKeyOffsets[node] <= memberAccession.repKey; node++) {
            nodeOffsets[node] = writeBuffer.size();
        }
        writeBuffer.append(reinterpret_cast<const char *>(&memberAccession.repKey), sizeof(DBKeyType));
        writeBuffer.append(reinterpret_cast<const char *>(&memberAccession.memberKey), sizeof(DBKeyType));
        writeBuffer.append(reinterpret_cast<const char *>(&memberAccession.accessionLength), sizeof(unsigned int));
        writeBuffer.append(accessions, memberAccession.accessionOffset, memberAccession.accessionLength);
    }
    for (; node <= MMseqsNode::count; node++) {
        nodeOffsets[node] = writeBuffer.size();
    }
    if (fwrite(writeBuffer.data(), sizeof(char), writeBuffer.size(), file) != writeBuffer.size()) {
        Debug(Debug::ERROR) << "Cannot write " << fileName << "\n";
        EXIT(EXIT_FAILURE);
    }
    MMseqsNode::addExchangeFlush(indexData, writtenBytes, nodeOffsets);
}

static void queueNextMemberAccession(MemberAccessionQueue &queue, std::vector<MMseqsNode::ExchangeBlock> &blocks, size_t block, const std::string &prefix) {
    const char *data = MMseqsNode::readExchange(blocks[block], MEMBER_ACCESSION_SIZE, prefix);
    if (data == NULL) {
        return;
    }
    QueuedMemberAccession queuedMemberAccession;
    memcpy(&queuedMemberAccession.repKey, data, sizeof(DBKeyType));
    memcpy(&queuedMemberAccession.memberKey, data + sizeof(DBKeyType), sizeof(DBKeyType));
    memcpy(&queuedMemberAccession.accessionLength, data + 2 * sizeof(DBKeyType), sizeof(unsigned int));
    queuedMemberAccession.accession = MMseqsNode::readExchange(blocks[block], queuedMemberAccession.accessionLength, prefix);
    queuedMemberAccession.block = block;
    queue.push(queuedMemberAccession);
}

static std::string readMemberAccessions(const std::string &sortedByMemberIndex, const std::string &sortedByMember, DBReader<DBKeyType> &headers, const std::vector<size_t> &nodeKeyOffsets,
                                        const std::string &headerOffsetsFile, size_t memoryLimit, const std::string &sortedByRep) {
    const CompareQueuedPair compareQueuedByMember = { ClusterPair::compareByMember };
    PairQueue memberQueue(compareQueuedByMember);
    std::vector<MMseqsNode::ExchangeBlock> memberBlocks;
    std::vector<FILE *> memberFiles;
    MMseqsNode::openPairs(memberQueue, memberBlocks, memberFiles, sortedByMemberIndex, sortedByMember);

    // the last header offset of a node is the size of its names
    const size_t localFirstKey = nodeKeyOffsets[MMseqsNode::rank];
    const size_t localEndKey = nodeKeyOffsets[MMseqsNode::rank + 1];
    uint64_t nameBytes = 0;
    FILE *headerOffsets = FileUtil::openFileOrDie(headerOffsetsFile.c_str(), "r", true);
    if (fseek(headerOffsets, -static_cast<long>(sizeof(uint64_t)), SEEK_END) != 0 || fread(&nameBytes, sizeof(uint64_t), 1, headerOffsets) != 1 || fclose(headerOffsets) != 0) {
        Debug(Debug::ERROR) << "Cannot read " << headerOffsetsFile << "\n";
        EXIT(EXIT_FAILURE);
    }
    const size_t keysPerLoad = std::max<size_t>(1, memoryLimit / 8 / (nameBytes / std::max<size_t>(1, localEndKey - localFirstKey) + 1));

    const std::string sortedByRepFile = sortedByRep + "." + SSTR(MMseqsNode::rank);
    FILE *sortedByRepHandle = FileUtil::openAndDelete(sortedByRepFile.c_str(), "w");
    std::string sortedByRepIndexData;
    size_t writtenBytes = 0;
    std::vector<MemberAccession> memberAccessions;
    std::string accessions;
    DBKeyType previousMemberKey = DB_KEY_INVALID;
    for (size_t keyBegin = localFirstKey; keyBegin < localEndKey; keyBegin += keysPerLoad) {
        const size_t keyEnd = std::min(localEndKey, keyBegin + keysPerLoad);
        headers.loadData(keyBegin, keyEnd);
        if (keyEnd < localEndKey) {
            headers.prefetchData(keyEnd, std::min(localEndKey, keyEnd + keysPerLoad));
        }
        while (memberQueue.empty() == false && memberQueue.top().pair.memberKey < keyEnd) {
            const QueuedPair queuedPair = memberQueue.top();
            memberQueue.pop();
            MMseqsNode::queueNextPair(memberQueue, memberBlocks, queuedPair.block, sortedByMember);
            if (queuedPair.pair.memberKey == previousMemberKey) {
                Debug(Debug::ERROR) << "Key " << queuedPair.pair.memberKey << " is in two clusters\n";
                EXIT(EXIT_FAILURE);
            }
            previousMemberKey = queuedPair.pair.memberKey;
            const std::string accession = Util::parseFastaHeader(headers.getData(queuedPair.pair.memberKey, 0));
            const MemberAccession memberAccession = { queuedPair.pair.repKey, queuedPair.pair.memberKey, accessions.size(), static_cast<unsigned int>(accession.size()) };
            accessions.append(accession);
            memberAccessions.push_back(memberAccession);
            // the flush is copied once more into the write buffer
            if (accessions.size() + memberAccessions.size() * sizeof(MemberAccession) >= memoryLimit / 8) {
                writeMemberAccessions(memberAccessions, accessions, nodeKeyOffsets, sortedByRepHandle, sortedByRepFile, sortedByRepIndexData, writtenBytes);
                memberAccessions.clear();
                accessions.clear();
            }
        }
    }
    if (memberQueue.empty() == false) {
        Debug(Debug::ERROR) << "Key " << memberQueue.top().pair.memberKey << " of the clustering is not in the sequence database\n";
        EXIT(EXIT_FAILURE);
    }
    writeMemberAccessions(memberAccessions, accessions, nodeKeyOffsets, sortedByRepHandle, sortedByRepFile, sortedByRepIndexData, writtenBytes);
    MMseqsNode::closeExchangeFiles(memberFiles, sortedByMember);
    return MMseqsNode::closeExchange(sortedByRepHandle, sortedByRepFile, sortedByRep + ".index", sortedByRepIndexData);
}

static void writeRepresentativeEntries(const std::string &sortedByRepIndex, const std::string &sortedByRep, const std::string &tsv, const Parameters &par) {
    MemberAccessionQueue repQueue;
    std::vector<MMseqsNode::ExchangeBlock> repBlocks;
    std::vector<FILE *> repFiles;
    MMseqsNode::openExchange(sortedByRepIndex, sortedByRep, repBlocks, repFiles);
    for (size_t block = 0; block < repBlocks.size(); block++) {
        queueNextMemberAccession(repQueue, repBlocks, block, sortedByRep);
    }
    const std::string dataFile = tsv + "." + SSTR(MMseqsNode::rank);
    const std::string indexFile = tsv + ".index." + SSTR(MMseqsNode::rank);
    const bool shouldCompress = par.dbOut == true && par.compressed == true;
    const int dbType = par.dbOut == true ? Parameters::DBTYPE_GENERIC_DB : Parameters::DBTYPE_OMIT_FILE;
    DBWriter writer(dataFile.c_str(), indexFile.c_str(), 1, shouldCompress, dbType);
    writer.open();
    std::string repAccession;
    DBKeyType currentRepKey = DB_KEY_INVALID;
    while (repQueue.empty() == false) {
        const QueuedMemberAccession queuedMemberAccession = repQueue.top();
        repQueue.pop();
        if (queuedMemberAccession.repKey != currentRepKey) {
            if (currentRepKey != DB_KEY_INVALID) {
                writer.writeEnd(currentRepKey, 0, par.dbOut);
            }
            // the representative itself comes first
            if (queuedMemberAccession.memberKey != queuedMemberAccession.repKey) {
                Debug(Debug::ERROR) << "The cluster of key " << queuedMemberAccession.repKey << " has no line for its representative\n";
                EXIT(EXIT_FAILURE);
            }
            currentRepKey = queuedMemberAccession.repKey;
            repAccession.assign(queuedMemberAccession.accession, queuedMemberAccession.accessionLength);
            writer.writeStart(0);
        }
        writer.writeAdd(repAccession.c_str(), repAccession.size(), 0);
        writer.writeAdd("\t", 1, 0);
        writer.writeAdd(queuedMemberAccession.accession, queuedMemberAccession.accessionLength, 0);
        writer.writeAdd("\n", 1, 0);
        queueNextMemberAccession(repQueue, repBlocks, queuedMemberAccession.block, sortedByRep);
    }
    if (currentRepKey != DB_KEY_INVALID) {
        writer.writeEnd(currentRepKey, 0, par.dbOut);
    }
    writer.close(par.dbOut == false);
    if (par.dbOut == false) {
        FileUtil::remove(indexFile.c_str());
    }
    MMseqsNode::closeExchangeFiles(repFiles, sortedByRep);
}

int createtsv(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    const bool isLinfinity = DBReader<DBKeyType>::getExtendedDbtype(FileUtil::parseDbType(par.db1.c_str())) & Parameters::DBTYPE_EXTENDED_LINFINITY;
    if (isLinfinity) {
        MMseqsNode::init();
        if (par.filenames.size() <= 3 || par.db2.compare(par.db1) != 0) {
            Debug(Debug::ERROR) << "createtsv needs the --createdb-mode " << Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH << " database as both query and target database: createtsv DB DB clu tsv\n";
            EXIT(EXIT_FAILURE);
        }
        if (par.fullHeader || par.targetTsvColumn != 1 || par.idxSeqSrc != 0 || par.dbOut) {
            Debug(Debug::ERROR) << "--full-header, --target-column, --idx-seq-src and --db-output are not supported for --createdb-mode " << Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH << " databases\n";
            EXIT(EXIT_FAILURE);
        }
        const int clusterDbtype = FileUtil::parseDbType(par.db3.c_str());
        if (Parameters::isEqualDbtype(clusterDbtype, Parameters::DBTYPE_CLUSTER_RES) == false
            || (DBReader<DBKeyType>::getExtendedDbtype(clusterDbtype) & Parameters::DBTYPE_EXTENDED_LINFINITY) == 0) {
            Debug(Debug::ERROR) << par.db3 << " is not a --createdb-mode " << Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH << " clustering\n";
            EXIT(EXIT_FAILURE);
        }
        DBReader<DBKeyType> headers(par.hdr1.c_str(), (par.db1 + ".nodeinfo").c_str(), par.threads, DBReader<DBKeyType>::USE_INDEX|DBReader<DBKeyType>::USE_DATA|DBReader<DBKeyType>::USE_LINFINITY);
        headers.open(DBReader<DBKeyType>::NOSORT);
        const std::vector<size_t> nodeKeyOffsets = headers.getFileKeyOffsets(MMseqsNode::count);
        const std::string sortedByMember = par.db4 + "_sortedbymem";
        const std::string sortedByMemberIndex = MMseqsNode::sendClusterPairs(par.db3, true, ClusterPair::compareByMember, nodeKeyOffsets, sortedByMember);
        const std::string sortedByRep = par.db4 + "_sortedbyrep";
        const std::string sortedByRepIndex = readMemberAccessions(sortedByMemberIndex, sortedByMember, headers, nodeKeyOffsets, par.hdr1 + ".index." + SSTR(MMseqsNode::rank), Util::computeMemory(par.splitMemoryLimit), sortedByRep);
        FileUtil::remove((sortedByMember + "." + SSTR(MMseqsNode::rank)).c_str());
        writeRepresentativeEntries(sortedByRepIndex, sortedByRep, par.db4, par);
        std::vector<std::string> exchanges;
        exchanges.push_back(sortedByMember);
        exchanges.push_back(sortedByRep);
        MMseqsNode::removeExchanges(sortedByRep, sortedByRep, exchanges);
        headers.close();
        return EXIT_SUCCESS;
    }

    bool queryNucs = Parameters::isEqualDbtype(FileUtil::parseDbType(par.db1.c_str()), Parameters::DBTYPE_NUCLEOTIDES);
    bool targetNucs = Parameters::isEqualDbtype(FileUtil::parseDbType(par.db2.c_str()), Parameters::DBTYPE_NUCLEOTIDES);
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    int queryHeaderType = (queryNucs) ? IndexReader::SRC_HEADERS : IndexReader::HEADERS;
    queryHeaderType = (par.idxSeqSrc == 0) ? queryHeaderType :  (par.idxSeqSrc == 1) ?  IndexReader::HEADERS : IndexReader::SRC_HEADERS;
    IndexReader qDbrHeader(par.db1, par.threads, queryHeaderType, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    IndexReader * tDbrHeader=NULL;
    DBReader<DBKeyType> * queryDB = qDbrHeader.sequenceReader;
    DBReader<DBKeyType> * targetDB = NULL;
    bool sameDB = (par.db2.compare(par.db1) == 0);
    const bool hasTargetDB = par.filenames.size() > 3;
    DBReader<DBKeyType>::Index * qHeaderIndex = qDbrHeader.sequenceReader->getIndex();
    DBReader<DBKeyType>::Index * tHeaderIndex = NULL;

    if (hasTargetDB) {
        if (sameDB) {
            tDbrHeader = &qDbrHeader;
            tHeaderIndex = qHeaderIndex;
            targetDB = queryDB;
        } else {

            int targetHeaderType = (targetNucs) ? IndexReader::SRC_HEADERS : IndexReader::HEADERS;
            targetHeaderType = (par.idxSeqSrc == 0) ? targetHeaderType :  (par.idxSeqSrc == 1) ?  IndexReader::HEADERS : IndexReader::SRC_HEADERS;

            tDbrHeader = new IndexReader(par.db2, par.threads, targetHeaderType, touch);
            tHeaderIndex = tDbrHeader->sequenceReader->getIndex();
            targetDB = tDbrHeader->sequenceReader;
        }
    }

    DBReader<DBKeyType> *reader;
    if (hasTargetDB) {

        reader = new DBReader<DBKeyType>(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<DBKeyType>::USE_INDEX|DBReader<DBKeyType>::USE_DATA);
    } else {

        reader = new DBReader<DBKeyType>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<DBKeyType>::USE_INDEX|DBReader<DBKeyType>::USE_DATA);
    }
    reader->open(DBReader<DBKeyType>::LINEAR_ACCCESS);

    uint16_t extended = DBReader<DBKeyType>::getExtendedDbtype(reader->getDbtype());
    bool needSET = false;
    std::map<DBKeyType, std::string> qSetToSource, tSetToSource;
    if (extended & Parameters::DBTYPE_EXTENDED_SET) {
        needSET = true;
        if (hasTargetDB) {
            qSetToSource = Util::readLookup((par.db1 + ".source"), 2);
            tSetToSource = Util::readLookup((par.db2 + ".source"), 2);
        } else {
            qSetToSource = Util::readLookup((par.db1 + ".source"), 2);
        }
    }

    const std::string& dataFile = hasTargetDB ? par.db4 : par.db3;
    const std::string& indexFile = hasTargetDB ? par.db4Index : par.db3Index;
    const bool shouldCompress = par.dbOut == true && par.compressed == true;
    const int dbType = par.dbOut == true ? Parameters::DBTYPE_GENERIC_DB : Parameters::DBTYPE_OMIT_FILE;
    DBWriter writer(dataFile.c_str(), indexFile.c_str(), par.threads, shouldCompress, dbType);
    writer.open();

    const size_t targetColumn = (par.targetTsvColumn == 0) ? SIZE_T_MAX :  par.targetTsvColumn - 1;
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

        const char *columnPointer[255];
        char *dbKey = new char[par.maxSeqLen + 1];

        std::string outputBuffer;
        outputBuffer.reserve(10 * 1024);

#pragma omp for schedule(dynamic, 1000)
        for (size_t i = 0; i < reader->getSize(); ++i) {
            DBKeyType queryKey = reader->getDbKey(i);
            size_t queryIndex;
            char *headerData;
            if(needSET == false) {
                queryIndex = queryDB->getId(queryKey);
                headerData = queryDB->getData(queryIndex, thread_idx);
                if (headerData == NULL) {
                    Debug(Debug::WARNING) << "Invalid header entry in query " << queryKey << "!\n";
                    continue;
                }
            }

            std::string queryHeader;
            if (needSET == true) {
                queryHeader = qSetToSource[queryKey];
            } else if (par.fullHeader) {
                queryHeader = "\"";
                queryHeader.append(headerData, qHeaderIndex[queryIndex].length - 2);
                queryHeader.append("\"");
            } else {
                queryHeader = Util::parseFastaHeader(headerData);
            }

            size_t entryIndex = 0;

            char *data = reader->getData(i, thread_idx);
            while (*data != '\0') {
                if(targetColumn != SIZE_T_MAX){
                    size_t foundElements = Util::getWordsOfLine(data, columnPointer, 255);
                    if (foundElements < targetColumn) {
                        Debug(Debug::WARNING) << "Not enough columns!" << "\n";
                        continue;
                    }
                    Util::parseKey(columnPointer[targetColumn], dbKey);
                }
                std::string targetAccession;
                if(targetColumn == SIZE_T_MAX){
                    targetAccession = "";
                } else if (hasTargetDB) {
                    DBKeyType targetKey = Util::fast_atoi<DBKeyType>(dbKey);
                    size_t targetIndex = targetDB->getId(targetKey);
                    char *targetData;
                    if(needSET == false) {
                        targetData = targetDB->getData(targetIndex, thread_idx);
                        if (targetData == NULL) {
                            Debug(Debug::WARNING) << "Invalid header entry in query " << queryKey << " and target " << targetKey << "!\n";
                            continue;
                        }
                        if (par.fullHeader) {
                            targetAccession = "\"";
                            targetAccession.append(targetData, tHeaderIndex[targetIndex].length - 2);
                            targetAccession.append("\"");
                        } else {
                            targetAccession = Util::parseFastaHeader(targetData);
                        }
                    } else {
                        targetAccession = tSetToSource[targetKey];
                    }
                } else {
                    targetAccession = dbKey;
                }

                if (par.firstSeqRepr && !entryIndex) {
                    queryHeader = targetAccession;
                }

                outputBuffer.append(queryHeader);
                outputBuffer.append("\t");
                outputBuffer.append(targetAccession);

                size_t offset = 0;
                if (targetColumn != 0) {
                    outputBuffer.append("\t");
                    offset = 0;
                } else {
                    offset = strlen(dbKey);
                }

                char *nextLine = Util::skipLine(data);
                outputBuffer.append(data + offset, (nextLine - (data + offset)) - 1);
                outputBuffer.append("\n");
                data = nextLine;
                entryIndex++;
            }
            writer.writeData(outputBuffer.c_str(), outputBuffer.length(), queryKey, thread_idx, par.dbOut);
            outputBuffer.clear();
        }
        delete[] dbKey;
    }
    writer.close(par.dbOut == false);

    if (par.dbOut == false) {
        if (hasTargetDB) {
            FileUtil::remove(par.db4Index.c_str());
        } else {
            FileUtil::remove(par.db3Index.c_str());
        }
    }

    reader->close();
    delete reader;
    if (hasTargetDB) {
        if (sameDB == false) {
            delete tDbrHeader;
        }
    }
    qSetToSource.clear();
    tSetToSource.clear();

    return EXIT_SUCCESS;
}
#undef SIZE_T_MAX
