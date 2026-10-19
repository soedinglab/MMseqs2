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

// the size of the exchanged member accession at record, 0 if bytes do not hold all of it
static size_t memberAccessionSize(const char *record, size_t bytes) {
    if (bytes < MEMBER_ACCESSION_SIZE) {
        return 0;
    }
    unsigned int accessionLength;
    memcpy(&accessionLength, record + 2 * sizeof(DBKeyType), sizeof(unsigned int));
    return (bytes >= MEMBER_ACCESSION_SIZE + accessionLength) ? MEMBER_ACCESSION_SIZE + accessionLength : 0;
}

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

static std::string readMemberAccessions(const std::string &sortedByMemberIndex, const std::string &sortedByMember, DBReader<DBKeyType> &headers, const std::vector<size_t> &nodeKeyOffsets,
                                        const std::string &headerOffsetsFile, size_t memoryLimit, const std::string &sortedByRep) {
    std::vector<MMseqsNode::ExchangeBlock> memberBlocks;
    std::vector<FILE *> memberFiles;
    MMseqsNode::openExchange(sortedByMemberIndex, sortedByMember, memberBlocks, memberFiles);

    // the last header offset of a node is the size of its names
    const size_t localFirstKey = nodeKeyOffsets[MMseqsNode::rank];
    const size_t localEndKey = nodeKeyOffsets[MMseqsNode::rank + 1];
    uint64_t nameBytes = 0;
    FILE *headerOffsets = FileUtil::openFileOrDie(headerOffsetsFile.c_str(), "r", true);
    if (fseek(headerOffsets, -static_cast<long>(sizeof(uint64_t)), SEEK_END) != 0 || fread(&nameBytes, sizeof(uint64_t), 1, headerOffsets) != 1 || fclose(headerOffsets) != 0) {
        Debug(Debug::ERROR) << "Cannot read " << headerOffsetsFile << "\n";
        EXIT(EXIT_FAILURE);
    }
    // a key also holds its name offset in the read buffer and in the load
    const size_t keyBytes = nameBytes / std::max<size_t>(1, localEndKey - localFirstKey) + sizeof(uint64_t) + sizeof(size_t);
    const size_t keysPerLoad = std::max<size_t>(1, memoryLimit / 8 / keyBytes);

    const std::string sortedByRepFile = sortedByRep + "." + SSTR(MMseqsNode::rank);
    FILE *sortedByRepHandle = FileUtil::openAndDelete(sortedByRepFile.c_str(), "w");
    std::string sortedByRepIndexData;
    size_t writtenBytes = 0;
    std::vector<MemberAccession> memberAccessions;
    std::string accessions;
    const int threads = Parameters::getInstance().threads;
    std::vector<std::string> threadAccessions(threads);
    std::vector<size_t> threadAccessionOffsets(threads + 1);
    std::vector<char> received;
    DBKeyType previousMemberKey = DB_KEY_INVALID;
    for (size_t keyBegin = localFirstKey; keyBegin < localEndKey; keyBegin += keysPerLoad) {
        const size_t keyEnd = std::min(localEndKey, keyBegin + keysPerLoad);
        headers.loadData(keyBegin, keyEnd);
        if (keyEnd < localEndKey) {
            headers.prefetchData(keyEnd, std::min(localEndKey, keyEnd + keysPerLoad));
        }
        const ClusterPair loadLast = { static_cast<DBKeyType>(keyEnd - 1), static_cast<DBKeyType>(keyEnd - 1) };
        size_t pairCount;
        while ((pairCount = MMseqsNode::readExchangeRange(memberBlocks, memoryLimit / 64, ClusterPair::recordSize, ClusterPair::compareByMember, loadLast, received, NULL, sortedByMember)) > 0) {
            ClusterPair *pairs = reinterpret_cast<ClusterPair *>(received.data());
            SORT_PARALLEL(pairs, pairs + pairCount, ClusterPair::compareByMember);
            const size_t firstAccession = memberAccessions.size();
            memberAccessions.resize(firstAccession + pairCount);
#pragma omp parallel
            {
                int thread_idx = 0;
#ifdef OPENMP
                thread_idx = omp_get_thread_num();
#endif
#pragma omp for schedule(dynamic, 1)
                for (int thread = 0; thread < threads; thread++) {
                    threadAccessions[thread].clear();
                    for (size_t i = pairCount * thread / threads; i < pairCount * (thread + 1) / threads; i++) {
                        const DBKeyType memberKey = pairs[i].memberKey;
                        if (memberKey == ((i > 0) ? pairs[i - 1].memberKey : previousMemberKey)) {
                            Debug(Debug::ERROR) << "Key " << memberKey << " is in two clusters\n";
                            EXIT(EXIT_FAILURE);
                        }
                        const std::string accession = Util::parseFastaHeader(headers.getData(memberKey, thread_idx));
                        const MemberAccession memberAccession = { pairs[i].repKey, memberKey, threadAccessions[thread].size(), static_cast<unsigned int>(accession.size()) };
                        memberAccessions[firstAccession + i] = memberAccession;
                        threadAccessions[thread].append(accession);
                    }
                }
            }
            threadAccessionOffsets[0] = accessions.size();
            for (int thread = 0; thread < threads; thread++) {
                threadAccessionOffsets[thread + 1] = threadAccessionOffsets[thread] + threadAccessions[thread].size();
            }
            accessions.resize(threadAccessionOffsets[threads]);
#pragma omp parallel for schedule(dynamic, 1)
            for (int thread = 0; thread < threads; thread++) {
                memcpy(&accessions[threadAccessionOffsets[thread]], threadAccessions[thread].data(), threadAccessions[thread].size());
                for (size_t i = pairCount * thread / threads; i < pairCount * (thread + 1) / threads; i++) {
                    memberAccessions[firstAccession + i].accessionOffset += threadAccessionOffsets[thread];
                }
            }
            previousMemberKey = pairs[pairCount - 1].memberKey;
            // the flush is copied once more into the write buffer
            if (accessions.size() + memberAccessions.size() * sizeof(MemberAccession) >= memoryLimit / 8) {
                writeMemberAccessions(memberAccessions, accessions, nodeKeyOffsets, sortedByRepHandle, sortedByRepFile, sortedByRepIndexData, writtenBytes);
                memberAccessions.clear();
                accessions.clear();
            }
        }
    }
    const ClusterPair lastPair = { DB_KEY_INVALID, DB_KEY_INVALID };
    if (MMseqsNode::readExchangeRange(memberBlocks, memoryLimit / 64, ClusterPair::recordSize, ClusterPair::compareByMember, lastPair, received, NULL, sortedByMember) > 0) {
        Debug(Debug::ERROR) << "Key " << reinterpret_cast<ClusterPair *>(received.data())->memberKey << " of the clustering is not in the sequence database\n";
        EXIT(EXIT_FAILURE);
    }
    writeMemberAccessions(memberAccessions, accessions, nodeKeyOffsets, sortedByRepHandle, sortedByRepFile, sortedByRepIndexData, writtenBytes);
    MMseqsNode::closeExchangeFiles(memberFiles, sortedByMember);
    return MMseqsNode::closeExchange(sortedByRepHandle, sortedByRepFile, sortedByRep + ".index", sortedByRepIndexData);
}

static void writeRepresentativeEntries(const std::string &sortedByRepIndex, const std::string &sortedByRep, const std::string &tsv, const Parameters &par) {
    std::vector<MMseqsNode::ExchangeBlock> repBlocks;
    std::vector<FILE *> repFiles;
    MMseqsNode::openExchange(sortedByRepIndex, sortedByRep, repBlocks, repFiles);
    const size_t memoryLimit = Util::computeMemory(par.splitMemoryLimit);
    const std::string dataFile = tsv + "." + SSTR(MMseqsNode::rank);
    const std::string indexFile = tsv + ".index." + SSTR(MMseqsNode::rank);
    DBWriter writer(dataFile.c_str(), indexFile.c_str(), 1, false, Parameters::DBTYPE_OMIT_FILE);
    writer.open();
    const ClusterPair lastPair = { DB_KEY_INVALID, DB_KEY_INVALID };
    std::vector<char> records;
    std::vector<size_t> recordOffsets;
    std::vector<MemberAccession> memberAccessions;
    std::vector<std::string> threadLines(par.threads);
    // the accession of the representative whose cluster the previous range ended in
    std::string repAccession;
    DBKeyType previousRepKey = DB_KEY_INVALID;
    size_t recordCount;
    while ((recordCount = MMseqsNode::readExchangeRange(repBlocks, memoryLimit / 8, memberAccessionSize, ClusterPair::compareByRep, lastPair, records, &recordOffsets, sortedByRep)) > 0) {
        memberAccessions.resize(recordCount);
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < recordCount; i++) {
            const char *record = records.data() + recordOffsets[i];
            MemberAccession &memberAccession = memberAccessions[i];
            memcpy(&memberAccession.repKey, record, sizeof(DBKeyType));
            memcpy(&memberAccession.memberKey, record + sizeof(DBKeyType), sizeof(DBKeyType));
            memcpy(&memberAccession.accessionLength, record + 2 * sizeof(DBKeyType), sizeof(unsigned int));
            memberAccession.accessionOffset = recordOffsets[i] + MEMBER_ACCESSION_SIZE;
        }
        SORT_PARALLEL(memberAccessions.begin(), memberAccessions.end(), MemberAccession::compareByRep);
#pragma omp parallel for schedule(dynamic, 1)
        for (int thread = 0; thread < par.threads; thread++) {
            std::string &lines = threadLines[thread];
            lines.clear();
            const size_t begin = recordCount * thread / par.threads;
            const size_t end = recordCount * (thread + 1) / par.threads;
            if (begin == end) {
                continue;
            }
            DBKeyType repKey = (begin > 0) ? memberAccessions[begin - 1].repKey : previousRepKey;
            const char *repName = repAccession.data();
            size_t repLength = repAccession.size();
            // a cluster that continues from the part before has its representative line earlier in this range or in a previous range
            const MemberAccession repLine = { repKey, repKey, 0, 0 };
            const size_t repLinePos = std::lower_bound(memberAccessions.begin(), memberAccessions.begin() + begin, repLine, MemberAccession::compareByRep) - memberAccessions.begin();
            if (repLinePos < begin && memberAccessions[repLinePos].memberKey == repKey) {
                repName = records.data() + memberAccessions[repLinePos].accessionOffset;
                repLength = memberAccessions[repLinePos].accessionLength;
            }
            for (size_t i = begin; i < end; i++) {
                const MemberAccession &memberAccession = memberAccessions[i];
                if (memberAccession.repKey != repKey) {
                    // the representative itself comes first
                    if (memberAccession.memberKey != memberAccession.repKey) {
                        Debug(Debug::ERROR) << "The cluster of key " << memberAccession.repKey << " has no line for its representative\n";
                        EXIT(EXIT_FAILURE);
                    }
                    repKey = memberAccession.repKey;
                    repName = records.data() + memberAccession.accessionOffset;
                    repLength = memberAccession.accessionLength;
                }
                lines.append(repName, repLength);
                lines.push_back('\t');
                lines.append(records.data() + memberAccession.accessionOffset, memberAccession.accessionLength);
                lines.push_back('\n');
            }
        }
        for (int thread = 0; thread < par.threads; thread++) {
            writer.writeAdd(threadLines[thread].data(), threadLines[thread].size(), 0);
        }
        previousRepKey = memberAccessions[recordCount - 1].repKey;
        const MemberAccession repLine = { previousRepKey, previousRepKey, 0, 0 };
        const size_t repLinePos = std::lower_bound(memberAccessions.begin(), memberAccessions.end(), repLine, MemberAccession::compareByRep) - memberAccessions.begin();
        if (memberAccessions[repLinePos].memberKey == previousRepKey) {
            repAccession.assign(records.data() + memberAccessions[repLinePos].accessionOffset, memberAccessions[repLinePos].accessionLength);
        }
    }
    writer.close(true, false);
    FileUtil::remove(indexFile.c_str());
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
