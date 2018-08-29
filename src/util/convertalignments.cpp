#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "PrefilteringIndexReader.h"
#include "FileUtil.h"

#ifdef OPENMP
#include <omp.h>
#endif

class HeaderIdReader {
public:
    HeaderIdReader(const std::string &dataName, bool noPreload)
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
                    reader = PrefilteringIndexReader::openNewHeaderReader(index, (dataName + "_h").c_str(), noPreload == false);
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


            if (noPreload == false) {
                reader->readMmapedDataInMemory();
                reader->mlock();
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

private:
    DBReader<unsigned int> *reader;
    DBReader<unsigned int> *index;
};

void printSeqBasedOnAln(std::string &out, const Sequence *seq, unsigned int offset, const std::string &bt, bool reverse) {
    unsigned int seqPos = 0;
    for (uint32_t i = 0; i < bt.size(); ++i) {
        switch (bt[i]) {
            case 'M':
                out.append(1, seq->subMat->int2aa[seq->int_sequence[offset + seqPos]]);
                seqPos++;
                break;
            case 'I':
                if (reverse == true) {
                    out.append(1, '-');
                } else {
                    out.append(1, seq->subMat->int2aa[seq->int_sequence[offset + seqPos]]);
                    seqPos++;
                }
                break;
            case 'D':
                if (reverse == true) {
                    out.append(1, seq->subMat->int2aa[seq->int_sequence[offset + seqPos]]);
                    seqPos++;
                } else {
                    out.append(1, '-');
                }
                break;
        }
    }
}

int convertalignments(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4);

    DBReader<unsigned int> *queryReader = NULL;
    DBReader<unsigned int> *targetReader = NULL;

    const bool sameDB = par.db1.compare(par.db2) == 0 ? true : false;
    const int format = par.formatAlignmentMode;
    const bool needSequenceDB =
            format == Parameters::FORMAT_ALIGNMENT_PAIRWISE
            || format == Parameters::FORMAT_ALIGNMENT_SAM;

    if (needSequenceDB) {
        targetReader = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str());
        targetReader->open(DBReader<unsigned int>::NOSORT);
        if (sameDB == true) {
            queryReader = targetReader;
        } else {
            queryReader = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str());
            queryReader->open(DBReader<unsigned int>::NOSORT);
        }
    }
    Debug(Debug::INFO) << "Query Header file: " << par.db1 << "_h\n";
    HeaderIdReader qHeaderDbr(par.db1.c_str(), par.noPreload);

    HeaderIdReader *tHeaderDbr;
    if(sameDB){
        tHeaderDbr = &qHeaderDbr;
    } else {
        Debug(Debug::INFO) << "Target Header file: " << par.db2 << "_h\n";
        tHeaderDbr = new HeaderIdReader(par.db2.c_str(), par.noPreload);
    }

    Debug(Debug::INFO) << "Alignment database: " << par.db3 << "\n";
    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str());
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

#ifdef OPENMP
    unsigned int totalThreads = par.threads;
#else
    unsigned int totalThreads = 1;
#endif

    unsigned int localThreads = totalThreads;
    if (alnDbr.getSize() <= totalThreads) {
        localThreads = alnDbr.getSize();
    }

    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), localThreads);
    resultWriter.open();
    Debug(Debug::INFO) << "Start writing file to " << par.db4 << "\n";
    bool isDb = par.dbOut;
    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0f, -0.2f);

#pragma omp parallel num_threads(localThreads)
    {
        Sequence *querySeq;
        Sequence *targetSeq;
        if (needSequenceDB) {
            querySeq = new Sequence(par.maxSeqLen, queryReader->getDbtype(), &subMat, 0, false, false, false);
            targetSeq = new Sequence(par.maxSeqLen, targetReader->getDbtype(), &subMat, 0, false, false, false);
        }

        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        std::string result;
        result.reserve(1024*1024);

        char buffer[1024];
        std::vector<Matcher::result_t> results;
        results.reserve(300);

#pragma omp  for schedule(dynamic, 10)
        for (size_t i = 0; i < alnDbr.getSize(); i++) {
            Debug::printProgress(i);

            unsigned int queryKey = alnDbr.getDbKey(i);
            char *data = alnDbr.getData(i);

            if (needSequenceDB) {
                querySeq->mapSequence(i, queryKey, queryReader->getDataByDBKey(queryKey));
            }

            std::string queryId = qHeaderDbr.getId(queryKey);
            Matcher::readAlignmentResults(results, data, true);
            unsigned int missMatchCount;
            for (size_t j = 0; j < results.size(); j++) {
                const Matcher::result_t &res = results[j];
                std::string targetId = tHeaderDbr->getId(res.dbKey);
                unsigned int gapOpenCount = 0;
                unsigned int alnLen = res.alnLength;
                if (res.backtrace.size() > 0) {
                    size_t matchCount = 0;
                    alnLen = 0;
                    for (size_t pos = 0; pos < res.backtrace.size(); pos++) {
                        int cnt = 0;
                        if (isdigit(res.backtrace[pos])) {
                            cnt += Util::fast_atoi<int>(res.backtrace.c_str() + pos);
                            while (isdigit(res.backtrace[pos])) {
                                pos++;
                            }
                        }
                        alnLen += cnt;

                        switch (res.backtrace[pos]) {
                            case 'M':
                                matchCount += cnt;
                                break;
                            case 'D':
                            case 'I':
                                gapOpenCount += 1;
                                break;
                        }
                    }
//                res.seqId = X / alnLen;
                    unsigned int identical = static_cast<unsigned int>( res.seqId * static_cast<float>(alnLen) + 0.5 );
                    missMatchCount = static_cast<unsigned int>( matchCount - identical);
                } else {
                    int adjustQstart = (res.qStartPos == -1) ? 0 : res.qStartPos;
                    int adjustDBstart = (res.dbStartPos == -1) ? 0 : res.dbStartPos;
                    float bestMatchEstimate = static_cast<float>(std::min(res.qEndPos - adjustQstart,
                                                                          res.dbEndPos - adjustDBstart));
                    missMatchCount = static_cast<unsigned int>( bestMatchEstimate * (1.0f - res.seqId) + 0.5 );
                }


                switch (format) {
                    case Parameters::FORMAT_ALIGNMENT_BLAST_TAB: {
                        int count = snprintf(buffer, sizeof(buffer),
                                             "%s\t%s\t%1.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2E\t%d\n",
                                             queryId.c_str(), targetId.c_str(), res.seqId, alnLen,
                                             missMatchCount, gapOpenCount,
                                             res.qStartPos + 1, res.qEndPos + 1,
                                             res.dbStartPos + 1, res.dbEndPos + 1,
                                             res.eval, res.score);

                        if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                            Debug(Debug::WARNING) << "Truncated line in entry" << j << "!\n";
                            continue;
                        }

                        result.append(buffer, count);
                        break;
                    }
                    case Parameters::FORMAT_ALIGNMENT_BLAST_WITH_LEN: {
                        int count = snprintf(buffer, sizeof(buffer),
                                             "%s\t%s\t%1.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2E\t%d\t%d\t%d\n",
                                             queryId.c_str(), targetId.c_str(), res.seqId, alnLen,
                                             missMatchCount, gapOpenCount,
                                             res.qStartPos + 1, res.qEndPos + 1,
                                             res.dbStartPos + 1, res.dbEndPos + 1,
                                             res.eval, res.score,
                                             res.qLen, res.dbLen);

                        if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                            Debug(Debug::WARNING) << "Truncated line in entry" << j << "!\n";
                            continue;
                        }

                        result.append(buffer, count);
                        break;
                    }
                    case Parameters::FORMAT_ALIGNMENT_PAIRWISE: {
                        int count = snprintf(buffer, sizeof(buffer),
                                             ">%s\t%s\t%1.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2E\t%d\n",
                                             queryId.c_str(), targetId.c_str(), res.seqId, alnLen, missMatchCount,
                                             gapOpenCount,
                                             res.qStartPos + 1, res.qEndPos + 1, res.dbStartPos + 1, res.dbEndPos + 1,
                                             res.eval, res.score);

                        if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                            Debug(Debug::WARNING) << "Truncated line in entry" << j << "!\n";
                            continue;
                        }

                        result.append(buffer, count);

                        const std::string &backtrace = Matcher::uncompressAlignment(res.backtrace);
                        printSeqBasedOnAln(result, querySeq, res.qStartPos, backtrace, false);
                        result.append(1, '\n');

                        targetSeq->mapSequence(i, res.dbKey, targetReader->getDataByDBKey(res.dbKey));

                        printSeqBasedOnAln(result, targetSeq, res.dbStartPos, backtrace, true);
                        result.append(1, '\n');
                        break;
                    }
                    case Parameters::FORMAT_ALIGNMENT_SAM:
                    default:
                        Debug(Debug::ERROR) << "Not implemented yet";
                        EXIT(EXIT_FAILURE);
                }
            }

            resultWriter.writeData(result.c_str(), result.size(), queryKey, thread_idx, isDb);
            results.clear();
            result.clear();
        }
        if (needSequenceDB) {
            delete querySeq;
            delete targetSeq;
        }
    }
    resultWriter.close();

    // tsv output
    if (isDb == false) {
        FileUtil::deleteFile(par.db4Index);
    }

    alnDbr.close();
    if (sameDB == false) {
        delete tHeaderDbr;
    }

    if (needSequenceDB) {
        queryReader->close();
        delete queryReader;
        if (sameDB == false) {
            targetReader->close();
            delete targetReader;
        }
    }

    return EXIT_SUCCESS;
}

