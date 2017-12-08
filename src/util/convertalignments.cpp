#include <string>
#include <vector>
#include <fstream>
#include <FileUtil.h>

#include "Alignment.h"
#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "PrefilteringIndexReader.h"

#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"

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

void printSeqBasedOnAln(std::ostream &out, const char *seq, unsigned int offset, const std::string &bt, bool reverse) {
    unsigned int seqPos = 0;
    for (uint32_t i = 0; i < bt.size(); ++i) {
        switch (bt[i]) {
            case 'M':
                out << seq[offset + seqPos];
                seqPos++;
                break;
            case 'I':
                if (reverse == true) {
                    out << '-';
                } else {
                    out << seq[offset + seqPos];
                    seqPos++;
                }
                break;
            case 'D':
                if (reverse == true) {
                    out << seq[offset + seqPos];
                    seqPos++;
                } else {
                    out << '-';
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
        std::string tHeaderName = (par.db2 + "_h");
        Debug(Debug::INFO) << "Target Header file: " << par.db2 << "_h\n";
        tHeaderDbr = new HeaderIdReader(par.db2.c_str(), par.noPreload);
    }

    Debug(Debug::INFO) << "Alignment database: " << par.db3 << "\n";
    DBReader<unsigned int> alnDbr(par.db3.c_str(), std::string(par.db3 + ".index").c_str());
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads);
    resultWriter.open();
    Debug(Debug::INFO) << "Start writing file to " << par.db4 << "\n";
    bool isDb = par.dbOut;

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < alnDbr.getSize(); i++) {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        unsigned int queryKey = alnDbr.getDbKey(i);
        char *data = alnDbr.getData(i);

        std::string querySeq;
        if (needSequenceDB) {
            querySeq = queryReader->getDataByDBKey(queryKey);
        }

        char buffer[1024];
        std::ostringstream ss;

        std::string queryId = qHeaderDbr.getId(queryKey);
        std::vector<Matcher::result_t> results = Matcher::readAlignmentResults(data, true);
        unsigned int missMatchCount;
        for (size_t j = 0; j < results.size(); j++) {
            const Matcher::result_t &res = results[j];
            std::string targetId = tHeaderDbr->getId(res.dbKey);
            unsigned int gapOpenCount = 0;
            unsigned int alnLen = res.alnLength;
            if(res.backtrace.size() > 0) {
                size_t matchCount = 0;

                std::vector<int> vec;
                alnLen = 0;
                for (size_t pos = 0; pos < res.backtrace.size(); pos++) {
                    int cnt=0;
                    if(isdigit(res.backtrace[pos])){
                        cnt += Util::fast_atoi<int>(res.backtrace.c_str()+pos);
                        while(isdigit(res.backtrace[pos])){
                            pos++;
                        }
                    }
                    alnLen+=cnt;

                    switch(res.backtrace[pos]){
                        case 'M':
                            matchCount+= cnt;
                            break;
                        case 'D':
                        case 'I':
                            gapOpenCount+=1;
                            break;
                    }
                }
//                res.seqId = X / alnLen;
                unsigned int identical = static_cast<unsigned int>( res.seqId * static_cast<float>(alnLen)  + 0.5 );
                missMatchCount = static_cast<unsigned int>( matchCount - identical);
            }else{
                int adjustQstart = (res.qStartPos==-1)? 0 : res.qStartPos;
                int adjustDBstart = (res.dbStartPos==-1)? 0 : res.dbStartPos;
                float bestMatchEstimate = static_cast<float>(std::min(res.qEndPos - adjustQstart,  res.dbEndPos-  adjustDBstart));
                missMatchCount = static_cast<unsigned int>( bestMatchEstimate *  (1.0f - res.seqId) + 0.5 );
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

                    ss << buffer;
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

                    ss << buffer;
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

                    ss << buffer;

                    const std::string &backtrace = Matcher::uncompressAlignment(res.backtrace);
                    printSeqBasedOnAln(ss, querySeq.c_str(), res.qStartPos, backtrace, false);
                    ss << '\n';

                    std::string targetSeq = targetReader->getDataByDBKey(res.dbKey);
                    printSeqBasedOnAln(ss, targetSeq.c_str(), res.dbStartPos, backtrace, true);
                    ss << '\n';
                    break;
                }
                case Parameters::FORMAT_ALIGNMENT_SAM:
                default:
                    Debug(Debug::ERROR) << "Not implemented yet";
                    EXIT(EXIT_FAILURE);
            }
        }

        std::string result = ss.str();
        if(isDb==false){
            resultWriter.writeData(result.c_str(), result.size(), queryKey, thread_idx, false);
        }else{
            resultWriter.writeData(result.c_str(), result.size(), queryKey, thread_idx);
        }
    }
    resultWriter.close();

    // remove NULL byte
    // \n\0 -> ' '\n
    if(isDb==false) {
        FileUtil::deleteFile(par.db4Index);
    }
    if (par.earlyExit) {
        Debug(Debug::INFO) << "Done. Exiting early now.\n";
        _Exit(EXIT_SUCCESS);
    }

    Debug(Debug::INFO) << "Done." << "\n";

    alnDbr.close();
    if(sameDB == false) {
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

