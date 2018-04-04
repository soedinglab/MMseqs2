//
// Created by mad on 11/28/17.
//
#include <string>
#include <vector>
#include <fstream>
#include <FileUtil.h>
#include <Orf.h>
#include <itoa.h>
#include <NucleotideMatrix.h>

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


int proteinaln2nucl(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4);

    Debug(Debug::INFO) << "Query  file: " << par.db1 << "\n";
    DBReader<unsigned int> * qdbr = new DBReader<unsigned int>(par.db1.c_str(), (par.db1 + ".index").c_str());
    qdbr->open(DBReader<unsigned int>::NOSORT);
    qdbr->readMmapedDataInMemory();

    DBReader<unsigned int> *tdbr = NULL;
//    BaseMatrix * subMat = new NucleotideMatrix(par.scoringMatrixFile.c_str(), 1.0, 0.0);

    Debug(Debug::INFO) << "Target  file: " << par.db2 << "\n";
    bool sameDB = false;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = qdbr;
    } else {
        tdbr = new DBReader<unsigned int>(par.db2.c_str(), (par.db2 + ".index").c_str());
        tdbr->open(DBReader<unsigned int>::NOSORT);
        tdbr->readMmapedDataInMemory();
    }


    Debug(Debug::INFO) << "Alignment database: " << par.db3 << "\n";
    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str());
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads);
    resultWriter.open();
    Debug(Debug::INFO) << "Start writing file to " << par.db4 << "\n";

#pragma omp parallel for schedule(dynamic, 10)
    for (size_t i = 0; i < alnDbr.getSize(); i++) {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        unsigned int alnKey = alnDbr.getDbKey(i);
        char *data = alnDbr.getData(i);

        char buffer[1024];
        std::string ss;
        ss.reserve(1024);

        unsigned int queryId = qdbr->getId(alnKey);
        char *querySeq = qdbr->getData(queryId);

        std::vector<Matcher::result_t> results = Matcher::readAlignmentResults(data, true);
        for (size_t j = 0; j < results.size(); j++) {
            Matcher::result_t &res = results[j];
            bool hasBacktrace = (res.backtrace.size() > 0);
            unsigned int targetId = tdbr->getId(results[j].dbKey);
            char * targetSeq = tdbr->getData(targetId);

            std::string newBacktrace;
            newBacktrace.reserve(1024);
            res.dbStartPos = res.dbStartPos*3;
            res.dbEndPos   = res.dbEndPos*3;
            res.dbLen      = res.dbLen*3;
            res.qStartPos  = res.qStartPos*3;
            res.qEndPos    = res.qEndPos*3;
            res.qLen       = res.qLen*3;
            size_t idCnt = 0;
            size_t alnLen = 0;

            int qPos=res.qStartPos;
            int tPos=res.dbStartPos;

            for (size_t pos = 0; pos < res.backtrace.size(); pos++) {
                int cnt=0;
                if(isdigit(res.backtrace[pos])){
                    cnt += Util::fast_atoi<int>(res.backtrace.c_str()+pos);
                    while(isdigit(res.backtrace[pos])){
                        pos++;
                    }
                }
                bool update = false;
                switch(res.backtrace[pos]){
                    case 'M':
                        for(int bt = 0; bt < cnt*3; bt++) {
                            idCnt += (querySeq[qPos] == targetSeq[tPos]);
                            tPos++;
                            qPos++;
                        }
                        update = true;
                        break;
                    case 'D':
                        for(int bt = 0; bt < cnt*3; bt++) {
                            tPos++;
                        }
                        update = true;
                        break;
                    case 'I':
                        for(int bt = 0; bt < cnt*3; bt++) {
                            qPos++;
                        }
                        update = true;
                        break;

                }
                if(update){
                    char *buffNext = Itoa::i32toa_sse2(cnt*3, buffer);
                    size_t len = buffNext - buffer;
                    alnLen += cnt*3;
                    newBacktrace.append(buffer, len-1);
                    newBacktrace.push_back(res.backtrace[pos]);
                }

            }
            res.backtrace = newBacktrace;
            res.seqId = static_cast<float>(idCnt)/ static_cast<float>(alnLen);
            // recompute alignment
            size_t len = Matcher::resultToBuffer(buffer, res, hasBacktrace, false);
            ss.append(buffer, len);
        }

        resultWriter.writeData(ss.c_str(), ss.length(), alnKey, thread_idx);
        ss.clear();
    }
    alnDbr.close();
    resultWriter.close();
    if (sameDB == false) {
        tdbr->close();
        delete tdbr;
    }
    Debug(Debug::INFO) << "Done." << "\n";

    return EXIT_SUCCESS;
}

