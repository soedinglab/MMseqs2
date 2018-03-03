//
// Created by mad on 11/28/17.
//
#include <string>
#include <vector>
#include <fstream>
#include <FileUtil.h>
#include <Orf.h>

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
    par.parseParameters(argc, argv, command, 2);

    Debug(Debug::INFO) << "Alignment database: " << par.db1 << "\n";
    DBReader<unsigned int> alnDbr(par.db1.c_str(), par.db1Index.c_str());
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter resultWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads);
    resultWriter.open();
    Debug(Debug::INFO) << "Start writing file to " << par.db2 << "\n";

#pragma omp parallel for schedule(dynamic, 10)
    for (size_t i = 0; i < alnDbr.getSize(); i++) {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        unsigned int alnKey = alnDbr.getDbKey(i);
        char *data = alnDbr.getData(i);

        std::string querySeq;

        char buffer[1024];
        std::string ss;
        ss.reserve(1024);


        std::vector<Matcher::result_t> results = Matcher::readAlignmentResults(data, true);
        for (size_t j = 0; j < results.size(); j++) {
            Matcher::result_t &res = results[j];
            bool hasBacktrace = (res.backtrace.size() > 0);
            std::string newBacktrace;
            newBacktrace.reserve(1024);
            res.dbStartPos = res.dbStartPos*3;
            res.dbEndPos   = res.dbEndPos*3;
            res.dbLen      = res.dbLen*3;
            res.qStartPos  = res.qStartPos*3;
            res.qEndPos    = res.qEndPos*3;
            res.qLen       = res.qLen*3;
            for (size_t pos = 0; pos < res.backtrace.size(); pos++) {
                int cnt=0;
                if(isdigit(res.backtrace[pos])){
                    cnt += Util::fast_atoi<int>(res.backtrace.c_str()+pos);
                    while(isdigit(res.backtrace[pos])){
                        pos++;
                    }
                }
                switch(res.backtrace[pos]){
                    case 'M':
                    case 'D':
                    case 'I':
                        newBacktrace.append(std::string(cnt*3, res.backtrace[pos]));
                        break;

                }

            }
            res.backtrace = newBacktrace;

            size_t len = Matcher::resultToBuffer(buffer, res, hasBacktrace, true);
            ss.append(buffer, len);
        }

        resultWriter.writeData(ss.c_str(), ss.length(), alnKey, thread_idx);
        ss.clear();
    }
    alnDbr.close();
    resultWriter.close();
    Debug(Debug::INFO) << "Done." << "\n";

    return EXIT_SUCCESS;
}

