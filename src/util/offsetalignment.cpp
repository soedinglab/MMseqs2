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


int offsetalignment(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4);

    int queryDbType = DBReader<unsigned int>::parseDbType(par.db1.c_str());
    int targetDbType = DBReader<unsigned int>::parseDbType(par.db2.c_str());
    if(queryDbType == -1 || targetDbType == -1){
        Debug(Debug::ERROR) << "Please recreate your database or add a .dbtype file to your sequence/profile database.\n";
        EXIT(EXIT_FAILURE);
    }
    if(queryDbType == DBReader<unsigned int>::DBTYPE_NUC && targetDbType == DBReader<unsigned int>::DBTYPE_NUC ){
        Debug(Debug::ERROR) << "It is not supported that both dbs are nucleotide databases.\n";
        EXIT(EXIT_FAILURE);
    }

    Debug(Debug::INFO) << "Query Header file: " << par.db1 << "_h\n";
    std::string qHeaderName = (par.db1 + "_h");
    DBReader<unsigned int> qHeaderDbr(qHeaderName.c_str(), std::string(qHeaderName + ".index").c_str());
    qHeaderDbr.open(DBReader<unsigned int>::NOSORT);
    std::string tHeaderName = (par.db2 + "_h");
    Debug(Debug::INFO) << "Target Header file: " << par.db2 << "_h\n";
    DBReader<unsigned int> tHeaderDbr(tHeaderName.c_str(), std::string(tHeaderName + ".index").c_str());
    tHeaderDbr.open(DBReader<unsigned int>::NOSORT);

    Debug(Debug::INFO) << "Alignment database: " << par.db3 << "\n";
    DBReader<unsigned int> alnDbr(par.db3.c_str(), std::string(par.db3 + ".index").c_str());
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads);
    resultWriter.open();
    Debug(Debug::INFO) << "Start writing file to " << par.db4 << "\n";

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < alnDbr.getSize(); i++) {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        unsigned int queryKey = alnDbr.getDbKey(i);
        char *data = alnDbr.getData(i);

        std::string querySeq;


        char buffer[1024];

        size_t queryId = qHeaderDbr.getId(queryKey);
        Orf::SequenceLocation loc;
        if(queryDbType == DBReader<unsigned int>::DBTYPE_NUC){
            char * qheader = qHeaderDbr.getData(queryId);
            loc = Orf::parseOrfHeader(qheader);
        }
        std::string ss;
        ss.reserve(1024);

        std::vector<Matcher::result_t> results = Matcher::readAlignmentResults(data, true);
        for (size_t j = 0; j < results.size(); j++) {
            Matcher::result_t &res = results[j];
            size_t targetId = tHeaderDbr.getId(res.dbKey);
            bool hasBacktrace = (res.backtrace.size() > 0);
            if(targetDbType == DBReader<unsigned int>::DBTYPE_NUC){
                char * theader = tHeaderDbr.getData(targetId);
                loc = Orf::parseOrfHeader(theader);
                res.dbStartPos = loc.from + res.dbStartPos*3;
                res.dbEndPos   = loc.from + res.dbEndPos*3;
                res.dbLen      = res.dbLen*3;
            }else{
                res.qStartPos = loc.from + res.qStartPos*3;
                res.qEndPos = loc.from + res.qEndPos*3;
                res.qLen = res.qLen*3;
            }
            size_t len = Matcher::resultToBuffer(buffer, res, hasBacktrace, false);
            ss.append(buffer, len);
        }

        resultWriter.writeData(ss.c_str(), ss.length(), queryKey, thread_idx);
        ss.clear();
    }
    qHeaderDbr.close();
    tHeaderDbr.close();
    alnDbr.close();
    resultWriter.close();
    Debug(Debug::INFO) << "Done." << "\n";

    return EXIT_SUCCESS;
}

