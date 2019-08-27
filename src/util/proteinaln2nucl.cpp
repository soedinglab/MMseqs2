#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "itoa.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"

#ifdef OPENMP
#include <omp.h>
#endif


int proteinaln2nucl(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> *qdbr_nuc = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    qdbr_nuc->open(DBReader<unsigned int>::NOSORT);
    qdbr_nuc->readMmapedDataInMemory();

    DBReader<unsigned int> *qdbr_aa = new DBReader<unsigned int>(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    qdbr_aa->open(DBReader<unsigned int>::NOSORT);
    qdbr_aa->readMmapedDataInMemory();

    DBReader<unsigned int> *tdbr_nuc = NULL;
    DBReader<unsigned int> *tdbr_aa = NULL;
//    NucleotideMatrix subMat(par.scoringMatrixFile.nucleotides, 1.0, 0.0);

    bool sameDB = false;
    if (par.db1.compare(par.db2) == 0 && par.db3.compare(par.db4) == 0) {
        sameDB = true;
        tdbr_nuc = qdbr_nuc;
        tdbr_aa = qdbr_aa;

    } else if (par.db1.compare(par.db2) != 0 && par.db3.compare(par.db4) != 0) {
        tdbr_nuc = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        tdbr_nuc->open(DBReader<unsigned int>::NOSORT);
        tdbr_nuc->readMmapedDataInMemory();

        tdbr_aa =  new DBReader<unsigned int>(par.db4.c_str(), par.db4Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        tdbr_aa->open(DBReader<unsigned int>::NOSORT);
        tdbr_aa->readMmapedDataInMemory();
    }
    else{
        Debug(Debug::ERROR) << "Either query database == target database for nucleotide and amino acid or != for both .\n";
        EXIT(EXIT_FAILURE);
    }

    if(Parameters::isEqualDbtype(qdbr_nuc->getDbtype(), Parameters::DBTYPE_NUCLEOTIDES) == false ||
       Parameters::isEqualDbtype(tdbr_nuc->getDbtype(), Parameters::DBTYPE_NUCLEOTIDES) == false ||
       Parameters::isEqualDbtype(qdbr_aa->getDbtype(), Parameters::DBTYPE_AMINO_ACIDS) == false ||
       Parameters::isEqualDbtype(tdbr_aa->getDbtype(), Parameters::DBTYPE_AMINO_ACIDS) == false){
        Debug(Debug::ERROR) << "Wrong query and target database input.\n";
        EXIT(EXIT_FAILURE);
    }

    DBReader<unsigned int> alnDbr(par.db5.c_str(), par.db5Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter resultWriter(par.db6.c_str(), par.db6Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    resultWriter.open();
    Debug::Progress progress(alnDbr.getSize());

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        char buffer[1024];
        std::string ss;
        ss.reserve(1024);

        std::vector<Matcher::result_t> results;
        results.reserve(300);

        std::string newBacktrace;
        newBacktrace.reserve(1024);

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < alnDbr.getSize(); i++) {
            progress.updateProgress();

            unsigned int alnKey = alnDbr.getDbKey(i);
            char *data = alnDbr.getData(i, thread_idx);

            unsigned int queryId = qdbr_nuc->getId(alnKey);
            char *nuclQuerySeq = qdbr_nuc->getData(queryId, thread_idx);
            char *aaQuerySeq = qdbr_aa->getDataByDBKey(alnKey, thread_idx);
            unsigned int nuclQuerySeqLen = qdbr_nuc->getSeqLen(queryId);

            bool qStartCodon = false;
            if (aaQuerySeq[0] == '*' )
                qStartCodon = true;

            Matcher::readAlignmentResults(results, data, true);
            for (size_t j = 0; j < results.size(); j++) {
                Matcher::result_t &res = results[j];
                bool hasBacktrace = (res.backtrace.size() > 0);

                if(!hasBacktrace ){
                    Debug(Debug::ERROR) << "This module only supports database "\
                                           "input with backtrace string.\n";
                    EXIT(EXIT_FAILURE);
                }



                unsigned int targetId = tdbr_nuc->getId(results[j].dbKey);
                char *nuclTargetSeq = tdbr_nuc->getData(targetId, thread_idx);
                char *aaTargetSeq = tdbr_aa->getDataByDBKey(results[j].dbKey, thread_idx);
                unsigned int nuclTargetSeqLen = tdbr_nuc->getSeqLen(targetId);

                bool tStartCodon = false;
                if (aaTargetSeq[0] == '*' )
                    tStartCodon = true;

                if (tStartCodon && res.dbStartPos == 0){
                    Debug(Debug::ERROR) << "Alignment contains unalignable character.\n";
                    EXIT(EXIT_FAILURE);
                }
                if (qStartCodon && res.qStartPos == 0){
                    Debug(Debug::ERROR) << "Alignment contains unalignable character.\n";
                    EXIT(EXIT_FAILURE);
                }


                res.dbStartPos = res.dbStartPos*3 + (tStartCodon? -3:0);
                res.dbEndPos   = res.dbEndPos*3+2 + (tStartCodon? -3:0);
                res.dbLen      = nuclTargetSeqLen;
                res.qStartPos  = res.qStartPos*3 + (qStartCodon? -3:0);
                res.qEndPos    = res.qEndPos*3+2 + (qStartCodon? -3:0);
                res.qLen       = nuclQuerySeqLen;
                size_t idCnt = 0;
                size_t alnLen = 0;

                int qPos = res.qStartPos;
                int tPos = res.dbStartPos;

                for (size_t pos = 0; pos < res.backtrace.size(); pos++) {
                    int cnt =0;
                    if (isdigit(res.backtrace[pos])){
                        cnt += Util::fast_atoi<int>(res.backtrace.c_str()+pos);
                        while (isdigit(res.backtrace[pos])){
                            pos++;
                        }
                    }
                    bool update = false;
                    switch (res.backtrace[pos]) {
                        case 'M':
                            for (int bt = 0; bt < cnt*3; bt++) {
                                idCnt += (nuclQuerySeq[qPos] == nuclTargetSeq[tPos]);
                                tPos++;
                                qPos++;
                            }
                            update = true;
                            break;
                        case 'D':
                            for (int bt = 0; bt < cnt*3; bt++) {
                                tPos++;
                            }
                            update = true;
                            break;
                        case 'I':
                            for (int bt = 0; bt < cnt*3; bt++) {
                                qPos++;
                            }
                            update = true;
                            break;

                    }
                    if (update) {

                        alnLen += cnt*3;
                        newBacktrace.append(SSTR(cnt*3));
                        newBacktrace.push_back(res.backtrace[pos]);
                    }

                }
                res.backtrace = newBacktrace;
                res.seqId = static_cast<float>(idCnt)/ static_cast<float>(alnLen);
                // recompute alignment
                size_t len = Matcher::resultToBuffer(buffer, res, hasBacktrace, false);
                ss.append(buffer, len);
                newBacktrace.clear();
            }

            resultWriter.writeData(ss.c_str(), ss.length(), alnKey, thread_idx);
            ss.clear();
            results.clear();
        }
    }
    resultWriter.close();
    alnDbr.close();
    if (sameDB == false) {
        tdbr_nuc->close();
        tdbr_aa->close();
        delete tdbr_nuc;
        delete tdbr_aa;
    }
    return EXIT_SUCCESS;
}

