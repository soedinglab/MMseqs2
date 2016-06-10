//
// Created by mad on 6/10/16.
//
#include <sys/time.h>
#include "DBWriter.h"
#include "Debug.h"
#include "Parameters.h"
#include "DBReader.h"
#include "Util.h"
#include "Matcher.h"

#ifdef OPENMP
#include <omp.h>
#endif


int doExtractalignedregion(Parameters & par) {
    DBReader<unsigned int> * qdbr = NULL;
    DBReader<unsigned int> * tdbr = NULL;
    bool sameDB = false;
    Debug(Debug::WARNING) << "Query  file: " << par.db1 << "\n";
    qdbr = new DBReader<unsigned int>( par.db1.c_str(), (par.db1+".index").c_str());
    qdbr->open(DBReader<unsigned int>::NOSORT);
    qdbr->readMmapedDataInMemory();
    Debug(Debug::WARNING) << "Target  file: " << par.db2  << "\n";
    if(par.db1.compare(par.db2) == 0){
        sameDB = true;
        tdbr = qdbr;
    }else{
        tdbr = new DBReader<unsigned int>( par.db2.c_str(), (par.db2+".index").c_str());
        tdbr->open(DBReader<unsigned int>::NOSORT);
        tdbr->readMmapedDataInMemory();
    }
    std::string qffindexHeaderDB = (par.db1 + "_h");
    Debug(Debug::WARNING) << "Query Header file: " << qffindexHeaderDB << "\n";
    std::string  dbffindexHeaderDB = (par.db2 + "_h");
    Debug(Debug::WARNING) << "Target Header file: " << dbffindexHeaderDB << "\n";
    Debug(Debug::WARNING) << "Alignment database: " << par.db3 << "\n";
    DBReader<unsigned int> alndbr(par.db3.c_str(), std::string(par.db3+".index").c_str());
    alndbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), par.threads);
    dbw.open();
    Debug(Debug::WARNING) << "Start writing file to " << par.db4 << "\n";
#pragma omp for schedule(dynamic, 1000)
    for(size_t i = 0; i < alndbr.getSize(); i++){
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        unsigned int queryKey = alndbr.getDbKey(i);
        char * data = alndbr.getData(i);
        char * qseq = NULL;
        if(par.extractMode == Parameters::EXTRACT_QUERY){
            qseq = qdbr->getDataByDBKey(queryKey);
        }

        std::vector<Matcher::result_t> results = Matcher::readAlignmentResults(data);
        for(size_t j = 0; j < results.size(); j++){
            Matcher::result_t res = results[j];
            char *tseq = tdbr->getDataByDBKey(res.dbKey);
            if(par.extractMode == Parameters::EXTRACT_QUERY) {
                std::string data = std::string(qseq + res.qStartPos, res.qEndPos - res.qStartPos);
                data.append("\n");
                dbw.write(data.c_str(), data.size(), SSTR( res.dbKey ).c_str(), thread_idx);
            }else if(par.extractMode == Parameters::EXTRACT_TARGET) {
                std::string data = std::string(tseq + res.dbStartPos, res.dbEndPos - res.dbStartPos);
                data.append("\n");
                dbw.write(data.c_str(), data.size(), SSTR( res.dbKey ).c_str(), thread_idx);
            }
        }
    }
    Debug(Debug::WARNING) << "Done." << "\n";
    dbw.close();
    alndbr.close();
    qdbr->close();
    delete qdbr;
    if(sameDB == false){
        tdbr->close();
        delete tdbr;
    }
    return 0;
}

int extractalignedregion(int argc, const char **argv) {

    std::string usage("Extract aligned regions from a alignment result.\n");
    usage.append("USAGE: <queryDB> <targetDB> <resultDB> <outDB>\n");
    usage.append("\nDesigned and implemented by Martin Steinegger <martin.steinegger@mpibpc.mpg.de>\n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.extractalignedregion, 4);

// never allow deletions
    par.allowDeletion = false;
    Debug(Debug::WARNING) << "Compute profile.\n";
    struct timeval start, end;
    gettimeofday(&start, NULL);

    int retCode = doExtractalignedregion(par);

    gettimeofday(&end, NULL);
    time_t sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for processing: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " <<
    (sec % 60) << "s\n";
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return retCode;
}

