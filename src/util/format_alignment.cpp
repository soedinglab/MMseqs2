//
// Created by mad on 10/21/15.
//

#include <string>
#include <vector>
#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"

#include "Debug.h"
#include "DBReader.h"


int formatalignment (int argc, const char * argv[])
{

    std::string alignmentDB = "";
    std::string outFile = "";

    std::string usage;
    usage.append("Convert a ffindex alignment database to BLAST tab or SAM flat file.\n");
    usage.append("USAGE: <queryDb> <targetDb> <alignmentDB> <outFile>\n");
    usage.append("\nDesigned and implemented by Martin Steinegger <martin.steinegger@campus.lmu.de>.\n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.formatalignment, 4);


    std::string qffindexHeaderDB = (par.db1 + "_h");
    Debug(Debug::WARNING) << "Query Header file: " << qffindexHeaderDB << "\n";
    DBReader q_header( qffindexHeaderDB.c_str(), (qffindexHeaderDB+".index").c_str());
    q_header.open(DBReader::NOSORT);

    std::string  dbffindexHeaderDB = (par.db2 + "_h");
    Debug(Debug::WARNING) << "Target Header file: " << dbffindexHeaderDB << "\n";
    DBReader db_header( dbffindexHeaderDB.c_str(), (dbffindexHeaderDB+".index").c_str());
    db_header.open(DBReader::NOSORT);


    Debug(Debug::WARNING) << "Alignment database: " << par.db3 << "\n";
    DBReader dbr_aln(par.db3.c_str(), std::string(par.db3+".index").c_str());
    dbr_aln.open(DBReader::NOSORT);

    FILE *fastaFP =  fopen(par.db4.c_str(), "w");

    Debug(Debug::WARNING) << "Start writing file to " << par.db4 << "\n";
    for(size_t i = 0; i < dbr_aln.getSize(); i++){
        std::string queryKey = dbr_aln.getDbKey(i).c_str();
        char * header = q_header.getDataByDBKey(queryKey.c_str());
        char * data = dbr_aln.getData(i);
        std::string queryId = Util::parseFastaHeader(header);
        std::vector<Matcher::result_t> results = Matcher::readAlignmentResults(data);
        for(size_t j = 0; j < results.size(); j++){
            Matcher::result_t res = results[j];
            char *headerLine = db_header.getDataByDBKey(res.dbKey.c_str());
            std::string targetId = Util::parseFastaHeader(headerLine);
            size_t missMatchCount = res.seqId * std::min(res.qLen, res.dbLen);
            size_t gapOpenCount = 0;
            fprintf(fastaFP, "%s\t%s\t%1.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2E\t%d\n",
                    queryId.c_str(), targetId.c_str(),  res.seqId, res.alnLength, missMatchCount, gapOpenCount,
                    res.qStartPos,res.qEndPos, res.dbStartPos,res.dbEndPos, res.eval, res.score);
//            std::cout << buffer;
//            fastaFP.wri
//            std::cout << queryId << "\t";
//            std::cout << targetId << "\t";
//            std::cout << res.seqId << "\t";
//            std::cout << res.alnLength << "\t";
//            std::cout << missMatchCount << "\t";
//            std::cout << gapOpenCount << "\t";
//            std::cout << res.qStartPos << "\t";
//            std::cout << res.qEndPos << "\t";
//            std::cout << res.dbStartPos << "\t";
//            std::cout << res.dbEndPos << "\t";
//            std::cout << res.eval << "\t";
//            std::cout << res.score << "\t";
            //for(size_t i = 1; i < elmCount; i++){
            //    std::cout << strtod(entry[i],NULL) <<" ";
            //}
            //std::cout << std::endl;
        }

    }
    Debug(Debug::WARNING) << "Done." << "\n";

    fclose(fastaFP);

    dbr_aln.close();
    q_header.close();
    db_header.close();

    return 0;
}
