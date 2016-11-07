//
// Created by mad on 10/21/15.
//

#include <string>
#include <vector>
#include <Alignment.h>
#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"

#include "Debug.h"
#include "DBReader.h"


void printSeqBasedOnAln(FILE * out, char *seq, unsigned int offset, std::string bt, bool reverse) {
    unsigned int seqPos = 0;
    for (uint32_t i = 0; i < bt.size(); ++i){
        switch(bt[i]){
            case 'M':
                fprintf(out, "%c", seq[offset + seqPos]);
                seqPos++;
                break;
            case 'I':
                if(reverse == true){
                    fprintf(out, "-");
                }else{
                    fprintf(out, "%c", seq[offset + seqPos]);
                    seqPos++;
                }
                break;
            case 'D':
                if(reverse == true){
                    fprintf(out, "%c", seq[offset + seqPos]);
                    seqPos++;
                }else{
                    fprintf(out, "-");
                }
                break;
        }
    }

}

int convertalignments(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4);

    DBReader<unsigned int> * query = NULL;
    DBReader<unsigned int> * target = NULL;
    bool sameDB = false;
    if(par.formatAlignmentMode != Parameters::FORMAT_ALIGNMENT_BLAST_TAB){
        Debug(Debug::WARNING) << "Query  file: " << par.db1 << "\n";
        query = new DBReader<unsigned int>( par.db1.c_str(), (par.db1+".index").c_str());
        query->open(DBReader<unsigned int>::NOSORT);
        query->readMmapedDataInMemory();
        Debug(Debug::WARNING) << "Target  file: " << par.db2  << "\n";
        if(par.db1.compare(par.db2) == 0){
            sameDB = true;
            target = query;
        }else{
            target = new DBReader<unsigned int>( par.db2.c_str(), (par.db2+".index").c_str());
            target->open(DBReader<unsigned int>::NOSORT);
            target->readMmapedDataInMemory();
        }
    }
    std::string qffindexHeaderDB = (par.db1 + "_h");
    Debug(Debug::WARNING) << "Query Header file: " << qffindexHeaderDB << "\n";
    DBReader<unsigned int> q_header( qffindexHeaderDB.c_str(), (qffindexHeaderDB+".index").c_str());
    q_header.open(DBReader<unsigned int>::NOSORT);
    q_header.readMmapedDataInMemory();
    std::string  dbffindexHeaderDB = (par.db2 + "_h");
    Debug(Debug::WARNING) << "Target Header file: " << dbffindexHeaderDB << "\n";
    DBReader<unsigned int> db_header( dbffindexHeaderDB.c_str(), (dbffindexHeaderDB+".index").c_str());
    db_header.open(DBReader<unsigned int>::NOSORT);
    db_header.readMmapedDataInMemory();
    Debug(Debug::WARNING) << "Alignment database: " << par.db3 << "\n";
    DBReader<unsigned int> dbr_aln(par.db3.c_str(), std::string(par.db3+".index").c_str());
    dbr_aln.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    FILE *fastaFP =  fopen(par.db4.c_str(), "w");

    Debug(Debug::WARNING) << "Start writing file to " << par.db4 << "\n";
    for(size_t i = 0; i < dbr_aln.getSize(); i++){
        unsigned int queryKey = dbr_aln.getDbKey(i);
        char * header = q_header.getDataByDBKey(queryKey);
        char * data = dbr_aln.getData(i);
        char * querySeq = NULL;
        if(par.formatAlignmentMode != Parameters::FORMAT_ALIGNMENT_BLAST_TAB){
            querySeq = query->getDataByDBKey(queryKey);
        }

        std::string queryId = Util::parseFastaHeader(header);
        std::vector<Matcher::result_t> results = Matcher::readAlignmentResults(data);
        for(size_t j = 0; j < results.size(); j++){
            Matcher::result_t res = results[j];
            char *headerLine = db_header.getDataByDBKey(res.dbKey);
            std::string targetId = Util::parseFastaHeader(headerLine);
            unsigned int missMatchCount = (unsigned int)(res.seqId * std::min(res.qLen, res.dbLen));
            unsigned int gapOpenCount = 0;
            if(par.formatAlignmentMode == Parameters::FORMAT_ALIGNMENT_BLAST_TAB){
                fprintf(fastaFP, "%s\t%s\t%1.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2E\t%d\n",
                        queryId.c_str(), targetId.c_str(),  res.seqId, res.alnLength, missMatchCount, gapOpenCount,
                        res.qStartPos + 1,res.qEndPos + 1, res.dbStartPos + 1, res.dbEndPos + 1, res.eval, res.score);
            } else if(par.formatAlignmentMode == Parameters::FORMAT_ALIGNMENT_PAIRWISE) {
                fprintf(fastaFP, ">%s\t%s\t%1.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2E\t%d\n",
                        queryId.c_str(), targetId.c_str(),  res.seqId, res.alnLength, missMatchCount, gapOpenCount,
                        res.qStartPos + 1,res.qEndPos + 1, res.dbStartPos + 1,res.dbEndPos + 1, res.eval, res.score);
                std::string backtrace = res.backtrace;
                printSeqBasedOnAln(fastaFP, querySeq, res.qStartPos, backtrace, false);
                fprintf(fastaFP, "\n");
                char * targetSeq = target->getDataByDBKey(res.dbKey);
                printSeqBasedOnAln(fastaFP, targetSeq, res.dbStartPos, backtrace, true);
                fprintf(fastaFP, "\n");
            } else if(par.formatAlignmentMode == Parameters::FORMAT_ALIGNMENT_SAM){
                ;
                //TODO
            }
        }
    }
    Debug(Debug::WARNING) << "Done." << "\n";
    fclose(fastaFP);

    dbr_aln.close();
    q_header.close();
    db_header.close();
    if(par.formatAlignmentMode != Parameters::FORMAT_ALIGNMENT_BLAST_TAB) {
        query->close();
        delete query;
        if(sameDB == false){
            target->close();
            delete target;
        }
    }
    return 0;
}


