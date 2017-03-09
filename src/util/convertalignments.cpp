#include <string>
#include <vector>

#include "Alignment.h"
#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"

#include "Debug.h"
#include "DBReader.h"

class HeaderIdReader {
public:
    HeaderIdReader(const char* dataName, const char* indexName, bool noPreload) {
        reader = new DBReader<unsigned int>(dataName, indexName);
        reader->open(DBReader<unsigned int>::NOSORT);
        if(noPreload == false) {
            reader->readMmapedDataInMemory();
        }

        isLookup = false;
        if (Util::endsWith(".lookupdb", dataName)) {
            isLookup = true;
        }
    }

    std::string getId(unsigned int key) {
        size_t id = reader->getId(key);
        const char *data = reader->getData(id);
        if(isLookup) {
            return std::string(data, reader->getSeqLens(id) - 2);
        }

        return Util::parseFastaHeader(data);
    }

    ~HeaderIdReader() {
        reader->close();
        delete reader;
    }

private:
    DBReader<unsigned int> *reader;
    bool isLookup;
};


void printSeqBasedOnAln(FILE *out, const char *seq, unsigned int offset, const std::string &bt, bool reverse) {
    unsigned int seqPos = 0;
    for (uint32_t i = 0; i < bt.size(); ++i) {
        switch (bt[i]) {
            case 'M':
                fprintf(out, "%c", seq[offset + seqPos]);
                seqPos++;
                break;
            case 'I':
                if (reverse == true) {
                    fprintf(out, "-");
                } else {
                    fprintf(out, "%c", seq[offset + seqPos]);
                    seqPos++;
                }
                break;
            case 'D':
                if (reverse == true) {
                    fprintf(out, "%c", seq[offset + seqPos]);
                    seqPos++;
                } else {
                    fprintf(out, "-");
                }
                break;
        }
    }

}

int convertalignments(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4);

    IndexReader *queryReader = NULL;
    IndexReader *targetReader = NULL;

    bool sameDB = false;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
    }

    if (par.formatAlignmentMode != Parameters::FORMAT_ALIGNMENT_BLAST_TAB) {
        targetReader = new IndexReader(par.db2.c_str(), par.db2Index.c_str(), par.noPreload);
        if (sameDB == true) {
            queryReader = targetReader;
        } else {
            queryReader = new IndexReader(par.db1.c_str(), par.db1Index.c_str(), par.noPreload);
        }
    }
    std::string qHeaderName = (par.db1 + "_h");
    Debug(Debug::INFO) << "Query Header file: " << qHeaderName << "\n";
    HeaderIdReader qHeaderDbr(qHeaderName.c_str(), (qHeaderName + ".index").c_str(), par.noPreload);

    HeaderIdReader *tHeaderDbr;
    if(sameDB){
        tHeaderDbr = &qHeaderDbr;
    } else {
        std::string tHeaderName = (par.db2 + "_h");
        Debug(Debug::INFO) << "Target Header file: " << tHeaderName << "\n";
        tHeaderDbr = new HeaderIdReader(tHeaderName.c_str(), (tHeaderName + ".index").c_str(), par.noPreload);
    }

    Debug(Debug::INFO) << "Alignment database: " << par.db3 << "\n";
    DBReader<unsigned int> alnDbr(par.db3.c_str(), std::string(par.db3 + ".index").c_str());
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    FILE *fastaFP = fopen(par.db4.c_str(), "w");
    Debug(Debug::INFO) << "Start writing file to " << par.db4 << "\n";
    for (size_t i = 0; i < alnDbr.getSize(); i++) {
        unsigned int queryKey = alnDbr.getDbKey(i);
        char *data = alnDbr.getData(i);

        std::string querySeq;
        if (par.formatAlignmentMode != Parameters::FORMAT_ALIGNMENT_BLAST_TAB) {
            querySeq = queryReader->getSequenceString(queryKey);
        }


        std::string queryId = qHeaderDbr.getId(queryKey);
        std::vector<Matcher::result_t> results = Matcher::readAlignmentResults(data);
        for (size_t j = 0; j < results.size(); j++) {
            const Matcher::result_t &res = results[j];

            std::string targetId = tHeaderDbr->getId(res.dbKey);
            unsigned int missMatchCount = (unsigned int) (res.seqId * std::min(res.qLen, res.dbLen));
            unsigned int gapOpenCount = 0;
            if (par.formatAlignmentMode == Parameters::FORMAT_ALIGNMENT_BLAST_TAB) {
                fprintf(fastaFP, "%s\t%s\t%1.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2E\t%d\n",
                        queryId.c_str(), targetId.c_str(), res.seqId, res.alnLength, missMatchCount, gapOpenCount,
                        res.qStartPos + 1, res.qEndPos + 1, res.dbStartPos + 1, res.dbEndPos + 1, res.eval, res.score);
            } else if (par.formatAlignmentMode == Parameters::FORMAT_ALIGNMENT_PAIRWISE) {
                fprintf(fastaFP, ">%s\t%s\t%1.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2E\t%d\n",
                        queryId.c_str(), targetId.c_str(), res.seqId, res.alnLength, missMatchCount, gapOpenCount,
                        res.qStartPos + 1, res.qEndPos + 1, res.dbStartPos + 1, res.dbEndPos + 1, res.eval, res.score);

                const std::string &backtrace = res.backtrace;
                printSeqBasedOnAln(fastaFP, querySeq.c_str(), res.qStartPos, backtrace, false);
                fprintf(fastaFP, "\n");

                std::string targetSeq = targetReader->getSequenceString(res.dbKey);
                printSeqBasedOnAln(fastaFP, targetSeq.c_str(), res.dbStartPos, backtrace, true);
                fprintf(fastaFP, "\n");
            } else if (par.formatAlignmentMode == Parameters::FORMAT_ALIGNMENT_SAM) { ;
                //TODO
            }
        }
    }
    fclose(fastaFP);

    if (par.earlyExit) {
        Debug(Debug::INFO) << "Done. Exiting early now.\n";
        _Exit(EXIT_SUCCESS);
    }

    Debug(Debug::INFO) << "Done." << "\n";

    alnDbr.close();
    if(sameDB == false) {
        delete tHeaderDbr;
    }

    if (par.formatAlignmentMode != Parameters::FORMAT_ALIGNMENT_BLAST_TAB) {
        delete queryReader;
        if (sameDB == false) {
            delete targetReader;
        }
    }

    return EXIT_SUCCESS;
}


