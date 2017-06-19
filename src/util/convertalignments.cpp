#include <string>
#include <vector>
#include <fstream>

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

class OutputWriter {
public:
    OutputWriter(const std::string &outname, bool db, unsigned int threads = 1)
            : db(db), writer(NULL), fs(NULL) {
        if (db == true) {
            std::string index(outname);
            index.append(".index");

            writer = new DBWriter(outname.c_str(), index.c_str(), threads);
            writer->open();
        } else {
            fs = new std::ofstream(outname);
            if (fs->fail()) {
                Debug(Debug::ERROR) << "Cannot read " << outname << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
    }

    void write(const std::string &data, unsigned int key, unsigned int thread = 0) {
        if (db == true) {
            writer->writeData(data.c_str(), data.size(), key, thread);
        } else {
            (*fs) << data;
        }
    }

    ~OutputWriter() {
        if (db == true) {
            writer->close();
            delete writer;
        } else {
            fs->close();
            delete fs;
        }
    }

private:
    const bool db;
    DBWriter *writer;
    std::ofstream *fs;
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

    bool sameDB = false;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
    }

    if (par.formatAlignmentMode != Parameters::FORMAT_ALIGNMENT_BLAST_TAB) {
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

    Debug(Debug::INFO) << "Start writing file to " << par.db4 << "\n";
    bool isDb = par.dbOut;
    unsigned int threads = 1;
    if (isDb == true) {
        threads = static_cast<unsigned int>(par.threads);
    }

    OutputWriter *writer = new OutputWriter(par.db4, isDb, threads);
#pragma omp parallel for schedule(static) num_threads(threads)
    for (size_t i = 0; i < alnDbr.getSize(); i++) {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        unsigned int queryKey = alnDbr.getDbKey(i);
        char *data = alnDbr.getData(i);

        std::string querySeq;
        if (par.formatAlignmentMode != Parameters::FORMAT_ALIGNMENT_BLAST_TAB) {
            querySeq = queryReader->getDataByDBKey(queryKey);
        }

        char buffer[1024];
        std::ostringstream ss;

        std::string queryId = qHeaderDbr.getId(queryKey);
        std::vector<Matcher::result_t> results = Matcher::readAlignmentResults(data, true);
        for (size_t j = 0; j < results.size(); j++) {
            const Matcher::result_t &res = results[j];

            std::string targetId = tHeaderDbr->getId(res.dbKey);
            unsigned int missMatchCount = static_cast<unsigned int>((1.0f-res.seqId) * res.alnLength);
            unsigned int gapOpenCount = 0;

            if(res.backtrace.size() > 0){
                for(size_t pos = 0; pos < res.backtrace.size(); pos++){
                    gapOpenCount += (res.backtrace[pos]=='I'||res.backtrace[pos]=='D') ;
                }
            }

            if (par.formatAlignmentMode == Parameters::FORMAT_ALIGNMENT_BLAST_TAB) {
                int count = snprintf(buffer, sizeof(buffer), "%s\t%s\t%1.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2E\t%d\n",
                    queryId.c_str(), targetId.c_str(), res.seqId, res.alnLength, missMatchCount, gapOpenCount,
                    res.qStartPos + 1, res.qEndPos + 1, res.dbStartPos + 1, res.dbEndPos + 1, res.eval, res.score);

                if (count < 0 || count>=sizeof(buffer)) {
                    Debug(Debug::WARNING) << "Truncated line in entry" << j << "!\n";
                    continue;
                }

                ss << buffer;
            } else if (par.formatAlignmentMode == Parameters::FORMAT_ALIGNMENT_PAIRWISE) {
                int count = snprintf(buffer, sizeof(buffer), ">%s\t%s\t%1.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2E\t%d\n",
                        queryId.c_str(), targetId.c_str(), res.seqId, res.alnLength, missMatchCount, gapOpenCount,
                        res.qStartPos + 1, res.qEndPos + 1, res.dbStartPos + 1, res.dbEndPos + 1, res.eval, res.score);

                if (count < 0 || count>=sizeof(buffer)) {
                    Debug(Debug::WARNING) << "Truncated line in entry" << j << "!\n";
                    continue;
                }

                ss << buffer;

                const std::string &backtrace =  Matcher::uncompressAlignment(res.backtrace);
                printSeqBasedOnAln(ss, querySeq.c_str(), res.qStartPos, backtrace, false);
                ss << '\n';

                std::string targetSeq = targetReader->getDataByDBKey(res.dbKey);
                printSeqBasedOnAln(ss, targetSeq.c_str(), res.dbStartPos, backtrace, true);
                ss << '\n';

            } else if (par.formatAlignmentMode == Parameters::FORMAT_ALIGNMENT_SAM) { ;
                ;
                //TODO
            }
        }
        
        std::string result = ss.str();
        writer->write(result, queryKey, thread_idx);
    }

    delete writer;

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
        queryReader->close();
        delete queryReader;
        if (sameDB == false) {
            queryReader->close();
            delete targetReader;
        }
    }

    return EXIT_SUCCESS;
}

