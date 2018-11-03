#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "IndexReader.h"
#include "FileUtil.h"
#include "TranslateNucl.h"
#include "Sequence.h"
#include "Orf.h"

#ifdef OPENMP
#include <omp.h>
#endif

void printSeqBasedOnAln(std::string &out, const char *seq, unsigned int offset,
                        const std::string &bt, bool reverse, bool isReverseStrand,
                        bool translateSequence, const TranslateNucl &translateNucl) {
    unsigned int seqPos = 0;
    char codon[3];
    for (uint32_t i = 0; i < bt.size(); ++i) {
        char seqChar = (isReverseStrand == true) ? Orf::complement(seq[offset - seqPos]) : seq[offset + seqPos];
        if (translateSequence) {
            codon[0] = (isReverseStrand == true) ? Orf::complement(seq[offset - seqPos])     : seq[offset + seqPos];
            codon[1] = (isReverseStrand == true) ? Orf::complement(seq[offset - (seqPos+1)]) : seq[offset + (seqPos+1)];
            codon[2] = (isReverseStrand == true) ? Orf::complement(seq[offset - (seqPos+2)]) : seq[offset + (seqPos+2)];
            seqChar = translateNucl.translateSingleCodon(codon);
        }
        switch (bt[i]) {
            case 'M':
                out.append(1, seqChar);
                seqPos += (translateSequence) ?  3 : 1;
                break;
            case 'I':
                if (reverse == true) {
                    out.append(1, '-');
                } else {
                    out.append(1, seqChar);
                    seqPos += (translateSequence) ?  3 : 1;
                }
                break;
            case 'D':
                if (reverse == true) {
                    out.append(1, seqChar);
                    seqPos += (translateSequence) ?  3 : 1;
                } else {
                    out.append(1, '-');
                }
                break;
        }
    }
}

/*
query       Query sequence label
target      Target sequenc label
evalue      E-value
gapopen     Number of gap opens
pident      Percentage of identical matches
nident      Number of identical matches
qstart      1-based start position of alignment in query sequence
qend        1-based end position of alignment in query sequence
qlen        Query sequence length
tstart      1-based start position of alignment in target sequence
tend        1-based end position of alignment in target sequence
tlen        Target sequence length
alnlen      Number of alignment columns
raw         Raw alignment score
bits        Bit score
cigar       Alignment as string M=letter pair, D=delete (gap in query), I=insert (gap in target)
qseq        Full-length query sequence
tseq        Full-length target sequence
qheader     Header of Query sequence
theader     Header of Target sequence
qaln        Aligned query sequence with gaps
taln        Aligned target sequence with gaps
qframe      Query frame (-3 to +3)
tframe      Target frame (-3 to +3)
mismatch    Number of mismatches
qcov        Fraction of query sequence covered by alignment
tcov        Fraction of target sequence covered by alignment
 */

std::vector<int> getOutputFormat(std::string outformat, bool &needdatabase, bool &needbacktrace) {
    std::vector<std::string> outformatSplit = Util::split(outformat, ",");
    std::vector<int> formatCodes;
    int code = 0;
    for (size_t i = 0; i < outformatSplit.size(); ++i) {
        if(outformatSplit[i].compare("query") == 0){ code = Parameters::OUTFMT_QUERY;}
        else if (outformatSplit[i].compare("target") == 0){ code = Parameters::OUTFMT_TARGET;}
        else if (outformatSplit[i].compare("evalue") == 0){ code = Parameters::OUTFMT_EVALUE;}
        else if (outformatSplit[i].compare("gapopen") == 0){ code = Parameters::OUTFMT_GAPOPEN;}
        else if (outformatSplit[i].compare("pident") == 0){ code = Parameters::OUTFMT_PIDENT;}
        else if (outformatSplit[i].compare("nident") == 0){ code = Parameters::OUTFMT_NIDENT;}
        else if (outformatSplit[i].compare("qstart") == 0){ code = Parameters::OUTFMT_QSTART;}
        else if (outformatSplit[i].compare("qend") == 0){ code = Parameters::OUTFMT_QEND;}
        else if (outformatSplit[i].compare("qlen") == 0){ code = Parameters::OUTFMT_QLEN;}
        else if (outformatSplit[i].compare("tstart") == 0){ code = Parameters::OUTFMT_TSTART;}
        else if (outformatSplit[i].compare("tend") == 0){ code = Parameters::OUTFMT_TEND;}
        else if (outformatSplit[i].compare("tlen") == 0){ code = Parameters::OUTFMT_TLEN;}
        else if (outformatSplit[i].compare("alnlen") == 0){ code = Parameters::OUTFMT_ALNLEN;}
        else if (outformatSplit[i].compare("raw") == 0){ needdatabase = true; code = Parameters::OUTFMT_RAW;}
        else if (outformatSplit[i].compare("bits") == 0){ code = Parameters::OUTFMT_BITS;}
        else if (outformatSplit[i].compare("cigar") == 0){ needbacktrace = true; code = Parameters::OUTFMT_CIGAR;}
        else if (outformatSplit[i].compare("qseq") == 0){ needdatabase = true; code = Parameters::OUTFMT_QSEQ;}
        else if (outformatSplit[i].compare("tseq") == 0){ needdatabase = true; code = Parameters::OUTFMT_TSEQ;}
        else if (outformatSplit[i].compare("qheader") == 0){ code = Parameters::OUTFMT_QHEADER;}
        else if (outformatSplit[i].compare("theader") == 0){ code = Parameters::OUTFMT_THEADER;}
        else if (outformatSplit[i].compare("qaln") == 0){ needbacktrace = true; needdatabase = true; code = Parameters::OUTFMT_QALN;}
        else if (outformatSplit[i].compare("taln") == 0){ needbacktrace = true; needdatabase = true; code = Parameters::OUTFMT_TALN;}
        else if (outformatSplit[i].compare("qframe") == 0){ code = Parameters::OUTFMT_QFRAME;}
        else if (outformatSplit[i].compare("tframe") == 0){ code = Parameters::OUTFMT_TFRAME;}
        else if (outformatSplit[i].compare("mismatch") == 0){ code = Parameters::OUTFMT_MISMATCH;}
        else if (outformatSplit[i].compare("qcov") == 0){ code = Parameters::OUTFMT_QCOV;}
        else if (outformatSplit[i].compare("tcov") == 0){ code = Parameters::OUTFMT_TCOV;}
        else if (outformatSplit[i].compare("empty") == 0){ code = Parameters::OUTFMT_EMPTY;}
        else {
            Debug(Debug::ERROR) << "Format code " << outformatSplit[i] << " does not exist.";
            EXIT(EXIT_FAILURE);
        }
        formatCodes.push_back(code);
    }
    return formatCodes;
}

int convertalignments(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4);

    const bool sameDB = par.db1.compare(par.db2) == 0 ? true : false;
    const int format = par.formatAlignmentMode;
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);

    bool needSequenceDB = false;
    bool needbacktrace = false;
    const std::vector<int> outcodes = getOutputFormat(par.outfmt, needSequenceDB, needbacktrace);

    Debug(Debug::INFO) << "Query Header database: " << par.hdr1 << "\n";
    const int indexReaderMode = IndexReader::NEED_HEADERS | (needSequenceDB ? IndexReader::NEED_SEQUENCES : 0);
    IndexReader qDbr(par.db1.c_str(), indexReaderMode, touch);

    IndexReader *tDbr;
    if (sameDB) {
        tDbr = &qDbr;
    } else {
        Debug(Debug::INFO) << "Target database: " << par.db2 << "\n";
        tDbr = new IndexReader(par.db2.c_str(), indexReaderMode, touch);
    }

    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0f, -0.2f);
    EvalueComputation *evaluer = NULL;
    bool queryNucs = false;
    bool targetNucs = false;
    bool isTranslatedSearch = false;
    bool queryProfile = false;
    bool targetProfile = false;
    if (needSequenceDB) {
        queryNucs = qDbr.getDbtype() == Sequence::NUCLEOTIDES;
        targetNucs = tDbr->getDbtype() == Sequence::NUCLEOTIDES;
        if((targetNucs == true || queryNucs == true ) && !(queryNucs == true && targetNucs == true)){
            isTranslatedSearch = true;
        }
        queryProfile = qDbr.getDbtype() == Sequence::HMM_PROFILE;
        targetProfile = tDbr->getDbtype() == Sequence::HMM_PROFILE;
        evaluer = new EvalueComputation(tDbr->sequenceReader->getAminoAcidDBSize(), &subMat, par.gapOpen, par.gapExtend);
    }

    Debug(Debug::INFO) << "Alignment database: " << par.db3 << "\n";
    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str());
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

#ifdef OPENMP
    unsigned int totalThreads = par.threads;
#else
    unsigned int totalThreads = 1;
#endif

    unsigned int localThreads = totalThreads;
    if (alnDbr.getSize() <= totalThreads) {
        localThreads = alnDbr.getSize();
    }

    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), localThreads);
    resultWriter.open();
    Debug(Debug::INFO) << "Start writing to " << par.db4 << "\n";

    const bool isDb = par.dbOut;
    TranslateNucl translateNucl(static_cast<TranslateNucl::GenCode>(par.translationTable));

#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        char buffer[1024];

        std::string result;
        result.reserve(1024*1024);

        std::string queryProfData;
        queryProfData.reserve(1024);

        std::string targetProfData;
        targetProfData.reserve(1024);

#pragma omp  for schedule(dynamic, 10)
        for (size_t i = 0; i < alnDbr.getSize(); i++) {
            Debug::printProgress(i);

            const unsigned int queryKey = alnDbr.getDbKey(i);
            char *querySeqData = NULL;
            queryProfData.clear();
            if (needSequenceDB) {
                size_t qId = qDbr.sequenceReader->getId(queryKey);
                querySeqData = qDbr.sequenceReader->getData(qId);
                if (queryProfile) {
                    Sequence::extractProfileConsensus(querySeqData, subMat, queryProfData);
                }
            }

            size_t qHeaderId = qDbr.headerReader->getId(queryKey);
            const char *qHeader = qDbr.headerReader->getData(qHeaderId);
            size_t qHeaderLen = qDbr.headerReader->getSeqLens(qHeaderId);
            std::string queryId = Util::parseFastaHeader(qHeader);

            char *data = alnDbr.getData(i);
            while (*data != '\0') {
                Matcher::result_t res = Matcher::parseAlignmentRecord(data, true);
                data = Util::skipLine(data);

                if (res.backtrace.size() == 0 && needbacktrace == true) {
                    Debug(Debug::ERROR) << "Backtrace cigar is missing in the alignment result. Please recompute the alignment with the -a flag.\n"
                                           "Command: mmseqs align " << par.db1 << " " << par.db2 << " " << par.db3 << " " << "alnNew -a\n";
                    EXIT(EXIT_FAILURE);
                }

                size_t tHeaderId = tDbr->headerReader->getId(res.dbKey);
                const char *tHeader = tDbr->headerReader->getData(tHeaderId);
                size_t tHeaderLen = tDbr->headerReader->getSeqLens(tHeaderId);
                std::string targetId = Util::parseFastaHeader(tHeader);

                unsigned int gapOpenCount = 0;
                unsigned int alnLen = res.alnLength;
                unsigned int missMatchCount = 0;
                unsigned int identical = 0;
                if (res.backtrace.size() > 0) {
                    size_t matchCount = 0;
                    alnLen = 0;
                    for (size_t pos = 0; pos < res.backtrace.size(); pos++) {
                        int cnt = 0;
                        if (isdigit(res.backtrace[pos])) {
                            cnt += Util::fast_atoi<int>(res.backtrace.c_str() + pos);
                            while (isdigit(res.backtrace[pos])) {
                                pos++;
                            }
                        }
                        alnLen += cnt;

                        switch (res.backtrace[pos]) {
                            case 'M':
                                matchCount += cnt;
                                break;
                            case 'D':
                            case 'I':
                                gapOpenCount += 1;
                                break;
                        }
                    }
//                res.seqId = X / alnLen;
                    identical = static_cast<unsigned int>( res.seqId * static_cast<float>(alnLen) + 0.5 );
                    //res.alnLength = alnLen;
                    missMatchCount = static_cast<unsigned int>( matchCount - identical);
                } else {
                    const int adjustQstart = (res.qStartPos == -1) ? 0 : res.qStartPos;
                    const int adjustDBstart = (res.dbStartPos == -1) ? 0 : res.dbStartPos;
                    const float bestMatchEstimate = static_cast<float>(std::min(res.qEndPos - adjustQstart, res.dbEndPos - adjustDBstart));
                    missMatchCount = static_cast<unsigned int>(bestMatchEstimate * (1.0f - res.seqId) + 0.5);
                }

                switch (format) {
                    case Parameters::FORMAT_ALIGNMENT_BLAST_TAB: {
                        if (outcodes.empty()) {
                            int count = snprintf(buffer, sizeof(buffer),
                                                 "%s\t%s\t%1.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2E\t%d\n",
                                                 queryId.c_str(), targetId.c_str(), res.seqId, alnLen,
                                                 missMatchCount, gapOpenCount,
                                                 res.qStartPos + 1, res.qEndPos + 1,
                                                 res.dbStartPos + 1, res.dbEndPos + 1,
                                                 res.eval, res.score);
                            if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                                Debug(Debug::WARNING) << "Truncated line in entry" << i << "!\n";
                                continue;
                            }
                            result.append(buffer, count);
                        } else {
                            char *targetSeqData = NULL;
                            targetProfData.clear();
                            if (needSequenceDB) {
                                size_t tId = tDbr->sequenceReader->getId(res.dbKey);
                                targetSeqData = tDbr->sequenceReader->getData(tId);
                                if (targetProfile) {
                                    Sequence::extractProfileConsensus(targetSeqData, subMat, targetProfData);
                                }
                            }
                            for(size_t i = 0; i < outcodes.size(); i++) {
                                switch (outcodes[i]) {
                                    case Parameters::OUTFMT_QUERY:
                                        result.append(queryId);
                                        break;
                                    case Parameters::OUTFMT_TARGET:
                                        result.append(targetId);
                                        break;
                                    case Parameters::OUTFMT_EVALUE:
                                        result.append(SSTR(res.eval));
                                        break;
                                    case Parameters::OUTFMT_GAPOPEN:
                                        result.append(SSTR(gapOpenCount));
                                        break;
                                    case Parameters::OUTFMT_PIDENT:
                                        result.append(SSTR(res.seqId));
                                        break;
                                    case Parameters::OUTFMT_NIDENT:
                                        result.append(SSTR(identical));
                                        break;
                                    case Parameters::OUTFMT_QSTART:
                                        result.append(SSTR(res.qStartPos + 1));
                                        break;
                                    case Parameters::OUTFMT_QEND:
                                        result.append(SSTR(res.qEndPos + 1));
                                        break;
                                    case Parameters::OUTFMT_QLEN:
                                        result.append(SSTR(res.qLen));
                                        break;
                                    case Parameters::OUTFMT_TSTART:
                                        result.append(SSTR(res.dbStartPos + 1));
                                        break;
                                    case Parameters::OUTFMT_TEND:
                                        result.append(SSTR(res.dbEndPos + 1));
                                        break;
                                    case Parameters::OUTFMT_TLEN:
                                        result.append(SSTR(res.dbLen));
                                        break;
                                    case Parameters::OUTFMT_ALNLEN:
                                        result.append(SSTR(alnLen));
                                        break;
                                    case Parameters::OUTFMT_RAW:
                                        result.append(SSTR(static_cast<int>(evaluer->computeRawScoreFromBitScore(res.score) + 0.5)));
                                        break;
                                    case Parameters::OUTFMT_BITS:
                                        result.append(SSTR(res.score));
                                        break;
                                    case Parameters::OUTFMT_CIGAR:
                                        result.append(SSTR(res.backtrace));
                                        break;
                                    case Parameters::OUTFMT_QSEQ:
                                        if (queryProfile) {
                                            result.append(queryProfData.c_str(), res.qLen);
                                        } else {
                                            result.append(querySeqData, res.qLen);
                                        }
                                        break;
                                    case Parameters::OUTFMT_TSEQ:
                                        if (targetProfile) {
                                            result.append(targetProfData.c_str(), res.dbLen);
                                        } else {
                                            result.append(targetSeqData, res.dbLen);
                                        }
                                        break;
                                    case Parameters::OUTFMT_QHEADER:
                                        result.append(qHeader, qHeaderLen-2);
                                        break;
                                    case Parameters::OUTFMT_THEADER:
                                        result.append(tHeader, tHeaderLen-2);
                                        break;
                                    case Parameters::OUTFMT_QALN:
                                        if (queryProfile) {
                                            printSeqBasedOnAln(result, queryProfData.c_str(), res.qStartPos,
                                                               Matcher::uncompressAlignment(res.backtrace), false, (res.qStartPos > res.qEndPos),
                                                               (isTranslatedSearch == true && queryNucs == true), translateNucl);
                                        } else {
                                            printSeqBasedOnAln(result, querySeqData, res.qStartPos,
                                                               Matcher::uncompressAlignment(res.backtrace), false, (res.qStartPos > res.qEndPos),
                                                               (isTranslatedSearch == true && queryNucs == true), translateNucl);
                                        }
                                        break;
                                    case Parameters::OUTFMT_TALN: {
                                        if (targetProfile) {
                                            printSeqBasedOnAln(result, targetProfData.c_str(), res.dbStartPos,
                                                               Matcher::uncompressAlignment(res.backtrace), true,
                                                               (res.dbStartPos > res.dbEndPos),
                                                               (isTranslatedSearch == true && targetNucs == true), translateNucl);
                                        } else {
                                            printSeqBasedOnAln(result, targetSeqData, res.dbStartPos,
                                                               Matcher::uncompressAlignment(res.backtrace), true,
                                                               (res.dbStartPos > res.dbEndPos),
                                                               (isTranslatedSearch == true && targetNucs == true), translateNucl);
                                        }
                                        break;
                                    }
                                    case Parameters::OUTFMT_MISMATCH:
                                        result.append(SSTR(missMatchCount));
                                        break;
                                    case Parameters::OUTFMT_QCOV:
                                        result.append(SSTR(res.qcov));
                                        break;
                                    case Parameters::OUTFMT_TCOV:
                                        result.append(SSTR(res.dbcov));
                                        break;
                                    case Parameters::OUTFMT_EMPTY:
                                        result.push_back('-');
                                        break;
                                }
                                if (i < outcodes.size() - 1) {
                                    result.push_back('\t');
                                }
                            }
                            result.push_back('\n');
                        }
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
                            Debug(Debug::WARNING) << "Truncated line in entry" << i << "!\n";
                            continue;
                        }

                        result.append(buffer, count);
                        break;
                    }
                    default:
                        Debug(Debug::ERROR) << "Not implemented yet";
                        EXIT(EXIT_FAILURE);
                }
            }

            resultWriter.writeData(result.c_str(), result.size(), queryKey, thread_idx, isDb);
            result.clear();
        }
    }
    resultWriter.close();

    // tsv output
    if (isDb == false) {
        FileUtil::deleteFile(par.db4Index);
    }

    alnDbr.close();
    if (needSequenceDB) {
        if (sameDB == false) {
            delete tDbr;
        }
        delete evaluer;
    }

    return EXIT_SUCCESS;
}
