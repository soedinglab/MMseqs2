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
#include "MemoryMapped.h"
#include "NcbiTaxonomy.h"
#include "MappingReader.h"
#include "ParquetDBWriter.h"

#define ZSTD_STATIC_LINKING_ONLY
#include <zstd.h>
#include "result_viz_prelude.html.zst.h"

#include <map>
#include <set>

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
qset        Query set
tset        Target set
 */

std::map<DBKeyType, DBKeyType> readKeyToSet(const std::string& file) {
    std::map<DBKeyType, DBKeyType> mapping;
    if (file.length() == 0) {
        return mapping;
    }

    MemoryMapped lookup(file, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    char* data = (char *) lookup.getData();
    char* end = data + lookup.mappedSize();
    const char* entry[255];
    while (data < end && *data != '\0') {
        const size_t columns = Util::getWordsOfLine(data, entry, 255);
        if (columns < 3) {
            Debug(Debug::WARNING) << "Not enough columns in lookup file " << file << "\n";
            continue;
        }
        mapping.emplace(Util::fast_atoi<DBKeyType>(entry[0]), Util::fast_atoi<DBKeyType>(entry[2]));
        data = Util::skipLine(data);
    }
    lookup.close();
    return mapping;
}


std::map<DBKeyType, std::string> readSetToSource(const std::string& file) {
    std::map<DBKeyType, std::string> mapping;
    if (file.length() == 0) {
        return mapping;
    }

    MemoryMapped source(file, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    char* data = (char *) source.getData();
    char* end = data + source.mappedSize();
    const char* entry[255];
    while (data < end && *data != '\0') {
        const size_t columns = Util::getWordsOfLine(data, entry, 255);
        if (columns < 2) {
            Debug(Debug::WARNING) << "Not enough columns in lookup file " << file << "\n";
            continue;
        }
        data = Util::skipLine(data);
        std::string source(entry[1], data - entry[1] - 1);
        mapping.emplace(Util::fast_atoi<DBKeyType>(entry[0]), source);
    }
    source.close();
    return mapping;
}

int convertalignments(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    const bool sameDB = par.db1.compare(par.db2) == 0 ? true : false;
    int format = par.formatAlignmentMode;
    bool addColumnHeaders = false;
    bool isParquet = false;
    if (format == Parameters::FORMAT_ALIGNMENT_BLAST_TAB_AS_PARQUET){
        isParquet = true;
        format = Parameters::FORMAT_ALIGNMENT_BLAST_TAB;
    }

    if (format == Parameters::FORMAT_ALIGNMENT_BLAST_TAB_WITH_HEADERS) {
        format = Parameters::FORMAT_ALIGNMENT_BLAST_TAB;
        addColumnHeaders = true;
    }

    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);

    bool needSequenceDB = false;
    bool needBacktrace = false;
    bool needFullHeaders = false;
    bool needLookup = false;
    bool needSource = false;
    bool needTaxonomy = false;
    bool needTaxonomyMapping = false;
    const std::vector<int> outcodes = Parameters::getOutputFormat(format, par.outfmt, needSequenceDB, needBacktrace, needFullHeaders,
                                                                  needLookup, needSource, needTaxonomyMapping, needTaxonomy);

    NcbiTaxonomy* t = NULL;
    if(needTaxonomy){
        std::string db2NoIndexName = PrefilteringIndexReader::dbPathWithoutIndex(par.db2);
        t = NcbiTaxonomy::openTaxonomy(db2NoIndexName);
    }
    MappingReader* mapping = NULL;
    if (needTaxonomy || needTaxonomyMapping) {
        std::string db2NoIndexName = PrefilteringIndexReader::dbPathWithoutIndex(par.db2);
        mapping = new MappingReader(db2NoIndexName);
    }

    bool isTranslatedSearch = false;

    int dbaccessMode = needSequenceDB ? (DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA) : (DBReader<DBKeyType>::USE_INDEX);

    std::map<DBKeyType, DBKeyType> qKeyToSet;
    std::map<DBKeyType, DBKeyType> tKeyToSet;
    if (needLookup) {
        std::string file1 = par.db1 + ".lookup";
        std::string file2 = par.db2 + ".lookup";
        qKeyToSet = readKeyToSet(file1);
        tKeyToSet = readKeyToSet(file2);
    }

    std::map<DBKeyType, std::string> qSetToSource;
    std::map<DBKeyType, std::string> tSetToSource;
    if (needSource) {
        std::string file1 = par.db1 + ".source";
        std::string file2 = par.db2 + ".source";
        qSetToSource = readSetToSource(file1);
        tSetToSource = readSetToSource(file2);
    }

    IndexReader qDbr(par.db1, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
    IndexReader qDbrHeader(par.db1, par.threads, IndexReader::SRC_HEADERS , (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);

    IndexReader *tDbr;
    IndexReader *tDbrHeader;
    if (sameDB) {
        tDbr = &qDbr;
        tDbrHeader= &qDbrHeader;
    } else {
        tDbr = new IndexReader(par.db2, par.threads, IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
        tDbrHeader = new IndexReader(par.db2, par.threads, IndexReader::SRC_HEADERS, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    }

    bool queryNucs = Parameters::isEqualDbtype(qDbr.sequenceReader->getDbtype(), Parameters::DBTYPE_NUCLEOTIDES);
    bool targetNucs = Parameters::isEqualDbtype(tDbr->sequenceReader->getDbtype(), Parameters::DBTYPE_NUCLEOTIDES);
    if (needSequenceDB) {
        // try to figure out if search was translated. This is can not be solved perfectly.
        bool seqtargetAA = false;
        if(Parameters::isEqualDbtype(tDbr->getDbtype(), Parameters::DBTYPE_INDEX_DB)){
            IndexReader tseqDbr(par.db2, par.threads, IndexReader::SEQUENCES, 0, IndexReader::PRELOAD_INDEX);
            seqtargetAA = Parameters::isEqualDbtype(tseqDbr.sequenceReader->getDbtype(), Parameters::DBTYPE_AMINO_ACIDS);
        } else if(targetNucs == true && queryNucs == true && par.searchType == Parameters::SEARCH_TYPE_AUTO){
            Debug(Debug::WARNING) << "It is unclear from the input if a translated or nucleotide search was performed\n "
                                     "Please provide the parameter --search-type 2 (translated) or 3 (nucleotide)\n";
            EXIT(EXIT_FAILURE);
        } else if(par.searchType == Parameters::SEARCH_TYPE_TRANSLATED){
            seqtargetAA = true;
        }

        if((targetNucs == true && queryNucs == false )  || (targetNucs == false && queryNucs == true ) || (targetNucs == true && seqtargetAA == true && queryNucs == true )  ){
            isTranslatedSearch = true;
        }
    }

    int gapOpen, gapExtend;
    SubstitutionMatrix * subMat= NULL;
    if (targetNucs == true && queryNucs == true && isTranslatedSearch == false) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
        gapOpen = par.gapOpen.values.nucleotide();
        gapExtend =  par.gapExtend.values.nucleotide();
    }else{
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);
        gapOpen = par.gapOpen.values.aminoacid();
        gapExtend = par.gapExtend.values.aminoacid();
    }
    EvalueComputation *evaluer = NULL;
    bool queryProfile = false;
    bool targetProfile = false;
    if (needSequenceDB) {
        queryProfile = Parameters::isEqualDbtype(qDbr.sequenceReader->getDbtype(), Parameters::DBTYPE_HMM_PROFILE);
        targetProfile = Parameters::isEqualDbtype(tDbr->sequenceReader->getDbtype(), Parameters::DBTYPE_HMM_PROFILE);
        evaluer = new EvalueComputation(tDbr->sequenceReader->getAminoAcidDBSize(), subMat, gapOpen, gapExtend);
    }

    DBReader<DBKeyType> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<DBKeyType>::USE_INDEX|DBReader<DBKeyType>::USE_DATA);
    alnDbr.open(isParquet ? DBReader<DBKeyType>::SORT_BY_ID : DBReader<DBKeyType>::LINEAR_ACCCESS);

    size_t localThreads = 1;
#ifdef OPENMP
    localThreads = std::max(std::min((size_t)par.threads, alnDbr.getSize()), (size_t)1);
#endif

    const bool shouldCompress = par.dbOut == true && par.compressed == true;
    const int dbType = par.dbOut == true ? Parameters::DBTYPE_GENERIC_DB : Parameters::DBTYPE_OMIT_FILE;
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), localThreads, shouldCompress, dbType);
    if (!isParquet){
        resultWriter.open();
    }

    const bool isDb = par.dbOut;
    TranslateNucl translateNucl(static_cast<TranslateNucl::GenCode>(par.translationTable));



    if (format == Parameters::FORMAT_ALIGNMENT_SAM) {
        char buffer[1024];
        std::set<DBKeyType> headerWritten;
        resultWriter.writeStart(0);
        std::string header = "@HD\tVN:1.4\tSO:queryname\n";
        resultWriter.writeAdd(header.c_str(), header.size(), 0);

        for (size_t i = 0; i < alnDbr.getSize(); i++) {
            char *data = alnDbr.getData(i, 0);
            while (*data != '\0') {
                char dbKeyBuffer[255 + 1];
                Util::parseKey(data, dbKeyBuffer);
                const DBKeyType dbKey = Util::fast_atoi<DBKeyType>(dbKeyBuffer);
                if (headerWritten.insert(dbKey).second) {
                    size_t tId = tDbr->sequenceReader->getId(dbKey);
                    if (tId == DB_ENTRY_NOT_FOUND) {
                        Debug(Debug::ERROR) << "Invalid target key " << dbKey << " in alignment entry " << alnDbr.getDbKey(i) << ".\n";
                        EXIT(EXIT_FAILURE);
                    }
                    unsigned int seqLen = tDbr->sequenceReader->getSeqLen(tId);
                    size_t tHeaderId = tDbrHeader->sequenceReader->getId(dbKey);
                    if (tHeaderId == DB_ENTRY_NOT_FOUND) {
                        Debug(Debug::ERROR) << "Invalid target header key " << dbKey << " in alignment entry " << alnDbr.getDbKey(i) << ".\n";
                        EXIT(EXIT_FAILURE);
                    }
                    const char *tHeader = tDbrHeader->sequenceReader->getData(tHeaderId, 0);
                    std::string targetId = Util::parseFastaHeader(tHeader);
                    int count = snprintf(buffer, sizeof(buffer), "@SQ\tSN:%s\tLN:%d\n", targetId.c_str(),
                                         (int32_t) seqLen);
                    if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                        Debug(Debug::WARNING) << "Truncated line in header " << i << "!\n";
                        continue;
                    }
                    resultWriter.writeAdd(buffer, count, 0);
                }
                resultWriter.writeEnd(0, 0, false, 0);
                data = Util::skipLine(data);
            }
        }
    } else if (format == Parameters::FORMAT_ALIGNMENT_HTML) {
        size_t dstSize = ZSTD_findDecompressedSize(result_viz_prelude_html_zst, result_viz_prelude_html_zst_len);
        char* dst = (char*)malloc(sizeof(char) * dstSize);
        size_t realSize = ZSTD_decompress(dst, dstSize, result_viz_prelude_html_zst, result_viz_prelude_html_zst_len);
        resultWriter.writeData(dst, realSize, 0, 0, false, false);
        const char* scriptBlock = "<script>render([";
        resultWriter.writeData(scriptBlock, strlen(scriptBlock), 0, 0, false, false);
        free(dst);
    } else if (addColumnHeaders == true && outcodes.empty() == false) {
        std::vector<std::string> outfmt = Util::split(par.outfmt, ",");
        std::string header(outfmt[0]);
        for(size_t i = 1; i < outfmt.size(); i++) {
            header.append(1, '\t');
            header.append(outfmt[i]);
        }
        header.append(1, '\n');
        resultWriter.writeData(header.c_str(), header.length(), 0, 0, false, false);
    }

    
    Debug::Progress progress(alnDbr.getSize());
    // per-thread buffers of processEntry
    struct ThreadBuffers {
        char buffer[1024];
        std::string result;
        std::string queryProfData;
        std::string queryHeaderBuffer;
        std::string targetProfData;
        std::string newBacktrace;
        std::string complementBuffer;
        ThreadBuffers() {
            result.reserve(1024*1024);
            queryProfData.reserve(1024);
            queryHeaderBuffer.reserve(1024);
            targetProfData.reserve(1024);
            newBacktrace.reserve(1024);
            complementBuffer.reserve(1024);
        }
    };
    std::vector<ThreadBuffers> threadBuffers(localThreads);

    // all output of one alignment entry: text into the result buffer, or rows into parquetWriter (parquet format)
    auto processEntry = [&](const DBKeyType queryKey, char *data, const unsigned int thread_idx, ParquetDBWriter *parquetWriter) {
            ThreadBuffers& tb = threadBuffers[thread_idx];
            char (&buffer)[1024] = tb.buffer;
            std::string& result = tb.result;
            std::string& queryProfData = tb.queryProfData;
            std::string& queryHeaderBuffer = tb.queryHeaderBuffer;
            std::string& targetProfData = tb.targetProfData;
            std::string& newBacktrace = tb.newBacktrace;
            std::string& complementBuffer = tb.complementBuffer;
            const TaxonNode * taxonNode = NULL;
            char *querySeqData = NULL;
            size_t querySeqLen = 0;
            queryProfData.clear();
            if (needSequenceDB) {
                size_t qId = qDbr.sequenceReader->getId(queryKey);
                querySeqData = qDbr.sequenceReader->getData(qId, thread_idx);
                querySeqLen = qDbr.sequenceReader->getSeqLen(qId);
                if (queryProfile) {
                    size_t queryEntryLen = qDbr.sequenceReader->getEntryLen(qId);
                    Sequence::extractProfileConsensus(querySeqData, queryEntryLen, *subMat, queryProfData);
                }
            }

            size_t qHeaderId = qDbrHeader.sequenceReader->getId(queryKey);
            const char *qHeader = qDbrHeader.sequenceReader->getData(qHeaderId, thread_idx);
            size_t qHeaderLen = qDbrHeader.sequenceReader->getSeqLen(qHeaderId);
            std::string queryId = Util::parseFastaHeader(qHeader);
            if (sameDB && needFullHeaders) {
                queryHeaderBuffer.assign(qHeader, qHeaderLen);
                qHeader = (char*) queryHeaderBuffer.c_str();
            }

            if (format == Parameters::FORMAT_ALIGNMENT_HTML) {
                const char* jsStart = "{\"query\": {\"accession\": \"%s\",\"sequence\": \"";
                int count = snprintf(buffer, sizeof(buffer), jsStart, queryId.c_str(), querySeqData);
                if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                    Debug(Debug::WARNING) << "Truncated line in entry " << queryKey << "!\n";
                    return;
                }
                result.append(buffer, count);
                if (queryProfile) {
                    result.append(queryProfData);
                } else {
                    result.append(querySeqData, querySeqLen);
                }
                result.append("\"}, \"alignments\": [\n");
            }

            int row_counter = 0;
            while (*data != '\0') {
                Matcher::result_t res = Matcher::parseAlignmentRecord(data, true);
                data = Util::skipLine(data);
                if (res.backtrace.empty() && needBacktrace == true) {
                    Debug(Debug::ERROR) << "Backtrace cigar is missing in the alignment result. Please recompute the alignment with the -a flag.\n"
                                           "Command: mmseqs align " << par.db1 << " " << par.db2 << " " << par.db3 << " " << "alnNew -a\n";
                    EXIT(EXIT_FAILURE);
                }

                size_t tHeaderId = tDbrHeader->sequenceReader->getId(res.dbKey);
                const char *tHeader = tDbrHeader->sequenceReader->getData(tHeaderId, thread_idx);
                size_t tHeaderLen = tDbrHeader->sequenceReader->getSeqLen(tHeaderId);
                std::string targetId = Util::parseFastaHeader(tHeader);

                unsigned int gapOpenCount = 0;
                unsigned int alnLen = res.alnLength;
                unsigned int missMatchCount = 0;
                unsigned int identical = 0;
                if (res.backtrace.empty() == false) {
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
                    identical = static_cast<unsigned int>(res.seqId * static_cast<float>(alnLen) + 0.5);
                    //res.alnLength = alnLen;
                    missMatchCount = static_cast<unsigned int>( matchCount - identical);
                } else {
                    const int adjustQstart = (res.qStartPos == -1) ? 0 : res.qStartPos;
                    const int adjustDBstart = (res.dbStartPos == -1) ? 0 : res.dbStartPos;
                    const float bestMatchEstimate = static_cast<float>(std::min(abs(res.qEndPos - adjustQstart), abs(res.dbEndPos - adjustDBstart)));
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
                                Debug(Debug::WARNING) << "Truncated line in entry " << queryKey << "!\n";
                                continue;
                            }
                            result.append(buffer, count);
                        } else {
                            char *targetSeqData = NULL;
                            targetProfData.clear();
                            unsigned int taxon = 0;
                            if (needTaxonomy || needTaxonomyMapping) {
                                taxon = mapping->lookup(res.dbKey);
                                if (taxon == 0) {
                                    taxonNode = NULL;
                                } else if (needTaxonomy) {
                                    taxonNode = t->taxonNode(taxon, false);
                                }
                            }
                            if (needSequenceDB) {
                                size_t tId = tDbr->sequenceReader->getId(res.dbKey);
                                targetSeqData = tDbr->sequenceReader->getData(tId, thread_idx);
                                if (targetProfile) {
                                    size_t targetEntryLen = tDbr->sequenceReader->getEntryLen(tId);
                                    Sequence::extractProfileConsensus(targetSeqData, targetEntryLen, *subMat, targetProfData);
                                }
                            }
                            for(size_t i = 0; i < outcodes.size(); i++) {
                                switch (outcodes[i]) {
                                    case Parameters::OUTFMT_QUERY:
                                        if (isParquet){
                                            parquetWriter->writeCell(queryId, i);
                                        } else {
                                            result.append(queryId);
                                        }
                                        break;
                                    case Parameters::OUTFMT_TARGET:
                                        if (isParquet){
                                           parquetWriter->writeCell(targetId, i);
                                        } else {
                                            result.append(targetId);
                                        }
                                        break;
                                    case Parameters::OUTFMT_EVALUE:
                                        if (isParquet){
                                            parquetWriter->writeCell(res.eval, i);
                                        }
                                        else {
                                            result.append(SSTR(res.eval));
                                        }
                                        break;
                                    case Parameters::OUTFMT_GAPOPEN:
                                        if (isParquet){
                                            parquetWriter->writeCell(static_cast<int>(gapOpenCount), i);
                                        } else {
                                            result.append(SSTR(gapOpenCount));
                                        }
                                        break;
                                    case Parameters::OUTFMT_FIDENT:
                                        if (isParquet){
                                            parquetWriter->writeCell(res.seqId, i);
                                        } else {
                                            result.append(SSTR(res.seqId));
                                        }
                                        break;
                                    case Parameters::OUTFMT_PIDENT:
                                        if (isParquet){
                                            parquetWriter->writeCell(res.seqId*100, i);
                                        } else {
                                            result.append(SSTR(res.seqId*100));
                                        }
                                        break;
                                    case Parameters::OUTFMT_NIDENT:
                                        if (isParquet){
                                            parquetWriter->writeCell(static_cast<int>(identical), i);
                                        } else {
                                        result.append(SSTR(identical));
                                        }
                                        break;
                                    case Parameters::OUTFMT_QSTART:
                                        if (isParquet){
                                            parquetWriter->writeCell(res.qStartPos + 1, i);
                                        } else {
                                            result.append(SSTR(res.qStartPos + 1));
                                        }
                                        break;
                                    case Parameters::OUTFMT_QEND:
                                        if (isParquet){
                                            parquetWriter->writeCell(res.qEndPos + 1, i);
                                        } else {
                                        result.append(SSTR(res.qEndPos + 1));
                                        }
                                        break;
                                    case Parameters::OUTFMT_QLEN:
                                        if (isParquet){
                                            parquetWriter->writeCell(static_cast<int>(res.qLen), i);
                                        } else {
                                        result.append(SSTR(res.qLen));
                                        }
                                        break;
                                    case Parameters::OUTFMT_TSTART:
                                        if (isParquet){
                                            parquetWriter->writeCell(res.dbStartPos + 1, i);
                                        } else {
                                        result.append(SSTR(res.dbStartPos + 1));
                                        }
                                        break;
                                    case Parameters::OUTFMT_TEND:
                                        if (isParquet){
                                            parquetWriter->writeCell(res.dbEndPos + 1, i);
                                        } else {
                                        result.append(SSTR(res.dbEndPos + 1));
                                        }
                                        break;
                                    case Parameters::OUTFMT_TLEN:
                                        if (isParquet){
                                            parquetWriter->writeCell(static_cast<int>(res.dbLen), i);
                                        } else {
                                        result.append(SSTR(res.dbLen));
                                        }
                                        break;
                                    case Parameters::OUTFMT_ALNLEN:
                                        if (isParquet){
                                            parquetWriter->writeCell(static_cast<int>(alnLen), i);
                                        } else {
                                        result.append(SSTR(alnLen));
                                        }
                                        break;
                                    case Parameters::OUTFMT_RAW:
                                        if (isParquet){
                                            parquetWriter->writeCell(static_cast<int>(evaluer->computeRawScoreFromBitScore(res.score) + 0.5), i);
                                        } else {
                                        result.append(SSTR(static_cast<int>(evaluer->computeRawScoreFromBitScore(res.score) + 0.5)));
                                        }
                                        break;
                                    case Parameters::OUTFMT_BITS:
                                        if (isParquet){
                                            parquetWriter->writeCell(res.score, i);
                                        } else {
                                        result.append(SSTR(res.score));
                                        }
                                        break;
                                    case Parameters::OUTFMT_CIGAR:
                                        if(isTranslatedSearch == true && targetNucs == true && queryNucs == true ){
                                            Matcher::result_t::protein2nucl(res.backtrace, newBacktrace);
                                            res.backtrace = newBacktrace;
                                        }
                                        if (isParquet){
                                            parquetWriter->writeCell(res.backtrace, i);
                                            } else {
                                            result.append(SSTR(res.backtrace));
                                        }
                                        newBacktrace.clear();
                                        break;
                                    case Parameters::OUTFMT_QSEQ:
                                        if (queryProfile) {
                                            if (isParquet){
                                                parquetWriter->writeCell(queryProfData, i);
                                            } else {
                                                result.append(queryProfData.c_str(), res.qLen);
                                            }
                                        } else {
                                            if (isParquet){
                                                parquetWriter->writeCell(std::string(querySeqData,res.qLen), i);
                                            } else {
                                                result.append(querySeqData, res.qLen);
                                            }
                                        }
                                        break;
                                    case Parameters::OUTFMT_TSEQ:
                                        if (targetProfile) {
                                            if (isParquet){
                                                parquetWriter->writeCell(targetProfData, i);
                                            } else {
                                                result.append(targetProfData.c_str(), res.dbLen);
                                            }                                        
                                        } else {
                                            if (isParquet){
                                                parquetWriter->writeCell(std::string(targetSeqData,  res.dbLen), i);
                                            } else {                                            
                                                result.append(targetSeqData, res.dbLen);
                                            }
                                        }
                                        break;
                                    case Parameters::OUTFMT_QHEADER:
                                        if (isParquet){
                                            parquetWriter->writeCell(std::string(qHeader, qHeaderLen),i);
                                        } else {  
                                            result.append(qHeader, qHeaderLen);
                                        }
                                        break;
                                    case Parameters::OUTFMT_THEADER:
                                        if (isParquet){
                                            parquetWriter->writeCell(std::string(tHeader, tHeaderLen), i);
                                        } else { 
                                            result.append(tHeader, tHeaderLen);
                                        }
                                        break;
                                    case Parameters::OUTFMT_QALN:
                                        if (queryProfile) {
                                           if (isParquet){
                                                parquetWriter->addSeqBasedOnAln(i, queryProfData.c_str(), res.qStartPos,
                                                               Matcher::uncompressAlignment(res.backtrace), false, (res.qStartPos > res.qEndPos),
                                                               (isTranslatedSearch == true && queryNucs == true), translateNucl);
                                           } else {
                                                printSeqBasedOnAln(result, queryProfData.c_str(), res.qStartPos,
                                                               Matcher::uncompressAlignment(res.backtrace), false, (res.qStartPos > res.qEndPos),
                                                               (isTranslatedSearch == true && queryNucs == true), translateNucl);
                                           }
                                        } else {
                                            if (isParquet){
                                                parquetWriter->addSeqBasedOnAln(i, querySeqData, res.qStartPos,
                                                               Matcher::uncompressAlignment(res.backtrace), false, (res.qStartPos > res.qEndPos),
                                                               (isTranslatedSearch == true && queryNucs == true), translateNucl);
                                            } else {
                                                printSeqBasedOnAln(result, querySeqData, res.qStartPos,
                                                               Matcher::uncompressAlignment(res.backtrace), false, (res.qStartPos > res.qEndPos),
                                                               (isTranslatedSearch == true && queryNucs == true), translateNucl);
                                            }    
                                        }
                                        break;
                                    case Parameters::OUTFMT_TALN: {
                                        if (targetProfile) {
                                            if (isParquet){
                                                parquetWriter->addSeqBasedOnAln(i, targetProfData.c_str(), res.dbStartPos,
                                                               Matcher::uncompressAlignment(res.backtrace), true,
                                                               (res.dbStartPos > res.dbEndPos),
                                                               (isTranslatedSearch == true && targetNucs == true), translateNucl);
                                            } else {
                                                printSeqBasedOnAln(result, targetProfData.c_str(), res.dbStartPos,
                                                               Matcher::uncompressAlignment(res.backtrace), true,
                                                               (res.dbStartPos > res.dbEndPos),
                                                               (isTranslatedSearch == true && targetNucs == true), translateNucl);
                                            }
                                        } else {
                                            if (isParquet){
                                                parquetWriter->addSeqBasedOnAln(i, targetSeqData, res.dbStartPos,
                                                               Matcher::uncompressAlignment(res.backtrace), true,
                                                               (res.dbStartPos > res.dbEndPos),
                                                               (isTranslatedSearch == true && targetNucs == true), translateNucl);
                                            } else {
                                                printSeqBasedOnAln(result, targetSeqData, res.dbStartPos,
                                                               Matcher::uncompressAlignment(res.backtrace), true,
                                                               (res.dbStartPos > res.dbEndPos),
                                                               (isTranslatedSearch == true && targetNucs == true), translateNucl);
                                            }
                                        }
                                        break;
                                    }
                                    case Parameters::OUTFMT_MISMATCH:
                                        if (isParquet){
                                            parquetWriter->writeCell(static_cast<int>(missMatchCount), i);
                                        } else { 
                                            result.append(SSTR(missMatchCount));
                                        }
                                        break;
                                    case Parameters::OUTFMT_QCOV:
                                        if (isParquet){
                                            parquetWriter->writeCell(res.qcov, i);
                                        } else { 
                                            result.append(SSTR(res.qcov));
                                        }
                                        break;
                                    case Parameters::OUTFMT_TCOV:
                                        if (isParquet){
                                            parquetWriter->writeCell(res.dbcov, i);
                                        } else { 
                                            result.append(SSTR(res.dbcov));
                                        }
                                        break;
                                    case Parameters::OUTFMT_QSET:
                                        if (isParquet){
                                            parquetWriter->writeCell(qSetToSource[qKeyToSet[queryKey]], i);
                                        } else { 
                                            result.append(SSTR(qSetToSource[qKeyToSet[queryKey]]));
                                        }
                                        break;
                                    case Parameters::OUTFMT_QSETID:
                                        if (isParquet){
                                            parquetWriter->writeCell(static_cast<int64_t>(qKeyToSet[queryKey]), i);
                                        } else {
                                            result.append(SSTR(qKeyToSet[queryKey]));
                                        }
                                        break;
                                    case Parameters::OUTFMT_TSET:
                                        if (isParquet){
                                            parquetWriter->writeCell(tSetToSource[tKeyToSet[res.dbKey]], i);
                                        } else {
                                            result.append(SSTR(tSetToSource[tKeyToSet[res.dbKey]]));
                                        }
                                        break;
                                    case Parameters::OUTFMT_TSETID:
                                        if (isParquet){
                                            parquetWriter->writeCell(static_cast<int64_t>(tKeyToSet[res.dbKey]), i);
                                        } else {
                                            result.append(SSTR(tKeyToSet[res.dbKey]));
                                        }
                                        break;
                                    case Parameters::OUTFMT_TAXID:
                                        if (isParquet){
                                            parquetWriter->writeCell(static_cast<int>(taxon), i);
                                        } else {
                                            result.append(SSTR(taxon));
                                        }
                                        break;
                                    case Parameters::OUTFMT_TAXNAME:
                                        if (isParquet){
                                            parquetWriter->writeCell((taxonNode != NULL) ? t->getString(taxonNode->nameIdx) : "unclassified", i);
                                        } else {
                                            result.append((taxonNode != NULL) ? t->getString(taxonNode->nameIdx) : "unclassified");
                                        }
                                        break;
                                    case Parameters::OUTFMT_TAXLIN:
                                        if (isParquet){
                                            parquetWriter->writeCell((taxonNode != NULL) ? t->taxLineage(taxonNode, true) : "unclassified", i);
                                        } else {
                                            result.append((taxonNode != NULL) ? t->taxLineage(taxonNode, true) : "unclassified");
                                        }
                                        break;
                                    case Parameters::OUTFMT_EMPTY:
                                        if (isParquet){
                                            parquetWriter->writeCell(std::string("-"), i);
                                        } else {
                                            result.push_back('-');
                                        }
                                        break;
                                    case Parameters::OUTFMT_QORFSTART:
                                        if (isParquet){
                                            parquetWriter->writeCell(res.queryOrfStartPos, i);
                                        } else {
                                            result.append(SSTR(res.queryOrfStartPos));
                                        }
                                        break;
                                    case Parameters::OUTFMT_QORFEND:
                                        if (isParquet){
                                            parquetWriter->writeCell(res.queryOrfEndPos, i);
                                        } else {
                                        result.append(SSTR(res.queryOrfEndPos));
                                        }
                                        break;
                                    case Parameters::OUTFMT_TORFSTART:
                                        if (isParquet){
                                            parquetWriter->writeCell(res.dbOrfStartPos, i);
                                        } else {
                                        result.append(SSTR(res.dbOrfStartPos));
                                        }
                                        break;
                                    case Parameters::OUTFMT_TORFEND:
                                        if (isParquet){
                                            parquetWriter->writeCell(res.dbOrfEndPos, i);
                                        } else {
                                            result.append(SSTR(res.dbOrfEndPos));
                                        }
                                        break;
                                    case Parameters::OUTFMT_PPOS: {
                                        float pPositive = 0;
                                        int matchCount = 0;
                                        if (res.backtrace.empty() == false) {
                                            int qPos = res.qStartPos;
                                            int tPos = res.dbStartPos;
                                            std::string unpackedBt = Matcher::uncompressAlignment(res.backtrace);
                                            for (size_t pos = 0; pos < unpackedBt.size(); pos++) {
                                                switch (unpackedBt[pos]) {
                                                    case 'M': {
                                                        char qRes = queryProfile ? queryProfData[qPos] : querySeqData[qPos];
                                                        char tRes = targetProfile ? targetProfData[tPos] : targetSeqData[tPos];
                                                        pPositive += (subMat->subMatrix[subMat->aa2num[(int)qRes]][subMat->aa2num[(int)tRes]] > 0);
                                                        matchCount += 1;
                                                        qPos++;
                                                        tPos++;
                                                        break;
                                                    }
                                                    case 'D':
                                                        tPos++;
                                                        break;
                                                    case 'I':
                                                        qPos++;
                                                        break;
                                                }
                                            }
                                            pPositive /= static_cast<float>(matchCount);
                                        }
                                        if (isParquet){
                                            parquetWriter->writeCell(pPositive, i);
                                        } else {
                                            result.append(SSTR(pPositive));
                                        }
                                        break;
                                    }
                                    case Parameters::OUTFMT_QFRAME: {
                                        int frame;
                                        if (res.qStartPos <= res.qEndPos) {
                                            frame = (res.qStartPos - 1) % 3 + 1;
                                        } else {
                                            frame = -1 * ((res.qLen - res.qStartPos) % 3 + 1);
                                        }
                                        if (isParquet){
                                            parquetWriter->writeCell(frame, i);
                                        } else {
                                            result.append(SSTR(frame));
                                        }
                                        break;
                                    }
                                    case Parameters::OUTFMT_TFRAME: {
                                        int frame;
                                        if (res.dbStartPos <= res.dbEndPos) {
                                            frame = (res.dbStartPos - 1) % 3 + 1;
                                        } else {
                                            frame = -1 * ((res.dbLen - res.dbStartPos) % 3 + 1);
                                        }
                                        if (isParquet){
                                            parquetWriter->writeCell(frame, i);
                                        } else {
                                            result.append(SSTR(frame));
                                        }
                                        break;
                                    } 
                                }
                                if (i < outcodes.size() - 1 && !isParquet) {
                                    result.push_back('\t');
                                }
                            }
                            if (!isParquet){
                            result.push_back('\n');
                            }
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
                            Debug(Debug::WARNING) << "Truncated line in entry " << queryKey << "!\n";
                            continue;
                        }

                        result.append(buffer, count);
                        break;
                    }
                    case Parameters::FORMAT_ALIGNMENT_SAM: {
                        bool forward = res.qEndPos > res.qStartPos;
                        int rawScore = static_cast<int>(evaluer->computeRawScoreFromBitScore(res.score) + 0.5);
                        uint32_t mapq = -4.343 * log(exp(static_cast<double>(-rawScore)));
                        mapq = (uint32_t) (mapq + 4.99);
                        mapq = mapq < 254 ? mapq : 254;
                        int count = snprintf(buffer, sizeof(buffer), "%s\t%d\t%s\t%d\t%d\t",  queryId.c_str(), (forward) ? 0 : 16, targetId.c_str(), std::min(res.dbStartPos + 1, res.dbEndPos + 1), mapq);
                        if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                            Debug(Debug::WARNING) << "Truncated line in entry " << queryKey << "!\n";
                            continue;
                        }
                        result.append(buffer, count);
                        if (isTranslatedSearch == true && queryNucs == true) {
                            Matcher::result_t::protein2nucl(res.backtrace, newBacktrace);
                            result.append(newBacktrace);
                            newBacktrace.clear();
                        } else {
                            result.append(res.backtrace);
                        }
                        result.append("\t*\t0\t0\t");
                        int start = std::min(res.qStartPos, res.qEndPos);
                        int end   = std::max(res.qStartPos, res.qEndPos);
                        if (forward == false && queryNucs == true) {
                            if (queryProfile) {
                                complementBuffer.assign(queryProfData.c_str() + start, (end + 1) - start);
                            } else {
                                complementBuffer.assign(querySeqData + start, (end + 1) - start);
                            }
                            std::reverse(complementBuffer.begin(), complementBuffer.end());
                            for (size_t i = 0; i < complementBuffer.size(); i++) {
                                complementBuffer[i] = Orf::complement(complementBuffer[i]);
                            }
                            result.append(complementBuffer);
                            complementBuffer.clear();
                        } else {
                            if (queryProfile) {
                                result.append(queryProfData.c_str() + start, (end + 1) - start);
                            } else {
                                result.append(querySeqData + start, (end + 1) - start);
                            }
                        }
                        count = snprintf(buffer, sizeof(buffer), "\t*\tAS:i:%d\tNM:i:%d\n", rawScore, missMatchCount);
                        if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                            Debug(Debug::WARNING) << "Truncated line in entry " << queryKey << "!\n";
                            continue;
                        }
                        result.append(buffer, count);
                        break;
                    }
                    case Parameters::FORMAT_ALIGNMENT_HTML: {
                        const char* jsAln = "{\"target\": \"%s\", \"seqId\": %1.3f, \"alnLen\": %d, \"mismatch\": %d, \"gapopen\": %d, \"qStartPos\": %d, \"qEndPos\": %d, \"dbStartPos\": %d, \"dbEndPos\": %d, \"eval\": %.2E, \"score\": %d, \"qLen\": %d, \"dbLen\": %d, \"qAln\": \"";
                        int count = snprintf(buffer, sizeof(buffer), jsAln,
                                             targetId.c_str(), res.seqId, alnLen,
                                             missMatchCount, gapOpenCount,
                                             res.qStartPos + 1, res.qEndPos + 1,
                                             res.dbStartPos + 1, res.dbEndPos + 1,
                                             res.eval, res.score,
                                             res.qLen, res.dbLen);

                        if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                            Debug(Debug::WARNING) << "Truncated line in entry " << queryKey << "!\n";
                            continue;
                        }
                        result.append(buffer, count);
                        if (queryProfile) {
                            printSeqBasedOnAln(result, queryProfData.c_str(), res.qStartPos,
                                               Matcher::uncompressAlignment(res.backtrace), false, (res.qStartPos > res.qEndPos),
                                               (isTranslatedSearch == true && queryNucs == true), translateNucl);
                        } else {
                            printSeqBasedOnAln(result, querySeqData, res.qStartPos,
                                               Matcher::uncompressAlignment(res.backtrace), false, (res.qStartPos > res.qEndPos),
                                               (isTranslatedSearch == true && queryNucs == true), translateNucl);
                        }
                        result.append("\", \"dbAln\": \"");
                        size_t tId = tDbr->sequenceReader->getId(res.dbKey);
                        char* targetSeqData = tDbr->sequenceReader->getData(tId, thread_idx);
                        if (targetProfile) {
                            size_t targetEntryLen = tDbr->sequenceReader->getEntryLen(tId);
                            Sequence::extractProfileConsensus(targetSeqData, targetEntryLen, *subMat, targetProfData);
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
                        result.append("\" },\n");
                        break;
                    }

//                    case Parameters::FORMAT_ALIGNMENT_GFF:{
//                        // for TBLASTX
//                        bool strand = res.qEndPos > res.qStartPos;
//                        int currStart = std::min(res.qStartPos, res.qEndPos);
//                        int currEnd = std::max(res.qStartPos, res.qEndPos);
//                        int currLen = currEnd - currStart;
//                        result.append(queryId);
//                        result.append("\tconserve\tprotein_match\t");
//                        result.append(SSTR(currStart+1));
//                        result.push_back('\t');
//                        result.append(SSTR(currEnd+1));
//                        result.push_back('\t');
//                        result.append(SSTR(currLen));
//                        result.push_back('\t');
//                        result.push_back((strand) ? '-' : '+');
//                        result.append("\t.\t");
//                        result.append("ID=");
//                        result.append(queryId);
//                        result.append(":hsp:");
//                        result.append(SSTR(counter));
//                        result.append(";");
//                        break;
//                    }
                    default:
                        Debug(Debug::ERROR) << "Not implemented yet";
                        EXIT(EXIT_FAILURE);
                }
                if (isParquet) {
                    // row complete, the writer flushes to carquet once enough data is pending
                    parquetWriter->endRow();
                }
                row_counter++;
            }

            if (format == Parameters::FORMAT_ALIGNMENT_HTML) {
                result.append("]},\n");
            }

            if (!isParquet) {
                resultWriter.writeData(result.c_str(), result.size(), queryKey, thread_idx, isDb);
                result.clear();
            }
    };

    if (isParquet) {
        // file layout, sort modes and writer options of convert2parquet; columns (init(outcodes)) and
        // file names (<output dir>/<output name>_<n>.parquet) stay those of convertalis
        const std::string dir = FileUtil::dirName(par.db4);
        const std::string base = FileUtil::baseName(par.db4);
        const size_t lastdot = base.find_last_of(".");
        const std::string noExtension = (lastdot == std::string::npos) ? base : base.substr(0, lastdot);
        ParquetDBWriter::ConvertHooks hooks;
        hooks.partNames = [&](size_t files) {
            ParquetDBWriter::removeNumberedParquet(dir, noExtension);
            std::vector<std::string> names;
            for (size_t f = 0; f < files; f++) {
                names.push_back(dir + "/" + noExtension + "_" + std::to_string(f) + ".parquet");
            }
            return names;
        };
        hooks.openPart = [&](ParquetDBWriter& writer, const std::string& file) {
            writer.init(file.c_str(), outcodes);
        };
        hooks.writeEntry = [&](ParquetDBWriter& writer, DBKeyType key, char* data, size_t, unsigned int thread_idx) {
            processEntry(key, data, thread_idx, &writer);
        };
        // one file per thread at least, as before
        hooks.minFiles = localThreads;
        ParquetDBWriter::writeDbAsParquet(alnDbr, par.db4, par.parquetSortMode, localThreads,
                                          ParquetDBWriter::writerOptionsFromParameters(par), NULL,
                                          std::vector<ParquetDBWriter::ColumnSpec>(), hooks);
    } else {
#pragma omp parallel num_threads(localThreads)
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for schedule(dynamic, 10)
            for (size_t i = 0; i < alnDbr.getSize(); i++) {
                progress.updateProgress();
                processEntry(alnDbr.getDbKey(i), alnDbr.getData(i, thread_idx), thread_idx, NULL);
            }
        }
    }
    if (format == Parameters::FORMAT_ALIGNMENT_HTML) {
        const char* endBlock = "]);</script>";
        resultWriter.writeData(endBlock, strlen(endBlock), 0, localThreads - 1, false, false);
    }
    // tsv output
    if (!isParquet){
        resultWriter.close(true);
        if (isDb == false) {
            FileUtil::remove(par.db4Index.c_str());
        }
    }
    if (needTaxonomy) {
        delete t;
    }
    if (mapping != NULL) {
        delete mapping;
    }
    
    alnDbr.close();

    if (sameDB == false) {
        delete tDbr;
        delete tDbrHeader;
    }
    if (needSequenceDB) {
        delete evaluer;
    }
    delete subMat;

    return EXIT_SUCCESS;
}
