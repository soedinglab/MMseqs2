#include "Debug.h"
#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Matcher.h"
#include "Util.h"
#include "itoa.h"

#include "Orf.h"

#include <unistd.h>
#include <climits>
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif

int splitsequence(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.maxSeqLen = 10000;
    par.sequenceOverlap = 300;
    par.parseParameters(argc, argv, command, true, 0, 0);
    int mode = DBReader<unsigned int>::USE_INDEX;
    if(par.sequenceSplitMode == Parameters::SEQUENCE_SPLIT_MODE_HARD){
        mode |= DBReader<unsigned int>::USE_DATA;
    }
    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, mode);
    reader.open(DBReader<unsigned int>::NOSORT);
    bool sizeLarger = false;
    for(size_t i = 0; i < reader.getSize(); i++){
        sizeLarger |= (reader.getSeqLen(i) > par.maxSeqLen);
    }

    // if no sequence needs to be splitted
    if (sizeLarger == false) {
        DBReader<unsigned int>::softlinkDb(par.db1, par.db2, DBFiles::SEQUENCE_DB);
        reader.close();
        return EXIT_SUCCESS;
    }

    DBReader<unsigned int> headerReader(par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    headerReader.open(DBReader<unsigned int>::NOSORT);

    if(par.sequenceSplitMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT && par.compressed == true) {
        Debug(Debug::WARNING) << "Sequence split mode (--sequence-split-mode 0) and compressed (--compressed 1) can not be combined.\nTurn compressed to 0";
        par.compressed = 0;
    }

    DBWriter sequenceWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, reader.getDbtype());
    sequenceWriter.open();

    DBWriter headerWriter(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, false, Parameters::DBTYPE_GENERIC_DB);
    headerWriter.open();

    size_t sequenceOverlap = par.sequenceOverlap;
    Debug::Progress progress(reader.getSize());
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        size_t querySize = 0;
        size_t queryFrom = 0;
        reader.decomposeDomainByAminoAcid(thread_idx, par.threads, &queryFrom, &querySize);
        if (querySize == 0) {
            queryFrom = 0;
        }
        char buffer[1024];

        for (unsigned int i = queryFrom; i < (queryFrom + querySize); ++i){
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            const char* data=NULL;
            if(par.sequenceSplitMode == Parameters::SEQUENCE_SPLIT_MODE_HARD){
                data = reader.getData(i, thread_idx);
            }
            size_t seqLen = reader.getSeqLen(i);
            char* header = headerReader.getData(i, thread_idx);
            Orf::SequenceLocation loc = Orf::parseOrfHeader(header);
            size_t from = 0;
            if(loc.id != UINT_MAX) {
                from = (loc.strand==Orf::STRAND_MINUS)? loc.to : loc.from;
            }
            unsigned int dbKey = (loc.id != UINT_MAX) ? loc.id : key;
            size_t splitCnt = (size_t) ceilf(static_cast<float>(seqLen) / static_cast<float>(par.maxSeqLen - sequenceOverlap));

            for (size_t split = 0; split < splitCnt; split++) {
                size_t len = std::min(par.maxSeqLen, seqLen - (split * par.maxSeqLen - split*sequenceOverlap));
                size_t startPos = split * par.maxSeqLen - split*sequenceOverlap;
                if(par.sequenceSplitMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT){
                    // +2 to emulate the \n\0
                    sequenceWriter.writeIndexEntry(key, reader.getOffset(i) + startPos, len+2, thread_idx);
                }else{
                    sequenceWriter.writeStart(thread_idx);
                    sequenceWriter.writeAdd(data + startPos, len, thread_idx);
                    char newLine = '\n';
                    sequenceWriter.writeAdd(&newLine, 1, thread_idx);
                    sequenceWriter.writeEnd(key, thread_idx, true);
                }
                size_t fromPos = from + startPos;
                size_t toPos = (from + startPos) + (len - 1);
                if(loc.id != UINT_MAX && loc.strand == Orf::STRAND_MINUS){
                    fromPos = (seqLen - 1) - (from + startPos);
                    toPos   = fromPos - std::min(fromPos, len);
                }

                size_t bufferLen = Orf::writeOrfHeader(buffer, dbKey, fromPos, toPos, 0, 0);
                headerWriter.writeData(buffer, bufferLen, key, thread_idx);
            }
        }
    }
    headerWriter.close(true);
    sequenceWriter.close(true);
    headerReader.close();
    reader.close();
    if(par.sequenceSplitMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT){
        DBReader<unsigned int>::softlinkDb(par.db1, par.db2, DBFiles::DATA);
    }
    // make identifiers stable
#pragma omp parallel
    {
#pragma omp single
        {
#pragma omp task
            {
                DBWriter::createRenumberedDB(par.hdr2, par.hdr2Index, "", "");
            }

#pragma omp task
            {
                DBWriter::createRenumberedDB(par.db2, par.db2Index, par.createLookup ? par.db1 : "", par.createLookup ? par.db1Index : "");
            }
        }
    }
    DBReader<unsigned int>::softlinkDb(par.db1, par.db2, DBFiles::SOURCE);

    return EXIT_SUCCESS;
}

