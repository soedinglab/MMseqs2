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
    par.parseParameters(argc, argv, command, 2);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> headerReader(par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    headerReader.open(DBReader<unsigned int>::NOSORT);

    DBWriter sequenceWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, reader.getDbtype());
    sequenceWriter.open();

    DBWriter headerWriter(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, false, Parameters::DBTYPE_OFFSETDB);
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
        Util::decomposeDomainByAminoAcid(reader.getDataSize(), reader.getSeqLens(), reader.getSize(),
                                         thread_idx, par.threads, &queryFrom, &querySize);
        if (querySize == 0) {
            queryFrom = 0;
        }
        char buffer[LINE_MAX];

        for (unsigned int i = queryFrom; i < (queryFrom + querySize); ++i){
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            const char* data = reader.getData(i, thread_idx);
            size_t dataLength = reader.getSeqLens(i);
            size_t seqLen = dataLength -2;
            char* header = headerReader.getData(i, thread_idx);
            Orf::SequenceLocation loc = Orf::parseOrfHeader(header);
            size_t from = 0;
            if(loc.id != UINT_MAX) {
                from = (loc.strand==Orf::STRAND_MINUS)? loc.to : loc.from;
            }
            unsigned int dbKey = (loc.id != UINT_MAX) ? loc.id : key;
            size_t splitCnt = (size_t) ceilf(static_cast<float>(seqLen) / static_cast<float>(par.maxSeqLen - sequenceOverlap));
            std::string headerAccession = Util::parseFastaHeader(header);

            for (size_t split = 0; split < splitCnt; split++) {
                size_t len = std::min(par.maxSeqLen, seqLen - (split * par.maxSeqLen - split*sequenceOverlap));
                sequenceWriter.writeStart(thread_idx);
                size_t startPos = split * par.maxSeqLen - split*sequenceOverlap;
                sequenceWriter.writeAdd(data + startPos, len, thread_idx);
                char newLine = '\n';
                sequenceWriter.writeAdd(&newLine, 1, thread_idx);
                sequenceWriter.writeEnd(key, thread_idx, true);
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

    // make identifiers stable
#pragma omp parallel
    {
#pragma omp single
        {
#pragma omp task
            {
                DBReader<unsigned int> frameHeaderReader(par.hdr2.c_str(), par.hdr2Index.c_str(),
                                                       par.threads,
                                                       DBReader<unsigned int>::USE_INDEX);
                frameHeaderReader.open(DBReader<unsigned int>::SORT_BY_ID_OFFSET);
                FILE *hIndex = fopen((par.hdr2Index + "_tmp").c_str(), "w");
                if (hIndex == NULL) {
                    Debug(Debug::ERROR) << "Could not open " << par.hdr2Index << "_tmp for writing!\n";
                    EXIT(EXIT_FAILURE);
                }
                for (size_t i = 0; i < frameHeaderReader.getSize(); i++) {
                    DBReader<unsigned int>::Index *idx = frameHeaderReader.getIndex(i);
                    char buffer[1024];
                    size_t len = DBWriter::indexToBuffer(buffer, i, idx->offset, frameHeaderReader.getSeqLens(i));
                    int written = fwrite(buffer, sizeof(char), len, hIndex);
                    if (written != (int) len) {
                        Debug(Debug::ERROR) << "Could not write to data file " << par.hdr2Index << "_tmp\n";
                        EXIT(EXIT_FAILURE);
                    }
                }
                fclose(hIndex);
                frameHeaderReader.close();
                std::rename((par.hdr2Index + "_tmp").c_str(), par.hdr2Index.c_str());
            }

#pragma omp task
            {
                DBReader<unsigned int> frameSequenceReader(par.db2.c_str(), par.db2Index.c_str(),
                                                         par.threads,
                                                         DBReader<unsigned int>::USE_INDEX);
                frameSequenceReader.open(DBReader<unsigned int>::SORT_BY_ID_OFFSET);

                FILE *sIndex = fopen((par.db2Index + "_tmp").c_str(), "w");
                if (sIndex == NULL) {
                    Debug(Debug::ERROR) << "Could not open " << par.db2Index << "_tmp for writing!\n";
                    EXIT(EXIT_FAILURE);
                }

                for (size_t i = 0; i < frameSequenceReader.getSize(); i++) {
                    DBReader<unsigned int>::Index *idx = (frameSequenceReader.getIndex(i));
                    char buffer[1024];
                    size_t len = DBWriter::indexToBuffer(buffer, i, idx->offset, frameSequenceReader.getSeqLens(i));
                    int written = fwrite(buffer, sizeof(char), len, sIndex);
                    if (written != (int) len) {
                        Debug(Debug::ERROR) << "Could not write to data file " << par.db2Index << "_tmp\n";
                        EXIT(EXIT_FAILURE);
                    }
                }
                fclose(sIndex);
                frameSequenceReader.close();
                std::rename((par.db2Index + "_tmp").c_str(), par.db2Index.c_str());
            }
        }
    }
    return EXIT_SUCCESS;
}

