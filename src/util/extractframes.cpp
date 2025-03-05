#include "Debug.h"
#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Matcher.h"
#include "Util.h"
#include "TranslateNucl.h"
#include "Orf.h"

#include <unistd.h>
#include <climits>
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif

void handleSingleFrame(TranslateNucl& translateNucl, DBWriter& sequenceWriter, DBWriter& headerWriter, unsigned int key, char* headerBuffer, const char* data, size_t seqLen, int frame, bool reverse, bool translate, char*& aaBuffer, size_t& aaBufferSize, int thread_idx) {
    data = data + frame;
    seqLen = seqLen - frame;
    if (translate == true) {
        if (seqLen < 3) {
            return;
        }
        size_t codonLength = (seqLen / 3) * 3;
        if ((codonLength + 1) > aaBufferSize) {
            aaBufferSize = codonLength * 1.5 + 1;
            aaBuffer = (char*)realloc(aaBuffer, aaBufferSize * sizeof(char));
        }
        translateNucl.translate(aaBuffer, data, codonLength);
        aaBuffer[codonLength / 3] = '\n';
        sequenceWriter.writeData(aaBuffer, (codonLength / 3) + 1, key, thread_idx);
        size_t bufferLen;
        if (reverse) {
            bufferLen = Orf::writeOrfHeader(headerBuffer, key, frame + codonLength, static_cast<size_t>(frame), 0, 0);
        } else {
            bufferLen = Orf::writeOrfHeader(headerBuffer, key, static_cast<size_t>(frame), frame + codonLength, 0, 0);
        }
        headerWriter.writeData(headerBuffer, bufferLen, key, thread_idx);
    } else {
        // +1: add newline, but remove it from the end pos
        sequenceWriter.writeData(data, seqLen + 1, key, thread_idx);
        size_t bufferLen;
        if (reverse) {
            bufferLen = Orf::writeOrfHeader(headerBuffer, key, seqLen - 1, static_cast<size_t>(frame), 0, 0);
        } else {
            bufferLen = Orf::writeOrfHeader(headerBuffer, key, static_cast<size_t>(frame), seqLen - 1, 0, 0);
        }
        headerWriter.writeData(headerBuffer, bufferLen, key, thread_idx);
    }
}

int extractframes(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::NOSORT);

    int outputDbtype = reader.getDbtype();
    if (par.translate) {
        outputDbtype = Parameters::DBTYPE_AMINO_ACIDS;
    }
    DBWriter sequenceWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, outputDbtype);
    sequenceWriter.open();

    DBWriter headerWriter(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, false, Parameters::DBTYPE_GENERIC_DB);
    headerWriter.open();

    unsigned int forwardFrames = Orf::getFrames(par.forwardFrames);
    unsigned int reverseFrames = Orf::getFrames(par.reverseFrames);

    Debug::Progress progress(reader.getSize());
    TranslateNucl translateNucl(static_cast<TranslateNucl::GenCode>(par.translationTable));
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

        size_t aaBufferSize = par.maxSeqLen + 3 + 1;
        char* aa = NULL;
        if (par.translate == true) {
            aa = (char*)malloc(aaBufferSize * sizeof(char));
        }

        char buffer[1024];

        std::string reverseComplementStr;
        reverseComplementStr.reserve(32000);

        for (unsigned int i = queryFrom; i < (queryFrom + querySize); ++i){
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            const char* data = reader.getData(i, thread_idx);
            size_t seqLen = reader.getSeqLen(i);

            if (forwardFrames & Orf::FRAME_1) {
                handleSingleFrame(translateNucl, sequenceWriter, headerWriter, key, buffer, data, seqLen, 0, false, par.translate, aa, aaBufferSize, thread_idx);
            }
            if (forwardFrames & Orf::FRAME_2) {
                handleSingleFrame(translateNucl, sequenceWriter, headerWriter, key, buffer, data, seqLen, 1, false, par.translate, aa, aaBufferSize, thread_idx);
            }
            if (forwardFrames & Orf::FRAME_3) {
                handleSingleFrame(translateNucl, sequenceWriter, headerWriter, key, buffer, data, seqLen, 2, false, par.translate, aa, aaBufferSize, thread_idx);
            }

            if (reverseFrames != 0) {
                // bool hasWrongChar = false;
                for (size_t pos = 0; pos < seqLen; ++pos) {
                    char reverseComplement = Orf::complement(data[seqLen - pos - 1]);
                    reverseComplement = (reverseComplement == '.') ? 'N' : reverseComplement;
                    reverseComplementStr.push_back(reverseComplement);
                    // hasWrongChar |= (reverseComplement == '.');
                }
                // if (hasWrongChar == true) {
                //     continue;
                // }
                reverseComplementStr.push_back('\n');
                data = reverseComplementStr.c_str();
            }

            if (reverseFrames & Orf::FRAME_1) {
                handleSingleFrame(translateNucl, sequenceWriter, headerWriter, key, buffer, data, seqLen, 0, true, par.translate, aa, aaBufferSize, thread_idx);
            }
                
            if (reverseFrames & Orf::FRAME_2) {
                handleSingleFrame(translateNucl, sequenceWriter, headerWriter, key, buffer, data, seqLen, 1, true, par.translate, aa, aaBufferSize, thread_idx);
            }
                
            if (reverseFrames & Orf::FRAME_3) {
                handleSingleFrame(translateNucl, sequenceWriter, headerWriter, key, buffer, data, seqLen, 2, true, par.translate, aa, aaBufferSize, thread_idx);
            }
            reverseComplementStr.clear();
        }
        if (aa != NULL) {
            free(aa);
        }
    }
    headerWriter.close(true);
    sequenceWriter.close(true);
    reader.close();

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

