#include "Debug.h"
#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Matcher.h"
#include "Util.h"
#include "TranslateNucl.h"
#include "itoa.h"

#include "Orf.h"

#include <unistd.h>
#include <climits>
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif

int extractframes(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::NOSORT);

    int outputDbtype = reader.getDbtype();
    if(par.translate) {
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

        char* aa = new char[par.maxSeqLen + 3 + 1];
        char buffer[1024];
        std::string reverseComplementStr;
        reverseComplementStr.reserve(32000);
        for (unsigned int i = queryFrom; i < (queryFrom + querySize); ++i){
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            const char* data = reader.getData(i, thread_idx);
            size_t seqLen = reader.getSeqLen(i);

            size_t bufferLen;
            if (forwardFrames & Orf::FRAME_1) {
                if (par.translate) {
                    size_t currSeqLen = seqLen + 1;
                    if (currSeqLen >= 3) {
                        if (currSeqLen > (3 * par.maxSeqLen)) {
                            currSeqLen = (3 * par.maxSeqLen);
                        }
                        size_t condonLength = currSeqLen / 3 * 3;
                        translateNucl.translate(aa, data, condonLength);
                        sequenceWriter.writeData(aa, (condonLength / 3), key, thread_idx);    
                        bufferLen = Orf::writeOrfHeader(buffer, key, static_cast<size_t >(0), seqLen - 1, 0, 0);
                        headerWriter.writeData(buffer, bufferLen, key, thread_idx);        
                    }
                } else {
                    sequenceWriter.writeData(data, seqLen + 1, key, thread_idx);
                    bufferLen = Orf::writeOrfHeader(buffer, key, static_cast<size_t >(0), seqLen - 1, 0, 0);
                    headerWriter.writeData(buffer, bufferLen, key, thread_idx);
                }
            }
            if (forwardFrames & Orf::FRAME_2) {
                if (par.translate) {
                    size_t currSeqLen = seqLen;
                    if (currSeqLen >= 3) {
                        if (currSeqLen > (3 * par.maxSeqLen)) {
                            currSeqLen = (3 * par.maxSeqLen);
                        }
                        size_t condonLength = currSeqLen / 3 * 3;
                        translateNucl.translate(aa, data + 1, condonLength);
                        sequenceWriter.writeData(aa, (condonLength / 3), key, thread_idx);    
                        bufferLen = Orf::writeOrfHeader(buffer, key, static_cast<size_t >(1), seqLen - 2, 0, 0);
                        headerWriter.writeData(buffer, bufferLen, key, thread_idx);        
                    }
                } else {
                    sequenceWriter.writeData(data + 1, seqLen, key, thread_idx);
                    bufferLen = Orf::writeOrfHeader(buffer, key, static_cast<size_t >(1), seqLen - 2, 0, 0);
                    headerWriter.writeData(buffer, bufferLen, key, thread_idx);
                }
            }
            if (forwardFrames & Orf::FRAME_3) {
                if (par.translate) {
                    size_t currSeqLen = seqLen - 1;
                    if (currSeqLen >= 3) {
                        if (currSeqLen > (3 * par.maxSeqLen)) {
                            currSeqLen = (3 * par.maxSeqLen);
                        }
                        size_t condonLength = currSeqLen / 3 * 3;
                        translateNucl.translate(aa, data + 2, condonLength);
                        sequenceWriter.writeData(aa, (condonLength / 3), key, thread_idx);    
                        bufferLen = Orf::writeOrfHeader(buffer, key, static_cast<size_t >(2), seqLen - 3, 0, 0);
                        headerWriter.writeData(buffer, bufferLen, key, thread_idx);        
                    }
                } else {
                    sequenceWriter.writeData(data + 2, seqLen - 1, key, thread_idx);
                    bufferLen = Orf::writeOrfHeader(buffer, key, static_cast<size_t >(2), seqLen - 3, 0, 0);
                    headerWriter.writeData(buffer, bufferLen, key, thread_idx);
                }
            }


            if(reverseFrames != 0){
                size_t sequenceLength =  seqLen;
                // bool hasWrongChar = false;
                for(size_t pos = 0; pos < sequenceLength; ++pos) {
                    char reverseComplement = Orf::complement(data[sequenceLength - pos - 1]);
                    reverseComplement = (reverseComplement == '.') ? 'N' : reverseComplement;
                    reverseComplementStr.push_back(reverseComplement);
                    // hasWrongChar |= (reverseComplement == '.');
                }
//                if(hasWrongChar == true){
//                    continue;
//                }
                reverseComplementStr.push_back('\n');
    
                seqLen = reverseComplementStr.size();
                data = reverseComplementStr.c_str();
            }

            if (reverseFrames & Orf::FRAME_1) {
                if (par.translate) {
                    size_t currSeqLen = seqLen;
                    if (currSeqLen >= 3) {
                        if (currSeqLen > (3 * par.maxSeqLen)) {
                            currSeqLen = (3 * par.maxSeqLen);
                        }
                        size_t condonLength = currSeqLen / 3 * 3;
                        translateNucl.translate(aa, data, condonLength);
                        sequenceWriter.writeData(aa, (condonLength / 3), key, thread_idx);    
                        bufferLen = Orf::writeOrfHeader(buffer, key, seqLen - 2, static_cast<size_t >(0), 0, 0);
                        headerWriter.writeData(buffer, bufferLen, key, thread_idx);        
                    }
                } else {
                    sequenceWriter.writeData(data, seqLen, key, thread_idx);
                    bufferLen = Orf::writeOrfHeader(buffer, key, seqLen - 2, static_cast<size_t >(0), 0, 0);
                    headerWriter.writeData(buffer, bufferLen, key, thread_idx);
                }
            }
                
            if (reverseFrames & Orf::FRAME_2) {
                if (par.translate) {
                    size_t currSeqLen = seqLen - 1;
                    if (currSeqLen >= 3) {
                        if (currSeqLen > (3 * par.maxSeqLen)) {
                            currSeqLen = (3 * par.maxSeqLen);
                        }
                        size_t condonLength = currSeqLen / 3 * 3;
                        translateNucl.translate(aa, data + 1, condonLength);
                        sequenceWriter.writeData(aa, (condonLength / 3), key, thread_idx);    
                        bufferLen = Orf::writeOrfHeader(buffer, key, seqLen - 3, static_cast<size_t >(1), 0, 0);
                        headerWriter.writeData(buffer, bufferLen, key, thread_idx);        
                    }
                } else {
                    sequenceWriter.writeData(data + 1, seqLen - 1, key, thread_idx);
                    bufferLen = Orf::writeOrfHeader(buffer, key, seqLen - 3, static_cast<size_t >(1), 0, 0);
                    headerWriter.writeData(buffer, bufferLen, key, thread_idx);
                }
            }
                
            if (reverseFrames & Orf::FRAME_3) {
                if (par.translate) {
                    size_t currSeqLen = seqLen - 2;
                    if (currSeqLen >= 3) {
                        if (currSeqLen > (3 * par.maxSeqLen)) {
                            currSeqLen = (3 * par.maxSeqLen);
                        }
                        size_t condonLength = currSeqLen / 3 * 3;
                        translateNucl.translate(aa, data + 2, condonLength);
                        sequenceWriter.writeData(aa, (condonLength / 3), key, thread_idx);    
                        bufferLen = Orf::writeOrfHeader(buffer, key, seqLen - 4, static_cast<size_t >(2), 0, 0);
                        headerWriter.writeData(buffer, bufferLen, key, thread_idx);        
                    }
                } else {
                    sequenceWriter.writeData(data + 2, seqLen - 2, key, thread_idx);
                    bufferLen = Orf::writeOrfHeader(buffer, key, seqLen - 4, static_cast<size_t >(2), 0, 0);
                    headerWriter.writeData(buffer, bufferLen, key, thread_idx);
                }                    
            }
            reverseComplementStr.clear();
        }
        delete[] aa;
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

