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

int extractframes(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::NOSORT);

    DBWriter sequenceWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, reader.getDbtype());
    sequenceWriter.open();

    DBWriter headerWriter(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, false, Parameters::DBTYPE_GENERIC_DB);
    headerWriter.open();

    unsigned int forwardFrames = Orf::getFrames(par.forwardFrames);
    unsigned int reverseFrames = Orf::getFrames(par.reverseFrames);
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
        std::string reverseComplementStr;
        reverseComplementStr.reserve(32000);
        for (unsigned int i = queryFrom; i < (queryFrom + querySize); ++i){
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            const char* data = reader.getData(i, thread_idx);
            size_t dataLength = reader.getEntryLen(i);

            size_t bufferLen;
            switch (forwardFrames){
                case Orf::FRAME_1:
                    // -1 to ignore the null byte copy the new line
                    sequenceWriter.writeData(data, dataLength - 1, key, thread_idx);
                    bufferLen = Orf::writeOrfHeader(buffer, key, static_cast<size_t >(0), dataLength - 3, 0, 0);
                    headerWriter.writeData(buffer, bufferLen, key, thread_idx);
                    break;
                case Orf::FRAME_2:
                    sequenceWriter.writeData(data + 1, dataLength - 2, key, thread_idx);
                    bufferLen = Orf::writeOrfHeader(buffer, key, static_cast<size_t >(1), dataLength - 4, 0, 0);
                    headerWriter.writeData(buffer, bufferLen, key, thread_idx);
                    break;
                case Orf::FRAME_3:
                    sequenceWriter.writeData(data + 2, dataLength - 3, key, thread_idx);
                    bufferLen = Orf::writeOrfHeader(buffer, key, static_cast<size_t >(2), dataLength - 5, 0, 0);
                    headerWriter.writeData(buffer, bufferLen, key, thread_idx);
                    break;
            }

            if(reverseFrames != 0){
                size_t sequenceLength =  dataLength -2;
                bool hasWrongChar = false;
                for(size_t pos = 0; pos < sequenceLength; ++pos) {
                    char reverseComplement = Orf::complement(data[sequenceLength - pos - 1]);
                    reverseComplement = (reverseComplement == '.') ? 'N' : reverseComplement;
                    reverseComplementStr.push_back(reverseComplement);
                    hasWrongChar |= (reverseComplement == '.');
                }
//                if(hasWrongChar == true){
//                    continue;
//                }
                reverseComplementStr.push_back('\n');
            }

            switch (reverseFrames){
                case Orf::FRAME_1:
                    sequenceWriter.writeData(reverseComplementStr.c_str(), reverseComplementStr.size(), key, thread_idx);
                    bufferLen = Orf::writeOrfHeader(buffer, key, reverseComplementStr.size() - 2, static_cast<size_t >(0), 0, 0);
                    headerWriter.writeData(buffer, bufferLen, key, thread_idx);
                    break;
                case Orf::FRAME_2:
                    sequenceWriter.writeData(reverseComplementStr.c_str()+1, reverseComplementStr.size()-1, key, thread_idx);
                    bufferLen = Orf::writeOrfHeader(buffer, key, reverseComplementStr.size() - 3, static_cast<size_t >(1), 0, 0);
                    headerWriter.writeData(buffer, bufferLen, key, thread_idx);
                    break;
                case Orf::FRAME_3:
                    sequenceWriter.writeData(reverseComplementStr.c_str()+2, reverseComplementStr.size()-2, key, thread_idx);
                    bufferLen = Orf::writeOrfHeader(buffer, key, reverseComplementStr.size() - 4, static_cast<size_t >(2), 0, 0);
                    headerWriter.writeData(buffer, bufferLen, key, thread_idx);
                    break;
            }
            reverseComplementStr.clear();
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

