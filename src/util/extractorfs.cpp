#include <unistd.h>
#include <climits>
#include <algorithm>
#include <vector>

#include "Debug.h"
#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"

#include "Orf.h"

#ifdef OPENMP
#include <omp.h>

#endif

unsigned int getFrames(std::string frames) {
    unsigned int result = 0;

    std::vector<std::string> frame = Util::split(frames, ",");

    if(std::find(frame.begin(), frame.end(), "1") != frame.end()) {
        result |= Orf::FRAME_1;
    }

    if(std::find(frame.begin(), frame.end(), "2") != frame.end()) {
        result |= Orf::FRAME_2;
    }

    if(std::find(frame.begin(), frame.end(), "3") != frame.end()) {
        result |= Orf::FRAME_3;
    }

    return result;
}

int extractorfs(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str());
    reader.open(DBReader<unsigned int>::NOSORT);

    std::string headerIn(par.db1);
    headerIn.append("_h");

    std::string headerInIndex(par.db1);
    headerInIndex.append("_h.index");

    DBReader<unsigned int> headerReader(headerIn.c_str(), headerInIndex.c_str());
    headerReader.open(DBReader<unsigned int>::NOSORT);

    DBWriter sequenceWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads);
    sequenceWriter.open();

    std::string headerOut(par.db2);
    headerOut.append("_h");

    std::string headerIndexOut(par.db2);
    headerIndexOut.append("_h.index");

    DBWriter headerWriter(headerOut.c_str(), headerIndexOut.c_str(), par.threads);
    headerWriter.open();

    unsigned int forwardFrames = getFrames(par.forwardFrames);
    unsigned int reverseFrames = getFrames(par.reverseFrames);

    unsigned int extendMode = 0;
    if(par.orfLongest)
        extendMode |= Orf::EXTEND_START;

    if(par.orfExtendMin)
        extendMode |= Orf::EXTEND_END;

    size_t total = 0;
    #pragma omp parallel
    {
        Orf orf;

        #pragma omp for schedule(static)
        for (unsigned int i = 0; i < reader.getSize(); ++i){
            Debug::printProgress(i);
            int thread_idx = 0;
            #ifdef OPENMP
            thread_idx = omp_get_thread_num();
            #endif

            std::string data(reader.getData(i));
            // remove newline in sequence
            data.erase(std::remove(data.begin(), data.end(), '\n'), data.end());

            if(!orf.setSequence(data.c_str(), data.length())) {
                Debug(Debug::WARNING) << "Invalid sequence with index " << i << "!\n";
                continue;
            }

            std::string header(headerReader.getData(i));
            // remove newline in header
            header.erase(std::remove(header.begin(), header.end(), '\n'), header.end());

            std::string headerTemplate;
            if (par.useHeader) {
                headerTemplate.assign(Util::parseFastaHeader(header));
                headerTemplate.append("_");
            }

            std::vector<Orf::SequenceLocation> res;
            orf.findAll(res, par.orfMinLength, par.orfMaxLength, par.orfMaxGaps, forwardFrames, reverseFrames, extendMode);
            size_t orfNum = 0;
            for (std::vector<Orf::SequenceLocation>::const_iterator it = res.begin(); it != res.end(); ++it) {
                Orf::SequenceLocation loc = *it;

                std::string id;
                if (par.useHeader) {
                    id = headerTemplate;
                    id.append(SSTR(orfNum));
                    orfNum++;
                } else {
                    #pragma omp critical
                    {
                        orfNum = total++;
                    }
                    id = SSTR(orfNum + par.identifierOffset);
                }

                if (id.length() >= 31) {
                    Debug(Debug::ERROR) << "Id: " << id << " is too long. Maximum of 32 characters are allowed.\n";
                    EXIT(EXIT_FAILURE);
                }

                if (par.orfSkipIncomplete && (loc.hasIncompleteStart || loc.hasIncompleteEnd))
                    continue;

                char buffer[LINE_MAX];
                snprintf(buffer, LINE_MAX, "%s [Orf: %zu, %zu, %d, %d, %d]\n", header.c_str(), loc.from, loc.to, loc.strand, loc.hasIncompleteStart, loc.hasIncompleteEnd);

                headerWriter.writeData(buffer, strlen(buffer), id.c_str(), thread_idx);

                std::string sequence = orf.view(loc);
                sequence.append("\n");
                sequenceWriter.writeData(sequence.c_str(), sequence.length(), id.c_str(), thread_idx);
            }
        }
    }

    headerWriter.close();
    sequenceWriter.close();
    headerReader.close();
    reader.close();
    
    return EXIT_SUCCESS;
}
