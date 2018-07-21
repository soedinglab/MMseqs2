#include <unistd.h>
#include <climits>
#include <algorithm>

#include "Debug.h"
#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Matcher.h"
#include "Util.h"
#include "itoa.h"

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

    DBReader<unsigned int> headerReader(par.hdr1.c_str(), par.hdr1Index.c_str());
    headerReader.open(DBReader<unsigned int>::NOSORT);

    DBWriter sequenceWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads);
    sequenceWriter.open();

    DBWriter headerWriter(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads);
    headerWriter.open();


    std::string dbSetToOrf = par.db2 + "_member_lookup";
    std::string dbIndexSetToOrf = par.db2 + "_member_lookup.index";
    DBWriter writerSetToOrf(dbSetToOrf.c_str(), dbIndexSetToOrf.c_str(), par.threads);
    writerSetToOrf.open();

    std::string dbOrfToSet = par.db2 + "_set_lookup";
    std::string dbIndexOrfToSet = par.db2 + "_set_lookup.index";
    DBWriter writerOrfToSet(dbOrfToSet.c_str(), dbIndexOrfToSet.c_str(), par.threads);
    writerOrfToSet.open();

    unsigned int forwardFrames = getFrames(par.forwardFrames);
    unsigned int reverseFrames = getFrames(par.reverseFrames);

    unsigned int total = 0;
#pragma omp parallel shared(total)
    {
        Orf orf(par.translationTable, par.useAllTableStarts);

#pragma omp for schedule(dynamic, 10)
        for (unsigned int i = 0; i < reader.getSize(); ++i){
            unsigned int id;
            Debug::printProgress(i);
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif

            std::string orfsBuffer;
            orfsBuffer.reserve(10000);

            unsigned int key = reader.getDbKey(i);
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

            std::vector<Orf::SequenceLocation> res;
            orf.findAll(res, par.orfMinLength, par.orfMaxLength, par.orfMaxGaps, forwardFrames, reverseFrames, par.orfStartMode);
            for (std::vector<Orf::SequenceLocation>::const_iterator it = res.begin(); it != res.end(); ++it) {
                Orf::SequenceLocation loc = *it;

                size_t offset = __sync_fetch_and_add(&total, 1);
                id = offset + par.identifierOffset;
                if (par.contigStartMode < 2 && (loc.hasIncompleteStart == par.contigStartMode)) {
                    continue;
                }
                if (par.contigEndMode   < 2 && (loc.hasIncompleteEnd   == par.contigEndMode)) {
                    continue;
                }

                char buffer[LINE_MAX];
                snprintf(buffer, LINE_MAX, "%s [Orf: %d, %zu, %zu, %d, %d, %d]\n", header.c_str(), key, loc.from, loc.to, loc.strand, loc.hasIncompleteStart, loc.hasIncompleteEnd);

                headerWriter.writeData(buffer, strlen(buffer), id, thread_idx);

                std::string sequence = orf.view(loc);
                sequence.append("\n");
//            if(loc.hasIncompleteStart == false){
//                sequence.insert (0, "TAG");
//            }
                sequenceWriter.writeData(sequence.c_str(), sequence.length(), id, thread_idx);


                // write an alignemnt-like record
                // Matcher::result_t(targetId, score, qCov, dbCov, seqId, eval, alnLength, qStart, qEnd, qLen, dbStart, dbEnd, dbLen, "");

                // the orf length is like its string minus 1 ('\n'):
                size_t orfLen = sequence.length() - 1;
                size_t setLen = data.length();

                // orf search returns the position right after the orf, this keeps consitency with alignemnt format
                // if strand == -1 (reverse), we need to recompute the coordinates with respect to the positive strand and swap them
                size_t setFromWithStrand = (loc.strand > 0) ? loc.from : (setLen - loc.from - 1);
                size_t setToWithStrand = (loc.strand > 0) ? (loc.to - 1) : (setLen - loc.to);

                Matcher::result_t orfToSetResult(key, 1, 1, 0, 1, 0, orfLen, 0, (orfLen - 1), orfLen, setFromWithStrand, setToWithStrand, setLen, "");
                char orfToSetBuffer[LINE_MAX];
                size_t len = Matcher::resultToBuffer(orfToSetBuffer, orfToSetResult, false);
                writerOrfToSet.writeData(orfToSetBuffer, len, id, thread_idx);

                Itoa::u32toa_sse2(static_cast<uint32_t>(id), buffer);
                orfsBuffer.append(buffer);
                orfsBuffer.append("\n");
            }
            writerSetToOrf.writeData(orfsBuffer.c_str(), orfsBuffer.length(), key, thread_idx);
        }
    }

    writerSetToOrf.close();
    writerOrfToSet.close();

    headerWriter.close();
    sequenceWriter.close(Sequence::NUCLEOTIDES);
    headerReader.close();
    reader.close();

    return EXIT_SUCCESS;
}

