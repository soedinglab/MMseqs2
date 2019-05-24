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


int extractorfs(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> headerReader(par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    headerReader.open(DBReader<unsigned int>::NOSORT);

    DBWriter sequenceWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_NUCLEOTIDES);
    sequenceWriter.open();

    DBWriter headerWriter(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, false, Parameters::DBTYPE_OFFSETDB);
    headerWriter.open();

    if ((par.orfStartMode == 1) && (par.contigStartMode < 2)) {
        Debug(Debug::ERROR) << "Parameter combination is illegal, orf-start-mode 1 can only go with contig-start-mode 2\n";
        EXIT(EXIT_FAILURE);
    }

    unsigned int forwardFrames = Orf::getFrames(par.forwardFrames);
    unsigned int reverseFrames = Orf::getFrames(par.reverseFrames);
    const char newline = '\n';
    Debug::Progress progress(reader.getSize());

#pragma omp parallel
    {
        Orf orf(par.translationTable, par.useAllTableStarts);
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

        std::vector<Orf::SequenceLocation> res;
        res.reserve(1000);
        for (unsigned int i = queryFrom; i < (queryFrom + querySize); ++i){
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            const char* data = reader.getData(i, thread_idx);
            size_t dataLength = reader.getSeqLens(i);
            size_t sequenceLength = dataLength - 2;
            if(!orf.setSequence(data, dataLength - 2)) {
                Debug(Debug::WARNING) << "Invalid sequence with index " << i << "!\n";
                continue;
            }

            const char* header = headerReader.getData(i, thread_idx);
            std::string headerAccession = Util::parseFastaHeader(header);
            orf.findAll(res, par.orfMinLength, par.orfMaxLength, par.orfMaxGaps, forwardFrames, reverseFrames, par.orfStartMode);
            for (std::vector<Orf::SequenceLocation>::const_iterator it = res.begin(); it != res.end(); ++it) {
                Orf::SequenceLocation loc = *it;

                if (par.contigStartMode < 2 && (loc.hasIncompleteStart == par.contigStartMode)) {
                    continue;
                }
                if (par.contigEndMode   < 2 && (loc.hasIncompleteEnd   == par.contigEndMode)) {
                    continue;
                }

                char buffer[LINE_MAX];

                std::pair<const char*, size_t> sequence = orf.getSequence(loc);
                size_t fromPos = loc.from;
                size_t toPos = loc.to;
                if(loc.strand == Orf::STRAND_MINUS){
                    fromPos = (sequenceLength - 1) - loc.from;
                    toPos   = (sequenceLength - 1) - loc.to;
                }
                Orf::writeOrfHeader(buffer, key, fromPos, toPos, loc.hasIncompleteStart, loc.hasIncompleteEnd);

//                snprintf(buffer, LINE_MAX, "%.*s [Orf: %d, %zu, %zu, %d, %d]\n", (unsigned int)(headerAccession.size()), headerAccession.c_str(),
//                          toPos, loc.hasIncompleteStart, loc.hasIncompleteEnd);
                sequenceWriter.writeStart(thread_idx);
                sequenceWriter.writeAdd(sequence.first, sequence.second, thread_idx);
                sequenceWriter.writeAdd(&newline, 1, thread_idx);
                sequenceWriter.writeEnd(key, thread_idx);

                headerWriter.writeData(buffer, strlen(buffer), key, thread_idx);


            }
            res.clear();
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
                DBReader<unsigned int> orfHeaderReader(par.hdr2.c_str(), par.hdr2Index.c_str(),
                                                       par.threads, DBReader<unsigned int>::USE_INDEX);
                orfHeaderReader.open(DBReader<unsigned int>::SORT_BY_ID_OFFSET);
                FILE *hIndex = fopen((par.hdr2Index + "_tmp").c_str(), "w");
                if (hIndex == NULL) {
                    Debug(Debug::ERROR) << "Could not open " << par.hdr2Index << "_tmp for writing!\n";
                    EXIT(EXIT_FAILURE);
                }
                for (size_t i = 0; i < orfHeaderReader.getSize(); i++) {
                    DBReader<unsigned int>::Index *idx = orfHeaderReader.getIndex(i);
                    char buffer[1024];
                    size_t len = DBWriter::indexToBuffer(buffer, i, idx->offset, orfHeaderReader.getSeqLens(i));
                    int written = fwrite(buffer, sizeof(char), len, hIndex);
                    if (written != (int) len) {
                        Debug(Debug::ERROR) << "Could not write to data file " << par.hdr2Index << "_tmp\n";
                        EXIT(EXIT_FAILURE);
                    }
                }
                fclose(hIndex);
                orfHeaderReader.close();
                std::rename((par.hdr2Index + "_tmp").c_str(), par.hdr2Index.c_str());
            }

#pragma omp task
            {
                DBReader<unsigned int> orfSequenceReader(par.db2.c_str(), par.db2Index.c_str(),
                                                         par.threads, DBReader<unsigned int>::USE_INDEX);
                orfSequenceReader.open(DBReader<unsigned int>::SORT_BY_ID_OFFSET);

                FILE *sIndex = fopen((par.db2Index + "_tmp").c_str(), "w");
                if (sIndex == NULL) {
                    Debug(Debug::ERROR) << "Could not open " << par.db2Index << "_tmp for writing!\n";
                    EXIT(EXIT_FAILURE);
                }

                for (size_t i = 0; i < orfSequenceReader.getSize(); i++) {
                    DBReader<unsigned int>::Index *idx = (orfSequenceReader.getIndex(i));
                    char buffer[1024];
                    size_t len = DBWriter::indexToBuffer(buffer, i, idx->offset, orfSequenceReader.getSeqLens(i));
                    int written = fwrite(buffer, sizeof(char), len, sIndex);
                    if (written != (int) len) {
                        Debug(Debug::ERROR) << "Could not write to data file " << par.db2Index << "_tmp\n";
                        EXIT(EXIT_FAILURE);
                    }
                }
                fclose(sIndex);
                orfSequenceReader.close();
                std::rename((par.db2Index + "_tmp").c_str(), par.db2Index.c_str());
            }
        }
    }
    return EXIT_SUCCESS;
}

