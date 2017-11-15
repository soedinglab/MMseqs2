#include "MsaFilter.h"
#include "Parameters.h"
#include "PSSMCalculator.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "CompressedA3M.h"

#include "kseq.h"
#include "kseq_buffer_reader.h"

KSEQ_INIT(kseq_buffer_t*, kseq_buffer_reader)

#ifdef OPENMP
#include <omp.h>
#endif

#include <libgen.h>

void setMsa2ProfileDefaults(Parameters *p) {
    p->msaType = 1;
}

int msa2profile(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setMsa2ProfileDefaults(&par);
    par.parseParameters(argc, argv, command, 2);

    struct timeval start, end;
    gettimeofday(&start, NULL);

    std::string msaData = par.db1;
    std::string msaIndex = par.db1Index;
    DBReader<unsigned int> *headerReader = NULL, *sequenceReader = NULL;
    if (par.msaType == 0) {
        msaData = par.db1 + "_ca3m.ffdata";
        msaIndex = par.db1 + "_ca3m.ffindex";

        std::string msaHeaderData = par.db1 + "_header.ffdata";
        std::string msaHeaderIndex = par.db1 + "_header.ffindex";
        std::string msaSequenceData = par.db1 + "_sequence.ffdata";
        std::string msaSequenceIndex = par.db1 + "_sequence.ffindex";

        headerReader = new DBReader<unsigned int>(msaHeaderData.c_str(), msaHeaderIndex.c_str());
        headerReader->open(DBReader<unsigned int>::SORT_BY_LINE);

        sequenceReader = new DBReader<unsigned int>(msaSequenceData.c_str(), msaSequenceIndex.c_str());
        sequenceReader->open(DBReader<unsigned int>::SORT_BY_LINE);
    }

    DBReader<unsigned int> qDbr(msaData.c_str(), msaIndex.c_str());
    qDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    Debug(Debug::INFO) << "Finding maximum sequence length and set size.\n";
    unsigned int maxSeqLength = 0;
    unsigned int maxSetSize = 0;
    unsigned int *msaSizes = qDbr.getSeqLens();

#pragma omp parallel for schedule(static) reduction(max:maxSeqLength, maxSetSize)  
    for (size_t id = 0; id < qDbr.getSize(); id++) {
        bool inHeader = false;
        unsigned int setSize = 0;
        unsigned int seqLength = 0;

        char *entryData = qDbr.getData(id);
        for (size_t i = 0; i < msaSizes[id]; ++i) {
            // state machine to get the max sequence length and set size from MSA
            switch (entryData[i]) {
                case '>':
                    if (seqLength > maxSeqLength) {
                        maxSeqLength = seqLength;
                    }
                    seqLength = 0;
                    inHeader = true;
                    setSize++;
                    break;
                case '\n':
                    if (inHeader) {
                        inHeader = false;
                    }
                    break;
                default:
                    if (!inHeader) {
                        seqLength++;
                    }
                    break;
            }
        }

        // don't forget the last entry in an MSA
        if (!inHeader && seqLength > 0) {
            if (seqLength > maxSeqLength) {
                maxSeqLength = seqLength;
            }
            setSize++;
        }

        if (setSize > maxSetSize) {
            maxSetSize = setSize;
        }
    }

    // for SIMD memory alignment
    maxSeqLength =  maxSeqLength    / (VECSIZE_INT*4);
    maxSeqLength = (maxSeqLength+1) * (VECSIZE_INT*4);

    unsigned int threads = (unsigned int) par.threads;
    DBWriter resultWriter(par.db2.c_str(), par.db2Index.c_str(), threads, DBWriter::BINARY_MODE);
    resultWriter.open();

    DBWriter sequenceWriter(std::string(par.db2 + "_seq").c_str(),
                            std::string(par.db2 + "_seq.index").c_str(),
                            threads);
    sequenceWriter.open();

    DBWriter headerWriter(std::string(par.db2 + "_seq_h").c_str(),
                          std::string(par.db2 + "_seq_h.index").c_str(),
                          threads);
    headerWriter.open();

    DBWriter consensusWriter(std::string(par.db2 + "_consensus").c_str(),
                             std::string(par.db2 + "_consensus.index").c_str(),
                             threads);
    consensusWriter.open();

    Debug(Debug::INFO) << "Compute profiles from MSAs.\n";
#pragma omp parallel
    {
        SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0f, -0.2f);
        PSSMCalculator calculator(&subMat, maxSeqLength, maxSetSize, par.pca, par.pcb);
        Sequence sequence(maxSeqLength + 1, subMat.aa2int, subMat.int2aa,
                          Sequence::AMINO_ACIDS, 0, false, par.compBiasCorrection != 0);

        char *msaContent = (char*) mem_align(ALIGN_INT, sizeof(char) * maxSeqLength * maxSetSize);

        char *seqBuffer = new char[maxSeqLength + 1];
        bool *maskedColumns = new bool[maxSeqLength];

        kseq_buffer_t d;
        kseq_t *seq = kseq_init(&d);

        char **msaSequences = (char**) mem_align(ALIGN_INT, sizeof(char*) * maxSetSize);

        MsaFilter filter(maxSeqLength, maxSetSize, &subMat);
#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < qDbr.getSize(); ++id) {
            Debug::printProgress(id);

            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif

            unsigned int queryKey = qDbr.getDbKey(id);

            size_t msaPos = 0;

            unsigned int setSize = 0;
            unsigned int centerLengthWithGaps = 0;
            unsigned int maskedCount = 0;

            bool fastaError = false;

            char *entryData = qDbr.getData(id);
            size_t entryLength = qDbr.getSeqLens(id);

            std::string msa;
            if (par.msaType == 0) {
                msa = CompressedA3M::extractA3M(entryData, entryLength, *sequenceReader, *headerReader);
                d.buffer = const_cast<char*>(msa.c_str());
                d.length = msa.length();
            } else {
                d.buffer = entryData;
                d.length = entryLength;
            }
            d.position = 0;

            while (kseq_read(seq) >= 0) {
                if (seq->name.l == 0 || seq->seq.l == 0) {
                    Debug(Debug::WARNING) << "Invalid fasta sequence "
                                          << setSize << " in entry " << id << "!\n";
                    fastaError = true;
                    break;
                }

                if (seq->seq.l > maxSeqLength) {
                    Debug(Debug::WARNING) << "Member sequence "
                                          << setSize << " in entry " << id << " too long!\n";
                    fastaError = true;
                    break;
                }

                // first sequence is always the query
                if (setSize == 0) {
                    centerLengthWithGaps = static_cast<unsigned int>(strlen(seq->seq.s));
                    for (size_t i = 0; i < centerLengthWithGaps; ++i) {
                        if (seq->seq.s[i] == '-') {
                            maskedColumns[i] = true;
                            maskedCount++;
                        } else {
                            maskedColumns[i] = false;
                        }
                    }

                    std::string header(seq->name.s);
                    if (seq->comment.l > 0) {
                        header.append(seq->comment.s);
                    }
                    header.append("\n");
                    headerWriter.writeData(header.c_str(), header.size(), queryKey, thread_idx);
                }

                sequence.mapSequence(0, 0, seq->seq.s);
                msaSequences[setSize] = msaContent + msaPos;

                size_t seqPos = 0;
                for (size_t i = 0; i < centerLengthWithGaps; ++i) {
                    if (maskedColumns[i] == true) {
                        continue;
                    }

                    // skip a3m lower letters
                    if (par.msaType == 1 && islower(seq->seq.s[i])) {
                        continue;
                    }

                    if (seq->seq.s[i] == '-'){
                        msaContent[msaPos++] = MultipleAlignment::GAP;
                    } else {
                        int aa = sequence.int_sequence[i];
                        msaContent[msaPos++] = aa;

                        if (setSize == 0) {
                            seqBuffer[seqPos++] = subMat.int2aa[aa];
                        }
                    }
                }

                // fill up the sequence buffer for the SIMD profile calculation
                size_t rowSize = msaPos / (VECSIZE_INT*4);
                rowSize = (rowSize+1) * (VECSIZE_INT*4);
                while(msaPos < rowSize) {
                    msaContent[msaPos++] = MultipleAlignment::GAP;
                }

                if (setSize == 0) {
                    seqBuffer[seqPos++] = '\n';
                    sequenceWriter.writeData(seqBuffer, seqPos, queryKey, thread_idx);
                }

                setSize++;
            }
            kseq_rewind(seq);


            if (fastaError == true) {
                Debug(Debug::WARNING) << "Invalid msa " << id << "! Skipping entry.\n";
                continue;
            }

            if (setSize == 0) {
                Debug(Debug::WARNING) << "Empty msa " << id << "! Skipping entry.\n";
                continue;
            }

            unsigned int centerLength = centerLengthWithGaps - maskedCount;

            size_t filteredSetSize = setSize;
            if (par.filterMsa == 1) {
                filter.filter(setSize, centerLength, static_cast<int>(par.cov * 100),
                              static_cast<int>(par.qid * 100), par.qsc,
                              static_cast<int>(par.filterMaxSeqId * 100), par.Ndiff,
                              (const char **) msaSequences, &filteredSetSize);
                filter.shuffleSequences((const char **) msaSequences, setSize);
            }

            std::pair<const char *, std::string> pssmRes =
                    calculator.computePSSMFromMSA(filteredSetSize, centerLength,
                                                  (const char **) msaSequences, par.wg);

            char *data = (char *) pssmRes.first;
            size_t dataSize = centerLength * Sequence::PROFILE_AA_SIZE * sizeof(char);
            for (size_t i = 0; i < dataSize; i++) {
                // Avoid a null byte result
                data[i] = data[i] ^ 0x80;
            }

            resultWriter.writeData(data, dataSize, queryKey, thread_idx);

            std::string consensusStr = pssmRes.second;
            consensusStr.push_back('\n');
            consensusWriter.writeData(consensusStr.c_str(), consensusStr.length(), queryKey, thread_idx);
        }
        kseq_destroy(seq);
        free(msaSequences);
        free(msaContent);

        delete[] seqBuffer;
        delete[] maskedColumns;
    }

    consensusWriter.close();
    headerWriter.close();
    sequenceWriter.close();
    resultWriter.close();

    char* path = strdup((par.db2 + "_h").c_str());
    std::string base = basename(path);
    FileUtil::symlinkAlias(par.db2 + "_seq_h", base);
    FileUtil::symlinkAlias(par.db2 + "_seq_h.index", base + ".index");
    free(path);

    qDbr.close();
    if (sequenceReader != NULL) {
        sequenceReader->close();
        delete sequenceReader;
    }

    if (headerReader != NULL) {
        headerReader->close();
        delete headerReader;
    }

    gettimeofday(&end, NULL);
    time_t sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "Time for processing: "
                       << (sec / 3600) << " h "
                       << (sec % 3600 / 60) << " m "
                       << (sec % 60) << "s\n";

    return EXIT_SUCCESS;
}
