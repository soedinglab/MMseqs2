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

void setMsa2ProfileDefaults(Parameters *p) {
    p->msaType = 1;

}

int msa2profile(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setMsa2ProfileDefaults(&par);
    par.parseParameters(argc, argv, command, 2);

    Debug(Debug::INFO) << "Compute profiles from MSAs.\n";
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

    unsigned int maxMsaSize = 0;
    unsigned int *msaSizes = qDbr.getSeqLens();
#pragma omp parallel for schedule(static) reduction(max:maxMsaSize)
    for (size_t id = 0; id < qDbr.getSize(); id++) {
        if (msaSizes[id] > maxMsaSize) {
            maxMsaSize = msaSizes[id];
        }
    }

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

#pragma omp parallel
    {
        SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0f, -0.2f);
        PSSMCalculator calculator(&subMat, par.maxSeqLen, par.pca, par.pcb);
        Sequence sequence(par.maxSeqLen, subMat.aa2int, subMat.int2aa,
                          Sequence::AMINO_ACIDS, 0, false, par.compBiasCorrection != 0);

        char *msaContent = (char*) mem_align(ALIGN_INT, sizeof(char) * maxMsaSize);
        size_t msaPos = 0;

        char *seqBuffer = new char[par.maxSeqLen + 1];
        bool *maskedColumns = new bool[par.maxSeqLen];

#pragma omp for schedule(static)
        for (size_t id = 0; id < qDbr.getSize(); id++) {
            Debug::printProgress(id);
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif

            unsigned int queryKey = qDbr.getDbKey(id);

            size_t setSize = 0;
            size_t centerLengthWithGaps = 0;
            size_t maskedCount = 0;

            size_t entryLength = qDbr.getSeqLens(id);

            bool fastaError = false;
            std::vector<char *> table;

            char *entryData = qDbr.getData(id);
            std::string msa;
            if (par.msaType == 0) {
                msa = CompressedA3M::extractA3M(entryData, entryLength, *sequenceReader, *headerReader);
            } else if (par.msaType == 1 || par.msaType == 2) {
                msa = std::string(entryData, entryLength);
            }

            kseq_buffer_t d((char*)msa.c_str(), msa.length());
            kseq_t *seq = kseq_init(&d);
            while (kseq_read(seq) >= 0) {
                if (seq->name.l == 0 || seq->seq.l == 0) {
                    Debug(Debug::WARNING) << "Invalid fasta entry " << setSize << "!\n";
                    fastaError = true;
                    break;
                }

                if (setSize == 0) {
                    size_t length = strlen(seq->seq.s);
                    for (size_t i = 0; i < length; ++i) {
                        if (seq->seq.s[i] == '-') {
                            maskedColumns[i] = true;
                            maskedCount++;
                        } else {
                            maskedColumns[i] = false;
                        }
                    }

                    centerLengthWithGaps = length;

                    if (centerLengthWithGaps > par.maxSeqLen) {
                            Debug(Debug::WARNING)
                                    << "MSA " << id << " is too long and will be cut short!"
                                    << "Please adjust --max-seq-len.\n";
                        centerLengthWithGaps = par.maxSeqLen;
                    }

                    std::string header(seq->name.s);
                    if (seq->comment.l > 0) {
                        header.append(seq->comment.s);
                    }
                    header.append("\n");
                    headerWriter.writeData(header.c_str(), header.size(), queryKey, thread_idx);
                }

                sequence.mapSequence(0, 0, seq->seq.s);
                table.push_back(msaContent + msaPos);

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

                while(msaPos % ALIGN_INT != 0) {
                    msaContent[msaPos++] = MultipleAlignment::GAP;
                }

                if (setSize == 0) {
                    seqBuffer[seqPos++] = '\n';
                    sequenceWriter.writeData(seqBuffer, seqPos, queryKey, thread_idx);
                }

                setSize++;
            }
            kseq_destroy(seq);

            if (fastaError == true || setSize == 0) {
                Debug(Debug::WARNING) << "Invalid msa " << id << "!\n";
                continue;
            }

            char **msaSequences = (char**) mem_align(ALIGN_INT, sizeof(char*) * table.size());
            size_t cnt = 0;
            for (std::vector<char*>::const_iterator it = table.begin(); it != table.end(); ++it) {
                msaSequences[cnt++] = (*it);
            }

            size_t centerLength = centerLengthWithGaps - maskedCount;

            if (par.filterMsa == true) {
                MsaFilter filter(centerLength, setSize, &subMat);
                MsaFilter::MsaFilterResult filterRes
                        = filter.filter((const char **) msaSequences,
                                        setSize,
                                        centerLength,
                                        static_cast<int>(par.cov * 100),
                                        static_cast<int>(par.qid * 100),
                                        par.qsc,
                                        static_cast<int>(par.filterMaxSeqId * 100),
                                        par.Ndiff);

                setSize = filterRes.setSize;
                for (size_t i = 0; i < filterRes.setSize; i++) {
                    msaSequences[i] = (char *) filterRes.filteredMsaSequence[i];
                }
            }

            std::pair<const char *, std::string> pssmRes =
                    calculator.computePSSMFromMSA(setSize, centerLength,
                                                  (const char **) msaSequences, par.wg);

            free(msaSequences);

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
        free(msaContent);
    }

    consensusWriter.close();
    headerWriter.close();
    sequenceWriter.close();
    resultWriter.close();

    FileUtil::symlinkAlias(par.db2 + "_seq_h", par.db2 + "_h");
    FileUtil::symlinkAlias(par.db2 + "_seq_h.index", par.db2 + "_h.index");

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
