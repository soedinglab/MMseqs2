#include "MsaFilter.h"
#include "Parameters.h"
#include "PSSMCalculator.h"
#include "DBReader.h"
#include "DBWriter.h"

#include "kseq.h"
#include "kseq_buffer_reader.h"

KSEQ_INIT(kseq_buffer_t*, kseq_buffer_reader)

#ifdef OPENMP
#include <omp.h>
#endif

int msa2profile(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);

    Debug(Debug::INFO) << "Compute profiles from MSAs.\n";
    struct timeval start, end;
    gettimeofday(&start, NULL);

    DBReader<unsigned int> qDbr(par.db1.c_str(), par.db1Index.c_str());
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

        char *msaContent = new char[maxMsaSize];
        size_t msaPos = 0;

        char *seqBuffer = new char[par.maxSeqLen];
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
            size_t centerLength = 0;
            size_t maskedCount = 0;

            size_t entrySize = qDbr.getSeqLens(id);

            bool fastaError = false;
            std::vector<char *> table;
            kseq_buffer_t d(qDbr.getData(id), entrySize);
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

                    centerLength = length;

                    if (centerLength > par.maxSeqLen) {
                            Debug(Debug::WARNING)
                                    << "MSA " << id << " is too long and will be cut short!"
                                    << "Please adjust --max-seq-len.\n";
                        centerLength = par.maxSeqLen;
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
                for (size_t i = 0; i < centerLength; ++i) {
                    if (maskedColumns[i] == true) {
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
                msaContent[msaPos++] = 0;

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

#if __cplusplus <= 199711L
            // Before C++11 this memory might be in any order
            // So we actually have to copy it around
            char **msaSequences = new char *[table.size()];
            size_t cnt = 0;
            for (std::vector<char*>::const_iterator it = table.begin(); it != table.end(); ++it) {
                msaSequences[cnt++] = (*it);
            }
#else
            // C++11 guarantees that the internal memory is continuous
            char **msaSequences = table.data();
#endif

            std::pair<const char *, std::string> pssmRes =
                    calculator.computePSSMFromMSA(setSize, (centerLength - maskedCount),
                                                  (const char **) msaSequences, par.wg);

            char *data = (char *) pssmRes.first;
            size_t dataSize = (centerLength - maskedCount) * Sequence::PROFILE_AA_SIZE * sizeof(char);
            for (size_t i = 0; i < dataSize; i++) {
                // Avoid a null byte result
                data[i] = data[i] ^ 0x80;
            }

            resultWriter.writeData(data, dataSize, queryKey, thread_idx);

            std::string consensusStr = pssmRes.second;
            consensusStr.push_back('\n');
            consensusWriter.writeData(consensusStr.c_str(), consensusStr.length(), queryKey, thread_idx);
        }
        delete msaContent;
    }

    consensusWriter.close();
    headerWriter.close();
    sequenceWriter.close();
    resultWriter.close();

    qDbr.close();
    gettimeofday(&end, NULL);
    time_t sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "Time for processing: "
                       << (sec / 3600) << " h "
                       << (sec % 3600 / 60) << " m "
                       << (sec % 60) << "s\n";

    return EXIT_SUCCESS;
}
