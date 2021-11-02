#include "Matcher.h"
#include "MsaFilter.h"
#include "Parameters.h"
#include "PSSMCalculator.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "CompressedA3M.h"
#include "MathUtil.h"

#include "kseq.h"
#include "KSeqBufferReader.h"

KSEQ_INIT(kseq_buffer_t*, kseq_buffer_reader)

#ifdef OPENMP
#include <omp.h>
#endif

enum {
    MSA_CA3M = 0,
    MSA_A3M  = 1,
    MSA_STOCKHOLM = 2
};

void setMsa2ProfileDefaults(Parameters *p) {
    p->msaType = MSA_STOCKHOLM;
}

int msa2profile(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setMsa2ProfileDefaults(&par);
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_PROFILE);

    std::vector<std::string> qid_str_vec = Util::split(par.qid, ",");
    std::vector<int> qid_vec;
    for (size_t qid_idx = 0; qid_idx < qid_str_vec.size(); qid_idx++) {
        float qid_float = strtod(qid_str_vec[qid_idx].c_str(), NULL);
        qid_vec.push_back(static_cast<int>(qid_float*100));
    }
    std::sort(qid_vec.begin(), qid_vec.end());

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

        headerReader = new DBReader<unsigned int>(msaHeaderData.c_str(), msaHeaderIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        headerReader->open(DBReader<unsigned int>::SORT_BY_LINE);

        sequenceReader = new DBReader<unsigned int>(msaSequenceData.c_str(), msaSequenceIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        sequenceReader->open(DBReader<unsigned int>::SORT_BY_LINE);
    }

    unsigned int mode = DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA;
    std::string lookupFile = msaData + ".lookup";
    if (FileUtil::fileExists(lookupFile.c_str())) {
        mode |= DBReader<unsigned int>::USE_LOOKUP;
    }
    DBReader<unsigned int> qDbr(msaData.c_str(), msaIndex.c_str(), par.threads, mode);
    qDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    Debug(Debug::INFO) << "Finding maximum sequence length and set size.\n";
    unsigned int maxSeqLength = 0;
    unsigned int maxSetSize = 0;
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

#pragma omp for schedule(dynamic, 10) reduction(max:maxSeqLength, maxSetSize)
        for (size_t id = 0; id < qDbr.getSize(); id++) {
            bool inHeader = false;
            unsigned int setSize = 0;
            unsigned int seqLength = 0;

            char *entryData = qDbr.getData(id, thread_idx);
            for (size_t i = 0; i < qDbr.getEntryLen(id); ++i) {
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
    }

    // for SIMD memory alignment
    maxSeqLength = (maxSeqLength) / (VECSIZE_INT * 4) + 2;
    maxSeqLength *= (VECSIZE_INT * 4);

    unsigned int threads = (unsigned int) par.threads;
    int type = Parameters::DBTYPE_HMM_PROFILE;
    if (par.pcmode == Parameters::PCMODE_CONTEXT_SPECIFIC) {
        type = DBReader<unsigned int>::setExtendedDbtype(type, Parameters::DBTYPE_EXTENDED_CONTEXT_PSEUDO_COUNTS);
    }
    DBWriter resultWriter(par.db2.c_str(), par.db2Index.c_str(), threads, par.compressed, type);
    resultWriter.open();

    DBWriter headerWriter(par.hdr2.c_str(), par.hdr2Index.c_str(), threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    headerWriter.open();

    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0f, -0.2f);

    Debug::Progress progress(qDbr.getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

        PSSMCalculator calculator(&subMat, maxSeqLength + 1, maxSetSize, par.pcmode, par.pca, par.pcb, par.gapOpen.values.aminoacid(), par.gapPseudoCount);

        Sequence sequence(maxSeqLength + 1, Parameters::DBTYPE_AMINO_ACIDS, &subMat, 0, false, par.compBiasCorrection != 0);

        char *msaContent = (char*) mem_align(ALIGN_INT, sizeof(char) * (maxSeqLength + 1) * maxSetSize);

        float *seqWeight = new float[maxSetSize];
        float *pNullBuffer = new float[maxSeqLength + 1];
        bool *maskedColumns = new bool[maxSeqLength + 1];
        std::string result;
        result.reserve((par.maxSeqLen + 1) * Sequence::PROFILE_READIN_SIZE * sizeof(char));

        kseq_buffer_t d;
        kseq_t *seq = kseq_init(&d);

        char **msaSequences = (char**) mem_align(ALIGN_INT, sizeof(char*) * maxSetSize);
        std::vector<Matcher::result_t> alnResults;
        alnResults.reserve(maxSetSize);

        const bool maskByFirst = par.matchMode == 0;
        const float matchRatio = par.matchRatio;
        MsaFilter filter(maxSeqLength + 1, maxSetSize, &subMat, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid());

#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < qDbr.getSize(); ++id) {
            progress.updateProgress();

            unsigned int queryKey = qDbr.getDbKey(id);

            size_t msaPos = 0;

            unsigned int setSize = 0;
            unsigned int centerLengthWithGaps = 0;
            unsigned int maskedCount = 0;

            bool fastaError = false;

            char *entryData = qDbr.getData(id, thread_idx);
            size_t entryLength = qDbr.getEntryLen(id);

            std::string msa;
            std::string backtrace;
            if (par.msaType == MSA_CA3M) {
                msa = CompressedA3M::extractA3M(entryData, entryLength - 2, *sequenceReader, *headerReader, thread_idx);
                d.buffer = const_cast<char*>(msa.c_str());
                d.length = msa.length();
                par.msaType = MSA_A3M;
            } else {
                d.buffer = entryData;
                d.length = entryLength - 1;
            }
            d.position = 0;

            // remove comment line that makes kseq_read fail
            if (d.length) {
                if (d.buffer[0] == '#') {
                    size_t pos = 0;
                    while (pos < d.length && d.buffer[pos] != '\n') {
                        pos++;
                    }

                    if (pos < d.length) {
                        pos++;
                        d.buffer += pos;
                        d.length -= pos;
                    } else {
                        d.buffer += pos;
                        d.length = 0;
                    }
                }
            }

            // allow skipping first sequence in case of consensus, etc
            if (par.skipQuery == true) {
                kseq_read(seq);
            }

            while (kseq_read(seq) >= 0) {
                if (seq->name.l == 0 || seq->seq.l == 0) {
                    Debug(Debug::WARNING) << "Invalid fasta sequence " << setSize << " in entry " << queryKey << "\n";
                    fastaError = true;
                    break;
                }

                if (seq->seq.l > maxSeqLength) {
                    Debug(Debug::WARNING) << "Member sequence " << setSize << " in entry " << queryKey << " too long\n";
                    fastaError = true;
                    break;
                }

                if ((par.msaType == 0 || par.msaType == 1) && strncmp("ss_", seq->name.s, strlen("ss_")) == 0) {
                    continue;
                }

                // first sequence is always the query
                if (setSize == 0) {
                    centerLengthWithGaps = seq->seq.l;
                    backtrace.reserve(centerLengthWithGaps);
                    if (maskByFirst == true) {
                        for (size_t i = 0; i < centerLengthWithGaps; ++i) {
                            if (seq->seq.s[i] == '-') {
                                maskedColumns[i] = true;
                                maskedCount++;
                            } else {
                                maskedColumns[i] = false;
                            }
                        }
                    }
                    if ((mode & DBReader<unsigned int>::USE_LOOKUP) == 0) {
                        std::string header(seq->name.s);
                        if (seq->comment.l > 0) {
                            header.append(" ");
                            header.append(seq->comment.s);
                        }
                        header.append("\n");
                        headerWriter.writeData(header.c_str(), header.size(), queryKey, thread_idx);
                    }
                }

                sequence.mapSequence(0, 0, seq->seq.s, seq->seq.l);
                msaSequences[setSize] = msaContent + msaPos;

                for (size_t i = 0; i < centerLengthWithGaps; ++i) {
                    if (maskByFirst == true && maskedColumns[i] == true) {
                        continue;
                    }

                    // skip a3m lower letters
                    if (par.msaType == MSA_A3M && islower(seq->seq.s[i])) {
                        continue;
                    }

                    msaContent[msaPos++] = (seq->seq.s[i] == '-') ? (int)MultipleAlignment::GAP : sequence.numSequence[i];
                }

                // construct backtrace for all but the query sequence
                if (false && setSize > 0) {
                    backtrace.clear();
                    for (size_t i = 0; i < centerLengthWithGaps; ++i) {
                        bool isMaskedColumn = (maskByFirst && maskedColumns[i]);
                        if (seq->seq.s[i] == '-' && isMaskedColumn) {
                            continue;
                        }
                        if (seq->seq.s[i] == '-') {
                            backtrace.push_back('I');
                        }
                        else if (isMaskedColumn || (!maskByFirst && msaSequences[0][i] == MultipleAlignment::GAP)
                                                || (par.msaType == MSA_A3M && islower(seq->seq.s[i]))) {
                            backtrace.push_back('D');
                        }
                        else {
                            backtrace.push_back('M');
                        }
                    }
                    alnResults.emplace_back(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, backtrace);
                }

                // fill up the sequence buffer for the SIMD profile calculation
                size_t rowSize = msaPos / (VECSIZE_INT*4);
                rowSize = (rowSize+1) * (VECSIZE_INT*4);
                while(msaPos < rowSize) {
                    msaContent[msaPos++] = MultipleAlignment::GAP;
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

            if (maskByFirst == false) {
                PSSMCalculator::computeSequenceWeights(seqWeight, centerLengthWithGaps,
                                                       setSize, const_cast<const char**>(msaSequences));

                // Replace GAP with ENDGAP for all end gaps
                // ENDGAPs are ignored for counting percentage (multi-domain proteins)
                for (unsigned int k = 0; k < setSize; ++k) {
                    for (unsigned int i = 0; i < centerLengthWithGaps && msaSequences[k][i] == MultipleAlignment::GAP; ++i)
                        msaSequences[k][i] = MultipleAlignment::ENDGAP;
                    for (unsigned int i = centerLengthWithGaps - 1; msaSequences[k][i] == MultipleAlignment::GAP; i--)
                        msaSequences[k][i] = MultipleAlignment::ENDGAP;
                }

                for (unsigned int l = 0; l < centerLengthWithGaps; l++) {
                    float res = 0;
                    float gap = 0;
                    // Add up percentage of gaps
                    for (unsigned int k = 0; k < setSize; ++k) {
                        if (msaSequences[k][l] < MultipleAlignment::GAP) {
                            res += seqWeight[k];
                        } else if (msaSequences[k][l] != MultipleAlignment::ENDGAP) {
                            gap += seqWeight[k];
                        } else if (msaSequences[k][l] == MultipleAlignment::ENDGAP) {
                            msaSequences[k][l] = MultipleAlignment::GAP;
                        }
                    }

                    maskedColumns[l] =  (gap / (res + gap)) > matchRatio;
                    maskedCount += maskedColumns[l] ? 1 : 0;
                }

                for (unsigned int k = 0; k < setSize; ++k) {
                    unsigned int currentCol = 0;
                    for (unsigned int l = 0; l < centerLengthWithGaps; ++l) {
                        if (maskedColumns[l] == false) {
                            msaSequences[k][currentCol++] = msaSequences[k][l];
                        }
                    }

                    for (unsigned int l = currentCol; l < centerLengthWithGaps; ++l) {
                        msaSequences[k][l] = MultipleAlignment::GAP;
                    }
                }

                /* update backtraces
                // TODO: check if this works for a3m as well (probably not...)
                for (unsigned int k = 0; k < setSize - 1; ++k) {
                    std::string::iterator readIt = alnResults[k].backtrace.begin();
                    std::string::iterator writeIt = readIt;
                    for (unsigned int l = 0; l < centerLengthWithGaps; ++l) {
                        if (!maskedColumns[l]) {
                            *writeIt = *readIt;
                            ++readIt;
                            ++writeIt;
                        } else {
                            if (*readIt == 'D') {
                                *writeIt = 'D';
                                ++writeIt;
                            }
                            ++readIt;
                        }
                    }
                    alnResults[k].backtrace.erase(writeIt, alnResults[k].backtrace.end());
                } */
            }
            unsigned int centerLength = centerLengthWithGaps - maskedCount;

            MultipleAlignment::MSAResult msaResult(centerLength, centerLength, setSize, msaSequences);
            size_t filteredSetSize = setSize;
            if (par.filterMsa == 1) {
                filteredSetSize = filter.filter(setSize, centerLength, static_cast<int>(par.covMSAThr * 100),
                                                qid_vec, par.qsc,
                                                static_cast<int>(par.filterMaxSeqId * 100), par.Ndiff, par.filterMinEnable,
                                                (const char **) msaSequences, true);
            }

            PSSMCalculator::Profile pssmRes =
                    calculator.computePSSMFromMSA(filteredSetSize, msaResult.centerLength,
                                                  (const char **) msaResult.msaSequence, alnResults, par.wg);
            if (par.compBiasCorrection == true) {
                SubstitutionMatrix::calcGlobalAaBiasCorrection(&subMat, pssmRes.pssm, pNullBuffer,
                                                               Sequence::PROFILE_AA_SIZE,
                                                               centerLength);
            }
            pssmRes.toBuffer((const unsigned char*)msaSequences[0], centerLength, subMat, result);

            if (mode & DBReader<unsigned int>::USE_LOOKUP) {
                size_t lookupId = qDbr.getLookupIdByKey(queryKey);
                std::string header = qDbr.getLookupEntryName(lookupId);
                header.append(1, '\n');
                headerWriter.writeData(header.c_str(), header.length(), queryKey, thread_idx);
            }
            resultWriter.writeData(result.c_str(), result.length(), queryKey, thread_idx);
            result.clear();
        }
        kseq_destroy(seq);
        free(msaSequences);
        free(msaContent);
        delete[] pNullBuffer;
        delete[] maskedColumns;
        delete[] seqWeight;
    }
    headerWriter.close(true);
    resultWriter.close(true);
    qDbr.close();

    DBReader<unsigned int>::copyDb(par.db1, par.db2, (DBFiles::Files)(DBFiles::LOOKUP | DBFiles::SOURCE));

    if (sequenceReader != NULL) {
        sequenceReader->close();
        delete sequenceReader;
    }

    if (headerReader != NULL) {
        headerReader->close();
        delete headerReader;
    }

    return EXIT_SUCCESS;
}
