#include "MsaFilter.h"
#include "Parameters.h"
#include "PSSMCalculator.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "CompressedA3M.h"
#include "MathUtil.h"
#include "AlignmentSymmetry.h"
#include "Matcher.h"

#include "kseq.h"
#include "KSeqBufferReader.h"

KSEQ_INIT(kseq_buffer_t*, kseq_buffer_reader)

#ifdef OPENMP
#include <omp.h>
#endif

void setMsa2ResultDefaults(Parameters *p) {
    p->msaType = 2;
}

int msa2result(int argc, const char **argv, const Command &command) {
//    return EXIT_FAILURE;
//}
    Parameters &par = Parameters::getInstance();
    setMsa2ResultDefaults(&par);
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
    DBReader<unsigned int> msaReader(msaData.c_str(), msaIndex.c_str(), par.threads, mode);
    msaReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    size_t maxMsaArea = 0;
    unsigned int maxSetSize = 0;
    unsigned int maxSeqLength = 0;
    unsigned int* setSizes = (unsigned int*)calloc((msaReader.getSize() + 1), sizeof(unsigned int));
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

#pragma omp for schedule(dynamic, 10) reduction(max:maxSeqLength, maxSetSize)
        for (size_t id = 0; id < msaReader.getSize(); ++id) {
            bool inHeader = false;
            unsigned int setSize = 0;
            unsigned int seqLength = 0;

            char *entryData = msaReader.getData(id, thread_idx);
            for (size_t i = 0; i < msaReader.getEntryLen(id); ++i) {
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
            setSizes[id] = setSize;
            if (setSize > maxSetSize) {
                maxSetSize = setSize;
            }
            size_t area = setSize * (seqLength + 1);
            if (area > maxMsaArea) {
                maxMsaArea = area;
            }
        }
    }
    AlignmentSymmetry::computeOffsetFromCounts(setSizes, msaReader.getSize());

    // for SIMD memory alignment
    maxSeqLength = (maxSeqLength) / (VECSIZE_INT * 4) + 2;
    maxSeqLength *= (VECSIZE_INT * 4);

    unsigned int threads = (unsigned int) par.threads;
    DBWriter sequenceWriter(par.db2.c_str(), par.db2Index.c_str(), threads, par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
    sequenceWriter.open();

    DBWriter headerWriter(par.hdr2.c_str(), par.hdr2Index.c_str(), threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    headerWriter.open();

    DBWriter resultWriter(par.db3.c_str(), par.db3Index.c_str(), threads, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    resultWriter.open();

    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0f, -0.2f);
    SubstitutionMatrix::FastMatrix fastMatrix = SubstitutionMatrix::createAsciiSubMat(subMat);
    // we need an evaluer to compute a normalized bitScore.
    // the validity of the evalue is questionable since we give the number of entries in the MSA db
    // and not the total number of residues
    // also, since not all against all comparisons were carried out, it is not straight-forward to
    // decide what to correct for.
    EvalueComputation evaluer(msaReader.getSize(), &subMat, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid());

    Debug::Progress progress(msaReader.getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif


        PSSMCalculator calculator(
            &subMat, maxSeqLength + 1, maxSetSize, par.pcmode, par.pca, par.pcb
#ifdef GAP_POS_SCORING
            , par.gapOpen.values.aminoacid()
            , par.gapPseudoCount
#endif
        );
        Sequence sequence(maxSeqLength + 1, Parameters::DBTYPE_AMINO_ACIDS, &subMat, 0, false, par.compBiasCorrection != 0);

        char *msaContent = (char*) mem_align(ALIGN_INT, sizeof(char) * (maxSeqLength + 1) * maxSetSize);

        float *seqWeight = new float[maxSetSize];
        char *maskedColumns = new char[maxSeqLength + 1];
        char *seqBuffer = new char[maxSeqLength + 1];

        std::string result;
        result.reserve((par.maxSeqLen + 1) * Sequence::PROFILE_READIN_SIZE * sizeof(char));

        kseq_buffer_t d;
        kseq_t *seq = kseq_init(&d);

        char **msaSequences = (char**) mem_align(ALIGN_INT, sizeof(char*) * maxSetSize);

        const bool maskByFirst = par.matchMode == 0;
        const float matchRatio = par.matchRatio;
        MsaFilter filter(maxSeqLength + 1, maxSetSize, &subMat, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid());

        char buffer[1024 + 32768*4];

#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < msaReader.getSize(); ++id) {
            progress.updateProgress();

            unsigned int queryKey = msaReader.getDbKey(id);

            size_t msaPos = 0;

            unsigned int setSize = 0;
            unsigned int centerLengthWithGaps = 0;
            unsigned int maskedCount = 0;

            bool fastaError = false;

            char *entryData = msaReader.getData(id, thread_idx);
            size_t entryLength = msaReader.getEntryLen(id);

            std::string msa;
            if (par.msaType == 0) {
                msa = CompressedA3M::extractA3M(entryData, entryLength - 2, *sequenceReader, *headerReader, thread_idx);
                d.buffer = const_cast<char*>(msa.c_str());
                d.length = msa.length();
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

            unsigned int startKey = setSizes[id];
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

                const char newline = '\n';
                headerWriter.writeStart(thread_idx);
                headerWriter.writeAdd(seq->name.s, seq->name.l, thread_idx);
                if (seq->comment.l > 0) {
                    const char space = ' ';
                    headerWriter.writeAdd(&space, 1, thread_idx);
                    headerWriter.writeAdd(seq->comment.s, seq->comment.l, thread_idx);
                }
                headerWriter.writeAdd(&newline, 1, thread_idx);
                headerWriter.writeEnd(startKey, thread_idx);
                sequenceWriter.writeStart(thread_idx);
                for (size_t i = 0; i < seq->seq.l; ++i) {
                    if (seq->seq.s[i] == '-') {
                        continue;
                    }
                    sequenceWriter.writeAdd(&(seq->seq.s[i]), 1, thread_idx);
                }
                sequenceWriter.writeAdd(&newline, 1, thread_idx);
                sequenceWriter.writeEnd(startKey, thread_idx);
                startKey++;

                // first sequence is always the query
                if (setSize == 0) {
                    centerLengthWithGaps = seq->seq.l;
//                    if (maskByFirst == true) {
//                        for (size_t i = 0; i < centerLengthWithGaps; ++i) {
//                            if (seq->seq.s[i] == '-') {
//                                maskedColumns[i] = true;
//                                maskedCount++;
//                            } else {
//                                maskedColumns[i] = false;
//                            }
//                        }
//                    }
                }

                sequence.mapSequence(0, 0, seq->seq.s, seq->seq.l);
                msaSequences[setSize] = msaContent + msaPos;

                for (size_t i = 0; i < centerLengthWithGaps; ++i) {
                    if (maskByFirst == true && maskedColumns[i] == 1) {
                        continue;
                    }

                    // skip a3m lower letters
                    if (par.msaType == 1 && islower(seq->seq.s[i])) {
                        continue;
                    }

                    msaContent[msaPos++] = (seq->seq.s[i] == '-') ? (int)MultipleAlignment::GAP : sequence.numSequence[i];
                }

                // fill up the sequence buffer for the SIMD profile calculation
                size_t rowSize = msaPos / (VECSIZE_INT*4);
                rowSize = (rowSize+1) * (VECSIZE_INT*4);
                while (msaPos < rowSize) {
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
                PSSMCalculator::computeSequenceWeights(seqWeight, centerLengthWithGaps, setSize, const_cast<const char**>(msaSequences));

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

                    maskedColumns[l] = ((gap / (res + gap)) > matchRatio) ? 1 : 0;
                    maskedCount += maskedColumns[l];
                }

                for (unsigned int k = 0; k < setSize; ++k) {
                    unsigned int currentCol = 0;
                    unsigned int currentMask = 0;
                    for (unsigned int l = 0; l < centerLengthWithGaps; ++l) {
                        if (maskedColumns[l] == 0) {
                            msaSequences[k][currentCol++] = msaSequences[k][l];
                        } else {
                            seqBuffer[currentMask++] = msaSequences[k][l];
                        }
                    }

                    for (unsigned int l = currentCol; l < centerLengthWithGaps; ++l) {
                        msaSequences[k][l] = seqBuffer[l - currentCol];
                    }
                }
            }

            unsigned int centerLength = centerLengthWithGaps - maskedCount;
            size_t filteredSetSize = setSize;
            if (par.filterMsa == 1)  {
                filteredSetSize = filter.filter(setSize, centerLength, static_cast<int>(par.covMSAThr * 100),
                                                qid_vec, par.qsc,
                                                static_cast<int>(par.filterMaxSeqId * 100), par.Ndiff, par.filterMinEnable,
                                                (const char **) msaSequences, true);
            }
            PSSMCalculator::Profile pssmRes = calculator.computePSSMFromMSA(filteredSetSize, centerLength, (const char **) msaSequences, par.wg);

            resultWriter.writeStart(thread_idx);
            for (size_t i = 0; i < setSize; ++i) {
                const char* currSeq = msaSequences[i];
                unsigned int currentCol = 0;
                unsigned int currentMask = 0;
                std::string bt;
                std::string currSeqNoGaps;
                std::string consSeqNoGaps;
                size_t numIdentical = 0;
                for (size_t j = 0; j < centerLengthWithGaps; ++j) {
                    bool takeFromEnd = false;
                    if (maskedColumns[j] == 1) {
                        takeFromEnd = true;
                        currentMask++;
                    } else {
                        currentCol++;
                    }

                    char conRes = takeFromEnd ? '-' : pssmRes.consensus[currentCol - 1];
                    char seqRes = currSeq[takeFromEnd ? (centerLength + (currentMask - 1)) : currentCol - 1];
                    seqRes = seqRes == MultipleAlignment::GAP ? '-' : subMat.num2aa[(int)seqRes];

                    if (conRes == '-' && seqRes == '-') {
                        continue;
                    } else if (conRes != '-' && seqRes == '-') {
                        bt.append(1, 'I');
                        consSeqNoGaps.append(1, conRes);
                    } else if (conRes == '-' && seqRes != '-') {
                        bt.append(1, 'D');
                        currSeqNoGaps.append(1, seqRes);
                    } else if (conRes != '-' && seqRes != '-') {
                        bt.append(1, 'M');
                        currSeqNoGaps.append(1, seqRes);
                        consSeqNoGaps.append(1, conRes);
                    }

                    if (conRes == seqRes) {
                        numIdentical++;
                    }
                }
//                std::cout << gapCons << '\n';
//                std::cout << bt << '\n';
//                std::cout << gapSeq << '\n';
//                for(size_t pos = 0; pos < centerLength; pos++){
//                    int aa = currSeq[pos];
//                    std::cout << (char)((aa != MultipleAlignment::GAP) ? subMat.num2aa[(int)aa] : '-' );
//                }
//                std::cout << '\n';

/*
 *         result_t(unsigned int dbkey,int score,
                 float qcov, float dbcov,
                 float seqId, double eval,
                 unsigned int alnLength,
                 int qStartPos,
                 int qEndPos,
                 unsigned int qLen,
                 int dbStartPos,
                 int dbEndPos,
                 unsigned int dbLen,
                 std::string backtrace)
 */

                float seqId = (float)numIdentical / bt.length();

                unsigned int key = setSizes[id] + i;
                // initialize res with some values
                Matcher::result_t res(key, 0, 1.0, 1.0, seqId, 0, bt.length(), 0, consSeqNoGaps.size() - 1, consSeqNoGaps.size(), 0, currSeqNoGaps.size() - 1, currSeqNoGaps.size(), bt);

                // and update them and compute the score
                Matcher::updateResultByRescoringBacktrace(consSeqNoGaps.c_str(), currSeqNoGaps.c_str(), fastMatrix.matrix, evaluer, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid(), res);

                unsigned int len = Matcher::resultToBuffer(buffer, res, true, true);
                resultWriter.writeAdd(buffer, len, thread_idx);
            }
            resultWriter.writeEnd(queryKey, thread_idx);
        }

        kseq_destroy(seq);
        free(msaSequences);
        free(msaContent);

        delete[] maskedColumns;
        delete[] seqWeight;
    }
    resultWriter.close();
    headerWriter.close(true);
    sequenceWriter.close(true);
    msaReader.close();

    DBReader<unsigned int>::softlinkDb(par.db1, par.db2, (DBFiles::Files)(DBFiles::LOOKUP | DBFiles::SOURCE));

    if (sequenceReader != NULL) {
        sequenceReader->close();
        delete sequenceReader;
    }

    if (headerReader != NULL) {
        headerReader->close();
        delete headerReader;
    }
    delete[] fastMatrix.matrix;
    delete[] fastMatrix.matrixData;

    return EXIT_SUCCESS;
}
