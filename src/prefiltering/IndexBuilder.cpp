#include "IndexBuilder.h"
#include "tantan.h"

char* getScoreLookup(BaseMatrix &matrix) {
    char *idScoreLookup = NULL;
    idScoreLookup = new char[matrix.alphabetSize];
    for (int aa = 0; aa < matrix.alphabetSize; aa++){
        short score = matrix.subMatrix[aa][aa];
        if (score > CHAR_MAX || score < CHAR_MIN) {
            Debug(Debug::WARNING) << "Truncating substitution matrix diagonal score!";
        }
        idScoreLookup[aa] = (char) score;
    }
    return idScoreLookup;
}

class ProbabilityMatrix {
public:
    ProbabilityMatrix(BaseMatrix &matrix) : alphabetSize(matrix.alphabetSize) {
        probMatrix = new double*[matrix.alphabetSize];
        probMatrixPointers = new const double*[matrix.alphabetSize];
        std::fill_n(hardMaskTable, 256, matrix.aa2int[(int) 'X']);
        for (int i = 0; i < matrix.alphabetSize; ++i) {
            probMatrix[i] = new double[matrix.alphabetSize];
            probMatrixPointers[i] = probMatrix[i];
            for (int j = 0; j < matrix.alphabetSize; ++j) {
                probMatrix[i][j] = matrix.probMatrix[i][j] / (matrix.pBack[i] * matrix.pBack[j]);
            }
        }
    }
    ~ProbabilityMatrix() {
        for (int i = 0; i < alphabetSize; ++i) {
            delete[] probMatrix[i];
        }
        delete[] probMatrix;
        delete[] probMatrixPointers;
    }

    char hardMaskTable[256];
    const double **probMatrixPointers;

private:
    const int alphabetSize;
    double **probMatrix;

};

class DbInfo {
public:
    DbInfo(size_t dbFrom, size_t dbTo, unsigned int effectiveKmerSize, bool profile, unsigned int* seqLengths) {
        tableSize = 0;
        aaDbSize = 0;
        size_t dbSize = dbTo - dbFrom;
        sequenceOffsets = new size_t[dbSize];
        sequenceOffsets[0] = 0;
        for (size_t id = dbFrom; id < dbTo; id++) {
            int seqLen;
            if (profile) {
                // remove /0 and convert to profile length
                seqLen = std::max(static_cast<int>(seqLengths[id]) - 1, 0) / static_cast<int>(Sequence::PROFILE_AA_SIZE);
            } else {
                // remove /n and /0
                seqLen = std::max(static_cast<int>(seqLengths[id]) - 2, 0);
            }
            aaDbSize += seqLen;
            size_t idFromNull = (id - dbFrom);
            if (id < dbTo - 1) {
                sequenceOffsets[idFromNull + 1] = sequenceOffsets[idFromNull] + seqLen;
            }
            if (Util::overlappingKmers(seqLen, effectiveKmerSize > 0)) {
                tableSize += 1;
            }
        }
    }

    ~DbInfo() {
        delete[] sequenceOffsets;
    }

    size_t tableSize;
    size_t aaDbSize;
    size_t *sequenceOffsets;
};


void IndexBuilder::fillDatabase(IndexTable *indexTable, SequenceLookup **maskedLookup, SequenceLookup **unmaskedLookup,
                                BaseMatrix &subMat, Sequence *seq,
                                DBReader<unsigned int> *dbr, size_t dbFrom, size_t dbTo, int kmerThr) {
    Debug(Debug::INFO) << "Index table: counting k-mers...\n";

    const bool isProfile = seq->getSeqType() == Sequence::HMM_PROFILE;
    const bool isProfileStateSeq = seq->getSeqType() == Sequence::PROFILE_STATE_SEQ;

    dbTo = std::min(dbTo, dbr->getSize());
    size_t dbSize = dbTo - dbFrom;
    DbInfo info(dbFrom, dbTo, seq->getEffectiveKmerSize(), isProfile, dbr->getSeqLens());

    SequenceLookup *sequenceLookup;
    if (unmaskedLookup != NULL && maskedLookup == NULL) {
        *unmaskedLookup = new SequenceLookup(dbSize, info.aaDbSize);
        sequenceLookup = *unmaskedLookup;
    } else if (unmaskedLookup == NULL && maskedLookup != NULL) {
        *maskedLookup = new SequenceLookup(dbSize, info.aaDbSize);
        sequenceLookup = *maskedLookup;
    } else if (unmaskedLookup != NULL && maskedLookup != NULL) {
        *unmaskedLookup = new SequenceLookup(dbSize, info.aaDbSize);
        *maskedLookup = new SequenceLookup(dbSize, info.aaDbSize);
        sequenceLookup = *maskedLookup;
    }

    // need to prune low scoring k-mers through masking
    ProbabilityMatrix *probMatrix = NULL;
    if (maskedLookup != NULL) {
        probMatrix = new ProbabilityMatrix(subMat);
    }

    // identical scores for memory reduction code
    char *idScoreLookup = getScoreLookup(subMat);

    size_t maskedResidues = 0;
    size_t totalKmerCount = 0;
    #pragma omp parallel
    {
        Indexer idxer(static_cast<unsigned int>(indexTable->getAlphabetSize()), seq->getKmerSize());
        Sequence s(seq->getMaxLen(), seq->getSeqType(), &subMat, seq->getKmerSize(), seq->isSpaced(), false);

        KmerGenerator *generator = NULL;
        if (isProfile) {
            generator = new KmerGenerator(seq->getKmerSize(), indexTable->getAlphabetSize(), kmerThr);
            generator->setDivideStrategy(s.profile_matrix);
        }

        unsigned int *buffer = new unsigned int[seq->getMaxLen()];
        char *charSequence = new char[seq->getMaxLen()];

        #pragma omp for schedule(dynamic, 100) reduction(+:totalKmerCount, maskedResidues)
        for (size_t id = dbFrom; id < dbTo; id++) {
            Debug::printProgress(id - dbFrom);

            s.resetCurrPos();
            char *seqData = dbr->getData(id);
            unsigned int qKey = dbr->getDbKey(id);
            s.mapSequence(id - dbFrom, qKey, seqData);

            // count similar or exact k-mers based on sequence type
            if (isProfile) {
                // Find out if we should also mask profiles
                totalKmerCount += indexTable->addSimilarKmerCount(&s, generator, &idxer, kmerThr, idScoreLookup);
                (*unmaskedLookup)->addSequence(s.int_consensus_sequence, s.L, id - dbFrom, info.sequenceOffsets[id - dbFrom]);
            } else {
                // Do not mask if column state sequences are used
                if (unmaskedLookup != NULL) {
                    (*unmaskedLookup)->addSequence(s.int_sequence, s.L, id - dbFrom, info.sequenceOffsets[id - dbFrom]);
                }
                if (maskedLookup != NULL) {
                    if(isProfileStateSeq == false) {
                        for (int i = 0; i < s.L; ++i) {
                            charSequence[i] = (char) s.int_sequence[i];
                        }
                        // s.print();
                        maskedResidues += tantan::maskSequences(charSequence,
                                                                charSequence + s.L,
                                                                50 /*options.maxCycleLength*/,
                                                                probMatrix->probMatrixPointers,
                                                                0.005 /*options.repeatProb*/,
                                                                0.05 /*options.repeatEndProb*/,
                                                                0.9 /*options.repeatOffsetProbDecay*/,
                                                                0, 0,
                                                                0.9 /*options.minMaskProb*/,
                                                                probMatrix->hardMaskTable);

                        for (int i = 0; i < s.L; i++) {
                            s.int_sequence[i] = charSequence[i];
                        }
                    }
                    (*maskedLookup)->addSequence(s.int_sequence, s.L, id - dbFrom, info.sequenceOffsets[id - dbFrom]);
                }

                totalKmerCount += indexTable->addKmerCount(&s, &idxer, buffer, kmerThr, idScoreLookup);
            }
        }

        delete[] charSequence;
        delete[] buffer;

        if (generator != NULL) {
            delete generator;
        }
    }

    if(probMatrix != NULL) {
        delete probMatrix;
    }

    Debug(Debug::INFO) << "\nIndex table: Masked residues: " << maskedResidues << "\n";
    if(totalKmerCount == 0) {
        Debug(Debug::ERROR) << "No k-mer could be extracted for the database " << dbr->getDataFileName() << ".\n"
                            << "Maybe the sequences length is less than 14 residues.\n";
        if (maskedResidues == true){
            Debug(Debug::ERROR) << " or contains only low complexity regions.";
            Debug(Debug::ERROR) << "Use --mask-mode 0 to deactivate the low complexity filter.\n";
        }
        EXIT(EXIT_FAILURE);
    }

    dbr->remapData();

    //TODO find smart way to remove extrem k-mers without harming huge protein families
//    size_t lowSelectiveResidues = 0;
//    const float dbSize = static_cast<float>(dbTo - dbFrom);
//    for(size_t kmerIdx = 0; kmerIdx < indexTable->getTableSize(); kmerIdx++){
//        size_t res = (size_t) indexTable->getOffset(kmerIdx);
//        float selectivityOfKmer = (static_cast<float>(res)/dbSize);
//        if(selectivityOfKmer > 0.005){
//            indexTable->getOffset()[kmerIdx] = 0;
//            lowSelectiveResidues += res;
//        }
//    }
//    Debug(Debug::INFO) << "Index table: Remove "<< lowSelectiveResidues <<" none selective residues\n";
//    Debug(Debug::INFO) << "Index table: init... from "<< dbFrom << " to "<< dbTo << "\n";

    indexTable->initMemory(info.tableSize);
    indexTable->init();

    Debug(Debug::INFO) << "Index table: fill...\n";
    #pragma omp parallel
    {
        Sequence s(seq->getMaxLen(), seq->getSeqType(), &subMat, seq->getKmerSize(), seq->isSpaced(), false);
        Indexer idxer(static_cast<unsigned int>(indexTable->getAlphabetSize()), seq->getKmerSize());
        IndexEntryLocalTmp *buffer = new IndexEntryLocalTmp[seq->getMaxLen()];

        KmerGenerator *generator = NULL;
        if (isProfile) {
            generator = new KmerGenerator(seq->getKmerSize(), indexTable->getAlphabetSize(), kmerThr);
            generator->setDivideStrategy(s.profile_matrix);
        }

        #pragma omp for schedule(dynamic, 100)
        for (size_t id = dbFrom; id < dbTo; id++) {
            s.resetCurrPos();
            Debug::printProgress(id - dbFrom);

            unsigned int qKey = dbr->getDbKey(id);
            if (isProfile) {
                s.mapSequence(id - dbFrom, qKey, dbr->getData(id));
                indexTable->addSimilarSequence(&s, generator, &idxer, kmerThr, idScoreLookup);
            } else {
                s.mapSequence(id - dbFrom, qKey, sequenceLookup->getSequence(id - dbFrom));
                indexTable->addSequence(&s, &idxer, buffer, kmerThr, idScoreLookup);
            }
        }

        if (generator != NULL) {
            delete generator;
        }

        delete [] buffer;
    }
    delete[] idScoreLookup;

    indexTable->sortDBSeqLists();
    Debug(Debug::INFO) << "\nIndex table: removing duplicate entries...\n";
    indexTable->revertPointer();
    Debug(Debug::INFO) << "Index table init done.\n\n";
}
