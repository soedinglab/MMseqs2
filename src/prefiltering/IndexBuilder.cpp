#include "IndexBuilder.h"
#include "tantan.h"
#include "ExtendedSubstitutionMatrix.h"
#include "Masker.h"

#ifdef OPENMP
#include <omp.h>
#endif

char* getScoreLookup(BaseMatrix &matrix) {
    char *idScoreLookup = NULL;
    idScoreLookup = new char[matrix.alphabetSize];
    for (int aa = 0; aa < matrix.alphabetSize; aa++){
        short score = matrix.subMatrix[aa][aa];
        if (score > SCHAR_MAX || score < SCHAR_MIN) {
            Debug(Debug::WARNING) << "Truncating substitution matrix diagonal score!";
        }
        idScoreLookup[aa] = (char) score;
    }
    return idScoreLookup;
}


class DbInfo {
public:
    DbInfo(size_t dbFrom, size_t dbTo, unsigned int effectiveKmerSize, DBReader<unsigned int> & reader) {
        tableSize = 0;
        aaDbSize = 0;
        size_t dbSize = dbTo - dbFrom;
        sequenceOffsets = new size_t[dbSize];
        sequenceOffsets[0] = 0;
        for (size_t id = dbFrom; id < dbTo; id++) {
            const int seqLen = reader.getSeqLen(id);
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


void IndexBuilder::fillDatabase(IndexTable *indexTable, SequenceLookup ** externalLookup, BaseMatrix &subMat,
                                ScoreMatrix & three, ScoreMatrix & two, Sequence *seq,
                                DBReader<unsigned int> *dbr, size_t dbFrom, size_t dbTo, int kmerThr,
                                bool mask, bool maskLowerCaseMode, float maskProb, int maskNrepeats, int targetSearchMode) {
    Debug(Debug::INFO) << "Index table: counting k-mers\n";

    const bool isProfile = Parameters::isEqualDbtype(seq->getSeqType(), Parameters::DBTYPE_HMM_PROFILE);
    const bool isTargetSimiliarKmerSearch = isProfile || targetSearchMode;
    dbTo = std::min(dbTo, dbr->getSize());
    size_t dbSize = dbTo - dbFrom;
    DbInfo* info = new DbInfo(dbFrom, dbTo, seq->getEffectiveKmerSize(), *dbr);

    *externalLookup = new SequenceLookup(dbSize, info->aaDbSize);
    SequenceLookup *sequenceLookup = *externalLookup;


    // identical scores for memory reduction code
    char *idScoreLookup = getScoreLookup(subMat);
    Debug::Progress progress(dbTo-dbFrom);
    bool needMasking = (mask == 1 || maskNrepeats > 0  || maskLowerCaseMode == 1);
    size_t maskedResidues = 0;
    size_t totalKmerCount = 0;
    #pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        // need to prune low scoring k-mers through masking
        Masker *masker = NULL;
        if (needMasking) {
            masker = new Masker(subMat);
        }

        unsigned int alphabetSize = (indexTable != NULL) ? static_cast<unsigned int>(indexTable->getAlphabetSize())
                                                         : static_cast<unsigned int>(subMat.alphabetSize);
        Indexer idxer(alphabetSize, seq->getKmerSize());
        Sequence s(seq->getMaxLen(), seq->getSeqType(), &subMat, seq->getKmerSize(), seq->isSpaced(), false, true, seq->getUserSpacedKmerPattern());

        KmerGenerator *generator = NULL;
        if (isTargetSimiliarKmerSearch && indexTable != NULL) {
            generator = new KmerGenerator(seq->getKmerSize(), indexTable->getAlphabetSize(), kmerThr);
            if(isProfile){
                generator->setDivideStrategy(s.profile_matrix);
            }else{
                generator->setDivideStrategy(&three, &two);
            }
        }

        unsigned int *buffer = static_cast<unsigned int*>(malloc(seq->getMaxLen() * sizeof(unsigned int)));
        unsigned int bufferSize = seq->getMaxLen();
        #pragma omp for schedule(dynamic, 100) reduction(+:totalKmerCount, maskedResidues)
        for (size_t id = dbFrom; id < dbTo; id++) {
            progress.updateProgress();

            s.resetCurrPos();
            char *seqData = dbr->getData(id, thread_idx);
            unsigned int qKey = dbr->getDbKey(id);

            s.mapSequence(id - dbFrom, qKey, seqData, dbr->getSeqLen(id));
            if(s.getMaxLen() >= bufferSize ){
                buffer = static_cast<unsigned int*>(realloc(buffer, s.getMaxLen() * sizeof(unsigned int)));
                bufferSize = seq->getMaxLen();
            }
            // count similar or exact k-mers based on sequence type
            if (isTargetSimiliarKmerSearch) {
                // Find out if we should also mask profiles
                if(indexTable != NULL){
                    totalKmerCount += indexTable->addSimilarKmerCount(&s, generator);
                }
                unsigned char * seq = (isProfile) ? s.numConsensusSequence : s.numSequence;

                sequenceLookup->addSequence(seq, s.L, id - dbFrom, info->sequenceOffsets[id - dbFrom]);

            } else {
                // Do not mask if column state sequences are used
                maskedResidues += masker->maskSequence(s, mask, maskProb, maskLowerCaseMode, maskNrepeats);
                sequenceLookup->addSequence(s.numSequence, s.L, id - dbFrom, info->sequenceOffsets[id - dbFrom]);

                if(indexTable != NULL){
                    totalKmerCount += indexTable->addKmerCount(&s, &idxer, buffer, kmerThr, idScoreLookup);
                }
            }
        }

        free(buffer);

        if (generator != NULL) {
            delete generator;
        }
        if(masker != NULL) {
            delete masker;
        }
    }



    Debug(Debug::INFO) << "Index table: Masked residues: " << maskedResidues << "\n";
    if(indexTable != NULL && totalKmerCount == 0) {
        Debug(Debug::WARNING) << "No k-mer could be extracted for the database " << dbr->getDataFileName() << ".\n"
                            << "Maybe the sequences length is less than 14 residues.\n";
        if (maskedResidues == true){
            Debug(Debug::WARNING) << " or contains only low complexity regions.";
            Debug(Debug::WARNING) << "Use --mask 0 to deactivate the low complexity filter.\n";
        }
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
    if(indexTable != NULL){
        indexTable->initMemory(info->tableSize);
        indexTable->init();
    }

    delete info;
    if(indexTable != NULL) {
        Debug::Progress progress2(dbTo - dbFrom);
        Debug(Debug::INFO) << "Index table: fill\n";
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            Sequence s(seq->getMaxLen(), seq->getSeqType(), &subMat, seq->getKmerSize(), seq->isSpaced(), false, true,
                       seq->getUserSpacedKmerPattern());
            unsigned int alphabetSize = (indexTable != NULL) ? static_cast<unsigned int>(indexTable->getAlphabetSize())
                                                             : static_cast<unsigned int>(subMat.alphabetSize);
            Indexer idxer(alphabetSize, seq->getKmerSize());
            IndexEntryLocalTmp *buffer = static_cast<IndexEntryLocalTmp *>(malloc(
                    seq->getMaxLen() * sizeof(IndexEntryLocalTmp)));
            size_t bufferSize = seq->getMaxLen();
            KmerGenerator *generator = NULL;
            if (isTargetSimiliarKmerSearch) {
                generator = new KmerGenerator(seq->getKmerSize(), indexTable->getAlphabetSize(), kmerThr);
                if (isProfile) {
                    generator->setDivideStrategy(s.profile_matrix);
                } else {
                    generator->setDivideStrategy(&three, &two);
                }
            }

#pragma omp for schedule(dynamic, 100)
            for (size_t id = dbFrom; id < dbTo; id++) {
                s.resetCurrPos();
                progress2.updateProgress();

                unsigned int qKey = dbr->getDbKey(id);
                if (isTargetSimiliarKmerSearch) {
                    s.mapSequence(id - dbFrom, qKey, dbr->getData(id, thread_idx), dbr->getSeqLen(id));
                    indexTable->addSimilarSequence(&s, generator, &buffer, bufferSize, &idxer);
                } else {
                    s.mapSequence(id - dbFrom, qKey, sequenceLookup->getSequence(id - dbFrom));
                    indexTable->addSequence(&s, &idxer, &buffer, bufferSize, kmerThr, idScoreLookup);
                }
            }

            if (generator != NULL) {
                delete generator;
            }

            free(buffer);
        }
    }
    if(idScoreLookup!=NULL){
        delete[] idScoreLookup;
    }
    if(indexTable != NULL){
        indexTable->revertPointer();
        indexTable->sortDBSeqLists();
    }
}
