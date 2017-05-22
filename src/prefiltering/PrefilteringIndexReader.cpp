#include "PrefilteringIndexReader.h"
#include "DBWriter.h"
#include "Prefiltering.h"
#include "ExtendedSubstitutionMatrix.h"
#include "FileUtil.h"
#include "tantan.h"

const char*  PrefilteringIndexReader::CURRENT_VERSION="3.3.2";
unsigned int PrefilteringIndexReader::VERSION = 0;
unsigned int PrefilteringIndexReader::META = 1;
unsigned int PrefilteringIndexReader::SCOREMATRIXNAME = 2;
unsigned int PrefilteringIndexReader::SCOREMATRIX2MER = 3;
unsigned int PrefilteringIndexReader::SCOREMATRIX3MER = 4;
unsigned int PrefilteringIndexReader::DBRINDEX = 5;
unsigned int PrefilteringIndexReader::HDRINDEX = 6;

unsigned int PrefilteringIndexReader::ENTRIES = 91;
unsigned int PrefilteringIndexReader::ENTRIESOFFSETS = 92;
unsigned int PrefilteringIndexReader::ENTRIESNUM = 93;
unsigned int PrefilteringIndexReader::SEQCOUNT = 94;
unsigned int PrefilteringIndexReader::SEQINDEXDATA = 95;
unsigned int PrefilteringIndexReader::SEQINDEXDATASIZE = 96;
unsigned int PrefilteringIndexReader::SEQINDEXSEQOFFSET = 97;
unsigned int PrefilteringIndexReader::UNMASKEDSEQINDEXDATA = 98;

bool PrefilteringIndexReader::checkIfIndexFile(DBReader<unsigned int>* reader) {
    char * version = reader->getDataByDBKey(VERSION);
    if(version == NULL){
        return false;
    }
    return (strncmp(version, CURRENT_VERSION, strlen(CURRENT_VERSION)) == 0 ) ? true : false;
}

void PrefilteringIndexReader::createIndexFile(std::string outDB, DBReader<unsigned int> *dbr, DBReader<unsigned int> *hdbr,
                                              BaseMatrix * subMat, int maxSeqLen, bool hasSpacedKmer,
                                              bool compBiasCorrection, const int split, int alphabetSize, int kmerSize,
                                              bool diagonalScoring, int maskMode, int seqType, int kmerThr, int threads) {
    std::string outIndexName(outDB); // db.sk6
    std::string spaced = (hasSpacedKmer == true) ? "s" : "";
    outIndexName.append(".").append(spaced).append("k").append(SSTR(kmerSize));

    DBWriter writer(outIndexName.c_str(), std::string(outIndexName).append(".index").c_str(), 1, DBWriter::BINARY_MODE);
    writer.open();

    if (seqType != Sequence::HMM_PROFILE) {
        ScoreMatrix *s3 = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 3);
        char* serialized3mer = ScoreMatrix::serialize(*s3);
        Debug(Debug::INFO) << "Write SCOREMATRIX3MER (" << SCOREMATRIX3MER << ")\n";
        writer.writeData(serialized3mer, ScoreMatrix::size(*s3), SCOREMATRIX3MER, 0);
        free(serialized3mer);
        ScoreMatrix::cleanup(s3);

        writer.alignToPageSize();

        ScoreMatrix *s2 = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 2);
        char* serialized2mer = ScoreMatrix::serialize(*s2);
        Debug(Debug::INFO) << "Write SCOREMATRIX2MER (" << SCOREMATRIX2MER << ")\n";
        writer.writeData(serialized2mer, ScoreMatrix::size(*s2), SCOREMATRIX2MER, 0);
        free(serialized2mer);
        ScoreMatrix::cleanup(s2);
    }

    Sequence seq(maxSeqLen, subMat->aa2int, subMat->int2aa, seqType, kmerSize, hasSpacedKmer, compBiasCorrection);

    for (int step = 0; step < split; step++) {
        size_t splitStart = 0;
        size_t splitSize  = 0;
        Util::decomposeDomainByAminoAcid(dbr->getAminoAcidDBSize(), dbr->getSeqLens(), dbr->getSize(),
                                         step, split, &splitStart, &splitSize);
        IndexTable *indexTable = new IndexTable(alphabetSize, kmerSize, false);

        SequenceLookup *unmaskedLookup = NULL;
        PrefilteringIndexReader::fillDatabase(dbr, &seq, indexTable, subMat,
                                              splitStart, splitStart + splitSize, diagonalScoring,
                                              maskMode, &unmaskedLookup, kmerThr, threads);

        indexTable->printStatistics(subMat->int2aa);


        // save the entries
        unsigned int entries_key = MathUtil::concatenate(ENTRIES, step);
        Debug(Debug::INFO) << "Write ENTRIES (" << entries_key << ")\n";
        char *entries = (char *) indexTable->getEntries();
        size_t entriesSize = indexTable->getTableEntriesNum() * indexTable->getSizeOfEntry();
        writer.writeData(entries, entriesSize, entries_key, 0);

        // save the size
        unsigned int entriesoffsets_key = MathUtil::concatenate(ENTRIESOFFSETS, step);
        Debug(Debug::INFO) << "Write ENTRIESOFFSETS (" << entriesoffsets_key << ")\n";

        char *offsets = (char*)indexTable->getOffsets();
        size_t offsetsSize = (indexTable->getTableSize() + 1) * sizeof(size_t);
        writer.writeData(offsets, offsetsSize, entriesoffsets_key, 0);
        indexTable->deleteEntries();

        SequenceLookup *lookup = indexTable->getSequenceLookup();
        unsigned int seqindexdata_key = MathUtil::concatenate(SEQINDEXDATA, step);
        Debug(Debug::INFO) << "Write SEQINDEXDATA (" << seqindexdata_key << ")\n";
        writer.writeData(lookup->getData(), lookup->getDataSize(), seqindexdata_key, 0);

        if (unmaskedLookup != NULL) {
            unsigned int unmaskedSeqindexdata_key = MathUtil::concatenate(UNMASKEDSEQINDEXDATA, step);
            Debug(Debug::INFO) << "Write UNMASKEDSEQINDEXDATA (" << unmaskedSeqindexdata_key << ")\n";
            writer.writeData(unmaskedLookup->getData(), unmaskedLookup->getDataSize(), unmaskedSeqindexdata_key, 0);
            delete unmaskedLookup;
        }

        unsigned int seqindex_datasize_key = MathUtil::concatenate(SEQINDEXDATASIZE, step);
        Debug(Debug::INFO) << "Write SEQINDEXDATASIZE (" << seqindex_datasize_key << ")\n";
        int64_t seqindexDataSize = lookup->getDataSize();
        char *seqindexDataSizePtr = (char *) &seqindexDataSize;
        writer.writeData(seqindexDataSizePtr, 1 * sizeof(int64_t), seqindex_datasize_key, 0);

        size_t *sequenceOffsets = lookup->getOffsets();
        size_t sequenceCount = lookup->getSequenceCount();
        unsigned int seqindex_seqoffset = MathUtil::concatenate(SEQINDEXSEQOFFSET, step);
        Debug(Debug::INFO) << "Write SEQINDEXSEQOFFSET (" << seqindex_seqoffset << ")\n";
        writer.writeData((char *) sequenceOffsets, (sequenceCount + 1) * sizeof(size_t), seqindex_seqoffset, 0);

        // meta data
        // ENTRIESNUM
        unsigned int entriesnum_key = MathUtil::concatenate(ENTRIESNUM, step);
        Debug(Debug::INFO) << "Write ENTRIESNUM (" << entriesnum_key << ")\n";
        uint64_t entriesNum = indexTable->getTableEntriesNum();
        char *entriesNumPtr = (char *) &entriesNum;
        writer.writeData(entriesNumPtr, 1 * sizeof(uint64_t), entriesnum_key, 0);
        // SEQCOUNT
        unsigned int tablesize_key = MathUtil::concatenate(SEQCOUNT, step);
        Debug(Debug::INFO) << "Write SEQCOUNT (" << tablesize_key << ")\n";
        size_t tablesize = {indexTable->getSize()};
        char *tablesizePtr = (char *) &tablesize;
        writer.writeData(tablesizePtr, 1 * sizeof(size_t), tablesize_key, 0);

        delete indexTable;
    }
    Debug(Debug::INFO) << "Write META (" << META << ")\n";
    int local = 1;
    int spacedKmer = (hasSpacedKmer) ? 1 : 0;
    int headers = (hdbr != NULL) ? 1 : 0;
    int metadata[] = {kmerSize, alphabetSize, maskMode, split, local, spacedKmer, kmerThr, seqType, headers};
    char *metadataptr = (char *) &metadata;
    writer.writeData(metadataptr, sizeof(metadata), META, 0);

    Debug(Debug::INFO) << "Write SCOREMATRIXNAME (" << SCOREMATRIXNAME << ")\n";
    writer.writeData(subMat->getMatrixName().c_str(), subMat->getMatrixName().length(), SCOREMATRIXNAME, 0);

    Debug(Debug::INFO) << "Write VERSION (" << VERSION << ")\n";
    writer.writeData((char *) CURRENT_VERSION, strlen(CURRENT_VERSION) * sizeof(char), VERSION, 0);

    Debug(Debug::INFO) << "Write DBRINDEX (" << DBRINDEX << ")\n";
    char* data = DBReader<unsigned int>::serialize(*dbr);
    writer.writeData(data, DBReader<unsigned int>::indexMemorySize(*dbr), DBRINDEX, 0);
    free(data);

    if (hdbr != NULL) {
        Debug(Debug::INFO) << "Write HDRINDEX (" << HDRINDEX << ")\n";
        data = DBReader<unsigned int>::serialize(*hdbr);
        writer.writeData(data, DBReader<unsigned int>::indexMemorySize(*hdbr), HDRINDEX, 0);
        free(data);
    }

    writer.close();
    Debug(Debug::INFO) << "Done. \n";
}

DBReader<unsigned int> *PrefilteringIndexReader::openNewReader(DBReader<unsigned int>*dbr) {
    char * data = dbr->getDataByDBKey(DBRINDEX);
    return DBReader<unsigned int>::unserialize(data);
}

SequenceLookup *PrefilteringIndexReader::getSequenceLookup(DBReader<unsigned int>*dbr, int split) {
    size_t sequenceCount    = *((size_t *)dbr->getDataByDBKey(MathUtil::concatenate(SEQCOUNT, split)));

    char * seqData = dbr->getDataByDBKey(MathUtil::concatenate(SEQINDEXDATA, split));
    size_t seqOffsetsId = dbr->getId(MathUtil::concatenate(SEQINDEXSEQOFFSET, split));
    char * seqOffsetsData = dbr->getData(seqOffsetsId);
    size_t seqOffsetLength = dbr->getSeqLens(seqOffsetsId);

    SequenceLookup *sequenceLookup = new SequenceLookup(sequenceCount);
    sequenceLookup->initLookupByExternalData(seqData, seqOffsetLength, (size_t *) seqOffsetsData);

    return sequenceLookup;
}

SequenceLookup *PrefilteringIndexReader::getUnmaskedSequenceLookup(DBReader<unsigned int>*dbr, int split) {
    unsigned int key = MathUtil::concatenate(UNMASKEDSEQINDEXDATA, split);
    size_t id;
    if ((id = dbr->getId(key)) == UINT_MAX) {
        return NULL;
    }

    char * seqData = dbr->getData(id);
    size_t seqOffsetsId = dbr->getId(MathUtil::concatenate(SEQINDEXSEQOFFSET, split));
    char * seqOffsetsData = dbr->getData(seqOffsetsId);
    size_t seqOffsetLength = dbr->getSeqLens(seqOffsetsId);

    size_t sequenceCount = *((size_t *)dbr->getDataByDBKey(MathUtil::concatenate(SEQCOUNT, split)));
    SequenceLookup *sequenceLookup = new SequenceLookup(sequenceCount);
    sequenceLookup->initLookupByExternalData(seqData, seqOffsetLength, (size_t *) seqOffsetsData);

    return sequenceLookup;
}

IndexTable *PrefilteringIndexReader::generateIndexTable(DBReader<unsigned int> *dbr, int split, bool diagonalScoring) {
    PrefilteringIndexData data = getMetadata(dbr);
    IndexTable *retTable;
    if (data.local) {
        retTable = new IndexTable(data.alphabetSize, data.kmerSize, true);
    }else {
        Debug(Debug::ERROR) << "Search mode is not valid.\n";
        EXIT(EXIT_FAILURE);
    }

    SequenceLookup * sequenceLookup = NULL;
    if(diagonalScoring == true) {
        sequenceLookup = getSequenceLookup(dbr, split);
    }

    int64_t entriesNum = *((int64_t *)dbr->getDataByDBKey(MathUtil::concatenate(ENTRIESNUM, split)));
    size_t sequenceCount = *((size_t *)dbr->getDataByDBKey(MathUtil::concatenate(SEQCOUNT, split)));

    char *entriesData = dbr->getDataByDBKey(MathUtil::concatenate(ENTRIES, split));
    char *entriesOffsetsData = dbr->getDataByDBKey(MathUtil::concatenate(ENTRIESOFFSETS, split));
    retTable->initTableByExternalData(sequenceCount, entriesNum,
                                      (IndexEntryLocal*) entriesData, (size_t *)entriesOffsetsData, sequenceLookup);
    return retTable;
}

void PrefilteringIndexReader::printSummary(DBReader<unsigned int> *dbr) {
    Debug(Debug::INFO) << "Index version: " << dbr->getDataByDBKey(VERSION) << "\n";
    int *metadata_tmp = (int *) dbr->getDataByDBKey(META);

    Debug(Debug::INFO) << "KmerSize:     " << metadata_tmp[0] << "\n";
    Debug(Debug::INFO) << "AlphabetSize: " << metadata_tmp[1] << "\n";
    Debug(Debug::INFO) << "Mask:         " << metadata_tmp[2] << "\n";
    Debug(Debug::INFO) << "Split:        " << metadata_tmp[3] << "\n";
    Debug(Debug::INFO) << "Type:         " << metadata_tmp[4] << "\n";
    Debug(Debug::INFO) << "Spaced:       " << metadata_tmp[5] << "\n";
    Debug(Debug::INFO) << "KmerScore:    " << metadata_tmp[6] << "\n";
    Debug(Debug::INFO) << "SequenceType: " << metadata_tmp[7] << "\n";
    Debug(Debug::INFO) << "Headers:      " << metadata_tmp[8] << "\n";

    Debug(Debug::INFO) << "ScoreMatrix:  " << dbr->getDataByDBKey(SCOREMATRIXNAME) << "\n";
}

PrefilteringIndexData PrefilteringIndexReader::getMetadata(DBReader<unsigned int> *dbr) {
    PrefilteringIndexData prefData;

    int *metadata_tmp = (int *) dbr->getDataByDBKey(META);

    prefData.kmerSize = metadata_tmp[0];
    prefData.alphabetSize = metadata_tmp[1];
    prefData.maskMode = metadata_tmp[2];
    prefData.split = metadata_tmp[3];
    prefData.local = metadata_tmp[4];
    prefData.spacedKmer = metadata_tmp[5];
    prefData.kmerThr = metadata_tmp[6];
    prefData.seqType = metadata_tmp[7];
    prefData.headers = metadata_tmp[8];

    return prefData;
}

std::string PrefilteringIndexReader::getSubstitutionMatrixName(DBReader<unsigned int> *dbr) {
    return std::string(dbr->getDataByDBKey(SCOREMATRIXNAME));
}
//
ScoreMatrix *PrefilteringIndexReader::get2MerScoreMatrix(DBReader<unsigned int> *dbr) {
    PrefilteringIndexData meta = getMetadata(dbr);
    return ScoreMatrix::unserialize(dbr->getDataByDBKey(SCOREMATRIX2MER), meta.alphabetSize, 2);
}

ScoreMatrix *PrefilteringIndexReader::get3MerScoreMatrix(DBReader<unsigned int> *dbr) {
    PrefilteringIndexData meta = getMetadata(dbr);
    return ScoreMatrix::unserialize(dbr->getDataByDBKey(SCOREMATRIX3MER), meta.alphabetSize, 3);
}

std::string PrefilteringIndexReader::searchForIndex(const std::string &pathToDB) {
    for (size_t spaced = 0; spaced < 2; spaced++) {
        for (size_t k = 5; k <= 7; k++) {
            std::string outIndexName(pathToDB); // db.sk6
            std::string s = (spaced == true) ? "s" : "";
            outIndexName.append(".").append(s).append("k").append(SSTR(k));
            if (FileUtil::fileExists(outIndexName.c_str()) == true) {
                return outIndexName;
            }
        }
    }

    return "";
}


void PrefilteringIndexReader::fillDatabase(DBReader<unsigned int> *dbr, Sequence *seq, IndexTable *indexTable,
                                           BaseMatrix *subMat, size_t dbFrom, size_t dbTo,
                                           bool diagonalScoring, int maskMode,
                                           SequenceLookup **unmaskedLookup /* for maskMode 2 */,
                                           int kmerThr, unsigned int threads) {
    Debug(Debug::INFO) << "Index table: counting k-mers...\n";
    // fill and init the index table
    size_t aaCount = 0;
    dbTo = std::min(dbTo, dbr->getSize());
    size_t maskedResidues = 0;
    size_t totalKmerCount = 0;
    size_t tableSize = 0;

    size_t dbSize = dbTo - dbFrom;
    size_t *sequenceOffSet = new size_t[dbSize];
    // identical scores for memory reduction code
    char *idScoreLookup = new char[subMat->alphabetSize];
    for (int aa = 0; aa < subMat->alphabetSize; aa++){
        short score = subMat->subMatrix[aa][aa];
        if (score > CHAR_MAX || score < CHAR_MIN) {
            Debug(Debug::WARNING) << "Truncating substitution matrix diagonal score!";
        }
        idScoreLookup[aa] = (char) score;
    }

    // need to prune low scoring k-mers
    // masking code
    SubstitutionMatrix mat("blosum62.out", 2.0, 0.0);
    double probMatrix[subMat->alphabetSize][subMat->alphabetSize];
    const double **probMatrixPointers = new const double*[subMat->alphabetSize];
    char hardMaskTable[256];
    std::fill_n(hardMaskTable, 256, subMat->aa2int[(int) 'X']);
    for (int i = 0; i < subMat->alphabetSize; ++i){
        probMatrixPointers[i] = probMatrix[i];
        for (int j = 0; j < subMat->alphabetSize; ++j){
            probMatrix[i][j]  = std::exp(0.324032 * mat.subMatrix[i][j]);
        }
    }

    size_t aaDbSize = 0;
    sequenceOffSet[0] = 0;
    for (size_t id = dbFrom; id < dbTo; id++) {
        int seqLen = std::max(static_cast<int>(dbr->getSeqLens(id)) - 2, 0);
        aaDbSize += seqLen; // remove /n and /0
        size_t idFromNull = (id - dbFrom);
        if (id < dbTo - 1) {
            sequenceOffSet[idFromNull + 1] = sequenceOffSet[idFromNull] + seqLen;
        }
        if (Util::overlappingKmers(seqLen, seq->getEffectiveKmerSize() > 0)) {
            tableSize += 1;
        }
    }

    SequenceLookup *sequenceLookup = new SequenceLookup(dbSize, aaDbSize);

    if (maskMode == 2) {
        *unmaskedLookup = new SequenceLookup(dbSize, aaDbSize);
    }

#pragma omp parallel
    {
        Indexer idxer(static_cast<unsigned int>(subMat->alphabetSize), seq->getKmerSize());
        Sequence s(seq->getMaxLen(), seq->aa2int, seq->int2aa,
                   seq->getSeqType(), seq->getKmerSize(), seq->isSpaced(), false);

        KmerGenerator *generator = NULL;
        if (seq->getSeqType() == Sequence::HMM_PROFILE) {
            generator = new KmerGenerator(seq->getKmerSize(), subMat->alphabetSize, kmerThr);
            generator->setDivideStrategy(s.profile_matrix);
        }

        unsigned int * buffer = new unsigned int[seq->getMaxLen()];
        char * charSequence = new char[seq->getMaxLen()];
#pragma omp for schedule(dynamic, 100) reduction(+:aaCount, totalKmerCount, maskedResidues)
        for (size_t id = dbFrom; id < dbTo; id++) {
            s.resetCurrPos();
            Debug::printProgress(id - dbFrom);
            char *seqData = dbr->getData(id);
            unsigned int qKey = dbr->getDbKey(id);
            s.mapSequence(id - dbFrom, qKey, seqData);
            // count similar or exact k-mers based on sequence type
            if (seq->getSeqType() == Sequence::HMM_PROFILE) {
                totalKmerCount += indexTable->addSimilarKmerCount(&s, generator, &idxer, kmerThr, idScoreLookup);
            } else {
                if (maskMode == 2) {
                    (*unmaskedLookup)->addSequence(&s, id - dbFrom, sequenceOffSet[id - dbFrom]);
                }

                if (maskMode >= 1) {
                    for (int i = 0; i < s.L; i++) {
                        charSequence[i] = (char) s.int_sequence[i];
                    }
                    maskedResidues += tantan::maskSequences(charSequence,
                                                            charSequence + s.L,
                                                            50 /*options.maxCycleLength*/,
                                                            probMatrixPointers,
                                                            0.005 /*options.repeatProb*/,
                                                            0.05 /*options.repeatEndProb*/,
                                                            0.9 /*options.repeatOffsetProbDecay*/,
                                                            0, 0,
                                                            0.9 /*options.minMaskProb*/,
                                                            hardMaskTable);
                    for (int i = 0; i < s.L; i++) {
                        s.int_sequence[i] = charSequence[i];
                    }
//                maskedResidues += Util::maskLowComplexity(subMat, &s, s.L, 12, 3,
//                                                          indexTable->getAlphabetSize(), seq->aa2int[(unsigned char) 'X'], true, true, true, true);
                }
                aaCount += s.L;
                totalKmerCount += indexTable->addKmerCount(&s, &idxer, buffer, kmerThr, idScoreLookup);
            }
            sequenceLookup->addSequence(&s, id - dbFrom, sequenceOffSet[id - dbFrom]);
        }

        delete [] charSequence;
        delete [] buffer;

        if (generator) {
            delete generator;
        }
    }
    delete[] sequenceOffSet;
    dbr->remapData();
    Debug(Debug::INFO) << "\n";
    Debug(Debug::INFO) << "Index table: Masked residues: " << maskedResidues << "\n";
    if(totalKmerCount == 0) {
        Debug(Debug::ERROR) << "No k-mer could be extracted for the database " << dbr->getDataFileName() << ".\n"
                            << "Maybe the sequences length is less than 14 residues";
        Debug(Debug::ERROR) << "\n";
        if (maskedResidues == true){
            Debug(Debug::ERROR) <<  " or contains only low complexity regions.";
            Debug(Debug::ERROR) << "Use --do-not-mask to deactivate low complexity filter.\n";
        }
        EXIT(EXIT_FAILURE);
    }
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
    size_t tableEntriesNum = 0;
    for (size_t i = 0; i < indexTable->getTableSize(); i++) {
        tableEntriesNum += (size_t) indexTable->getOffset(i);
    }

    if (diagonalScoring == true) {
        indexTable->initMemory(tableEntriesNum, sequenceLookup, tableSize);
    } else {
        indexTable->initMemory(tableEntriesNum, NULL, tableSize);
    }
    indexTable->init();

    Debug(Debug::INFO) << "Index table: fill...\n";

#pragma omp parallel
    {
        Sequence s(seq->getMaxLen(), seq->aa2int, seq->int2aa,
                   seq->getSeqType(), seq->getKmerSize(), seq->isSpaced(), false);
        Indexer idxer(static_cast<unsigned int>(subMat->alphabetSize), seq->getKmerSize());
        IndexEntryLocalTmp * buffer = new IndexEntryLocalTmp[seq->getMaxLen()];

        KmerGenerator *generator = NULL;
        if (seq->getSeqType() == Sequence::HMM_PROFILE) {
            generator = new KmerGenerator(seq->getKmerSize(), subMat->alphabetSize, kmerThr);
            generator->setDivideStrategy(s.profile_matrix);
        }

        #pragma omp for schedule(dynamic, 100)
        for (size_t id = dbFrom; id < dbTo; id++) {
            s.resetCurrPos();
            Debug::printProgress(id - dbFrom);

            unsigned int qKey = dbr->getDbKey(id);
            if (seq->getSeqType() == Sequence::HMM_PROFILE) {
                char *seqData = dbr->getData(id);
                s.mapSequence(id - dbFrom, qKey, seqData);
                indexTable->addSimilarSequence(&s, generator, &idxer, kmerThr, idScoreLookup);
            } else {
                s.mapSequence(id - dbFrom, qKey, sequenceLookup->getSequence(id - dbFrom));
                indexTable->addSequence(&s, &idxer, buffer, kmerThr, idScoreLookup);
            }
        }
        delete [] buffer;

        if (generator) {
            delete generator;
        }
    }
    if (diagonalScoring == false) {
        delete sequenceLookup;
    }
    delete [] idScoreLookup;
    if ((dbTo-dbFrom) > 10000) {
        Debug(Debug::INFO) << "\n";
    }
    Debug(Debug::INFO) << "Index table: removing duplicate entries...\n";
    indexTable->revertPointer();
    Debug(Debug::INFO) << "Index table init done.\n\n";

}

IndexTable* PrefilteringIndexReader::generateIndexTable(DBReader<unsigned int> *dbr, Sequence *seq, BaseMatrix *subMat,
                                                        int alphabetSize, int kmerSize, size_t dbFrom, size_t dbTo,
                                                        bool diagonalScoring, int maskMode, int kmerThr,
                                                        unsigned int threads) {

    IndexTable *indexTable = new IndexTable(alphabetSize, kmerSize, false);
    fillDatabase(dbr, seq, indexTable, subMat, dbFrom, dbTo, diagonalScoring, maskMode, NULL, kmerThr, threads);
    return indexTable;
}
