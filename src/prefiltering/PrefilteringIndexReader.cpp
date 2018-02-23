#include "PrefilteringIndexReader.h"
#include "DBWriter.h"
#include "Prefiltering.h"
#include "ExtendedSubstitutionMatrix.h"
#include "FileUtil.h"
#include "tantan.h"

const char*  PrefilteringIndexReader::CURRENT_VERSION="5.0.0";
unsigned int PrefilteringIndexReader::VERSION = 0;
unsigned int PrefilteringIndexReader::META = 1;
unsigned int PrefilteringIndexReader::SCOREMATRIXNAME = 2;
unsigned int PrefilteringIndexReader::SCOREMATRIX2MER = 3;
unsigned int PrefilteringIndexReader::SCOREMATRIX3MER = 4;
unsigned int PrefilteringIndexReader::DBRINDEX = 5;
unsigned int PrefilteringIndexReader::HDRINDEX = 6;

unsigned int PrefilteringIndexReader::ENTRIES = 7;
unsigned int PrefilteringIndexReader::ENTRIESOFFSETS = 8;
unsigned int PrefilteringIndexReader::ENTRIESNUM = 9;
unsigned int PrefilteringIndexReader::SEQCOUNT = 10;
unsigned int PrefilteringIndexReader::SEQINDEXDATA = 11;
unsigned int PrefilteringIndexReader::SEQINDEXDATASIZE = 12;
unsigned int PrefilteringIndexReader::SEQINDEXSEQOFFSET = 13;
unsigned int PrefilteringIndexReader::UNMASKEDSEQINDEXDATA = 14;
unsigned int PrefilteringIndexReader::GENERATOR = 15;

extern const char* version;

bool PrefilteringIndexReader::checkIfIndexFile(DBReader<unsigned int>* reader) {
    char * version = reader->getDataByDBKey(VERSION);
    if(version == NULL){
        return false;
    }
    return (strncmp(version, CURRENT_VERSION, strlen(CURRENT_VERSION)) == 0 ) ? true : false;
}

void PrefilteringIndexReader::createIndexFile(std::string outDB, DBReader<unsigned int> *dbr, DBReader<unsigned int> *hdbr,
                                              BaseMatrix * subMat, int maxSeqLen, bool hasSpacedKmer,
                                              bool compBiasCorrection, int alphabetSize, int kmerSize,
                                              bool diagonalScoring, int maskMode, int seqType, int kmerThr, int threads) {
    std::string outIndexName(outDB); // db.sk6
    std::string spaced = (hasSpacedKmer == true) ? "s" : "";
    outIndexName.append(".").append(spaced).append("k").append(SSTR(kmerSize));

    DBWriter writer(outIndexName.c_str(), std::string(outIndexName).append(".index").c_str(), 1, DBWriter::BINARY_MODE);
    writer.open();

    if (seqType != Sequence::HMM_PROFILE||seqType != Sequence::PROFILE_STATE_SEQ ) {
        int alphabetSize = subMat->alphabetSize;
        subMat->alphabetSize = subMat->alphabetSize-1;
        ScoreMatrix *s3 = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 3);
        ScoreMatrix *s2 = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 2);
        subMat->alphabetSize = alphabetSize;
        char* serialized3mer = ScoreMatrix::serialize(*s3);
        Debug(Debug::INFO) << "Write SCOREMATRIX3MER (" << SCOREMATRIX3MER << ")\n";
        writer.writeData(serialized3mer, ScoreMatrix::size(*s3), SCOREMATRIX3MER, 0);
        writer.alignToPageSize();
        free(serialized3mer);
        ScoreMatrix::cleanup(s3);

        char* serialized2mer = ScoreMatrix::serialize(*s2);
        Debug(Debug::INFO) << "Write SCOREMATRIX2MER (" << SCOREMATRIX2MER << ")\n";
        writer.writeData(serialized2mer, ScoreMatrix::size(*s2), SCOREMATRIX2MER, 0);
        writer.alignToPageSize();
        free(serialized2mer);
        ScoreMatrix::cleanup(s2);
    }

    Sequence seq(maxSeqLen, seqType, subMat, kmerSize, hasSpacedKmer, compBiasCorrection);
    // remove x (not needed in index)
    int adjustAlphabetSize = (seqType == Sequence::NUCLEOTIDES || seqType == Sequence::AMINO_ACIDS)
                             ? alphabetSize -1: alphabetSize;
    IndexTable *indexTable = new IndexTable(adjustAlphabetSize, kmerSize, false);

    SequenceLookup *unmaskedLookup = NULL;
    PrefilteringIndexReader::fillDatabase(dbr, &seq, indexTable, subMat,
                                          0, dbr->getSize(), diagonalScoring,
                                          maskMode, &unmaskedLookup, kmerThr, threads);

    indexTable->printStatistics(subMat->int2aa);

    // save the entries
    Debug(Debug::INFO) << "Write ENTRIES (" << ENTRIES << ")\n";
    char *entries = (char *) indexTable->getEntries();
    size_t entriesSize = indexTable->getTableEntriesNum() * indexTable->getSizeOfEntry();
    writer.writeData(entries, entriesSize, ENTRIES, 0);
    writer.alignToPageSize();

    // save the size
    Debug(Debug::INFO) << "Write ENTRIESOFFSETS (" << ENTRIESOFFSETS << ")\n";

    char *offsets = (char*)indexTable->getOffsets();
    size_t offsetsSize = (indexTable->getTableSize() + 1) * sizeof(size_t);
    writer.writeData(offsets, offsetsSize, ENTRIESOFFSETS, 0);
    writer.alignToPageSize();
    indexTable->deleteEntries();

    SequenceLookup *lookup = indexTable->getSequenceLookup();
    Debug(Debug::INFO) << "Write SEQINDEXDATA (" << SEQINDEXDATA << ")\n";
    writer.writeData(lookup->getData(), (lookup->getDataSize() + 1) * sizeof(char), SEQINDEXDATA, 0);
    writer.alignToPageSize();

    if (unmaskedLookup != NULL) {
        Debug(Debug::INFO) << "Write UNMASKEDSEQINDEXDATA (" << UNMASKEDSEQINDEXDATA << ")\n";
        writer.writeData(unmaskedLookup->getData(), (unmaskedLookup->getDataSize() + 1) * sizeof(char), UNMASKEDSEQINDEXDATA, 0);
        writer.alignToPageSize();
        delete unmaskedLookup;
    }

    Debug(Debug::INFO) << "Write SEQINDEXDATASIZE (" << SEQINDEXDATASIZE << ")\n";
    int64_t seqindexDataSize = lookup->getDataSize();
    char *seqindexDataSizePtr = (char *) &seqindexDataSize;
    writer.writeData(seqindexDataSizePtr, 1 * sizeof(int64_t), SEQINDEXDATASIZE, 0);
    writer.alignToPageSize();

    size_t *sequenceOffsets = lookup->getOffsets();
    size_t sequenceCount = lookup->getSequenceCount();
    Debug(Debug::INFO) << "Write SEQINDEXSEQOFFSET (" << SEQINDEXSEQOFFSET << ")\n";
    writer.writeData((char *) sequenceOffsets, (sequenceCount + 1) * sizeof(size_t), SEQINDEXSEQOFFSET, 0);
    writer.alignToPageSize();

    // meta data
    // ENTRIESNUM
    Debug(Debug::INFO) << "Write ENTRIESNUM (" << ENTRIESNUM << ")\n";
    uint64_t entriesNum = indexTable->getTableEntriesNum();
    char *entriesNumPtr = (char *) &entriesNum;
    writer.writeData(entriesNumPtr, 1 * sizeof(uint64_t), ENTRIESNUM, 0);
    writer.alignToPageSize();
    // SEQCOUNT
    Debug(Debug::INFO) << "Write SEQCOUNT (" << SEQCOUNT << ")\n";
    size_t tablesize = {indexTable->getSize()};
    char *tablesizePtr = (char *) &tablesize;
    writer.writeData(tablesizePtr, 1 * sizeof(size_t), SEQCOUNT, 0);
    writer.alignToPageSize();

    delete indexTable;

    Debug(Debug::INFO) << "Write META (" << META << ")\n";
    int local = 1;
    int spacedKmer = (hasSpacedKmer) ? 1 : 0;
    int headers = (hdbr != NULL) ? 1 : 0;
    int metadata[] = {kmerSize, alphabetSize, maskMode, local, spacedKmer, kmerThr, seqType, headers};
    char *metadataptr = (char *) &metadata;
    writer.writeData(metadataptr, sizeof(metadata), META, 0);
    writer.alignToPageSize();
    printMeta(metadata);

    Debug(Debug::INFO) << "Write SCOREMATRIXNAME (" << SCOREMATRIXNAME << ")\n";
    writer.writeData(subMat->getMatrixName().c_str(), subMat->getMatrixName().length(), SCOREMATRIXNAME, 0);
    writer.alignToPageSize();

    Debug(Debug::INFO) << "Write VERSION (" << VERSION << ")\n";
    writer.writeData((char *) CURRENT_VERSION, strlen(CURRENT_VERSION) * sizeof(char), VERSION, 0);
    writer.alignToPageSize();

    Debug(Debug::INFO) << "Write DBRINDEX (" << DBRINDEX << ")\n";
    char* data = DBReader<unsigned int>::serialize(*dbr);
    writer.writeData(data, DBReader<unsigned int>::indexMemorySize(*dbr), DBRINDEX, 0);
    writer.alignToPageSize();
    free(data);

    if (hdbr != NULL) {
        Debug(Debug::INFO) << "Write HDRINDEX (" << HDRINDEX << ")\n";
        data = DBReader<unsigned int>::serialize(*hdbr);
        writer.writeData(data, DBReader<unsigned int>::indexMemorySize(*hdbr), HDRINDEX, 0);
        writer.alignToPageSize();
        free(data);
    }

    Debug(Debug::INFO) << "Write GENERATOR (" << GENERATOR << ")\n";
    writer.writeData(version, strlen(version), GENERATOR, 0);
    writer.alignToPageSize();

    writer.close();
    Debug(Debug::INFO) << "Done. \n";
}

DBReader<unsigned int> *PrefilteringIndexReader::openNewHeaderReader(DBReader<unsigned int>*dbr, const char* dataFileName, bool touch) {
    size_t id = dbr->getId(HDRINDEX);
    char *data = dbr->getData(id);
    if (touch) {
        dbr->touchData(id);
    }

    DBReader<unsigned int> *reader = DBReader<unsigned int>::unserialize(data);
    reader->setDataFile(dataFileName);
    reader->open(DBReader<unsigned int>::NOSORT);

    return reader;
}

DBReader<unsigned int> *PrefilteringIndexReader::openNewReader(DBReader<unsigned int>*dbr, bool touch) {
    size_t id = dbr->getId(DBRINDEX);
    char *data = dbr->getData(id);
    if (touch) {
        dbr->touchData(id);
    }

    DBReader<unsigned int> *reader = DBReader<unsigned int>::unserialize(data);
    reader->open(DBReader<unsigned int>::NOSORT);

    return reader;
}

SequenceLookup *PrefilteringIndexReader::getSequenceLookup(DBReader<unsigned int>*dbr, bool touch) {
    size_t id;
    if ((id = dbr->getId(SEQINDEXDATA)) == UINT_MAX) {
        return NULL;
    }

    char * seqData = dbr->getData(id);

    size_t seqOffsetsId = dbr->getId(SEQINDEXSEQOFFSET);
    char * seqOffsetsData = dbr->getData(seqOffsetsId);

    size_t seqDataSizeId = dbr->getId(SEQINDEXDATASIZE);
    int64_t seqDataSize = *((int64_t *)dbr->getData(seqOffsetsId));

    size_t sequenceCountId = dbr->getId(SEQCOUNT);
    size_t sequenceCount = *((size_t *)dbr->getData(sequenceCountId));

    if (touch) {
        dbr->touchData(id);
        dbr->touchData(seqOffsetsId);
    }

    SequenceLookup *sequenceLookup = new SequenceLookup(sequenceCount);
    sequenceLookup->initLookupByExternalData(seqData, seqDataSize, (size_t *) seqOffsetsData);

    return sequenceLookup;
}

SequenceLookup *PrefilteringIndexReader::getUnmaskedSequenceLookup(DBReader<unsigned int>*dbr, bool touch) {
    size_t id;
    if ((id = dbr->getId(UNMASKEDSEQINDEXDATA)) == UINT_MAX) {
        return NULL;
    }

    char * seqData = dbr->getData(id);

    size_t seqOffsetsId = dbr->getId(SEQINDEXSEQOFFSET);
    char * seqOffsetsData = dbr->getData(seqOffsetsId);

    size_t seqDataSizeId = dbr->getId(SEQINDEXDATASIZE);
    int64_t seqDataSize = *((int64_t *)dbr->getData(seqOffsetsId));

    size_t sequenceCountId = dbr->getId(SEQCOUNT);
    size_t sequenceCount = *((size_t *)dbr->getData(sequenceCountId));

    if (touch) {
        dbr->touchData(id);
        dbr->touchData(seqOffsetsId);
    }

    SequenceLookup *sequenceLookup = new SequenceLookup(sequenceCount);
    sequenceLookup->initLookupByExternalData(seqData, seqDataSize, (size_t *) seqOffsetsData);

    return sequenceLookup;
}

IndexTable *PrefilteringIndexReader::generateIndexTable(DBReader<unsigned int> *dbr, bool diagonalScoring, bool touch) {
    PrefilteringIndexData data = getMetadata(dbr);
    IndexTable *retTable;
    if (data.local) {
        int adjustAlphabetSize = (data.seqType == Sequence::NUCLEOTIDES || data.seqType == Sequence::AMINO_ACIDS)
                                 ? (data.alphabetSize - 1) : data.alphabetSize;
        retTable = new IndexTable(adjustAlphabetSize, data.kmerSize, true);
    }else {
        Debug(Debug::ERROR) << "Search mode is not valid.\n";
        EXIT(EXIT_FAILURE);
    }

    SequenceLookup * sequenceLookup = NULL;
    if(diagonalScoring == true) {
        sequenceLookup = getSequenceLookup(dbr, touch);
    }

    size_t entriesNumId = dbr->getId(ENTRIESNUM);
    int64_t entriesNum = *((int64_t *)dbr->getData(entriesNumId));
    size_t sequenceCountId = dbr->getId(SEQCOUNT);
    size_t sequenceCount = *((size_t *)dbr->getData(sequenceCountId));

    size_t entriesDataId = dbr->getId(ENTRIES);
    char *entriesData = dbr->getData(entriesDataId);

    size_t entriesOffsetsDataId = dbr->getId(ENTRIESOFFSETS);
    char *entriesOffsetsData = dbr->getData(entriesOffsetsDataId);

    if (touch) {
        dbr->touchData(entriesNumId);
        dbr->touchData(sequenceCountId);
        dbr->touchData(entriesDataId);
        dbr->touchData(entriesOffsetsDataId);
    }

    retTable->initTableByExternalData(sequenceCount, entriesNum,
                                      (IndexEntryLocal*) entriesData, (size_t *)entriesOffsetsData, sequenceLookup);
    return retTable;
}

void PrefilteringIndexReader::printMeta(int *metadata_tmp) {
    Debug(Debug::INFO) << "KmerSize:     " << metadata_tmp[0] << "\n";
    Debug(Debug::INFO) << "AlphabetSize: " << metadata_tmp[1] << "\n";
    Debug(Debug::INFO) << "Mask:         " << metadata_tmp[2] << "\n";
    Debug(Debug::INFO) << "Type:         " << metadata_tmp[3] << "\n";
    Debug(Debug::INFO) << "Spaced:       " << metadata_tmp[4] << "\n";
    Debug(Debug::INFO) << "KmerScore:    " << metadata_tmp[5] << "\n";
    Debug(Debug::INFO) << "SequenceType: " << metadata_tmp[6] << "\n";
    Debug(Debug::INFO) << "Headers:      " << metadata_tmp[7] << "\n";
}

void PrefilteringIndexReader::printSummary(DBReader<unsigned int> *dbr) {
    Debug(Debug::INFO) << "Index version: " << dbr->getDataByDBKey(VERSION) << "\n";

    size_t id;
    if ((id = dbr->getId(GENERATOR)) != UINT_MAX) {
        Debug(Debug::INFO) << "Generated by:  " << dbr->getData(id) << "\n";
    }

    int *metadata_tmp = (int *) dbr->getDataByDBKey(META);

    printMeta(metadata_tmp);

    Debug(Debug::INFO) << "ScoreMatrix:  " << dbr->getDataByDBKey(SCOREMATRIXNAME) << "\n";
}

PrefilteringIndexData PrefilteringIndexReader::getMetadata(DBReader<unsigned int> *dbr) {
    PrefilteringIndexData prefData;

    int *metadata_tmp = (int *) dbr->getDataByDBKey(META);

    prefData.kmerSize = metadata_tmp[0];
    prefData.alphabetSize = metadata_tmp[1];
    prefData.maskMode = metadata_tmp[2];
    prefData.local = metadata_tmp[3];
    prefData.spacedKmer = metadata_tmp[4];
    prefData.kmerThr = metadata_tmp[5];
    prefData.seqType = metadata_tmp[6];
    prefData.headers = metadata_tmp[7];

    return prefData;
}

std::string PrefilteringIndexReader::getSubstitutionMatrixName(DBReader<unsigned int> *dbr) {
    return std::string(dbr->getDataByDBKey(SCOREMATRIXNAME));
}
//
ScoreMatrix *PrefilteringIndexReader::get2MerScoreMatrix(DBReader<unsigned int> *dbr, bool touch) {
    PrefilteringIndexData meta = getMetadata(dbr);
    size_t id = dbr->getId(SCOREMATRIX2MER);
    char *data = dbr->getData(id);

    if (touch) {
        dbr->touchData(id);
    }

    return ScoreMatrix::unserialize(data, meta.alphabetSize-1, 2);
}

ScoreMatrix *PrefilteringIndexReader::get3MerScoreMatrix(DBReader<unsigned int> *dbr, bool touch) {
    PrefilteringIndexData meta = getMetadata(dbr);
    size_t id = dbr->getId(SCOREMATRIX3MER);
    char *data = dbr->getData(id);

    if (touch) {
        dbr->touchData(id);
    }
    // remove x (not needed in index)
    return ScoreMatrix::unserialize(data, meta.alphabetSize-1, 3);
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
    dbTo = std::min(dbTo, dbr->getSize());
    size_t maskedResidues = 0;
    size_t totalKmerCount = 0;
    size_t tableSize = 0;

    size_t dbSize = dbTo - dbFrom;
    size_t *sequenceOffSet = new size_t[dbSize];

    const bool isProfileStateSeq = seq->getSeqType() == Sequence::PROFILE_STATE_SEQ;

    // identical scores for memory reduction code
    char *idScoreLookup = new char[subMat->alphabetSize];
    if(isProfileStateSeq == false){
        for (int aa = 0; aa < subMat->alphabetSize; aa++){
            short score = subMat->subMatrix[aa][aa];
            if (score > CHAR_MAX || score < CHAR_MIN) {
                Debug(Debug::WARNING) << "Truncating substitution matrix diagonal score!";
            }
            idScoreLookup[aa] = (char) score;
        }
    }
    // need to prune low scoring k-mers
    // masking code
    double probMatrix[subMat->alphabetSize][subMat->alphabetSize];
    const double **probMatrixPointers = new const double*[subMat->alphabetSize];
    char hardMaskTable[256];
    if(isProfileStateSeq == false) {
        std::fill_n(hardMaskTable, 256, subMat->aa2int[(int) 'X']);
        for (int i = 0; i < subMat->alphabetSize; ++i) {
            probMatrixPointers[i] = probMatrix[i];
            for (int j = 0; j < subMat->alphabetSize; ++j) {
                probMatrix[i][j] = subMat->probMatrix[i][j] / (subMat->pBack[i] * subMat->pBack[j]);
            }
        }
    }

    const bool isProfile = seq->getSeqType() == Sequence::HMM_PROFILE;
    size_t aaDbSize = 0;
    sequenceOffSet[0] = 0;
    for (size_t id = dbFrom; id < dbTo; id++) {
        int seqLen;
        if (isProfile) {
            // remove /0 and convert to profile length
            seqLen = std::max(static_cast<int>(dbr->getSeqLens(id)) - 1, 0) / static_cast<int>(Sequence::PROFILE_AA_SIZE);
        } else {
            // remove /n and /0
            seqLen = std::max(static_cast<int>(dbr->getSeqLens(id)) - 2, 0);
        }
        aaDbSize += seqLen;
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
        Indexer idxer(static_cast<unsigned int>(indexTable->getAlphabetSize()), seq->getKmerSize());
        Sequence s(seq->getMaxLen(), seq->getSeqType(), subMat, seq->getKmerSize(), seq->isSpaced(), false);

        KmerGenerator *generator = NULL;
        if (isProfile) {
            generator = new KmerGenerator(seq->getKmerSize(), indexTable->getAlphabetSize(), kmerThr);
            generator->setDivideStrategy(s.profile_matrix);
        }

        unsigned int * buffer = new unsigned int[seq->getMaxLen()];
        char * charSequence = new char[seq->getMaxLen()];

#pragma omp for schedule(dynamic, 100) reduction(+:totalKmerCount, maskedResidues)
        for (size_t id = dbFrom; id < dbTo; id++) {
            s.resetCurrPos();
            Debug::printProgress(id - dbFrom);
            char *seqData = dbr->getData(id);
            unsigned int qKey = dbr->getDbKey(id);
            s.mapSequence(id - dbFrom, qKey, seqData);

            if (maskMode == 2) {

                if (isProfile) {
                    (*unmaskedLookup)->addSequence(s.int_consensus_sequence, s.L, id - dbFrom, sequenceOffSet[id - dbFrom]);
                }else{
                    (*unmaskedLookup)->addSequence(s.int_sequence, s.L, id - dbFrom, sequenceOffSet[id - dbFrom]);
                }
            }

            // count similar or exact k-mers based on sequence type
            if (isProfile) {
                totalKmerCount += indexTable->addSimilarKmerCount(&s, generator, &idxer, kmerThr, idScoreLookup);
            } else { // Find out if we should also mask profiles
                if (maskMode == 1 || maskMode == 2) {
                    // Do not mask if column state sequences are used
                    if(isProfileStateSeq == false) {
                        for (int i = 0; i < s.L; ++i) {
                            charSequence[i] = (char) s.int_sequence[i];
                        }
//                        s.print();
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
                    }
                }


//                if (!isProfile) {
                totalKmerCount += indexTable->addKmerCount(&s, &idxer, buffer, kmerThr, idScoreLookup);
//                }
            }
            if (isProfile) {
                sequenceLookup->addSequence(s.int_consensus_sequence, s.L, id - dbFrom, sequenceOffSet[id - dbFrom]);
            }else{
                sequenceLookup->addSequence(s.int_sequence, s.L, id - dbFrom, sequenceOffSet[id - dbFrom]);
            }
        }

        delete [] charSequence;
        delete [] buffer;

        if (generator != NULL) {
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
        Sequence s(seq->getMaxLen(), seq->getSeqType(), subMat, seq->getKmerSize(), seq->isSpaced(), false);
        Indexer idxer(static_cast<unsigned int>(indexTable->getAlphabetSize()), seq->getKmerSize());
        IndexEntryLocalTmp * buffer = new IndexEntryLocalTmp[seq->getMaxLen()];

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
#pragma omp parallel for
    for (size_t i = 0; i < indexTable->getTableSize(); i++) {
        size_t entrySize;
        IndexEntryLocal * entries = indexTable->getDBSeqList(i, &entrySize);
        std::sort(entries, entries + entrySize, IndexEntryLocal::comapreByIdAndPos);
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
