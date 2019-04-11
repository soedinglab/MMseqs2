#include "LinsearchIndexReader.h"
#include <FileUtil.h>
#include <PrefilteringIndexReader.h>
#include "Debug.h"
#include "Timer.h"
#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"
#include "ReducedMatrix.h"
#include "KmerIndex.h"
#include "kmersearch.h"

extern const char* version;


int kmerindexdb(int argc, const char **argv, const Command &command) {
    MMseqsMPI::init(argc, argv);

    Parameters &par = Parameters::getInstance();
    setLinearFilterDefault(&par);
    par.parseParameters(argc, argv, command, 2, false, 0, MMseqsParameter::COMMAND_CLUSTLINEAR);

    DBReader<unsigned int> seqDbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    seqDbr.open(DBReader<unsigned int>::NOSORT);
    int querySeqType = seqDbr.getDbtype();

    setKmerLengthAndAlphabet(par, seqDbr.getAminoAcidDBSize(), querySeqType);
    std::vector<MMseqsParameter*>* params = command.params;
    par.printParameters(command.cmd, argc, argv, *params);
    Debug(Debug::INFO) << "Database size: "  << seqDbr.getSize() << " type: " << seqDbr.getDbTypeName() << "\n";
    std::string indexDB = LinsearchIndexReader::indexName(par.db2);
    if (par.checkCompatible && FileUtil::fileExists(indexDB.c_str())) {
        Debug(Debug::INFO) << "Check index " << indexDB << "\n";
        DBReader<unsigned int> index(indexDB.c_str(), (indexDB + ".index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        index.open(DBReader<unsigned int>::NOSORT);
        if (LinsearchIndexReader::checkIfIndexFile(&index) && LinsearchIndexReader::isIndexCompatible(index, par, seqDbr.getDbtype())) {
            Debug(Debug::INFO) << "Index is already up to date and compatible. Force recreation with --check-compatibility 0 parameter.\n";
            return EXIT_SUCCESS;
        } else {
            Debug(Debug::WARNING) << "Index is incompatible and will be recreated.\n";
        }
    }

    BaseMatrix *subMat;
    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.c_str(), 1.0, 0.0);
    }else {
        if (par.alphabetSize == 21) {
            subMat = new SubstitutionMatrix(par.scoringMatrixFile.c_str(), 2.0, 0.0);
        } else {
            SubstitutionMatrix sMat(par.scoringMatrixFile.c_str(), 2.0, 0.0);
            subMat = new ReducedMatrix(sMat.probMatrix, sMat.subMatrixPseudoCounts, sMat.aa2int, sMat.int2aa, sMat.alphabetSize, par.alphabetSize, 2.0);
        }
    }

    //seqDbr.readMmapedDataInMemory();
    const size_t KMER_SIZE = par.kmerSize;
    size_t chooseTopKmer = par.kmersPerSequence;

    size_t memoryLimit;
    if (par.splitMemoryLimit > 0) {
        memoryLimit = static_cast<size_t>(par.splitMemoryLimit) * 1024;
    } else {
        memoryLimit = static_cast<size_t>(Util::getTotalSystemMemory() * 0.9);
    }
    Debug(Debug::INFO) << "\n";
    size_t totalKmers = computeKmerCount(seqDbr, KMER_SIZE, chooseTopKmer);
    size_t totalSizeNeeded = computeMemoryNeededLinearfilter(totalKmers);
    Debug(Debug::INFO) << "Estimated memory consumption " << totalSizeNeeded/1024/1024 << " MB\n";
    // compute splits
    size_t splits = static_cast<size_t>(std::ceil(static_cast<float>(totalSizeNeeded) / memoryLimit));
//    size_t splits = 2;
    if (splits > 1) {
        // security buffer
        splits += 1;
    }

    Debug(Debug::INFO) << "Process file into " << splits << " parts\n";
    std::vector<std::string> splitFiles;
    KmerPosition *hashSeqPair = NULL;

    size_t writePos = 0;
    size_t mpiRank = 0;
    size_t adjustedKmerSize = KMER_SIZE;
#ifdef HAVE_MPI
    splits = std::max(static_cast<size_t>(MMseqsMPI::numProc), splits);
    size_t fromSplit = 0;
    size_t splitCount = 1;
    mpiRank = MMseqsMPI::rank;
    // if split size is great than nodes than we have to
    // distribute all splits equally over all nodes
    unsigned int * splitCntPerProc = new unsigned int[MMseqsMPI::numProc];
    memset(splitCntPerProc, 0, sizeof(unsigned int) * MMseqsMPI::numProc);
    for(size_t i = 0; i < splits; i++){
        splitCntPerProc[i % MMseqsMPI::numProc] += 1;
    }
    for(int i = 0; i < MMseqsMPI::rank; i++){
        fromSplit += splitCntPerProc[i];
    }
    splitCount = splitCntPerProc[MMseqsMPI::rank];
    delete[] splitCntPerProc;

    for(size_t split = fromSplit; split < fromSplit+splitCount; split++) {
        std::string splitFileName = par.db2 + "_split_" +SSTR(split);
        size_t splitKmerCount = (splits > 1) ? static_cast<size_t >(static_cast<double>(totalKmers/splits) * 1.2) : totalKmers;
        KmerSearch::ExtractKmerAndSortResult kmerRet  = KmerSearch::extractKmerAndSort(splitKmerCount, split, splits, seqDbr, par, subMat, KMER_SIZE, chooseTopKmer, 1, par.adjustKmerLength);
        hashSeqPair = kmerRet.kmers;
        // assign rep. sequence to same kmer members
        // The longest sequence is the first since we sorted by kmer, seq.Len and id
        if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
            writePos = LinsearchIndexReader::pickCenterKmer<Parameters::DBTYPE_NUCLEOTIDES>(hashSeqPair, splitKmerCount);
        }else{
            writePos = LinsearchIndexReader::pickCenterKmer<Parameters::DBTYPE_AMINO_ACIDS>(hashSeqPair, splitKmerCount);
        }

        LinsearchIndexReader::writeKmerIndexToDisk(splitFileName, hashSeqPair, writePos);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(mpiRank == 0){
        for(size_t split = 0; split < splits; split++) {
            std::string splitFileName = par.db2 + "_split_" +SSTR(split);
            splitFiles.push_back(splitFileName);
        }
    }
#else
    for(size_t split = 0; split < splits; split++) {
        std::string splitFileName = par.db2 + "_split_" +SSTR(split);
        size_t splitKmerCount = (splits > 1) ? static_cast<size_t >(static_cast<double>(totalKmers/splits) * 1.2) : totalKmers;
        KmerSearch::ExtractKmerAndSortResult kmerRet = KmerSearch::extractKmerAndSort(splitKmerCount, split, splits, seqDbr, par, subMat, KMER_SIZE, chooseTopKmer, 1, par.adjustKmerLength);
        hashSeqPair = kmerRet.kmers;
        adjustedKmerSize = std::max(adjustedKmerSize, kmerRet.adjustedKmer);
        // assign rep. sequence to same kmer members
        // The longest sequence is the first since we sorted by kmer, seq.Len and id
        if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
            writePos = LinsearchIndexReader::pickCenterKmer<Parameters::DBTYPE_NUCLEOTIDES>(hashSeqPair, splitKmerCount);
        }else{
            writePos = LinsearchIndexReader::pickCenterKmer<Parameters::DBTYPE_AMINO_ACIDS>(hashSeqPair, splitKmerCount);
        }

        if(splits > 1){
            LinsearchIndexReader::writeKmerIndexToDisk(splitFileName, hashSeqPair, writePos);
            delete [] hashSeqPair;
            hashSeqPair = NULL;
        }

        splitFiles.push_back(splitFileName);
    }
#endif
    if(mpiRank == 0){
        // write result
        DBWriter dbw(indexDB.c_str(), (indexDB+".index").c_str(), 1, par.compressed, Parameters::DBTYPE_INDEX_DB );
        dbw.open();

        Debug(Debug::INFO) << "Write VERSION (" << PrefilteringIndexReader::VERSION << ")\n";
        dbw.writeData((char *) PrefilteringIndexReader::CURRENT_VERSION, strlen(PrefilteringIndexReader::CURRENT_VERSION) * sizeof(char), PrefilteringIndexReader::VERSION, 0);
        dbw.alignToPageSize();

        Debug(Debug::INFO) << "Write META (" << PrefilteringIndexReader::META << ")\n";
        const int mask = par.maskMode > 0;
        const int spacedKmer = (par.spacedKmer) ? 1 : 0;
        const bool sameDB = (par.db1 == par.db2);
        const int headers1 =  1;
        const int headers2 = (sameDB) ? 1 : 0;
        const int seqType = seqDbr.getDbtype();
        const int srcSeqType = DBReader<unsigned int>::parseDbType(par.db2.c_str());
        // Reuse the compBiasCorr field to store the adjustedKmerSize, It is not needed in the linsearch
        int metadata[] = {static_cast<int>(par.maxSeqLen), static_cast<int>(KMER_SIZE), static_cast<int>(adjustedKmerSize), subMat->alphabetSize, mask, spacedKmer, 0, seqType, srcSeqType, headers1, headers2};
        char *metadataptr = (char *) &metadata;
        dbw.writeData(metadataptr, sizeof(metadata), PrefilteringIndexReader::META, 0);
        dbw.alignToPageSize();

        Timer timer;
        if(splits > 1) {
            seqDbr.unmapData();
            if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                LinsearchIndexReader::mergeAndWriteIndex<Parameters::DBTYPE_NUCLEOTIDES>(dbw, splitFiles, subMat->alphabetSize, adjustedKmerSize);
            }else{
                LinsearchIndexReader::mergeAndWriteIndex<Parameters::DBTYPE_AMINO_ACIDS>(dbw, splitFiles, subMat->alphabetSize, adjustedKmerSize);
            }
        } else {
            if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                LinsearchIndexReader::writeIndex<Parameters::DBTYPE_NUCLEOTIDES>(dbw, hashSeqPair, writePos, subMat->alphabetSize, adjustedKmerSize);
            }else{
                LinsearchIndexReader::writeIndex<Parameters::DBTYPE_AMINO_ACIDS>(dbw, hashSeqPair, writePos, subMat->alphabetSize, adjustedKmerSize);
            }
        }
        if(hashSeqPair){
            delete [] hashSeqPair;
            hashSeqPair = NULL;
        }
        // SEQCOUNT
        Debug(Debug::INFO) << "Write SEQCOUNT (" << PrefilteringIndexReader::SEQCOUNT << ")\n";
        size_t tablesize = {seqDbr.getSize()};
        char *tablesizePtr = (char *) &tablesize;
        dbw.writeData(tablesizePtr, 1 * sizeof(size_t), PrefilteringIndexReader::SEQCOUNT, 0);
        dbw.alignToPageSize();

        Debug(Debug::INFO) << "Write SCOREMATRIXNAME (" << PrefilteringIndexReader::SCOREMATRIXNAME << ")\n";
        char* subData = BaseMatrix::serialize(subMat);
        dbw.writeData(subData, BaseMatrix::memorySize(subMat), PrefilteringIndexReader::SCOREMATRIXNAME, 0);
        dbw.alignToPageSize();
        free(subData);

        if (par.spacedKmerPattern.empty() != false) {
            Debug(Debug::INFO) << "Write SPACEDPATTERN (" << PrefilteringIndexReader::SPACEDPATTERN << ")\n";
            dbw.writeData(par.spacedKmerPattern.c_str(), par.spacedKmerPattern.length(), PrefilteringIndexReader::SPACEDPATTERN, 0);
            dbw.alignToPageSize();
        }

        seqDbr.close();

        DBReader<unsigned int> dbr1(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        dbr1.open(DBReader<unsigned int>::NOSORT);
        Debug(Debug::INFO) << "Write DBR1INDEX (" << PrefilteringIndexReader::DBR1INDEX << ")\n";
        char* data = DBReader<unsigned int>::serialize(dbr1);
        size_t offsetIndex = dbw.getOffset(0);
        dbw.writeData(data, DBReader<unsigned int>::indexMemorySize(dbr1), PrefilteringIndexReader::DBR1INDEX, 0);
        dbw.alignToPageSize();

        Debug(Debug::INFO) << "Write DBR1DATA (" << PrefilteringIndexReader::DBR1DATA << ")\n";
        size_t offsetData = dbw.getOffset(0);
        dbw.writeStart(0);
        for(size_t fileIdx = 0; fileIdx < dbr1.getDataFileCnt(); fileIdx++) {
            dbw.writeAdd(dbr1.getDataForFile(fileIdx), dbr1.getDataSizeForFile(fileIdx), 0);
        }
        dbw.writeEnd( PrefilteringIndexReader::DBR1DATA, 0);
        dbw.alignToPageSize();
        free(data);

        if (sameDB == true) {
            dbw.writeIndexEntry(PrefilteringIndexReader::DBR2INDEX, offsetIndex, DBReader<unsigned int>::indexMemorySize(dbr1)+1, 0);
            dbw.writeIndexEntry(PrefilteringIndexReader::DBR2DATA,  offsetData,  dbr1.getTotalDataSize()+1, 0);
            dbr1.close();
        }else{
            dbr1.close();
            DBReader<unsigned int> dbr2(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
            dbr2.open(DBReader<unsigned int>::NOSORT);
            Debug(Debug::INFO) << "Write DBR2INDEX (" << PrefilteringIndexReader::DBR2INDEX << ")\n";
            data = DBReader<unsigned int>::serialize(dbr2);
            dbw.writeData(data, DBReader<unsigned int>::indexMemorySize(dbr2), PrefilteringIndexReader::DBR2INDEX, 0);
            dbw.alignToPageSize();
            Debug(Debug::INFO) << "Write DBR2DATA (" << PrefilteringIndexReader::DBR2DATA << ")\n";
            dbw.writeStart(0);
            for(size_t fileIdx = 0; fileIdx < dbr2.getDataFileCnt(); fileIdx++) {
                dbw.writeAdd(dbr2.getDataForFile(fileIdx), dbr2.getDataSizeForFile(fileIdx), 0);
            }
            dbw.writeEnd(PrefilteringIndexReader::DBR2DATA, 0);
            dbw.alignToPageSize();
            free(data);
            dbr2.close();
        }

        {
            Debug(Debug::INFO) << "Write HDR1INDEX (" << PrefilteringIndexReader::HDR1INDEX << ")\n";
            DBReader<unsigned int> hdbr1(par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
            hdbr1.open(DBReader<unsigned int>::NOSORT);

            data = DBReader<unsigned int>::serialize(hdbr1);
            size_t offsetIndex = dbw.getOffset(0);
            dbw.writeData(data, DBReader<unsigned int>::indexMemorySize(hdbr1), PrefilteringIndexReader::HDR1INDEX, 0);
            dbw.alignToPageSize();
            Debug(Debug::INFO) << "Write HDR1DATA (" << PrefilteringIndexReader::HDR1DATA << ")\n";
            size_t offsetData = dbw.getOffset(0);
            dbw.writeStart(0);
            for(size_t fileIdx = 0; fileIdx < hdbr1.getDataFileCnt(); fileIdx++) {
                dbw.writeAdd(hdbr1.getDataForFile(fileIdx), hdbr1.getDataSizeForFile(fileIdx), 0);
            }
            dbw.writeEnd(PrefilteringIndexReader::HDR1DATA, 0);
            dbw.alignToPageSize();
            free(data);
            if (sameDB == true) {
                dbw.writeIndexEntry(PrefilteringIndexReader::HDR2INDEX, offsetIndex, DBReader<unsigned int>::indexMemorySize(hdbr1)+1, 0);
                dbw.writeIndexEntry(PrefilteringIndexReader::HDR2DATA,  offsetData, hdbr1.getTotalDataSize()+1, 0);
                hdbr1.close();
            }else{
                hdbr1.close();
                DBReader<unsigned int> hdbr2(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
                hdbr2.open(DBReader<unsigned int>::NOSORT);
                Debug(Debug::INFO) << "Write HDR2INDEX (" <<PrefilteringIndexReader::HDR2INDEX << ")\n";
                data = DBReader<unsigned int>::serialize(hdbr2);
                dbw.writeData(data, DBReader<unsigned int>::indexMemorySize(hdbr2), PrefilteringIndexReader::HDR2INDEX, 0);
                dbw.alignToPageSize();
                Debug(Debug::INFO) << "Write HDR2DATA (" << PrefilteringIndexReader::HDR2DATA << ")\n";
                dbw.writeStart(0);
                for(size_t fileIdx = 0; fileIdx < hdbr2.getDataFileCnt(); fileIdx++) {
                    dbw.writeAdd(hdbr2.getDataForFile(fileIdx), hdbr2.getDataSizeForFile(fileIdx), 0);
                }
                dbw.writeEnd(PrefilteringIndexReader::HDR2DATA, 0);
                dbw.alignToPageSize();
                free(data);
            }
        }

        Debug(Debug::INFO) << "Write GENERATOR (" << PrefilteringIndexReader::GENERATOR << ")\n";
        dbw.writeData(version, strlen(version), PrefilteringIndexReader::GENERATOR, 0);
        dbw.alignToPageSize();

        Debug(Debug::INFO) << "Time for fill: " << timer.lap() << "\n";
        // add missing entries to the result (needed for clustering)
        dbw.close();
    }
    // free memory
    delete subMat;
    if(mpiRank != 0){
        if(hashSeqPair){
            delete [] hashSeqPair;
        }
        seqDbr.close();
    }

    return EXIT_SUCCESS;
}

