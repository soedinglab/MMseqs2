/*
 * createdb
 * written by Martin Steinegger <martin.steinegger@mpibpc.mpg.de>.
 * modified by Maria Hauser <mhauser@genzentrum.lmu.de> (splitting into sequences/headers databases)
 * modified by Milot Mirdita <milot@mirdita.de>
 */

#include "FileUtil.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "KSeqWrapper.h"
#include "itoa.h"

int createdb(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    std::vector<std::string> filenames(par.filenames);
    std::string dataFile = filenames.back();
    filenames.pop_back();
    for (size_t i = 0; i < filenames.size(); i++) {
        if (FileUtil::directoryExists(filenames[i].c_str()) == true) {
            Debug(Debug::ERROR) << "File " << filenames[i] << " is a directory.\n";
            EXIT(EXIT_FAILURE);
        }
    }

    int dbType = -1;
    if (par.dbType == 2) {
        dbType = Parameters::DBTYPE_NUCLEOTIDES;
    } else if(par.dbType == 1){
        dbType = Parameters::DBTYPE_AMINO_ACIDS;
    }

    std::string indexFile = dataFile + ".index";
    if (par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT && par.shuffleDatabase) {
        Debug(Debug::WARNING) << "Shuffle database can not be combined with --createdb-mode 0.\n";
        Debug(Debug::WARNING) << "We recompute with --shuffle 0.\n";
        par.shuffleDatabase = false;
    }

    if (par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT && par.filenames[0] == "stdin") {
        Debug(Debug::WARNING) << "Createdb-mode 0 can not be combined with stdin input.\n";
        Debug(Debug::WARNING) << "We recompute with --createdb-mode 1.\n";
        par.createdbMode = Parameters::SEQUENCE_SPLIT_MODE_HARD;
    }

    const unsigned int shuffleSplits = par.shuffleDatabase ? 32 : 1;
    if (par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT && par.compressed) {
        Debug(Debug::WARNING) << "Compressed database can not be combined with --createdb-mode 0.\n";
        Debug(Debug::WARNING) << "We recompute with --compressed 0.\n";
        par.compressed = 0;
    }

    std::string hdrDataFile = dataFile + "_h";
    std::string hdrIndexFile = dataFile + "_h.index";

    unsigned int entries_num = 0;
    size_t sampleCount = 0;

    const char newline = '\n';

    const size_t testForNucSequence = 100;
    size_t isNuclCnt = 0;
    Debug::Progress progress;
    std::vector<unsigned short>* sourceLookup = new std::vector<unsigned short>[shuffleSplits]();
    for (size_t i = 0; i < shuffleSplits; ++i) {
        sourceLookup[i].reserve(16384);
    }
    Debug(Debug::INFO) << "Converting sequences\n";

    std::string sourceFile = dataFile + ".source";

    redoComputation:
    FILE *source = fopen(sourceFile.c_str(), "w");
    if (source == NULL) {
        Debug(Debug::ERROR) << "Can not open " << sourceFile << " for writing!\n";
        EXIT(EXIT_FAILURE);
    }
    DBWriter hdrWriter(hdrDataFile.c_str(), hdrIndexFile.c_str(), shuffleSplits, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    hdrWriter.open();
    DBWriter seqWriter(dataFile.c_str(), indexFile.c_str(), shuffleSplits, par.compressed, (dbType == -1) ? Parameters::DBTYPE_OMIT_FILE : dbType );
    seqWriter.open();
    for (size_t fileIdx = 0; fileIdx < filenames.size(); fileIdx++) {
        unsigned int numEntriesInCurrFile = 0;
        std::string header;
        header.reserve(1024);

        char buffer[4096];
        size_t len = snprintf(buffer, sizeof(buffer), "%zu\t%s\n", fileIdx, FileUtil::baseName(filenames[fileIdx]).c_str());
        int written = fwrite(buffer, sizeof(char), len, source);
        if (written != (int) len) {
            Debug(Debug::ERROR) << "Cannot write to source file " << sourceFile << "\n";
            EXIT(EXIT_FAILURE);
        }

        KSeqWrapper* kseq = KSeqFactory(filenames[fileIdx].c_str());
        while (kseq->ReadEntry()) {
            progress.updateProgress();
            const KSeqWrapper::KSeqEntry &e = kseq->entry;
            if (e.name.l == 0) {
                Debug(Debug::ERROR) << "Fasta entry: " << entries_num << " is invalid.\n";
                EXIT(EXIT_FAILURE);
            }

            // header
            if(par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_HARD){
                header.append(e.name.s, e.name.l);
                if (e.comment.l > 0) {
                    header.append(" ", 1);
                    header.append(e.comment.s, e.comment.l);
                }

                std::string headerId = Util::parseFastaHeader(header.c_str());
                if (headerId.empty()) {
                    // An identifier is necessary for these two cases, so we should just give up
                    Debug(Debug::WARNING) << "Can not extract identifier from entry " << entries_num << ".\n";
                }
                header.push_back('\n');
            }
            unsigned int id = par.identifierOffset + entries_num;
            if (dbType == -1) {
                // check for the first 10 sequences if they are nucleotide sequences
                if (sampleCount < 10 || (sampleCount % 100) == 0) {
                    if (sampleCount < testForNucSequence) {
                        size_t cnt = 0;
                        for (size_t i = 0; i < e.sequence.l; i++) {
                            switch (toupper(e.sequence.s[i])) {
                                case 'T':
                                case 'A':
                                case 'G':
                                case 'C':
                                case 'U':
                                case 'N':
                                    cnt++;
                                    break;
                            }
                        }
                        const float nuclDNAFraction = static_cast<float>(cnt) / static_cast<float>(e.sequence.l);
                        if (nuclDNAFraction > 0.9) {
                            isNuclCnt += true;
                        }
                    }
                    sampleCount++;
                }
                if(par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT && e.multiline == true){
                    Debug(Debug::WARNING) << "Multiline fasta can not be combined with --createdb-mode 0.\n";
                    Debug(Debug::WARNING) << "We recompute with --createdb-mode 1.\n";
                    par.createdbMode = Parameters::SEQUENCE_SPLIT_MODE_HARD;
                    progress.reset(SIZE_MAX);
                    hdrWriter.close();
                    seqWriter.close();
                    delete kseq;
                    fclose(source);
                    for (size_t i = 0; i < shuffleSplits; ++i) {
                        sourceLookup[i].clear();
                    }
                    goto redoComputation;
                }
            }

            // Finally write down the entry
            unsigned int splitIdx = id % shuffleSplits;
            sourceLookup[splitIdx].emplace_back(fileIdx);
            if(par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT){
                // +2 to emulate the \n\0
                hdrWriter.writeIndexEntry(id, e.headerOffset, (e.sequenceOffset-e.headerOffset)+1, 0);
                seqWriter.writeIndexEntry(id, e.sequenceOffset, e.sequence.l+2, 0);
            }else{
                hdrWriter.writeData(header.c_str(), header.length(), id, splitIdx);
                seqWriter.writeStart(splitIdx);
                seqWriter.writeAdd(e.sequence.s, e.sequence.l, splitIdx);
                seqWriter.writeAdd(&newline, 1, splitIdx);
                seqWriter.writeEnd(id, splitIdx, true);
            }

            entries_num++;
            numEntriesInCurrFile++;
            header.clear();
        }
        delete kseq;
    }
    Debug(Debug::INFO) << "\n";
    fclose(source);
    hdrWriter.close(true);
    seqWriter.close(true);
    if(dbType == -1) {
        if (isNuclCnt == sampleCount) {
            dbType = Parameters::DBTYPE_NUCLEOTIDES;
        } else {
            dbType = Parameters::DBTYPE_AMINO_ACIDS;
        }
        seqWriter.writeDbtypeFile(seqWriter.getDataFileName(), dbType ,par.compressed);
    }
    Debug(Debug::INFO) << "Database type: " << Parameters::getDbTypeName(dbType) << "\n";



    if(entries_num == 0){
        Debug(Debug::ERROR) << "The input files have no entry: ";
        for (size_t fileIdx = 0; fileIdx < filenames.size(); fileIdx++) {
            Debug(Debug::ERROR) << " - " << filenames[fileIdx] << "\n";
        }
        Debug(Debug::ERROR) << "Please check your input files."
                               "Only files in fasta/fastq[.gz|bz2] are supported \n";

        EXIT(EXIT_FAILURE);
    }

    // fix ids
    DBWriter::createRenumberedDB(dataFile, indexFile, "", "", DBReader<unsigned int>::LINEAR_ACCCESS);
    DBWriter::createRenumberedDB(hdrDataFile, hdrIndexFile, "", "", DBReader<unsigned int>::LINEAR_ACCCESS);
    if(par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT) {
        if(filenames.size() == 1){
            FileUtil::symlinkAbs(filenames[0], dataFile);
            FileUtil::symlinkAbs(filenames[0], hdrDataFile);
        }else{
            for (size_t fileIdx = 0; fileIdx < filenames.size(); fileIdx++) {
                FileUtil::symlinkAbs(filenames[fileIdx], dataFile+"."+SSTR(fileIdx));
                FileUtil::symlinkAbs(filenames[fileIdx], hdrDataFile+"."+SSTR(fileIdx));
            }
        }
    }
    DBReader<unsigned int> readerHeader(hdrDataFile.c_str(), hdrIndexFile.c_str(), 1, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    readerHeader.open(DBReader<unsigned int>::NOSORT);

    // create lookup file
    std::string lookupDataFile = dataFile + ".lookup";
    std::string lookupIndexFile = lookupDataFile + ".index";
    DBWriter lookupFile(lookupDataFile.c_str(), lookupIndexFile.c_str(), 1, Parameters::WRITER_ASCII_MODE, Parameters::DBTYPE_OMIT_FILE);
    lookupFile.open();

    char lookupBuffer[32768];
    unsigned int splitIdx = 0;
    unsigned int splitCounter = 0;
    DBReader<unsigned int>::LookupEntry entry;
    for (unsigned int id = 0; id < readerHeader.getSize(); id++) {
        size_t splitSize = sourceLookup[splitIdx].size();
        if (splitSize == 0 || splitCounter > sourceLookup[splitIdx].size() - 1) {
            splitIdx++;
            splitCounter = 0;
        }
        char *header = readerHeader.getData(id, 0);
        entry.id = id;
        entry.entryName = Util::parseFastaHeader(header);
        entry.fileNumber = sourceLookup[splitIdx][splitCounter];
        if (entry.entryName.empty()) {
            // An identifier is necessary for these two cases, so we should just give up
            Debug(Debug::WARNING) << "Can not extract identifier from entry " << entries_num << ".\n";
        }
        size_t len = readerHeader.lookupEntryToBuffer(lookupBuffer, entry);
        lookupFile.writeData(lookupBuffer, len, 0, 0, false, false);
        splitCounter++;
    }
    lookupFile.close(true);
    FileUtil::remove(lookupIndexFile.c_str());
    readerHeader.close();
    delete[] sourceLookup;

    return EXIT_SUCCESS;
}
