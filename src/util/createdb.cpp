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
        if (FileUtil::fileExists(filenames[i].c_str()) == false) {
            Debug(Debug::ERROR) << "File " << filenames[i] << " does not exist.\n";
            EXIT(EXIT_FAILURE);
        }
        if (FileUtil::directoryExists(filenames[i].c_str()) == true) {
            Debug(Debug::ERROR) << "File " << filenames[i] << " is a directory.\n";
            EXIT(EXIT_FAILURE);
        }
    }

    KSeqWrapper *kseq = KSeqFactory(filenames[0].c_str());
    // check what kind of datbase it is
    bool isNuclDb = (par.dbType == 2) ? true : false;
    if (par.dbType == 0) {
        size_t isNuclCnt = 0;
        if (kseq->ReadEntry()) {
            const KSeqWrapper::KSeqEntry &e = kseq->entry;

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
            float nuclDNAFraction = static_cast<float>(cnt) / static_cast<float>(e.sequence.l);
            if (nuclDNAFraction > 0.9) {
                isNuclCnt += true;
            }
        }
        if (isNuclCnt) {
            isNuclDb = true;
        }
    }
    delete kseq;


    int dbType = Parameters::DBTYPE_AMINO_ACIDS;
    if (par.dbType == 2 || (par.dbType == 0 && isNuclDb == true)) {
        dbType = Parameters::DBTYPE_NUCLEOTIDES;
    }

    std::string indexFile = dataFile + ".index";
    if (par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT && par.shuffleDatabase) {
        Debug(Debug::WARNING) << "Shuffle database can not be combined with --createdb-mode 0.\n";
        Debug(Debug::WARNING) << "We recompute with --shuffle 0.\n";
        par.shuffleDatabase = false;
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
    size_t count = 0;
    size_t sampleCount = 0;

    const char newline = '\n';

    const size_t testForNucSequence = 100;
    size_t isNuclCnt = 0;
    std::string dataStr;
    dataStr.reserve(1000000);
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
    DBWriter seqWriter(dataFile.c_str(), indexFile.c_str(), shuffleSplits, par.compressed, dbType);
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

        kseq = KSeqFactory(filenames[fileIdx].c_str());
        while (kseq->ReadEntry()) {
            progress.updateProgress();
            const KSeqWrapper::KSeqEntry &e = kseq->entry;
            if (e.name.l == 0) {
                Debug(Debug::ERROR) << "Fasta entry: " << entries_num << " is invalid.\n";
                EXIT(EXIT_FAILURE);
            }

            // header
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
            unsigned int id = par.identifierOffset + entries_num;
            if (par.dbType == 0) {
                // check for the first 10 sequences if they are nucleotide sequences
                if (count < 10 || (count % 100) == 0) {
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
                bool redoComp = false;
                if (isNuclCnt == sampleCount || isNuclCnt == testForNucSequence) {
                    isNuclDb = true;
                } else if (isNuclDb == true && isNuclCnt != sampleCount) {
                    Debug(Debug::WARNING) << "Database does not look like a DNA database anymore.\n";
                    Debug(Debug::WARNING) << "We recompute as protein database.\n";
                    dbType = Parameters::DBTYPE_AMINO_ACIDS;
                    redoComp = true;
                }
                if(par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT && e.multiline == true){
                    Debug(Debug::WARNING) << "Multiline fasta can not be combined with --createdb-mode 0.\n";
                    Debug(Debug::WARNING) << "We recompute with --createdb-mode 1.\n";
                    par.createdbMode = Parameters::SEQUENCE_SPLIT_MODE_HARD;
                    redoComp = true;
                }
                if(redoComp){
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
            header.push_back('\n');
            if(par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT){
                // +2 to emulate the \n\0
                hdrWriter.writeIndexEntry(id, e.offset, header.size()+2, 0);
                seqWriter.writeIndexEntry(id, e.offset + header.size(), e.sequence.l+2, 0);
            }else{
                hdrWriter.writeData(header.c_str(), header.length(), id, splitIdx);
                seqWriter.writeStart(splitIdx);
                seqWriter.writeAdd(e.sequence.s, e.sequence.l, splitIdx);
                seqWriter.writeAdd(&newline, 1, splitIdx);
                seqWriter.writeEnd(id, splitIdx, true);
            }

            entries_num++;
            numEntriesInCurrFile++;
            count++;
            header.clear();
        }
        delete kseq;
    }
    Debug(Debug::INFO) << "\n";
    fclose(source);
    hdrWriter.close(true);
    seqWriter.close(true);


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
    DBWriter::createRenumberedDB(dataFile, indexFile, "", DBReader<unsigned int>::LINEAR_ACCCESS);
    DBWriter::createRenumberedDB(hdrDataFile, hdrIndexFile, "", DBReader<unsigned int>::LINEAR_ACCCESS);
    if(par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT) {
        for (size_t fileIdx = 0; fileIdx < filenames.size(); fileIdx++) {
            if(par.sequenceSplitMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT){
                FileUtil::symlinkAbs(filenames[0], dataFile+"."+SSTR(fileIdx));
                FileUtil::symlinkAbs(filenames[0], hdrDataFile+"."+SSTR(fileIdx));
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
