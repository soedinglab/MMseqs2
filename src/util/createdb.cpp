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

void renumberIdsInIndexByOffsetOrder(char * dataName, char * indexName) {
    DBReader<unsigned int> reader(dataName, indexName, 1, DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    std::string newSequenceIndex = std::string(indexName) + ".0";
    FILE *indexFile = fopen(newSequenceIndex.c_str(), "w");
    if (indexFile == NULL) {
        perror(newSequenceIndex.c_str());
        EXIT(EXIT_FAILURE);
    }

    char buffer[1024];
    for (unsigned int id = 0; id < reader.getSize(); id++) {
        DBReader<unsigned int>::Index *idx = reader.getIndex(id);
        idx->id = id;
        DBWriter::writeIndexEntryToFile(indexFile, buffer, *idx, reader.getSeqLens(id));
    }

    fclose(indexFile);
    reader.close();

    int result = rename(newSequenceIndex.c_str(), indexName);
    if (result != 0) {
        perror("Error renaming file");
    }
}

int createdb(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2, true, Parameters::PARSE_VARIADIC);

    if (par.maxSeqLen == Parameters::MAX_SEQ_LEN) {
        par.maxSeqLen = Parameters::MAX_SEQ_LEN - 1;
    }


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
            if (isNuclDb == false) {
                Debug(Debug::WARNING) << "Assuming DNA database, forcing parameter --dont-split-seq-by-len true\n";
                par.splitSeqByLen = false;
            }
            isNuclDb = true;
        }
    }
    delete kseq;


    int dbType = Parameters::DBTYPE_AMINO_ACIDS;
    if (par.dbType == 2 || (par.dbType == 0 && isNuclDb == true)) {
        dbType = Parameters::DBTYPE_NUCLEOTIDES;
    }

    std::string indexFile = dataFile + ".index";
    const unsigned int shuffleSplits = par.shuffleDatabase ? 32 : 1;
    DBWriter seqWriter(dataFile.c_str(), indexFile.c_str(), shuffleSplits, par.compressed, dbType);
    seqWriter.open();

    std::string hdrDataFile = dataFile + "_h";
    std::string hdrIndexFile = dataFile + "_h.index";
    DBWriter hdrWriter(hdrDataFile.c_str(), hdrIndexFile.c_str(), shuffleSplits, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    hdrWriter.open();

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

    for (size_t fileIdx = 0; fileIdx < filenames.size(); fileIdx++) {
        unsigned int numEntriesInCurrFile = 0;
        std::string splitHeader;
        splitHeader.reserve(1024);
        std::string header;
        header.reserve(1024);
        std::string splitId;
        splitId.reserve(1024);
        kseq = KSeqFactory(filenames[fileIdx].c_str());
        while (kseq->ReadEntry()) {
            progress.updateProgress();
            const KSeqWrapper::KSeqEntry &e = kseq->entry;
            if (e.name.l == 0) {
                Debug(Debug::ERROR) << "Fasta entry: " << entries_num << " is invalid.\n";
                EXIT(EXIT_FAILURE);
            }

            size_t splitCnt = 1;
            if (par.splitSeqByLen == true) {
                splitCnt = (size_t) ceilf(static_cast<float>(e.sequence.l) / static_cast<float>(par.maxSeqLen));
            }

            // header
            header.append(e.name.s, e.name.l);
            if (e.comment.l > 0) {
                header.append(" ", 1);
                header.append(e.comment.s, e.comment.l);
            }

            std::string headerId = Util::parseFastaHeader(header);
            if (headerId.empty()) {
                // An identifier is necessary for these two cases, so we should just give up
                Debug(Debug::WARNING) << "Can not extract identifier from entry " << entries_num << ".\n";

            }
            for (size_t split = 0; split < splitCnt; split++) {
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
                    if (isNuclCnt == sampleCount || isNuclCnt == testForNucSequence) {
                        if (isNuclDb == false) {
                            Debug(Debug::WARNING) << "Assume it is a DNA database.\n";
                            Debug(Debug::WARNING) << "Set parameter --dont-split-seq-by-len\n";
                            par.splitSeqByLen = false;
                            splitCnt = 1;
                        }
                        isNuclDb = true;
                    } else if (isNuclDb == true && isNuclCnt != sampleCount) {
                        Debug(Debug::ERROR) << "Database does not look like a DNA database anymore. Sorry our prediction went wrong.\n";
                        Debug(Debug::ERROR) << "Please recompute with --dbtype 1 flag.\n";
                        EXIT(EXIT_FAILURE);
                    }
                }

                splitId.append(headerId);
                if (splitCnt > 1) {
                    splitId.append("_");
                    splitId.append(SSTR(split));
                }
                // For split entries replace the found identifier by identifier_splitNumber
                // Also add another hint that it was split to the end of the header
                splitHeader.append(header);
                if (par.splitSeqByLen == true && splitCnt > 1) {
                    if (headerId.empty() == false) {
                        size_t pos = splitHeader.find(headerId);
                        if (pos != std::string::npos) {
                            splitHeader.erase(pos, headerId.length());
                            splitHeader.insert(pos, splitId);
                        }
                    }
                    splitHeader.append(" Split=");
                    splitHeader.append(SSTR(split));
                }

                // space is needed for later parsing
                splitHeader.append(" ", 1);
                splitHeader.append("\n");

                // Finally write down the entry
                unsigned int splitIdx = id % shuffleSplits;
                sourceLookup[splitIdx].emplace_back(fileIdx);

                hdrWriter.writeData(splitHeader.c_str(), splitHeader.length(), id, splitIdx);

                if (par.splitSeqByLen) {
                    size_t len = std::min(par.maxSeqLen, e.sequence.l - split * par.maxSeqLen);
                    seqWriter.writeStart(splitIdx);
                    seqWriter.writeAdd(e.sequence.s + split * par.maxSeqLen, len, splitIdx);
                    seqWriter.writeAdd(&newline, 1, splitIdx);
                    seqWriter.writeEnd(id, splitIdx, true);
                } else {
                    seqWriter.writeStart(splitIdx);
                    seqWriter.writeAdd(e.sequence.s, e.sequence.l, splitIdx);
                    seqWriter.writeAdd(&newline, 1, splitIdx);
                    seqWriter.writeEnd(id, splitIdx, true);
                }
                splitHeader.clear();
                splitId.clear();

                entries_num++;
                numEntriesInCurrFile++;
                count++;
            }
            header.clear();
        }
        delete kseq;
    }
    Debug(Debug::INFO) << "\n";
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
    renumberIdsInIndexByOffsetOrder(seqWriter.getDataFileName(), seqWriter.getIndexFileName());
    renumberIdsInIndexByOffsetOrder(hdrWriter.getDataFileName(), hdrWriter.getIndexFileName());

    DBReader<unsigned int> readerHeader(hdrWriter.getDataFileName(), hdrWriter.getIndexFileName(), 1, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    readerHeader.open(DBReader<unsigned int>::NOSORT);

    // create lookup file
    std::string lookupDataFile = dataFile + ".lookup";
    std::string lookupIndexFile = lookupDataFile + ".index";
    DBWriter lookupFile(lookupDataFile.c_str(), lookupIndexFile.c_str(), 1, Parameters::WRITER_ASCII_MODE, Parameters::DBTYPE_OMIT_FILE);
    lookupFile.open();

    char lookupBuffer[32768];
    const char tab = '\t';
    unsigned int splitIdx = 0;
    unsigned int splitCounter = 0;
    for (unsigned int id = 0; id < readerHeader.getSize(); id++) {
        size_t splitSize = sourceLookup[splitIdx].size();
        if (splitSize == 0 || splitCounter > sourceLookup[splitIdx].size() - 1) {
            splitIdx++;
            splitCounter = 0;
        }

        char *header = readerHeader.getData(id, 0);
        std::string headerId = Util::parseFastaHeader(header);
        if (headerId.empty()) {
            // An identifier is necessary for these two cases, so we should just give up
            Debug(Debug::WARNING) << "Can not extract identifier from entry " << entries_num << ".\n";
        }
        lookupFile.writeStart(0);
        char *tmpBuff = Itoa::u32toa_sse2(id, lookupBuffer);
        *(tmpBuff - 1) = '\t';
        lookupFile.writeAdd(lookupBuffer, tmpBuff - lookupBuffer, 0);
        lookupFile.writeAdd(headerId.c_str(), headerId.length(), 0);
        lookupFile.writeAdd(&tab, 1, 0);
        tmpBuff = Itoa::u32toa_sse2(sourceLookup[splitIdx][splitCounter], lookupBuffer);
        *(tmpBuff - 1) = '\n';
        lookupFile.writeAdd(lookupBuffer, tmpBuff - lookupBuffer, 0);
        lookupFile.writeEnd(id, 0, false);

        splitCounter++;
    }
    lookupFile.close(true);
    FileUtil::remove(lookupIndexFile.c_str());
    delete[] sourceLookup;

    return EXIT_SUCCESS;
}
