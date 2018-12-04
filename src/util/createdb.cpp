/*
 * createdb
 * written by Martin Steinegger <martin.steinegger@mpibpc.mpg.de>.
 * modified by Maria Hauser <mhauser@genzentrum.lmu.de> (splitting into sequences/headers databases)
 * modified by Milot Mirdita <milot@mirdita.de>
 */

#include <random>

#include "itoa.h"
#include "FileUtil.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Parameters.h"
#include "Util.h"
#include "KSeqWrapper.h"

int createdb(int argn, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argn, argv, command, 2, true, Parameters::PARSE_VARIADIC);

    if (par.maxSeqLen == Parameters::MAX_SEQ_LEN) {
        par.maxSeqLen = Parameters::MAX_SEQ_LEN - 1;
    }

    std::vector<std::string> filenames(par.filenames);

    std::string data_filename = filenames.back();
    filenames.pop_back();
    std::string index_filename = data_filename + ".index";

    std::string data_filename_hdr(data_filename);
    data_filename_hdr.append("_h");

    std::string index_filename_hdr(data_filename);
    index_filename_hdr.append("_h.index");

    std::string lookupFileName(data_filename);
    lookupFileName.append(".lookup");
    FILE* lookupFile = fopen(lookupFileName.c_str(), "w");
    if(lookupFile == NULL) {
        Debug(Debug::ERROR) << "Could not open " << lookupFileName << " for writing.";
        EXIT(EXIT_FAILURE);
    }

    for(size_t i = 0; i < filenames.size(); i++){
        if(FileUtil::fileExists(filenames[i].c_str())==false){
            Debug(Debug::ERROR) << "File " << filenames[i] << " does not exist.\n";
            EXIT(EXIT_FAILURE);
        }
        if(FileUtil::directoryExists(filenames[i].c_str())==true){
            Debug(Debug::ERROR) << "File " << filenames[i] << " is a directory.\n";
            EXIT(EXIT_FAILURE);
        }
    }

    DBWriter out_writer(data_filename.c_str(), index_filename.c_str());
    DBWriter out_hdr_writer(data_filename_hdr.c_str(), index_filename_hdr.c_str());
    out_writer.open();
    out_hdr_writer.open();

    unsigned int entries_num = 0;
    size_t count = 0;
    size_t sampleCount = 0;

    const size_t testForNucSequence = 100;
    size_t isNuclCnt = 0;
    bool assmeNuclDb = (par.dbType == 2) ? true : false;

    const char tab = '\t';
    char lookupBuffer[32768];

    // keep number of entries in each file
    unsigned int numEntriesInCurrFile = 0;
    unsigned int * fileToNumEntries =  new unsigned int[filenames.size()];
    for (size_t i = 0; i < filenames.size(); i++) {
        numEntriesInCurrFile = 0;
        std::string splitHeader;
        splitHeader.reserve(1024);
        std::string header;
        header.reserve(1024);
        std::string splitId;
        splitId.reserve(1024);
        KSeqWrapper *kseq = KSeqFactory(filenames[i].c_str());
        while (kseq->ReadEntry()) {
            Debug::printProgress(count);
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
                header.append(1, ' ');
                header.append(e.comment.s,e.comment.l);
            }

            std::string headerId = Util::parseFastaHeader(header);
            if (headerId == "") {
                // An identifier is necessary for these two cases, so we should just give up
                Debug(Debug::WARNING) << "Could not extract identifier from entry " << entries_num << ".\n";

            }
            for (size_t split = 0; split < splitCnt; split++) {
                unsigned int id = par.identifierOffset + entries_num;
                if(par.dbType == 0){
                    // check for the first 10 sequences if they are nucleotide sequences
                    if(count < 10 || (count % 100) == 0){
                        if(sampleCount < testForNucSequence){
                            size_t cnt=0;
                            for(size_t i = 0; i < e.sequence.l; i++){
                                switch(toupper(e.sequence.s[i]))
                                {
                                    case 'T':
                                    case 'A':
                                    case 'G':
                                    case 'C':
                                    case 'N': cnt++;
                                        break;
                                }
                            }
                            float nuclDNAFraction = static_cast<float>(cnt)/static_cast<float>(e.sequence.l);
                            if(nuclDNAFraction > 0.9){
                                isNuclCnt += true;
                            }
                        }
                        sampleCount++;
                    }
                    if (isNuclCnt == sampleCount || isNuclCnt == testForNucSequence) {
                        if(assmeNuclDb==false){
                            Debug(Debug::WARNING) << "Assume it is a DNA database.\n";
                            Debug(Debug::WARNING) << "Set parameter --dont-split-seq-by-len\n";
                            par.splitSeqByLen = false;
                            splitCnt = 1;
                        }
                        assmeNuclDb=true;
                    }else if(assmeNuclDb == true && isNuclCnt != sampleCount){
                        Debug(Debug::WARNING) << "Database does not look like a DNA database anymore. Sorry our prediction went wrong.\n";
                        Debug(Debug::WARNING) << "Please recompute with --dbtype 1 flag.\n";
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
                    if (headerId != "") {
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
                splitHeader.append(1, ' ');
                splitHeader.append(1, '\n');

                // Finally write down the entry
                out_hdr_writer.writeData(splitHeader.c_str(), splitHeader.length(), id);
                splitHeader.clear();
                splitId.clear();

                if (par.splitSeqByLen) {
                    size_t len = std::min(par.maxSeqLen, e.sequence.l - split * par.maxSeqLen);
                    out_writer.writeStart(0);
                    out_writer.writeAdd(e.sequence.s + split * par.maxSeqLen, len, 0);
                    char newLine = '\n';
                    out_writer.writeAdd(&newLine, 1, 0);
                    out_writer.writeEnd(id, 0, true);
                } else {
                    out_writer.writeStart(0);
                    out_writer.writeAdd(e.sequence.s, e.sequence.l, 0);
                    char newLine = '\n';
                    out_writer.writeAdd(&newLine, 1, 0);
                    out_writer.writeEnd(id, 0, true);
                }

                if (par.shuffleDatabase == false) {
                    char *tmpBuff = Itoa::i32toa_sse2(id, lookupBuffer);
                    *(tmpBuff-1) = '\t';
                    fwrite(lookupBuffer, sizeof(char), tmpBuff-lookupBuffer, lookupFile);
                    fwrite(headerId.c_str(), sizeof(char), headerId.length(), lookupFile);
                    fwrite(&tab, sizeof(char), 1, lookupFile);
                    tmpBuff = Itoa::i32toa_sse2(i, lookupBuffer);
                    *(tmpBuff-1) = '\n';
                    fwrite(lookupBuffer, sizeof(char), tmpBuff-lookupBuffer, lookupFile);
                }

                entries_num++;
                numEntriesInCurrFile++;
                count++;
            }
            header.clear();
        }
        fileToNumEntries[i] = numEntriesInCurrFile;
        delete kseq;
    }

    int dbType = Sequence::AMINO_ACIDS;
    if (par.dbType == 2 || (par.dbType == 0 && (isNuclCnt == sampleCount || isNuclCnt == testForNucSequence))) {
        dbType = Sequence::NUCLEOTIDES;
    }
    out_writer.close(dbType);
    out_hdr_writer.close();

    // shuffle data
    if (par.shuffleDatabase == true) {
        DBReader<unsigned int> readerSequence(out_writer.getDataFileName(), out_writer.getIndexFileName());
        readerSequence.open( DBReader<unsigned int>::NOSORT);
        readerSequence.readMmapedDataInMemory();
        DBReader<unsigned int>::Index * indexSequence  = readerSequence.getIndex();
        unsigned int * lengthSequence  = readerSequence.getSeqLens();

        DBReader<unsigned int> readerHeader(out_hdr_writer.getDataFileName(), out_hdr_writer.getIndexFileName());
        readerHeader.open( DBReader<unsigned int>::NOSORT);
        DBReader<unsigned int>::Index * indexHeader  = readerHeader.getIndex();
        unsigned int * lengthHeader  = readerHeader.getSeqLens();

        // Each file is identified by its first unshuffled entry
        // Intialize keyToFileAfterShuf to the sequential of the files
        unsigned int * keyToFileAfterShuf = new unsigned int[readerSequence.getSize()];
        unsigned int firstEntryOfFile = 0;
        for (size_t i = 0; i < filenames.size(); i++) {
            unsigned int numEntriesInCurrFile = fileToNumEntries[i];
            for (unsigned int j = firstEntryOfFile; j < (firstEntryOfFile + numEntriesInCurrFile); ++j ) {
                keyToFileAfterShuf[j] = i;
            }
            firstEntryOfFile += numEntriesInCurrFile;
        }
        delete [] fileToNumEntries;

        std::default_random_engine generator(0);
        std::uniform_int_distribution<unsigned int> distribution(0,readerSequence.getSize()-1);
        for (unsigned int n = 0; n < readerSequence.getSize(); n++) {
            unsigned int n_new = distribution(generator);
            std::swap(indexSequence[n_new], indexSequence[n]);
            std::swap(lengthSequence[n_new], lengthSequence[n]);
            std::swap(indexHeader[n_new], indexHeader[n]);
            std::swap(lengthHeader[n_new], lengthHeader[n]);
            std::swap(keyToFileAfterShuf[n_new], keyToFileAfterShuf[n]);
        }
        DBWriter out_writer_shuffled(data_filename.c_str(), index_filename.c_str());
        out_writer_shuffled.open();
        for (unsigned int n = 0; n < readerSequence.getSize(); n++) {
            unsigned int id = par.identifierOffset + n;
            const char * data = readerSequence.getData() + indexSequence[n].offset;
            out_writer_shuffled.writeData(data, lengthSequence[n]-1, id);
        }
        readerSequence.close();
        out_writer_shuffled.close(dbType);

        DBWriter out_hdr_writer_shuffled(data_filename_hdr.c_str(), index_filename_hdr.c_str());
        out_hdr_writer_shuffled.open();
        readerHeader.readMmapedDataInMemory();
        for (unsigned int n = 0; n < readerHeader.getSize(); n++) {
            unsigned int id = par.identifierOffset + n;
            const char *data = readerHeader.getData() + indexHeader[n].offset;
            std::string headerId = Util::parseFastaHeader(data);
            char *tmpBuff = Itoa::u32toa_sse2(id, lookupBuffer);
            *(tmpBuff-1) = '\t';
            fwrite(lookupBuffer, sizeof(char), tmpBuff-lookupBuffer, lookupFile);
            fwrite(headerId.c_str(), sizeof(char), headerId.length(), lookupFile);
            fwrite(&tab, sizeof(char), 1, lookupFile);
            tmpBuff = Itoa::u32toa_sse2(keyToFileAfterShuf[n], lookupBuffer);
            *(tmpBuff-1) = '\n';
            fwrite(lookupBuffer, sizeof(char), tmpBuff-lookupBuffer, lookupFile);
            out_hdr_writer_shuffled.writeData(data, lengthHeader[n]-1, id);
        }
        out_hdr_writer_shuffled.close();
        readerHeader.close();
        delete[] keyToFileAfterShuf;
    }
    fclose(lookupFile);

    return EXIT_SUCCESS;
}
