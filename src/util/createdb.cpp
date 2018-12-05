/*
 * createdb
 * written by Martin Steinegger <martin.steinegger@mpibpc.mpg.de>.
 * modified by Maria Hauser <mhauser@genzentrum.lmu.de> (splitting into sequences/headers databases)
 * modified by Milot Mirdita <milot@mirdita.de>
 */

#include <cstdio>

#include <map>
#include <fstream>
#include <unistd.h>
#include <math.h>
#include <itoa.h>
#include <random>

#include "FileUtil.h"
#include "DBWriter.h"
#include "Debug.h"
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
    std::string lookupFileNameIndex = lookupFileName+".index";

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

    KSeqWrapper *kseq = KSeqFactory(filenames[0].c_str());

    // check what kind of datbase it is
    bool isNuclDb = (par.dbType == 2) ? true : false;
    if(par.dbType == 0) {
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
            if(isNuclDb==false){
                Debug(Debug::WARNING) << "Assume it is a DNA database.\n";
                Debug(Debug::WARNING) << "Set parameter --dont-split-seq-by-len\n";
                par.splitSeqByLen = false;
            }
            isNuclDb=true;
        }
        delete kseq;
    }

    int dbType = Parameters::DBTYPE_AMINO_ACIDS;
    if (par.dbType == 2 ||  (par.dbType == 0 && isNuclDb == true) ) {
        dbType = Parameters::DBTYPE_NUCLEOTIDES;
    }
    int SPLITS_SHUFFEL = (par.shuffleDatabase) ? 32 : 1;
    DBWriter out_writer(data_filename.c_str(), index_filename.c_str(), SPLITS_SHUFFEL, par.compressed, dbType);
    DBWriter out_hdr_writer(data_filename_hdr.c_str(), index_filename_hdr.c_str(), SPLITS_SHUFFEL, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    DBWriter lookupFile(lookupFileName.c_str(), lookupFileNameIndex.c_str(), SPLITS_SHUFFEL, Parameters::WRITER_ASCII_MODE, Parameters::DBTYPE_GENERIC_DB);

    out_writer.open();
    out_hdr_writer.open();
    lookupFile.open();

    unsigned int entries_num = 0;
    size_t count = 0;
    size_t sampleCount = 0;

    const size_t testForNucSequence = 100;
    size_t isNuclCnt = 0;
    std::string dataStr;
    dataStr.reserve(1000000);
    // keep number of entries in each file
    unsigned int numEntriesInCurrFile = 0;
    unsigned int * fileToNumEntries =  new unsigned int[filenames.size()];
    for (size_t fileIdx = 0; fileIdx < filenames.size(); fileIdx++) {
        numEntriesInCurrFile = 0;
        std::string splitHeader;
        splitHeader.reserve(1024);
        std::string header;
        header.reserve(1024);
        std::string splitId;
        splitId.reserve(1024);
        char lookupBuffer[32768];
        KSeqWrapper *kseq = KSeqFactory(filenames[fileIdx].c_str());
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
                header.append(" ", 1);
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
                        if(isNuclDb==false){
                            Debug(Debug::WARNING) << "Assume it is a DNA database.\n";
                            Debug(Debug::WARNING) << "Set parameter --dont-split-seq-by-len\n";
                            par.splitSeqByLen = false;
                            splitCnt = 1;
                        }
                        isNuclDb=true;
                    }else if(isNuclDb == true && isNuclCnt != sampleCount){
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
                splitHeader.append(" ", 1);
                splitHeader.append("\n");

                // Finally write down the entry
                out_hdr_writer.writeData(splitHeader.c_str(), splitHeader.length(), id);

                size_t fileNo = count%SPLITS_SHUFFEL;
                lookupFile.writeStart(fileNo);

                char newline='\n';
                char tab='\t';
                char * tmpBuff = Itoa::u32toa_sse2(id, lookupBuffer);
                *(tmpBuff-1) = '\t';
                lookupFile.writeAdd(lookupBuffer, tmpBuff-lookupBuffer, fileNo);
                lookupFile.writeAdd(splitId.c_str(), splitId.length(), fileNo);
                lookupFile.writeAdd(&tab, 1, fileNo);
                tmpBuff = Itoa::u32toa_sse2(fileIdx, lookupBuffer);
                *(tmpBuff-1) = '\n';
                lookupFile.writeAdd(lookupBuffer, tmpBuff-lookupBuffer, fileNo);
                lookupFile.writeEnd(id, fileNo, false);


                if (par.splitSeqByLen) {
                    size_t len = std::min(par.maxSeqLen, e.sequence.l - split * par.maxSeqLen);
                    out_writer.writeStart(fileNo);
                    out_writer.writeAdd(e.sequence.s + split * par.maxSeqLen, len, fileNo);
                    out_writer.writeAdd(&newline, 1, fileNo);
                    out_writer.writeEnd(id, fileNo, true);
                } else {
                    out_writer.writeStart(fileNo);
                    out_writer.writeAdd(e.sequence.s, e.sequence.l, fileNo);
                    out_writer.writeAdd(&newline, 1, fileNo);
                    out_writer.writeEnd(id, fileNo, true);
                }
                splitHeader.clear();
                splitId.clear();

                entries_num++;
                numEntriesInCurrFile++;
                count++;
            }
            header.clear();
        }
        fileToNumEntries[fileIdx] = numEntriesInCurrFile;
        delete kseq;
    }

    out_hdr_writer.close();
    out_writer.close();
    lookupFile.close();
    FileUtil::deleteFile(lookupFileNameIndex);

    return EXIT_SUCCESS;
}
