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

    unsigned int entries_num = 1;
    size_t count = 0;
    const size_t testForNucSequence = 10;
    size_t isNuclCnt = 0;

    for (size_t i = 0; i < filenames.size(); i++) {
        std::string splitHeader;
        splitHeader.reserve(1024);
        std::string header;
        header.reserve(1024);
        std::string splitId;
        splitId.reserve(1024);
        char lookupBuffer[32768];
        KSeqWrapper *kseq = KSeqFactory(filenames[i].c_str());
        while (kseq->ReadEntry()) {
            Debug::printProgress(count);
            const KSeqWrapper::KSeqEntry &e = kseq->entry;
            if (e.name.length() == 0) {
                Debug(Debug::ERROR) << "Fasta entry: " << entries_num << " is invalid.\n";
                EXIT(EXIT_FAILURE);
            }

            size_t splitCnt = 1;
            if (par.splitSeqByLen == true) {
                splitCnt = (size_t) ceilf(static_cast<float>(e.sequence.length()) / static_cast<float>(par.maxSeqLen));
            }

            // header
            header.append(e.name);
            if (e.comment.length() > 0) {
                header.append(" ", 1);
                header.append(e.comment);
            }

            std::string headerId = Util::parseFastaHeader(header);
            if (headerId == "") {
                // An identifier is necessary for these two cases, so we should just give up
                Debug(Debug::WARNING) << "Could not extract identifier from entry " << entries_num << ".\n";

            }
            for (size_t split = 0; split < splitCnt; split++) {
                splitId.append(headerId);
                if (splitCnt > 1) {
                    splitId.append("_");
                    splitId.append(SSTR(split));
                }

                unsigned int id = par.identifierOffset + entries_num;
                char * tmpBuff = Itoa::i32toa_sse2(id, lookupBuffer);
                *(tmpBuff-1) = '\t';
                fwrite(lookupBuffer, sizeof(char), tmpBuff-lookupBuffer, lookupFile);
                fwrite(splitId.c_str(), sizeof(char), splitId.length(), lookupFile);
                char newline='\n';
                fwrite(&newline, sizeof(char), 1, lookupFile);

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
                splitHeader.clear();
                splitId.clear();

                // sequence
                const std::string &sequence = e.sequence;
                // check for the first 10 sequences if they are nucleotide sequences
                if(count < testForNucSequence){
                    size_t cnt=0;
                    for(size_t i = 0; i < sequence.size(); i++){
                        switch(toupper(sequence[i]))
                        {
                            case 'T':
                            case 'A':
                            case 'G':
                            case 'C':
                            case 'N': cnt++;
                            break;
                        }
                    }
                    if(cnt == sequence.size()){
                        isNuclCnt += true;
                    }
                }


                if(par.splitSeqByLen){
                    size_t len = std::min(par.maxSeqLen, sequence.length() - split * par.maxSeqLen);
                    std::string splitString(sequence.c_str() + split * par.maxSeqLen, len);
                    splitString.append("\n");
                    out_writer.writeData(splitString.c_str(), splitString.length(), id);
                }else{
                    out_writer.writeData(sequence.c_str(), sequence.length(), id);
                }

                entries_num++;
                count++;
            }
            header.clear();
        }
        delete kseq;
    }

    int dbType = DBReader<unsigned int>::DBTYPE_AA;
    if(isNuclCnt==count || isNuclCnt == testForNucSequence){
        dbType = DBReader<unsigned int>::DBTYPE_NUC;
    }
    fclose(lookupFile);
    out_hdr_writer.close();
    out_writer.close(dbType);

    return EXIT_SUCCESS;
}




