/*
 * createdb
 * written by Martin Steinegger <hauser@genzentrum.lmu.de>.
 * modified by Maria Hauser <mhauser@genzentrum.lmu.de> (splitting into sequences/headers databases)
 * modified by Milot Mirdita <milot@mirdita.de>
 */

#include <cstdio>

#include <map>
#include <fstream>
#include <unistd.h>
#include <math.h>

#include "DBWriter.h"
#include "Debug.h"
#include "Parameters.h"
#include "Util.h"
#include "FileUtil.h"

#include "kseq.h"

KSEQ_INIT(int, read)
int createdb(int argn, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argn, argv, command, 2);

    if(par.maxSeqLen == Parameters::MAX_SEQ_LEN){
        par.maxSeqLen = Parameters::MAX_SEQ_LEN - 1;
    }

    FILE *fasta_file = FileUtil::openFileOrDie(par.db1.c_str(), "r", true);

    std::string data_filename = par.db2;
    std::string index_filename = par.db2Index;

    std::string data_filename_hdr(data_filename);
    data_filename_hdr.append("_h");

    std::string index_filename_hdr(data_filename);
    index_filename_hdr.append("_h.index");

    std::ofstream lookupStream;
    if(!par.useHeader) {
        std::string lookupFile = par.db2;
        lookupFile.append(".lookup");
        lookupStream.open(lookupFile);
        if(lookupStream.fail()) {
            Debug(Debug::ERROR) << "Could not open " << lookupFile << " for writing.";
            EXIT(EXIT_FAILURE);
        }
    }

    DBWriter out_writer(data_filename.c_str(), index_filename.c_str());
    DBWriter out_hdr_writer(data_filename_hdr.c_str(), index_filename_hdr.c_str());
    out_writer.open();
    out_hdr_writer.open();

    bool doMapping = par.db3.length() > 0;
    std::map<std::string, size_t> mapping;
    if (doMapping) {
        const char *mapping_filename = par.db3.c_str();
        mapping = Util::readMapping(mapping_filename);
    }

    size_t entries_num = 1;

    kseq_t *seq = kseq_init(fileno(fasta_file));
    while (kseq_read(seq) >= 0) {
        if (seq->name.l == 0) {
            Debug(Debug::ERROR) << "Fasta entry: " << entries_num << " is invalid.\n";
            return EXIT_FAILURE;
        }

        size_t splitCnt = 1;
        if(par.splitSeqByLen == true){
            splitCnt = (size_t) ceilf(static_cast<float>(seq->seq.l) / static_cast<float>(par.maxSeqLen));
        }

        // header
        std::string header(seq->name.s, seq->name.l);
        if (seq->comment.l > 0) {
            header.append(" ", 1);
            header.append(seq->comment.s, seq->comment.l);
        }

        std::string headerId = Util::parseFastaHeader(header);
        if(headerId == "") {
            // An identifier is necessary for these two cases, so we should just give up
            if(par.useHeader || doMapping) {
                Debug(Debug::ERROR) << "Could not extract identifier from entry " << entries_num << "!.\n";
                return EXIT_FAILURE;
            } else {
                Debug(Debug::WARNING) << "Could not extract identifier from entry " << entries_num << ".\n";
            }
        }

        for(size_t split = 0; split < splitCnt; split++){
            std::string splitId(headerId);
            if (splitCnt > 1){
                splitId.append("_");
                splitId.append(SSTR(split));
            }

            std::string id;
            if (doMapping) {
                if (mapping.find(splitId) == mapping.end()) {
                    Debug(Debug::ERROR) << "Could not find entry: " << splitId << " in mapping file.\n";
                    return EXIT_FAILURE;
                }
                id = SSTR(mapping[splitId]);
            } else if(par.useHeader) {
                id = splitId;
            } else {
                id = SSTR(par.identifierOffset + entries_num);
            }

            if(par.useHeader == false) {
                lookupStream << id << "\t" << splitId << "\n";
            }

            // For split entries replace the found identifier by identifier_splitNumber
            // Also add another hint that it was split to the end of the header
            std::string splitHeader(header);
            if(par.splitSeqByLen == true && splitCnt > 1) {
                if(headerId != "") {
                    size_t pos = splitHeader.find(headerId);
                    if(pos != std::string::npos) {
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
            out_hdr_writer.writeData(splitHeader.c_str(), splitHeader.length(), id.c_str());

            // sequence
            std::string sequence = seq->seq.s;
            size_t len = std::min(par.maxSeqLen, sequence.length() - split * par.maxSeqLen);
            std::string splitString(sequence.c_str() + split*par.maxSeqLen, len);
            splitString.append("\n");
            out_writer.writeData(splitString.c_str(), splitString.length(), id.c_str());

            entries_num++;
        }
    }
    kseq_destroy(seq);

    if(par.useHeader == false) {
        lookupStream.close();
    }
    fclose(fasta_file);
    out_hdr_writer.close();
    out_writer.close();

    return EXIT_SUCCESS;
}




