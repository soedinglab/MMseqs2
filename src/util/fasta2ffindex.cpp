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
int createdb(int argn, const char **argv) {
    std::string usage("Converts a fasta database to ffindex.\n");
    usage.append("USAGE: <fastaDB>  <ffindexDB> [mappingFasta]\n");
    usage.append("\nDesigned and implemented by Martin Steinegger <martin.steinegger@mpibpc.mpg.de>.\n");

    Parameters par;
    par.parseParameters(argn, argv, usage, par.createdb, 2);

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

    std::map<std::string, size_t> mapping;
    if (par.db3.length() > 0) {
        const char *mapping_filename = par.db3.c_str();
        mapping = Util::readMapping(mapping_filename);
    }

    size_t entries_num = 1;

    std::string header_line;
    header_line.reserve(10000);

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
        for(size_t split = 0; split < splitCnt; split++){
            std::string id;
            if (par.db3.length() > 0) {
                std::string key = Util::parseFastaHeader(seq->name.s);
                if (mapping.find(key) == mapping.end()) {
                    Debug(Debug::ERROR) << "Could not find entry: " << key << " in mapping file.\n";
                    return EXIT_FAILURE;
                }
                id = SSTR(mapping[key]);
                std::string headerId = Util::parseFastaHeader(seq->name.s);
                lookupStream << id << "\t" << headerId << "\n";
            } else if(par.useHeader) {
                id = Util::parseFastaHeader(seq->name.s);
                if (par.splitSeqByLen){
                    id.append("_");
                    id.append(SSTR(split));
                };
            } else {
                id = SSTR(par.identifierOffset + entries_num);
                std::string headerId = Util::parseFastaHeader(seq->name.s);
                lookupStream << id << "\t" << headerId << "\n";
            }

            // header
            header_line.append(seq->name.s, seq->name.l);
            if (seq->comment.l) {
                header_line.append(" ", 1);
                header_line.append(seq->comment.s, seq->comment.l);
            }
            if(par.splitSeqByLen == true) {
                header_line.append(" Split=");
                header_line.append(SSTR(split));
            }
            // space is needed for later parsing
            header_line.append(" ", 1);
            header_line.append("\n");

            out_hdr_writer.write(header_line.c_str(), header_line.length(), id.c_str());
            header_line.clear();

            // sequence
            std::string sequence = seq->seq.s;
            sequence.append("\n");
            size_t len = std::min(par.maxSeqLen, sequence.length() - split * par.maxSeqLen);
            out_writer.write(sequence.c_str() + (split*par.maxSeqLen), len, id.c_str());
            entries_num++;
        }
    }
    kseq_destroy(seq);

    if(!par.useHeader) {
        lookupStream.close();
    }
    fclose(fasta_file);
    out_hdr_writer.close();
    out_writer.close();

    return EXIT_SUCCESS;
}




