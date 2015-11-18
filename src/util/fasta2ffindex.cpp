/*
 * createdb
 * written by Martin Steinegger <hauser@genzentrum.lmu.de>.
 * modified by Maria Hauser <mhauser@genzentrum.lmu.de> (splitting into sequences/headers databases)
 * modified by Milot Mirdita <milot@mirdita.de>
 */

#define _GNU_SOURCE 1
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include <cstdio> // fclose

#include <string>
#include <map>

#include "Debug.h"
#include "Util.h"
#include "DBWriter.h"

#include "kseq.h"
#include <unistd.h> // read
KSEQ_INIT(int, read)

void usage() {
    Debug(Debug::INFO) <<  "Converts a fasta database to ffindex.\n";
    Debug(Debug::INFO) <<  "USAGE: <fastaDB>  <ffindexDB> [mappingFasta]\n";
    Debug(Debug::INFO) <<  "\nDesigned and implemented by Martin Steinegger <martin.steinegger@campus.lmu.de>.\n";
}

int createdb(int argn, const char **argv) {
    if (argn < 2 || argn > 3) {
        usage();
        return EXIT_FAILURE;
    }

    FILE *fasta_file = Util::openFileOrDie(argv[0], "r");

    const char *data_filename = argv[1];

    std::string index_filename_str(data_filename);
    index_filename_str.append(".index");

    std::string data_filename_hdr_str(data_filename);
    data_filename_hdr_str.append("_h");

    std::string index_filename_hdr_str(data_filename);
    index_filename_hdr_str.append("_h.index");

    DBWriter out_writer(data_filename, index_filename_str.c_str());
    DBWriter out_hdr_writer(data_filename_hdr_str.c_str(), index_filename_hdr_str.c_str());

    out_writer.open();
    out_hdr_writer.open();

    const char *mapping_filename = NULL;
    if (argn == 3) {
        mapping_filename = argv[2];
    }

    std::map<std::string, size_t> mapping;
    if (mapping_filename != NULL) {
        mapping = Util::readMapping(mapping_filename);
    }

    size_t entries_num = 0;

    std::string header_line;
    header_line.reserve(10000);

    kseq_t *seq = kseq_init(fileno(fasta_file));
    while (kseq_read(seq) >= 0) {
        if (seq->name.l == 0) {
            Debug(Debug::ERROR) << "Fasta entry: " << entries_num << " is invalid.\n";
            return EXIT_FAILURE;
        }
        std::string id;
        if (mapping_filename != NULL) {
            std::string key = Util::parseFastaHeader(seq->name.s);
            if (mapping.find(key) == mapping.end()) {
                Debug(Debug::ERROR) << "Could not find entry: " << key << " in mapping file.\n";
                return EXIT_FAILURE;
            }
            id = SSTR(mapping[key]);
        } else {
            id = SSTR(entries_num);
        }

        // header
        header_line.append(seq->name.s, seq->name.l);
        if (seq->comment.l) {
            header_line.append(" ", 1);
            header_line.append(seq->comment.s, seq->comment.l);
        }
        header_line.append("\n");

        out_hdr_writer.write(header_line.c_str(), header_line.length(), id.c_str());
        header_line.clear();

        // sequence
        std::string sequence = seq->seq.s;
        sequence.append("\n");

        out_writer.write(sequence.c_str(), sequence.length(), id.c_str());

        entries_num++;
    }
    kseq_destroy(seq);

    fclose(fasta_file);
    out_hdr_writer.close();
    out_writer.close();

    return EXIT_SUCCESS;
}



