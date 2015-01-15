#define _GNU_SOURCE 1
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

#include "Orf.h"
#include "Util.h"

#include "kseq.h"
#define MAX_FILENAME_LIST_FILES 4096

KSEQ_INIT(int, read)

void printUsageExtractOrf()
{
    fprintf(stderr, "Converts a fasta database to ffindex.\n");
    fprintf(stderr, "USAGE: <fastaDB>  <ffindexDB>\n"
            "\nDesigned and implemented by Martin Steinegger <martin.steinegger@campus.lmu.de>.\n");
}

int extractorf(int argn, const char **argv)
{
    int err = EXIT_SUCCESS;

    if(argn  <  3)
    {
        printUsageExtractOrf();
        return EXIT_FAILURE;
    }

    char *fasta_filename = (char *)   argv[optind++];
    char *data_filename  = (char *)   argv[optind++];

    std::string index_filename_str(data_filename);
    index_filename_str.append(".index");
    char *index_filename = (char *) index_filename_str.c_str();

    std::string data_filename_hdr_str(data_filename);
    data_filename_hdr_str.append("_h");
    char *data_filename_hdr  = (char *)data_filename_hdr_str.c_str();

    std::string index_filename_hdr_str(data_filename);
    index_filename_hdr_str.append("_h.index");
    char *index_filename_hdr = (char *)index_filename_hdr_str.c_str();

    FILE *data_file, *index_file, *fasta_file, *data_file_hdr, *index_file_hdr;
    struct stat st;

    if(stat(data_filename, &st) == 0) { errno = EEXIST; perror(data_filename); return EXIT_FAILURE; }
    data_file  = fopen(data_filename, "w+");
    if( data_file == NULL) { perror(data_filename); return EXIT_FAILURE; }

    if(stat(index_filename, &st) == 0) { errno = EEXIST; perror(index_filename); return EXIT_FAILURE; }
    index_file = fopen(index_filename, "w+");
    if(index_file == NULL) { perror(index_filename); return EXIT_FAILURE; }

    fasta_file = fopen(fasta_filename, "r");
    if(fasta_file == NULL) { perror(fasta_filename); return EXIT_FAILURE; }

    if(stat(data_filename_hdr, &st) == 0) { errno = EEXIST; perror(data_filename_hdr); return EXIT_FAILURE; }
    data_file_hdr  = fopen(data_filename_hdr, "w+");
    if( data_file_hdr == NULL) { perror(data_filename_hdr); return EXIT_FAILURE; }

    if(stat(index_filename_hdr, &st) == 0) { errno = EEXIST; perror(index_filename_hdr); return EXIT_FAILURE; }
    index_file_hdr = fopen(index_filename_hdr, "w+");
    if(index_file_hdr == NULL) { perror(index_filename_hdr); return EXIT_FAILURE; }

    kseq_t *seq;
    seq = kseq_init(fileno(fasta_file));
    size_t offset_header = 0;
    size_t offset_sequence = 0;
    size_t entries_num = 0;
    std::string header_line;

    std::vector<std::string> start;
    start.push_back("ATG");
    start.push_back("AUG");

    header_line.reserve(10000);
    while (kseq_read(seq) >= 0) {
        std::string id = Util::parseFastaHeader(std::string(seq->name.s));
        header_line.append(seq->name.s, seq->name.l);
        if(seq->comment.l)  {
            header_line.append(" ",1);
            header_line.append(seq->comment.s,seq->comment.l);
        }
        header_line.append("\n");

        // header
        //ffindex_insert_memory(data_file_hdr, index_file_hdr, &offset_header,  (char *)  header_line.c_str(), header_line.length(),  (char *) id.c_str());
        // sequence
        char sequence_buffer[100000];
        strcpy(sequence_buffer, seq->seq.s);

        std::vector<Orf::SequenceLocation*> res;
        Orf::FindOrfs(sequence_buffer, res, 42);

        for(std::vector<Orf::SequenceLocation *>::const_iterator it = res.begin(); it != res.end(); it++) {
            Orf::SequenceLocation* loc = (Orf::SequenceLocation*)*it;
            printf("%zu, %zu, %d, %d, %d \n", loc->from, loc->to, loc->strand, loc->uncertainty_from, loc->uncertainty_to);
            //printf("%s\n", sequence.substr(loc->from, loc->to).c_str());
        }

        //sequence.append("\n");
        //ffindex_insert_memory(data_file,     index_file,     &offset_sequence, (char *) sequence.c_str(),  sequence.length() , (char *) id.c_str());
        entries_num++;
        header_line.clear();
    }
    kseq_destroy(seq);
    fclose(fasta_file);
    fclose(data_file);
    fclose(data_file_hdr);

    /* Sort the index entries and write back */
    rewind(index_file);
    ffindex_index_t* index = ffindex_index_parse(index_file, entries_num);
    if(index == NULL)
    {
        perror("ffindex_index_parse failed");
        return EXIT_FAILURE;
    }
    fclose(index_file);
    ffindex_sort_index_file(index);

    index_file = fopen(index_filename, "w");
    if(index_file == NULL)
    {
        perror(index_filename);
        return EXIT_FAILURE;
    }
    err += ffindex_write(index, index_file);

    rewind(index_file_hdr);
    ffindex_index_t* index_hdr = ffindex_index_parse(index_file_hdr, entries_num);
    if(index == NULL)
    {
        perror("ffindex_index_parse failed");
        return EXIT_FAILURE;
    }
    fclose(index_file_hdr);
    ffindex_sort_index_file(index_hdr);

    index_file_hdr = fopen(index_filename_hdr, "w");
    if(index_file_hdr == NULL)
    {
        perror(index_filename_hdr);
        return EXIT_FAILURE;
    }
    err += ffindex_write(index_hdr, index_file_hdr);

    return err;
}

