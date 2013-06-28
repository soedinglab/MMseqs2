/*
 * FFindex
 * written by Andy Hauser <hauser@genzentrum.lmu.de>.
 * Please add your name here if you distribute modified versions.
 * 
 * FFindex is provided under the Create Commons license "Attribution-ShareAlike
 * 3.0", which basically captures the spirit of the Gnu Public License (GPL).
 * 
 * See:
 * http://creativecommons.org/licenses/by-sa/3.0/
 */

#define _GNU_SOURCE 1
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


#include "ffindex.h"
#include "ffutil.h"

#define MAX_FILENAME_LIST_FILES 4096


void usage(char *program_name)
{
    fprintf(stderr, "USAGE: %s -v | [-s] [-e] data_filename index_filename fasta_filename [headers_filename headers_index_filename]\n"
            "\t-s\tsort index file\n"
            "\t-e\tcreate separate indexes for the headers and the sequences\n"
            "\nDesigned and implemented by Andreas W. Hauser <hauser@genzentrum.lmu.de>.\n", program_name);
}

int main(int argn, char **argv)
{
    int sort = 0, version = 0, headers = 0;
    int opt, err = EXIT_SUCCESS;
    while ((opt = getopt(argn, argv, "sve")) != -1)
    {
        switch (opt)
        {
            case 's':
                sort = 1;
                break;
            case 'v':
                version = 1;
                break;
            case 'e':
                headers = 1;
                break;
            default:
                usage(argv[0]);
                return EXIT_FAILURE;
        }
    }

    if(version == 1)
    {
        /* Don't you dare running it on a platform where byte != 8 bits */
        printf("%s version %.2f, off_t = %zd bits\n", argv[0], FFINDEX_VERSION, sizeof(off_t) * 8);
        return EXIT_SUCCESS;
    }

    if(argn - optind < 3)
    {
        usage(argv[0]);
        return EXIT_FAILURE;
    }


    char *data_filename  = argv[optind++];
    char *index_filename = argv[optind++];
    char *fasta_filename = argv[optind++];
    char *data_filename_hdr;
    char *index_filename_hdr;
    if(headers == 1){
        data_filename_hdr = argv[optind++];
        index_filename_hdr = argv[optind++];
    }
    FILE *data_file, *index_file, *fasta_file, *data_file_hdr, *index_file_hdr;

    struct stat st;

    if(stat(data_filename, &st) == 0) { errno = EEXIST; perror(data_filename); return EXIT_FAILURE; }
    data_file  = fopen(data_filename, "w");
    if( data_file == NULL) { perror(data_filename); return EXIT_FAILURE; }

    if(stat(index_filename, &st) == 0) { errno = EEXIST; perror(index_filename); return EXIT_FAILURE; }
    index_file = fopen(index_filename, "w+");
    if(index_file == NULL) { perror(index_filename); return EXIT_FAILURE; }

    fasta_file = fopen(fasta_filename, "r");
    if(fasta_file == NULL) { perror(fasta_filename); return EXIT_FAILURE; }

    if (headers == 1){
        if(stat(data_filename_hdr, &st) == 0) { errno = EEXIST; perror(data_filename_hdr); return EXIT_FAILURE; }
        data_file_hdr  = fopen(data_filename_hdr, "w");
        if( data_file_hdr == NULL) { perror(data_filename_hdr); return EXIT_FAILURE; }

        if(stat(index_filename_hdr, &st) == 0) { errno = EEXIST; perror(index_filename_hdr); return EXIT_FAILURE; }
        index_file_hdr = fopen(index_filename_hdr, "w+");
        if(index_file_hdr == NULL) { perror(index_filename_hdr); return EXIT_FAILURE; }
    }

    size_t fasta_size;
    char *fasta_data = ffindex_mmap_data(fasta_file, &fasta_size);
    size_t offset = 0;
    size_t hdr_offset = 0;
    size_t seq_offset = 0;
    size_t from_length = 0;
    char name[FFINDEX_MAX_ENTRY_NAME_LENTH];
    for(size_t fasta_offset = 1; fasta_offset < fasta_size; fasta_offset++) // position after first ">"
    {
        /* entry name is the UniProt ID or other ID until a blank space occurs*/
        size_t pos = 0;
        while (*(fasta_data + fasta_offset + pos) !=  ' '){
            name[pos] = *(fasta_data + fasta_offset + pos);
            pos++;
        }
        name[pos] = '\0';

        if (headers == 0){
            from_length = 1;
            while(fasta_offset < fasta_size && !(*(fasta_data + fasta_offset) == '>' && *(fasta_data + fasta_offset - 1) == '\n'))
            {
                fasta_offset++;
                from_length++;
            }
            ffindex_insert_memory(data_file, index_file, &offset, fasta_data + (fasta_offset - from_length), from_length, name);
        }
        else{
            // header
            from_length = 1;
            while(fasta_offset < fasta_size && !(*(fasta_data + fasta_offset - 1) == '\n'))
            {
                fasta_offset++;
                from_length++;
            }
            ffindex_insert_memory(data_file_hdr, index_file_hdr, &hdr_offset, fasta_data + (fasta_offset - from_length), from_length, name);
            
            // sequence
            from_length = 0;
            while(fasta_offset < fasta_size && !(*(fasta_data + fasta_offset) == '>' && *(fasta_data + fasta_offset - 1) == '\n'))
            {
                fasta_offset++;
                from_length++;
            }
            ffindex_insert_memory(data_file, index_file, &seq_offset, fasta_data + (fasta_offset - from_length), from_length, name);
        }
    }
    fclose(data_file);

    /* Sort the index entries and write back */
    if(sort)
    {
        rewind(index_file);
        ffindex_index_t* index = ffindex_index_parse(index_file, 0);
        if(index == NULL)
        {
            perror("ffindex_index_parse failed");
            exit(EXIT_FAILURE);
        }
        fclose(index_file);
        ffindex_sort_index_file(index);
        index_file = fopen(index_filename, "w");
        if(index_file == NULL) { perror(index_filename); return EXIT_FAILURE; }
        err += ffindex_write(index, index_file);

        if (headers == 1){
            rewind(index_file_hdr);
            ffindex_index_t* index_hdr = ffindex_index_parse(index_file_hdr, 0);
            if(index == NULL)
            { 
                perror("ffindex_index_parse failed");
                exit(EXIT_FAILURE);
            }
            fclose(index_file_hdr);
            ffindex_sort_index_file(index_hdr);
            index_file_hdr = fopen(index_filename_hdr, "w");
            if(index_file_hdr == NULL) { perror(index_filename_hdr); return EXIT_FAILURE; }
            err += ffindex_write(index_hdr, index_file_hdr);
        }
    }


    return err;
}

