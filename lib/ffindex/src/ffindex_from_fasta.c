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
    fprintf(stderr, "USAGE: %s -v | [-s] data_filename index_filename fasta_filename\n"
                    "\t-s\tsort index file\n"
                    "\nDesigned and implemented by Andreas W. Hauser <hauser@genzentrum.lmu.de>.\n", program_name);
}

int main(int argn, char **argv)
{
  int sort = 0, version = 0;
  int opt, err = EXIT_SUCCESS;
  while ((opt = getopt(argn, argv, "sv")) != -1)
  {
    switch (opt)
    {
      case 's':
        sort = 1;
        break;
      case 'v':
        version = 1;
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
  FILE *data_file, *index_file, *fasta_file;

  struct stat st;

  if(stat(data_filename, &st) == 0) { errno = EEXIST; perror(data_filename); return EXIT_FAILURE; }
  data_file  = fopen(data_filename, "w");
  if( data_file == NULL) { perror(data_filename); return EXIT_FAILURE; }

  if(stat(index_filename, &st) == 0) { errno = EEXIST; perror(index_filename); return EXIT_FAILURE; }
  index_file = fopen(index_filename, "w+");
  if(index_file == NULL) { perror(index_filename); return EXIT_FAILURE; }

  fasta_file = fopen(fasta_filename, "r");
  if(fasta_file == NULL) { perror(fasta_filename); return EXIT_FAILURE; }

  size_t fasta_size;
  char *fasta_data = ffindex_mmap_data(fasta_file, &fasta_size);
  size_t offset = 0;
  size_t from_length = 0;
  char name[FFINDEX_MAX_ENTRY_NAME_LENTH];
  int seq_id = 1;
  for(size_t fasta_offset = 1; fasta_offset < fasta_size; fasta_offset++) // position after first ">"
  {
    from_length = 1;
    while(fasta_offset < fasta_size && !(*(fasta_data + fasta_offset) == '>' && *(fasta_data + fasta_offset - 1) == '\n'))
    {
      fasta_offset++;
      from_length++;
    }
    sprintf(name, "%d", seq_id++);
    ffindex_insert_memory(data_file, index_file, &offset, fasta_data + (fasta_offset - from_length), from_length, name);
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
  }

  return err;
}

/* vim: ts=2 sw=2 et: */
