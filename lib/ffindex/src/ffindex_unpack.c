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
 *
 * ffindex_apply
 * apply a program to each FFindex entry
*/

#define _GNU_SOURCE 1
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>


#include "ffindex.h"
#include "ffutil.h"


int main(int argn, char **argv)
{
  if(argn < 4)
  {
    fprintf(stderr, "USAGE: %s DATA_FILENAME INDEX_FILENAME OUT_DIR\n"
                    "\nDesigned and implemented by Andy Hauser <hauser@genzentrum.lmu.de>.\n",
                    argv[0]);
    return -1;
  }
  char *data_filename  = argv[1];
  char *index_filename = argv[2];
  char *out_dir = argv[3];

  FILE *data_file  = fopen(data_filename,  "r");
  FILE *index_file = fopen(index_filename, "r");

  if( data_file == NULL) { fferror_print(__FILE__, __LINE__, argv[0], data_filename);  exit(EXIT_FAILURE); }
  if(index_file == NULL) { fferror_print(__FILE__, __LINE__, argv[0], index_filename);  exit(EXIT_FAILURE); }

  size_t data_size;
  char *data = ffindex_mmap_data(data_file, &data_size);

  ffindex_index_t* index = ffindex_index_parse(index_file, 0);
  if(index == NULL)
  {
    fferror_print(__FILE__, __LINE__, "ffindex_index_parse", index_filename);
    exit(EXIT_FAILURE);
  }

  if(chdir(out_dir) < 0){ fferror_print(__FILE__, __LINE__, argv[0], out_dir);  exit(EXIT_FAILURE); }

  size_t range_start = 0;
  size_t range_end = index->n_entries;

  // Foreach entry
  //#pragma omp parallel for
  for(size_t entry_index = range_start; entry_index < range_end; entry_index++)
  {
    //fprintf(stderr, "index %ld\n", entry_index);

    ffindex_entry_t* entry = ffindex_get_entry_by_index(index, entry_index);
    if(entry == NULL) { perror(entry->name); continue; }

    FILE *output_file = fopen(entry->name, "w");

    // Write file data to child's stdin.
    char *filedata = ffindex_get_data_by_entry(data, entry);
    size_t written = fwrite(filedata, entry->length - 1, 1, output_file);
    if(written < 1)   { perror(entry->name); break; }

    fclose(output_file);
  }

  return 0;
}

/* vim: ts=2 sw=2 et
 */
