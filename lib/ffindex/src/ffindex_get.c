/*
 * Ffindex
 * written by Andy Hauser <hauser@genzentrum.lmu.de>.
 * Please add your name here if you distribute modified versions.
 * 
 * Ffindex is provided under the Create Commons license "Attribution-ShareAlike
 * 3.0", which basically captures the spirit of the Gnu Public License (GPL).
 * 
 * See:
 * http://creativecommons.org/licenses/by-sa/3.0/
*/

#define _GNU_SOURCE 1
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "ffindex.h"
#include "ffutil.h"

void usage(char* program_name)
{
    fprintf(stderr, "USAGE: %s data_filename index_filename entry name(s)\n"
                    "-n\tuse index of entry instead of entry name\n"
                    "\nDesigned and implemented by Andy Hauser <hauser@genzentrum.lmu.de>.\n",
                    program_name);
}

int main(int argn, char **argv)
{
  int by_index = 0;
  int opt;
  while ((opt = getopt(argn, argv, "n")) != -1)
  {
    switch (opt)
    {
      case 'n':
        by_index = 1;
        break;
      default:
        usage(argv[0]);
        return EXIT_FAILURE;
    }
  }
  if(argn < 3)
  {
    usage(argv[0]);
    return EXIT_FAILURE;
  }
  char *data_filename  = argv[optind++];
  char *index_filename = argv[optind++];

  FILE *data_file  = fopen(data_filename,  "r");
  FILE *index_file = fopen(index_filename, "r");

  if( data_file == NULL) { fferror_print(__FILE__, __LINE__, "ffindex_get", data_filename);  exit(EXIT_FAILURE); }
  if(index_file == NULL) { fferror_print(__FILE__, __LINE__, "ffindex_get", index_filename);  exit(EXIT_FAILURE); }

  size_t data_size;
  char *data = ffindex_mmap_data(data_file, &data_size);

  ffindex_index_t* index = ffindex_index_parse(index_file, 0);
  if(index == NULL)
  {
    fferror_print(__FILE__, __LINE__, "ffindex_index_parse", index_filename);
    exit(EXIT_FAILURE);
  }

  if(by_index)
  {
    for(int i = optind; i < argn; i++)
    {
      size_t index_n = atol(argv[i]) - 1; // offset from 0 but specify from 1

      ffindex_entry_t* entry = ffindex_get_entry_by_index(index, index_n);
      if(entry == NULL)
      {
        errno = ENOENT; 
        fferror_print(__FILE__, __LINE__, "ffindex_get entry index out of range", argv[i]);
      }
      else
      {
        char *filedata = ffindex_get_data_by_entry(data, entry);
        if(filedata == NULL)
        {
          errno = ENOENT; 
          fferror_print(__FILE__, __LINE__, "ffindex_get entry index out of range", argv[i]);
        }
        else
          fwrite(filedata, entry->length - 1, 1, stdout);
      }
    }
  }
  else // by name
  {
    for(int i = optind; i < argn; i++)
    {
      char *filename = argv[i];

      ffindex_entry_t* entry = ffindex_get_entry_by_name(index, filename);
      if(entry == NULL)
      {
        errno = ENOENT; 
        fferror_print(__FILE__, __LINE__, "ffindex_get key not found in index", filename);
      }
      else
      {
        char *filedata = ffindex_get_data_by_entry(data, entry);
        if(filedata == NULL)
        {
          errno = ENOENT; 
          fferror_print(__FILE__, __LINE__, "ffindex_get key not found in index", filename);
        }
        else
          fwrite(filedata, entry->length - 1, 1, stdout);
      }
    }

      /* Alternative code using (slower) ffindex_fopen */
      /*
         FILE *file = ffindex_fopen(data, index, filename);
         if(file == NULL)
         {
         errno = ENOENT; 
         fferror_print(__FILE__, __LINE__, "ffindex_fopen file not found in index", filename);
         }
         else
         {
         char line[LINE_MAX];
         while(fgets(line, LINE_MAX, file) != NULL)
         printf("%s", line);
         }
         */
  }

  return 0;
}

/* vim: ts=2 sw=2 et
 */
