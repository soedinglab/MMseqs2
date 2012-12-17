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

#include <sys/types.h>
#include <sys/stat.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <unistd.h>


#include "ffindex.h"
#include "ffutil.h"

#define MAX_FILENAME_LIST_FILES 4096

void usage(char *program_name)
{
    fprintf(stderr, "USAGE: %s [-s|-u|-v] [-t] [-f file]* index_filename [filename]*\n"
                    "\t-f file\tfile each line containing a filename\n"
                    "\t\t-f can be specified up to %d times\n"
                    "\t-s\tsort index file\n"
                    "\t-u\tunlink entry (remove from index only)\n"
                    "\t-v\tprint version and other info then exit\n"
                    "\nDesigned and implemented by Andreas W. Hauser <hauser@genzentrum.lmu.de>.\n",
                    program_name, MAX_FILENAME_LIST_FILES);
}

int main(int argn, char **argv)
{
  int sort = 0, unlink = 0, version = 0, use_tree = 1;
  int opt, err = EXIT_SUCCESS;
  char* list_filenames[MAX_FILENAME_LIST_FILES];
  size_t list_filenames_index = 0;
  while ((opt = getopt(argn, argv, "stuvf:")) != -1)
  {
    switch (opt)
    {
      case 'f':
        list_filenames[list_filenames_index++] = optarg;
        break;
      case 's':
        sort = 1;
        break;
      case 't':
        use_tree = 1;
        break;
      case 'u':
        unlink = 1;
        break;
      case 'v':
        version = 1;
        break;
      default:
        fprintf(stderr, "Option %c not recognized\n", opt);
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

  if(optind >= argn)
  {
    usage(argv[0]);
    return EXIT_FAILURE;
  }

  char *index_filename = argv[optind++];
  FILE *index_file;

  index_file = fopen(index_filename, "r+");
  if(index_file == NULL) { perror(index_filename); return EXIT_FAILURE; }

  ffindex_index_t* index = ffindex_index_parse(index_file, 0);
  if(index == NULL) { perror("ffindex_index_parse failed"); return (EXIT_FAILURE); }

  fclose(index_file);

  /* Unlink entries */
  if(unlink)
  {
    if(use_tree)
    {
      /* Build tree */
      index = ffindex_index_as_tree(index);

      /* For each list_file unlink all entries */
      if(list_filenames_index > 0)
        for(int i = 0; i < list_filenames_index; i++)
        {
          printf("Unlinking entries from '%s'\n", list_filenames[i]);
          FILE *list_file = fopen(list_filenames[i], "r");
          if( list_file == NULL) { perror(list_filenames[i]); return EXIT_FAILURE; }

          /* unlink entries in file, one per line */
          char path[PATH_MAX];
          while(fgets(path, PATH_MAX, list_file) != NULL)
            index = ffindex_unlink(index, ffnchomp(path, strlen(path)));
        }

      /* unlink entries specified by args */
      for(int i = optind; i < argn; i++)
        index = ffindex_unlink(index, argv[i]);
    }
    else
    {
      char** sorted_names_to_unlink = malloc(FFINDEX_MAX_INDEX_ENTRIES_DEFAULT * sizeof(char *));
      if(sorted_names_to_unlink == NULL)
        fferror_print(__FILE__, __LINE__, __func__, "malloc failed");
      /* For each list_file unlink all entries */
      if(list_filenames_index > 0)
        for(int i = 0; i < list_filenames_index; i++)
        {
          printf("Unlinking entries from '%s'\n", list_filenames[i]);
          FILE *list_file = fopen(list_filenames[i], "r");
          if( list_file == NULL) { perror(list_filenames[i]); return EXIT_FAILURE; }

          /* unlink entries in file, one per line */
          char path[PATH_MAX];
          while(fgets(path, PATH_MAX, list_file) != NULL)
            sorted_names_to_unlink[i++] = ffnchomp(strdup(path), strlen(path));
          ffindex_unlink_entries(index, sorted_names_to_unlink, i);
        }

      /* unlink entries specified by args */
      int y = 0;
      for(int i = optind; i < argn; i++, y++)
        sorted_names_to_unlink[y] = argv[i];
      ffindex_unlink_entries(index, sorted_names_to_unlink, y);

      /* Sort the index entries and write back */
      if(sort)
      {
        ffindex_sort_index_file(index);
        index_file = fopen(index_filename, "w");
        if(index_file == NULL) { perror(index_filename); return EXIT_FAILURE; }
        err += ffindex_write(index, index_file);
      }
    }
  }

  /* Write index back */
  index_file = fopen(index_filename, "w");
  if(index_file == NULL) { perror(index_filename); return EXIT_FAILURE; }
  err += ffindex_write(index, index_file);
  return err;
}

/* vim: ts=2 sw=2 et
*/
