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
    fprintf(stderr, "USAGE: %s DATA_FILENAME INDEX_FILENAME PROGRAM [PROGRAM_ARGS]*\n"
                    "\nDesigned and implemented by Andy Hauser <hauser@genzentrum.lmu.de>.\n",
                    argv[0]);
    return -1;
  }
  char *data_filename  = argv[1];
  char *index_filename = argv[2];
  char *program_name   = argv[3];
  char **program_argv = argv + 3;

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
  
  // Ignore SIGPIPE
  struct sigaction handler;
  handler.sa_handler = SIG_IGN;
  sigemptyset(&handler.sa_mask);
  handler.sa_flags = 0;
  sigaction(SIGPIPE, &handler, NULL);

  size_t range_start = 0;
  size_t range_end = index->n_entries;

  // Foreach entry
//#pragma omp parallel for
  for(size_t entry_index = range_start; entry_index < range_end; entry_index++)
  {
    //fprintf(stderr, "index %ld\n", entry_index);
    int ret = 0;
    ffindex_entry_t* entry = ffindex_get_entry_by_index(index, entry_index);
    if(entry == NULL) { perror(entry->name); continue; }

    int pipefd[2];
    ret = pipe(pipefd);
    if(ret != 0) { perror(entry->name); continue; }

    pid_t child_pid = fork();
    if(child_pid == 0)
    {
      fclose(data_file);
      fclose(index_file);
      close(pipefd[1]);

      // Make pipe from parent our new stdin
      int newfd = dup2(pipefd[0], fileno(stdin));
      if(newfd < 0) { fprintf(stdout, "%d %d\n", pipefd[0], newfd); perror(entry->name); }
      close(pipefd[0]);

      // exec program with the pipe as stdin
      execvp(program_name, program_argv);
      // never reached
    }
    else if(child_pid > 0)
    {
      // Read end is for child only
      close(pipefd[0]);

      // Write file data to child's stdin.
      char *filedata = ffindex_get_data_by_entry(data, entry);
      ssize_t written = 0;
      while(written < entry->length)
      {
        int w = write(pipefd[1], filedata + written, entry->length - written);
        if(w < 0 && errno != EPIPE)   { perror(entry->name); break; }
        else if(w == 0 && errno != 0) { perror(entry->name); break; }
        else
          written += w;
      }

      close(pipefd[1]); // child gets EOF
      waitpid(child_pid, NULL, 0);
    }
    else
    {
      perror(entry->name);
      exit(errno);
    }
  }

  return 0;
}

/* vim: ts=2 sw=2 et
 */
