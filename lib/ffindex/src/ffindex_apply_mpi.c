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
#include <unistd.h>
#include <fcntl.h>


#include <mpi.h>


#include "ffindex.h"
#include "ffutil.h"

char* read_buffer;

int ffindex_apply_by_entry(char *data, ffindex_index_t* index, ffindex_entry_t* entry, char* program_name, char** program_argv, FILE* data_file_out, FILE* index_file_out, size_t *offset)
{
  int ret = 0;
  int capture_stdout = (data_file_out != NULL);

  int pipefd_stdin[2];
  int pipefd_stdout[2];

  ret = pipe(pipefd_stdin);
  if(ret != 0) { fprintf(stderr, "ERROR in pipe stdin!\n"); perror(entry->name); return errno; }

  if(capture_stdout)
  {
    ret = pipe(pipefd_stdout);
    if(ret != 0) { fprintf(stderr, "ERROR in pipe stdout!\n"); perror(entry->name); return errno; }
  }

  // Flush so child doesn't copy and also flushes, leading to duplicate output
  fflush(data_file_out);
  fflush(index_file_out);

  pid_t child_pid = fork();

  if(child_pid == 0)
  {
    close(pipefd_stdin[1]);
    if(capture_stdout)
    {
      fclose(data_file_out);
      fclose(index_file_out);
      close(pipefd_stdout[0]);
    }

    // Make pipe from parent our new stdin
    int newfd_in = dup2(pipefd_stdin[0], fileno(stdin));
    if(newfd_in < 0) { fprintf(stderr, "ERROR in dup2 in %d %d\n", pipefd_stdin[0], newfd_in); perror(entry->name); }
    close(pipefd_stdin[0]);

    if(capture_stdout)
    {
      int newfd_out = dup2(pipefd_stdout[1], fileno(stdout));
      if(newfd_out < 0) { fprintf(stderr, "ERROR in dup2 out %d %d\n", pipefd_stdout[1], newfd_out); perror(entry->name); }
      close(pipefd_stdout[1]);
    }

    // exec program with the pipe as stdin
    execvp(program_name, program_argv);
    // never reached
  }
  else if(child_pid > 0)
  {
    // parent writes to and possible reads from child

    int flags = 0;

    // Read end is for child only
    close(pipefd_stdin[0]);

    if(capture_stdout)
      close(pipefd_stdout[1]);

    char *filedata = ffindex_get_data_by_entry(data, entry);

    if(capture_stdout)
    {
      flags = fcntl(pipefd_stdout[0], F_GETFL, 0);
      fcntl(pipefd_stdout[0], F_SETFL, flags | O_NONBLOCK);
    }

    // Write file data to child's stdin.
    ssize_t written = 0;
    size_t to_write = entry->length - 1; // Don't write ffindex trailing '\0'
    char* b = read_buffer;
    while(written < to_write)
    {
      size_t rest = to_write - written;
      int batch_size = PIPE_BUF;
      if(rest < PIPE_BUF)
        batch_size = rest;

      ssize_t w = write(pipefd_stdin[1], filedata + written, batch_size);
      if(w < 0 && errno != EPIPE)
      { fprintf(stderr, "ERROR in child!\n"); perror(entry->name); break; }
      else
        written += w;

      if(capture_stdout)
      {
        // To avoid blocking try to read some data
        ssize_t r = read(pipefd_stdout[0], b, PIPE_BUF);
        if(r > 0)
          b += r;
      }
    }
    close(pipefd_stdin[1]); // child gets EOF

    if(capture_stdout)
    {
      // Read rest
      fcntl(pipefd_stdout[0], F_SETFL, flags); // Remove O_NONBLOCK
      ssize_t r;
      while((r = read(pipefd_stdout[0], b, PIPE_BUF)) > 0)
        b += r;
      close(pipefd_stdout[0]);

      ffindex_insert_memory(data_file_out, index_file_out, offset, read_buffer, b - read_buffer, entry->name);
    }

    waitpid(child_pid, NULL, 0);
  }
  else
  {
    fprintf(stderr, "ERROR in fork()\n");
    perror(entry->name);
    return errno;
  }
  return EXIT_SUCCESS;
}

int main(int argn, char **argv)
{
  int mpi_error,
      mpi_rank,
      mpi_num_procs;

  mpi_error = MPI_Init(&argn, &argv);
  mpi_error = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  mpi_error = MPI_Comm_size(MPI_COMM_WORLD, &mpi_num_procs);

  int opt;
  char *data_filename_out  = NULL,
       *index_filename_out = NULL;

  while ((opt = getopt(argn, argv, "d:i:")) != -1)
  {
    switch (opt)
    {
      case 'd':
        data_filename_out = optarg;
        break;
      case 'i':
        index_filename_out = optarg;
        break;
    }
  }

  if(argn - optind < 3)
  {
    fprintf(stderr, "Not enough arguments %d.\n", optind - argn);
    fprintf(stderr, "USAGE: %s -d DATA_FILENAME_OUT -i INDEX_FILENAME_OUT DATA_FILENAME INDEX_FILENAME -- PROGRAM [PROGRAM_ARGS]*\n"
                    "\nDesigned and implemented by Andy Hauser <hauser@genzentrum.lmu.de>.\n",
                    argv[0]);
    return -1;
  }
  read_buffer = malloc(400 * 1024 * 1024);
  char *data_filename  = argv[optind++];
  char *index_filename = argv[optind++];
  char *program_name   = argv[optind];
  char **program_argv = argv + optind;

  FILE *data_file  = fopen(data_filename,  "r");
  FILE *index_file = fopen(index_filename, "r");

  if( data_file == NULL) { fferror_print(__FILE__, __LINE__, argv[0], data_filename);  exit(EXIT_FAILURE); }
  if(index_file == NULL) { fferror_print(__FILE__, __LINE__, argv[0], index_filename);  exit(EXIT_FAILURE); }

  FILE *data_file_out = NULL, *index_file_out = NULL;
  // Setup one output FFindex for each MPI process
  if(data_filename_out != NULL && index_filename_out != NULL)
  {
    char* data_filename_out_rank  = malloc(FILENAME_MAX);
    char* index_filename_out_rank = malloc(FILENAME_MAX);
    snprintf( data_filename_out_rank, FILENAME_MAX, "%s.%d", data_filename_out,  mpi_rank);
    snprintf(index_filename_out_rank, FILENAME_MAX, "%s.%d", index_filename_out, mpi_rank);
    data_file_out  = fopen(data_filename_out_rank,  "w+");
    index_file_out = fopen(index_filename_out_rank, "w+");

    if( data_file_out == NULL) { fferror_print(__FILE__, __LINE__, argv[0], data_filename_out);  exit(EXIT_FAILURE); }
    if(index_file_out == NULL) { fferror_print(__FILE__, __LINE__, argv[0], index_filename_out);  exit(EXIT_FAILURE); }
  }

  int capture_stdout = (data_file_out != NULL);

  size_t data_size;
  char *data = ffindex_mmap_data(data_file, &data_size);

  ffindex_index_t* index = ffindex_index_parse(index_file, 0);
  if(index == NULL)
  {
    fferror_print(__FILE__, __LINE__, "ffindex_index_parse", index_filename);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  
  // Ignore SIGPIPE
  struct sigaction handler;
  handler.sa_handler = SIG_IGN;
  sigemptyset(&handler.sa_mask);
  handler.sa_flags = 0;
  sigaction(SIGPIPE, &handler, NULL);

  size_t batch_size, range_start, range_end;

  if(index->n_entries >= mpi_num_procs)
    batch_size = index->n_entries / mpi_num_procs;
  else
    batch_size = 0;
  range_start = mpi_rank * batch_size;
  range_end = range_start + batch_size;


  size_t offset = 0;
  // Foreach entry
  if(batch_size > 0)
    for(size_t entry_index = range_start; entry_index < range_end; entry_index++)
    {
      ffindex_entry_t* entry = ffindex_get_entry_by_index(index, entry_index);
      if(entry == NULL) { perror(entry->name); return errno; }
      int error = ffindex_apply_by_entry(data, index, entry, program_name, program_argv, data_file_out, index_file_out, &offset);
      if(error != 0)
        { perror(entry->name); break; }
    }
  ssize_t left_over = index->n_entries - (batch_size * mpi_num_procs);
  if(mpi_rank < left_over)
  {
    size_t left_over_entry_index = (batch_size * mpi_num_procs) + mpi_rank;
    ffindex_entry_t* entry = ffindex_get_entry_by_index(index, left_over_entry_index);
    if(entry == NULL) { perror(entry->name); return errno; }
    //fprintf(stderr, "handling left over: %ld\n", left_over_entry_index);
    int error = ffindex_apply_by_entry(data, index, entry, program_name, program_argv, data_file_out, index_file_out, &offset);
    if(error != 0)
      perror(entry->name);
  }

  if(capture_stdout)
    fclose(data_file_out);
  if(index_file_out != NULL)
    fclose(index_file_out);

  MPI_Barrier(MPI_COMM_WORLD);


  // merge FFindexes in master
  if(data_filename_out != NULL && mpi_rank == 0)
  {
    char* merge_command  = malloc(FILENAME_MAX * 5);
    for(int i = 0; i < mpi_num_procs; i++)
    {
      snprintf( merge_command, FILENAME_MAX, "ffindex_build -as %s %s -d %s.%d -i %s.%d",
                data_filename_out, index_filename_out, data_filename_out, i, index_filename_out, i);
      //puts(merge_command);
      system(merge_command);
    }
  }

  MPI_Finalize();

  return EXIT_SUCCESS;
}

/* vim: ts=2 sw=2 et
 */
