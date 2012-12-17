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

#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <limits.h>
#include <libgen.h>
#include <search.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "ext/fmemopen.h" /* For OS not yet implementing this new standard function */
#include "ffutil.h"
#include "ffindex.h"

/* XXX Use page size? */
#define FFINDEX_BUFFER_SIZE 4096

char* ffindex_copyright_text = "Designed and implemented by Andy Hauser <hauser@genzentrum.lmu.de>.";

char* ffindex_copyright()
{
  return ffindex_copyright_text;
}


/* Insert a memory chunk (string even without \0) into ffindex */
int ffindex_insert_memory(FILE *data_file, FILE *index_file, size_t *offset, char *from_start, size_t from_length, char *name)
{
    int myerrno = 0;
    size_t offset_before = *offset;
    size_t write_size = fwrite(from_start, sizeof(char), from_length, data_file);
    *offset += write_size;
    if(from_length != write_size)
      fferror_print(__FILE__, __LINE__, __func__, name);

    /* Seperate by '\0' and thus also make sure at least one byte is written */
    char buffer[1] = {'\0'};
    if(fwrite(buffer, sizeof(char), 1, data_file) != 1)
      perror("ffindex_insert_memory");
    *offset += 1;
    if(ferror(data_file) != 0)
      goto EXCEPTION_ffindex_insert_memory;

    /* write index entry */
    fprintf(index_file, "%s\t%zd\t%zd\n", name, offset_before, *offset - offset_before);

    return myerrno;

EXCEPTION_ffindex_insert_memory:
    {
      fferror_print(__FILE__, __LINE__, __func__, "");
      return myerrno;
    }
}


/* Insert all file from a list into ffindex */
int ffindex_insert_list_file(FILE *data_file, FILE *index_file, size_t *start_offset, FILE *list_file)
{
  size_t offset = *start_offset;
  char path[PATH_MAX];
  while(fgets(path, PATH_MAX, list_file) != NULL)
    ffindex_insert_file(data_file, index_file, &offset, ffnchomp(path, strlen(path)), basename(path));

  /* update return value */
  *start_offset = offset;
  return 0;
}


/* Insert all files from directory into ffindex */
int ffindex_insert_dir(FILE *data_file, FILE *index_file, size_t *start_offset, char *input_dir_name)
{
  DIR *dir = opendir(input_dir_name);
  if(dir == NULL)
  {
    fferror_print(__FILE__, __LINE__, __func__, input_dir_name);
    return -1;
  }

  size_t input_dir_name_len = strlen(input_dir_name);
  char path[PATH_MAX];
  strncpy(path, input_dir_name, NAME_MAX);
  if(input_dir_name[input_dir_name_len - 1] != '/')
  {
    path[input_dir_name_len] = '/';
    input_dir_name_len += 1;
  }

  size_t offset = *start_offset;
  struct dirent *entry;
  while((entry = readdir(dir)) != NULL)
  {
    if(entry->d_name[0] == '.')
      continue;
    strncpy(path + input_dir_name_len, entry->d_name, NAME_MAX);
    struct stat sb;
    if(stat(path, &sb) == -1)
      fferror_print(__FILE__, __LINE__, __func__, path);
    if(!S_ISREG(sb.st_mode))
      continue;
    ffindex_insert_file(data_file, index_file, &offset, path, entry->d_name);
  }
  closedir(dir);

  /* update return value */
  *start_offset = offset;

  return 0;
}


/* Insert one file by path into ffindex */
int ffindex_insert_file(FILE *data_file, FILE *index_file, size_t *offset, const char *path, char *name)
{
    FILE *file = fopen(path, "r");
    if(file == NULL)
      return errno;

    int ret = ffindex_insert_filestream(data_file, index_file, offset, file, name);
    fclose(file);
    return ret;
}

/* Insert one file by handle into ffindex */
int ffindex_insert_filestream(FILE *data_file, FILE *index_file, size_t *offset, FILE* file, char *name)
{
    int myerrno = 0;
    /* copy and paste file to data file */
    char buffer[FFINDEX_BUFFER_SIZE];
    size_t offset_before = *offset;
    size_t read_size;
    while((read_size = fread(buffer, sizeof(char), sizeof(buffer), file)) > 0)
    {
      size_t write_size = fwrite(buffer, sizeof(char), read_size, data_file);
      *offset += write_size;
      if(read_size != write_size)
        fferror_print(__FILE__, __LINE__, __func__, name);
    }
    if(ferror(file))
      warn("fread");

    /* Seperate by '\0' and thus also make sure at least one byte is written */
    buffer[0] = '\0';
    if(fwrite(buffer, sizeof(char), 1, data_file) != 1)
      perror("ffindex_insert_filestream");
    *offset += 1;
    if(ferror(data_file) != 0)
      goto EXCEPTION_ffindex_insert_file;

    /* write index entry */
    fprintf(index_file, "%s\t%zd\t%zd\n", name, offset_before, *offset - offset_before);

    if(ferror(file) != 0)
      goto EXCEPTION_ffindex_insert_file;

    return myerrno;

EXCEPTION_ffindex_insert_file:
    {
      fferror_print(__FILE__, __LINE__, __func__, "");
      return myerrno;
    }
}

/* XXX not implemented yet */
int ffindex_restore(FILE *data_file, FILE *index_file, char *input_dir_name)
{
  return -1;
}


char* ffindex_mmap_data(FILE *file, size_t* size)
{
  struct stat sb;
  fstat(fileno(file), &sb);
  *size = sb.st_size;
  int fd =  fileno(file);
  if(fd < 0)
  {
    fferror_print(__FILE__, __LINE__, __func__, "mmap failed");
    return MAP_FAILED;
  }
  return (char*)mmap(NULL, *size, PROT_READ, MAP_PRIVATE, fd, 0);
}


static int ffindex_compare_entries_by_name(const void *pentry1, const void *pentry2)
{   
  ffindex_entry_t* entry1 = (ffindex_entry_t*)pentry1;
  ffindex_entry_t* entry2 = (ffindex_entry_t*)pentry2;
  return strncmp(entry1->name, entry2->name, FFINDEX_MAX_ENTRY_NAME_LENTH);
}


ffindex_entry_t* ffindex_bsearch_get_entry(ffindex_index_t *index, char *name)
{
  ffindex_entry_t search;
  strncpy(search.name, name, FFINDEX_MAX_ENTRY_NAME_LENTH);
  return (ffindex_entry_t*)bsearch(&search, index->entries, index->n_entries, sizeof(ffindex_entry_t), ffindex_compare_entries_by_name);
}


ffindex_index_t* ffindex_index_parse(FILE *index_file, size_t num_max_entries)
{
  if(num_max_entries == 0)
    num_max_entries = FFINDEX_MAX_INDEX_ENTRIES_DEFAULT;
  size_t nbytes = sizeof(ffindex_index_t) + (sizeof(ffindex_entry_t) * num_max_entries);
  ffindex_index_t *index = (ffindex_index_t *)malloc(nbytes);
  if(index == NULL)
  {
    fprintf(stderr, "Failed to allocate %ld bytes\n", nbytes);
    fferror_print(__FILE__, __LINE__, __func__, "malloc failed");
    return NULL;
  }
  index->num_max_entries = num_max_entries;

  index->file = index_file;
  index->index_data = ffindex_mmap_data(index_file, &(index->index_data_size));
  if(index->index_data_size == 0)
    warn("Problem with data file. Is it empty or is another process readning it?");
  if(index->index_data == MAP_FAILED)
    return NULL;
  index->type = SORTED_ARRAY; /* XXX Assume a sorted file for now */
  int i = 0;
  char* d = index->index_data;
  char* end;
  /* Faster than scanf per line */
  for(i = 0; d < (index->index_data + index->index_data_size); i++)
  {
    int p;
    for(p = 0; *d != '\t'; d++)
      index->entries[i].name[p++] = *d;
    index->entries[i].name[p] = '\0';
    index->entries[i].offset = strtol(d, &end, 10);
    d = end;
    index->entries[i].length  = strtol(d, &end, 10);
    d = end + 1; /* +1 for newline */
  }

  index->n_entries = i;

  if(index->n_entries == 0)
    warn("index with 0 entries");

  return index;
}

ffindex_entry_t* ffindex_get_entry_by_index(ffindex_index_t *index, size_t entry_index)
{
  if(entry_index < index->n_entries)
    return &index->entries[entry_index];
  else
    return NULL;
}

/* Using a function for this looks like overhead. But a more advanced data format,
 * say a compressed one, can do it's magic here. 
 */
char* ffindex_get_data_by_offset(char* data, size_t offset)
{
  return data + offset;
}


char* ffindex_get_data_by_entry(char *data, ffindex_entry_t* entry)
{
  return ffindex_get_data_by_offset(data, entry->offset);
}


char* ffindex_get_data_by_name(char *data, ffindex_index_t *index, char *name)
{
  ffindex_entry_t* entry = ffindex_bsearch_get_entry(index, name);

  if(entry == NULL)
    return NULL;

  return ffindex_get_data_by_entry(data, entry);
}


char* ffindex_get_data_by_index(char *data, ffindex_index_t *index, size_t entry_index)
{
  ffindex_entry_t* entry = ffindex_get_entry_by_index(index, entry_index);

  if(entry == NULL)
    return NULL;

  return ffindex_get_data_by_entry(data, entry);
}


FILE* ffindex_fopen_by_entry(char *data, ffindex_entry_t* entry)
{
  char *filedata = ffindex_get_data_by_offset(data, entry->offset);
  return fmemopen(filedata, entry->length, "r");
}


FILE* ffindex_fopen_by_name(char *data, ffindex_index_t *index, char *filename)
{
  ffindex_entry_t* entry = ffindex_bsearch_get_entry(index, filename);

  if(entry == NULL)
    return NULL;

  return ffindex_fopen_by_entry(data, entry);
}


void ffindex_sort_index_file(ffindex_index_t *index)
{
  qsort(index->entries, index->n_entries, sizeof(ffindex_entry_t), ffindex_compare_entries_by_name);
}


int ffindex_write(ffindex_index_t* index, FILE* index_file)
{
  /* Use tree if available */
  if(index->type == TREE)
    return ffindex_tree_write(index, index_file);

  for(size_t i = 0; i < index->n_entries; i++)
  {
    ffindex_entry_t ffindex_entry = index->entries[i];
    if(fprintf(index_file, "%s\t%zd\t%zd\n", ffindex_entry.name, ffindex_entry.offset, ffindex_entry.length) < 0)
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}


ffindex_index_t* ffindex_unlink_entries(ffindex_index_t* index, char** sorted_names_to_unlink, int n_names)
{
  int i = index->n_entries - 1;
  /* walk list of names to delete */
  for(int n = n_names - 1; n >= 0;  n--)
  {
    char* name_to_unlink = sorted_names_to_unlink[n];
    /* walk index entries */
    for(; i >= 0; i--)
    {
      int cmp = strncmp(name_to_unlink, index->entries[i].name, FFINDEX_MAX_ENTRY_NAME_LENTH);
      if(cmp == 0) /* found entry */
      {
        /* Move entries after the unlinked ones to close the gap */
        size_t n_entries_to_move = index->n_entries - i - 1;
        if(n_entries_to_move > 0) /* not last element of array */
          memmove(index->entries + i, index->entries + i + 1, n_entries_to_move * sizeof(ffindex_entry_t));
        index->n_entries--;
        break;
      }
      else if(cmp > 0) /* not found */
        break;
    }
  }

  return index;
}


ffindex_index_t* ffindex_unlink(ffindex_index_t* index, char* name_to_unlink)
{
  /* Use tree if available */
  if(index->type == TREE)
    return ffindex_tree_unlink(index, name_to_unlink);

  ffindex_entry_t* entry = ffindex_bsearch_get_entry(index, name_to_unlink);
  if(entry == NULL)
  {
    fprintf(stderr, "Warning: could not find '%s'\n", name_to_unlink);
    return index;
  }
  /* Move entries after the unlinked one to close the gap */
  size_t n_entries_to_move = index->entries + index->n_entries - entry - 1;
  if(n_entries_to_move > 0) /* not last element of array */
    memmove(entry, entry + 1, n_entries_to_move * sizeof(ffindex_entry_t));
  index->n_entries--;
  return index;
}

/* tree version */

ffindex_entry_t *ffindex_tree_get_entry(ffindex_index_t* index, char* name)
{
  ffindex_entry_t search;
  strncpy(search.name, name, FFINDEX_MAX_ENTRY_NAME_LENTH);
  return (ffindex_entry_t *)tfind((const void *)&search, &index->tree_root, ffindex_compare_entries_by_name);
}


ffindex_index_t* ffindex_tree_unlink(ffindex_index_t* index, char* name_to_unlink)
{
  if(index->tree_root == NULL)
  {
    fferror_print(__FILE__, __LINE__, __func__, "tree is NULL");
    return NULL;
  }
  ffindex_entry_t search;
  strncpy(search.name, name_to_unlink, FFINDEX_MAX_ENTRY_NAME_LENTH);
  tdelete((const void *)&search, &index->tree_root, ffindex_compare_entries_by_name);
  return index;
}

ffindex_index_t* ffindex_index_as_tree(ffindex_index_t* index)
{
  index->tree_root = NULL;
  for(size_t i = 0; i < index->n_entries; i++)
  {
    ffindex_entry_t *entry = &index->entries[i];
    tsearch((const void *)entry, &index->tree_root, ffindex_compare_entries_by_name);
    //entry = *(ffindex_entry_t **)tsearch((const void *)entry, &index->tree_root, ffindex_compare_entries_by_name);
    //printf("entry find: %s\n", entry->name); 
  }
  index->type = TREE;
  return index;
}


int ffindex_tree_write(ffindex_index_t* index, FILE* index_file)
{
  int ret = EXIT_SUCCESS;
  void action(const void *node, const VISIT which, const int depth)
  {
    ffindex_entry_t *entry;
    switch (which)
    {
      case preorder:
        break;
      case endorder:
        break;
      case postorder:
      case leaf:
        entry = *(ffindex_entry_t **) node;
        if(fprintf(index_file, "%s\t%zd\t%zd\n", entry->name, entry->offset, entry->length) < 0)
          ret = EXIT_FAILURE;
        break;
    }                                        
  }
  twalk(index->tree_root, action);
  return ret;
}

/* vim: ts=2 sw=2 et
*/
