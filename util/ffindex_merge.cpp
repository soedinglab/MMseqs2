/* 
 * Written by Maria Hauser <mhauser@genzentrum.lmu.de>
 *
 */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <limits.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

#include <list>

typedef struct {
    char name[FFINDEX_MAX_ENTRY_NAME_LENTH];
    float pref_score;
    char* rest;
} pref_entry_t;


void usage (char *program_name)
{
    fprintf(stderr, "USAGE: %s data_filename_1 index_filename_1 data_filename_2 index_filename_2 data_filename_out index_filename_out\n"
            "\nWritten by Maria Hauser <mhauser@genzentrum.lmu.de>.\n", program_name);
}

std::list<pref_entry_t*>* parse_pref_list(char* pref_list){
    char name[FFINDEX_MAX_ENTRY_NAME_LENTH];
    char pref_score_str[1000];
    char* rest;



}

void store_data (char* data1, ffindex_index_t* index1, char* data2, ffindex_index_t* index2, std::list<pref_entry_t*>** db_img){
    size_t offset = 0;
    for (size_t i = 0; i < index1->n_entries; i++){
        ffindex_entry_t* e = ffindex_get_entry_by_index(index1, i);
        char* entry_data = ffindex_get_data_by_entry(data1, e);
        std::list<pref_entry_t*>* l = parse_pref_list(entry_data1);
    }

}


int main(int argn, char **argv)
{
    if (argn != 7){
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    char* data_fname_1 = argv[1];
    char* index_fname_1 = argv[2];
    char* data_fname_2 = argv[3];
    char* index_fname_2 = argv[4];
    char* data_fname_out = argv[5];
    char* index_fname_out = argv[6];

    FILE *data_file_1, *index_file_1, *data_file_2, *index_file_2, *data_file_out, *index_file_out;

    char *data1, *data2;

    ffindex_index_t *index1, *index2;

    size_t data_size_1, data_size_2;

    struct stat st;

    // open the first data file
    if(stat(data_fname_1, &st) == 0) { errno = EEXIST; perror(data_fname_1); return EXIT_FAILURE; }
    data_file_1  = fopen(data_fname_1, "r");
    if( data_file_1 == NULL) { perror(data_fname_1); return EXIT_FAILURE; }
    data1 = ffindex_mmap_data(data_file_1, &data_size_1);

    // open the second data file
    if(stat(data_fname_2, &st) == 0) { errno = EEXIST; perror(data_fname_2); return EXIT_FAILURE; }
    data_file_2  = fopen(data_fname_2, "r");
    if( data_file_2 == NULL) { perror(data_fname_2); return EXIT_FAILURE; }
    data2 = ffindex_mmap_data(data_file_2, &data_size_2);

    // open the first index
    if(stat(index_fname_1, &st) == 0) { errno = EEXIST; perror(index_fname_1); return EXIT_FAILURE; }
    index_file_1  = fopen(index_fname_1, "r");
    if( index_file_1 == NULL) { perror(index_fname_1); return EXIT_FAILURE; }
    index1 = ffindex_index_parse(index_file_1, 0);

    // open the second index
    if(stat(index_fname_2, &st) == 0) { errno = EEXIST; perror(index_fname_2); return EXIT_FAILURE; }
    index_file_2  = fopen(index_fname_2, "r");
    if( index_file_2 == NULL) { perror(index_fname_2); return EXIT_FAILURE; }
    index2 = ffindex_index_parse(index_file_2, 0);

    // create the structures to store the intermediate state of the database in the memory
    std::list<pref_entry_t*>** db_img = new std::list<pref_entry_t*>*[index1->n_entries + index2->n_entries];
    for (size_t i = 0; i < (index1->n_entries + index2->n_entries); i++){
        db_img[i] = new std::list<pref_entry_t*>;
    }

    // go through index 1 and merge all datasets 
    size_t offset = 0;
    for (size_t i = 0; i < index1->n_entries; i++){
        ffindex_entry_t* e1 = ffindex_get_entry_by_index(index1, i);
        char* entry_data1 = ffindex_get_data_by_entry(data1, e1);


        
        
        char* entry_data2 = ffindex_get_data_by_name(data2, index2, e1->name);
        if (entry_data2 != NULL){
            char* entry_data_merged = merge_data(entry_data1, entry_data2);
            ffindex_insert_memory(data_file_out, index_file_out, &offset, entry_data_merged, strlen(entry_data_merged), e1->name);
        }
        else
            ffindex_insert_memory(data_file_out, index_file_out, &offset, entry_data1, strlen(entry_data1), e1->name);
    }

    // go through index 2 and put all the remaining datasets into the output index
    for (size_t i = 0; i < index2->n_entries; i++){
        ffindex_entry_t* e2 = ffindex_get_entry_by_index(index2, i);
        char* entry_data2 = ffindex_get_data_by_entry(data2, e2);
        char* entry_data1 = ffindex_get_data_by_name(data1, index1, e2->name);
        if (entry_data1 == NULL)
            ffindex_insert_memory(data_file_out, index_file_out, &offset, entry_data2, strlen(entry_data2), e2->name);
    }

    // create the output ffindex database
    if(stat(data_fname_out, &st) == 0) { errno = EEXIST; perror(data_fname_out); exit(EXIT_FAILURE); }
    if(stat(index_fname_out, &st) == 0) { errno = EEXIST; perror(index_fname_out); exit(EXIT_FAILURE); }

    data_file_out = fopen(data_fname_out, "w");
    index_file_out = fopen(index_fname_out, "w");

    if( data_file_out == NULL){ perror(data_fname_out); exit(EXIT_FAILURE); } 
    if( index_file_out == NULL) { perror(index_fname_out); exit(EXIT_FAILURE); }


}

