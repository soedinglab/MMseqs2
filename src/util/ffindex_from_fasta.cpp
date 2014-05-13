/*
 * FFindex
 * written by Andy Hauser <hauser@genzentrum.lmu.de>.
 * modified by Maria Hauser <mhauser@genzentrum.lmu.de> (splitting into sequences/headers databases)
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
#include <string>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "CommandDeclarations.h"
#include "Util.h"



extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}


#define MAX_FILENAME_LIST_FILES 4096

// GenBank                           gi|gi-number|gb|accession|locus
// ENA Data Library                  gi|gi-number|emb|accession|locus
// DDBJ, DNA Database of Japan       gi|gi-number|dbj|accession|locus
// SWISS-PROT                        sp|accession|name
// Brookhaven Protein Data Bank (1)  pdb|entry|chain
// NBRF PIR                          pir||entry
// Protein Research Foundation       prf||name
// PDB EBI                           PDB:1ECY_A mol:protein length:142  ECOTIN
// TrEMBL                            tr|accession|name
// Patents                           pat|country|number
// GenInfo Backbone Id               bbs|number
// General database identifier       gnl|database|identifier
// NCBI Reference Sequence           ref|accession|locus
// Local Sequence identifier         lcl|identifier
// id only                           identifier
// id only                           identifier

//
//std::string parseFastaHeader(std::string header){
//    std::vector<std::string> arr = Util::split(header,"|");
//    for(int i = 0; i < arr.size(); i++)
//        std::cout << arr[i] << std::endl;
//	if(arr.size() > 1) {
//		if (Util::startWith("cl|", header)  ||
//		    Util::startWith("sp|", header)  ||
//            Util::startWith("tr|", header)  ||
//	   	    Util::startWith("ref|", header) ||
//            Util::startWith("pdb|", header)  ||
//            Util::startWith("bbs|", header)  ||
//            Util::startWith("lcl|", header)  ||
//            Util::startWith("pir||", header) ||
//            Util::startWith("prf||", header)  )
//            return arr[1];
//        else if (Util::startWith("gnl|", header) ||
//                 Util::startWith("pat|", header))
//            return arr[2];
//		else if (Util::startWith("gi|", header))
//            return arr[3];
//        
//	} else {
//		arr = Util::split(header," ");
//		return arr[0];
//	}
//}
void usage(const char *program_name)
{
    fprintf(stderr, "USAGE: %s -v | [-s] [-e] data_filename index_filename fasta_filename [headers_filename headers_index_filename]\n"
            "\t-s\tsort index file\n"
            "\t-e\tcreate separate indexes for the headers and the sequences\n"
            "\nDesigned and implemented by Andreas W. Hauser <hauser@genzentrum.lmu.de>.\n", program_name);
}

int createdb(int argn, const char **argv)
{
    int sort = 0, version = 0, headers = 0;
    int opt, err = EXIT_SUCCESS;
    while ((opt = getopt(argn,(char **) argv, "sve")) != -1)
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
    
    
    char *data_filename  = (char *)argv[optind++];
    char *index_filename = (char *)argv[optind++];
    char *fasta_filename = (char *)argv[optind++];
    char *data_filename_hdr;
    char *index_filename_hdr;
    if(headers == 1){
        data_filename_hdr = (char *)argv[optind++];
        index_filename_hdr = (char *)argv[optind++];
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
    size_t entries_num = 0;
    char name[FFINDEX_MAX_ENTRY_NAME_LENTH];
    for(size_t fasta_offset = 1; fasta_offset < fasta_size; fasta_offset++) // position after first ">"
    {
        /* entry name is the UniProt ID between pipe characters or other ID until a blank space occurs*/
        size_t local_offset = 0;
        int mode = 0;
        
        while(*(fasta_data + fasta_offset + local_offset) != '\n'){
            if(*(fasta_data + fasta_offset + local_offset) == '|'){
                mode = 1;
                break;
            }
            else if (*(fasta_data + fasta_offset + local_offset) == ' '){
                mode = 0;
                break;
            }
            local_offset++;
        }
        
        local_offset = 0;
        
        size_t pos = 0;
        if (mode == 0){
            while ((*(fasta_data + fasta_offset + local_offset) !=  ' ') && (*(fasta_data + fasta_offset + local_offset) != '\n')){
                name[pos] = *(fasta_data + fasta_offset + local_offset);
                local_offset++;
                pos++;
            }
        }
        else{
            while (*(fasta_data + fasta_offset + local_offset) != '|'){
                local_offset++;
            }
            local_offset++;
            while (*(fasta_data + fasta_offset + local_offset) !=  '|'){
                name[pos] = *(fasta_data + fasta_offset + local_offset);
                local_offset++;
                pos++;
            }
        }
        name[pos] = '\0';
        entries_num++;
        
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
            from_length-=1; // last \n should not be included
            ffindex_insert_memory(data_file_hdr, index_file_hdr, &hdr_offset, fasta_data + (fasta_offset - from_length), from_length, name);
            from_length+=1;
            // sequence
            from_length = 0;
            while(fasta_offset < fasta_size && !(*(fasta_data + fasta_offset) == '>' && *(fasta_data + fasta_offset - 1) == '\n'))
            {
                fasta_offset++;
                from_length++;
            }
            fasta_offset -= 1;
            from_length-=1; // last \n should not be included
            ffindex_insert_memory(data_file, index_file, &seq_offset, fasta_data + (fasta_offset - from_length), from_length, name);
            from_length+=1;
            fasta_offset += 1;
        }
    }
    fclose(data_file);
    
    /* Sort the index entries and write back */
    if(sort)
    {
        rewind(index_file);
        ffindex_index_t* index = ffindex_index_parse(index_file, entries_num);
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
            ffindex_index_t* index_hdr = ffindex_index_parse(index_file_hdr, entries_num);
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

