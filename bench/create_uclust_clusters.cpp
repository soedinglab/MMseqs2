#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include <sys/stat.h>

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}


void usage (const char* program_name){
    std::cout << "This program takes a UCLUST clustering results file and the corresponding DB in ffindex format and generates a ffindex database containing all clusters in FASTA format.\n";
    std::cout << "Usage: " << program_name << " -h | -c sorted_uclust_clustering_file -f ffindex_db_base -o ffindex_out_base\n";
}

int main (int argc, const char **argv){
    // parse options 
    std::string clu_filename;

    std::string db_data_filename;
    std::string db_index_filename;
    
    std::string out_data_filename;
    std::string out_index_filename;

    if (argc < 7){
        usage(argv[0]);
        exit(1);
    }
    for(int i=1; i<argc; ++i){
        std::string arg(argv[i]);
        if( arg == "-c" ){
            if (++i < argc){
                clu_filename = std::string(argv[i]);
            }
            else{
                std::cout << "No value provided for " << arg << "\n";
                return 0;
            }
            continue;
        }
        if( arg == "-f" ){
            if (++i < argc){
                db_data_filename = std::string(argv[i]);
                db_index_filename = db_data_filename + ".index";
            }
            else{
                std::cout << "No value provided for " << arg << "\n";
                return 0;
            }
            continue;
        }
        if( arg == "-o" ){
            if (++i < argc){
                out_data_filename = std::string(argv[i]);
                out_index_filename = out_data_filename + ".index";
            }
            else{
                std::cout << "No value provided for " << arg << "\n";
                return 0;
            }
            continue;
        }
        std::cout << "Unknown parameter found: " << arg << "\n";
        return 0;
    }
    
    // read the sequence database in ffindex format

    FILE *data_file  = fopen(db_data_filename.c_str(),  "r");
    FILE *index_file = fopen(db_index_filename.c_str(), "r");

    if( data_file == NULL) { fferror_print(__FILE__, __LINE__, argv[0], db_data_filename.c_str());  exit(1); }
    if(index_file == NULL) { fferror_print(__FILE__, __LINE__, argv[0], db_index_filename.c_str());  exit(1); }

    size_t data_size;
    char* data = ffindex_mmap_data(data_file, &data_size);

    ffindex_index_t* index = ffindex_index_parse(index_file, 0);
    if(index == NULL)
    {
        fferror_print(__FILE__, __LINE__, "ffindex_index_parse", (char*)db_index_filename.c_str());
        exit(EXIT_FAILURE);
    }

    // create a ffindex db where the clusters will be stored
    struct stat st;
    if(stat(out_data_filename.c_str(), &st) == 0) { errno = EEXIST; perror(out_data_filename.c_str()); exit(EXIT_FAILURE); }
    FILE* out_data_file  = fopen(out_data_filename.c_str(), "w");
    if( out_data_file == NULL) { perror(out_data_filename.c_str()); exit(EXIT_FAILURE); }

    if(stat(out_index_filename.c_str(), &st) == 0) { errno = EEXIST; perror(out_index_filename.c_str()); exit(EXIT_FAILURE); }
    FILE* out_index_file = fopen(out_index_filename.c_str(), "w+");
    if(out_index_file == NULL) { perror(out_index_filename.c_str()); exit(EXIT_FAILURE); }

    // buffer size for the file reading = 256 KB
    size_t buf_size = 262144;
    char* buf = new char[buf_size];
    size_t offset = 0;

    // cluster buffer size = 128 MB
    size_t clu_buf_size = 134217728;
    char* clu = new char[clu_buf_size];
    size_t clu_pos = 0;

    char clu_id[7];
    clu_id[6] = '\0';
    
    char seq_id[7];  
    seq_id[6] = '\0';

    int clu_num = 0;

    // read uclust clusters and store each cluster in the ffindex file
    std::ifstream in_clu(clu_filename.c_str());
    if (!in_clu.is_open()){ perror(out_data_filename.c_str()); exit(EXIT_FAILURE); }

    while(!in_clu.eof()){
        in_clu.getline(buf, buf_size);
        if (buf[0] == '#' || buf[0] == 'C')
            continue;
  
        //std::cout << buf << "\n";
       // parse the UCLUST clustering file
        char* pch = strtok(buf, "\t");
        
        size_t cnt = 1;
        int curr_clu_num = 0;
        while (pch != NULL) {
            pch = strtok(NULL, "\t");
            
            // read the current cluster number
            if (cnt == 1){
                curr_clu_num = atoi(pch);
            }
            // next cluster starts here, insert the last cluster into ffindex
            if (curr_clu_num > clu_num){
                if (clu_num % 1000 == 0)
                    std::cout << clu_num << " clusters done.\n";
                //std::cout << "Inserting cluster " << clu_id << "\n";
                ffindex_insert_memory(out_data_file, out_index_file, &offset, clu, clu_pos, clu_id);
                //std::cout << "----Cluster:----\n" << clu << "\n-------------------\n";
                clu_num = curr_clu_num;
                clu_pos = 0;
            }
            else if (curr_clu_num < clu_num){
                std::cerr << "Please sort your input file by the cluster number, e. g. with the command \nLC_COLLATE=C LANG=C sort -nk2 $this_file > $sorted_file\n";
                exit(EXIT_FAILURE);
            }
            // read the UniProt ID
            // assumed id length 6 characters! (UniProt IDs)
            if (cnt == 8){
                strncpy(seq_id, pch+3, 6);
                //std::cout << "seq_id: " << seq_id << "\n";
                if (buf[0] == 'S')
                    strncpy(clu_id, pch+3, 6);
                // append the sequence to the cluster
                ffindex_entry_t* seq_entry = ffindex_bsearch_get_entry(index, seq_id);
                strncpy(clu + clu_pos, ffindex_get_data_by_entry(data, seq_entry), seq_entry->length);
                //std::cout << "Got data: \n" << ffindex_get_data_by_entry(data, seq_entry) << "\n";
                //std::cout << "Written data at " << clu_pos << ": \n" << &clu[clu_pos] << "\n";
                // -1 removes \0 character
                clu_pos += seq_entry->length - 1;
                if (clu_pos > clu_buf_size)
                    std::cerr << "ERROR: BUFFER OVERFLOW (cluster " << clu_num << ")\n";
            }
            cnt++;
        }
    }
    in_clu.close();
    // insert the last cluster
    //std::cout << "Inserting cluster " << clu_id << "\n";
    ffindex_insert_memory(out_data_file, out_index_file, &offset, clu, clu_pos, clu_id);

    // sort the index
    std::cout << "Sorting the index\n";
    rewind(out_index_file);
    ffindex_index_t* out_index = ffindex_index_parse(out_index_file, 0);
    if(out_index == NULL)
    {  
        perror("ffindex_index_parse failed");
        exit(EXIT_FAILURE);
    }
    fclose(out_index_file);

    ffindex_sort_index_file(out_index);
    
    out_index_file = fopen(out_index_filename.c_str(), "w");
    if(out_index_file == NULL) { perror(out_index_filename.c_str()); return EXIT_FAILURE; }

    ffindex_write(out_index, out_index_file);

    fclose(out_index_file);
    delete[] buf;
    delete[] clu;
}
