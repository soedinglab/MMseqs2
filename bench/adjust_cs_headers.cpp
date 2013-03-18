extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>

int main (int argc, char** argv){

    char* dataf = argv[1];
    char* indexf = argv[2];
    char* outf = argv[3];

    FILE* dataFile = fopen(dataf, "r");
    FILE* indexFile = fopen(indexf, "r");

    if( indexFile == NULL) { fferror_print(__FILE__, __LINE__, "run_sw", indexf);  exit(EXIT_FAILURE); }
    if( dataFile == NULL) { fferror_print(__FILE__, __LINE__, "run_sw", dataf);  exit(EXIT_FAILURE); }

    size_t data_size = 0;
    char* data = ffindex_mmap_data(dataFile, &data_size);

    ffindex_index_t* index = ffindex_index_parse(indexFile, 0);
    if(index == NULL)
    {
        fferror_print(__FILE__, __LINE__, "ffindex_index_parse", indexf);
        exit(EXIT_FAILURE);
    }

    std::ofstream out;
    out.open(outf);

    for(int i = 0; i < index->n_entries; i++){
        if (i % 100000 == 0)
            std::cout << i << "\n";
        ffindex_entry_t* e = ffindex_get_entry_by_index(index, i);
        if(e == NULL) { perror(e->name); continue; }

        char* seq = ffindex_get_data_by_entry(data, e);
        char* name = strtok (e->name,".");

        // write a new header
        out << ">" << name << "\n";

        // skip the header from the sequence and write the rest
        int write = 0;
        for (int j = 0; j < e->length-1; j++){
            if (write == 0 && seq[j] == '\n'){
                write = 1;
            }
            else if (write == 1){
                out << seq[j];
            }
        }
    }
    out.close();

}
