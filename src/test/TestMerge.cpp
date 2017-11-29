#include <iostream>
#include <list>
#include <algorithm>
#include <math.h>
#include <IndexTable.h>
#include <Concat.h>


#include "Clustering.h"
#include "SetElement.h"

#include "DBReader.h"
#include "DBWriter.h"

const char* binary_name = "test_merge";

int main(int argc, char **argv)
{
    size_t fileCnt = std::max(2, argc - 1);
    FILE ** infiles = new FILE*[fileCnt];

    if(argc == 1){
        infiles[0]=fopen("/tmp/test1", "rb");
        infiles[1]=fopen("/tmp/test2", "rb");
    }
    for(size_t i = 0; i < fileCnt; i++){
        infiles[i] = fopen(argv[i + 1], "rb");
    }
    FILE * outfile = fopen("/tmp/test_out","wb");
    Concat::concatFiles(infiles, 2, outfile);
    for(size_t i = 0; i < fileCnt; i++){
        fclose(infiles[i]);
    }
    fclose(outfile);
}
