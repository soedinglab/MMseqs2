#include "Concat.h"
#include <vector>

const char* binary_name = "test_merge";

int main(int argc, char **argv) {
    std::vector<FILE*> files;
    for(size_t i = 1; i <= argc; i++){
        files.emplace_back(fopen(argv[i], "rb"));
    }
    
    FILE* outfile = fopen("/tmp/test_out","wb");
    Concat::concatFiles(files, outfile);
    fclose(outfile);

    for(size_t i = 0; i < files.size(); i++){
        fclose(files[i]);
    }
}
