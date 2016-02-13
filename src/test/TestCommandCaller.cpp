#include "CommandCaller.h"
#include "Parameters.h"

#include <string>

int main(int argc, const char **argv) {
    Parameters par;
    CommandCaller caller2;
    caller2.addVariable("PREFILTER1_PAR", "--max-seqs 10");
    const char *parameter[] = {"/Users/mad/Documents/databases/db_small/db_small", "/tmp/out", "/tmp/mmseqs", NULL};
    std::string program = par.mmdir;
    program.append("/src/workflow/cascaded_clustering.sh");
    caller2.callProgram(program.c_str(), 3, parameter);
}
