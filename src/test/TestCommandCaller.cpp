#include <iostream>
#include <list>
#include "CommandCaller.h"


int main(int argc, const char **argv)
{
/*    CommandCaller caller;
    caller.addVariable("VAR1", "BOOM BABY");
    caller.callProgram("/Users/mad/Documents/workspace/mmseqs/src/workflow/Search.sh",argv, argc);
*/

    CommandCaller caller2;
    caller2.addVariable("PREFILTER1_PAR", "--max-seqs 10");
    const char * parameter[] = {"/Users/mad/Documents/databases/db_small/db_small","/tmp/out","/tmp/mmseqs"};
    caller2.callProgram("/Users/mad/Documents/workspace/mmseqs/src/workflow/cascaded_clustering.sh",parameter, 3);
}
