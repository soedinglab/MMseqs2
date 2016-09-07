#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Parameters.h"
#include "Util.h"

#include <sys/time.h>


int mergeffindex(int argc, const char **argv) {
    std::string usage("\nMerge multiple ffindex files based on similar id into one file. \n");
    usage.append("Written by Martin Steinegger (martin.steinegger@mpibpc.mpg.de).\n\n");
    usage.append("USAGE: <QueryDB> <outDB> <ffindexDB1> ... <ffindexDBn>\n");
    Parameters par;
    par.parseParameters(argc, argv, usage, par.mergeffindex, 2, true, true);

    struct timeval start, end;
    gettimeofday(&start, NULL);
    std::vector<std::pair<std::string, std::string>> filenames;
    for (int i = 2; i < argc; i++) {
        std::string arg = argv[i];
        if(arg[0] == '-') {
            i++;
            continue;
        }

        filenames.emplace_back(arg, arg + ".index");
    }

    std::vector<std::string> prefixes = Util::split(par.mergePrefixes, ",");

    DBReader<unsigned int> qdbr(par.db1.c_str(), std::string(par.db1 + ".index").c_str());
    qdbr.open(DBReader<unsigned int>::NOSORT);

    DBWriter writer(par.db2.c_str(), std::string(par.db2 + ".index").c_str());
    writer.open();
    writer.mergeFiles(qdbr, filenames, prefixes);
    writer.close();

    qdbr.close();

    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for merging: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";

    return 0;
}
