#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Parameters.h"
#include "Util.h"

#include <sys/time.h>


int mergedbs(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2, true, Parameters::PARSE_VARIADIC);

    struct timeval start, end;
    gettimeofday(&start, NULL);

    if (par.filenames.size() <= 2) {
        Debug(Debug::ERROR) << "Not enough databases for merging passed!\n";
        EXIT(EXIT_FAILURE);
    }

    std::vector<std::pair<std::string, std::string>> filenames;
    for (size_t i = 2; i < par.filenames.size(); ++i) {
        filenames.emplace_back(par.filenames[i], par.filenames[i] + ".index");
    }

    std::vector<std::string> prefixes = Util::split(par.mergePrefixes, ",");

    DBReader<unsigned int> qdbr(par.db1.c_str(), par.db1Index.c_str(), DBReader<unsigned int>::USE_INDEX);
    qdbr.open(DBReader<unsigned int>::NOSORT);

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str());
    writer.open();
    writer.mergeFiles(qdbr, filenames, prefixes);
    writer.close();

    qdbr.close();

    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for merging: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";

    return 0;
}
