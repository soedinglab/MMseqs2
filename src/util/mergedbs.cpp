#include "IndexReader.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Parameters.h"
#include "Util.h"

int mergedbs(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    if (par.filenames.size() <= 2) {
        Debug(Debug::ERROR) << "Not enough databases for merging passed!\n";
        EXIT(EXIT_FAILURE);
    }

    std::vector<std::pair<std::string, std::string>> filenames;
    for (size_t i = 2; i < par.filenames.size(); ++i) {
        filenames.emplace_back(par.filenames[i], par.filenames[i] + ".index");
    }

    std::vector<std::string> prefixes = Util::split(par.mergePrefixes, ",");
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    IndexReader qDbr(par.db1, par.threads,  IndexReader::SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, DBReader<unsigned int>::USE_INDEX);
    int dbtype = FileUtil::parseDbType(filenames[0].first.c_str());
    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), 1, par.compressed, dbtype);
    writer.open();
    writer.mergeFiles(*qDbr.sequenceReader, filenames, prefixes);
    writer.close();


    return EXIT_SUCCESS;
}
