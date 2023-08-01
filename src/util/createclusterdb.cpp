#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "FastSort.h"
#include "Parameters.h"

#ifdef OPENMP
#include <omp.h>
#endif

int createclusearchdb(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);
    DBReader<unsigned int> clusterReader(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                         DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    clusterReader.open(DBReader<unsigned int>::NOSORT);
    std::vector<std::string> suffixes = Util::split(par.dbSuffixList, ",");
    suffixes.insert(suffixes.begin(), "");
    for(size_t prefix = 0; prefix < suffixes.size(); prefix++) {
        std::string db1 = par.db1 + suffixes[prefix];
        std::string db1Index = par.db1 + suffixes[prefix] + ".index";
        DBReader<unsigned int> reader(db1.c_str(), db1Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        reader.open(DBReader<unsigned int>::NOSORT);
        reader.readMmapedDataInMemory();

        std::string repDbSeq = par.db3 + suffixes[prefix];
        std::string repDbSeqIdx = par.db3 + suffixes[prefix] + ".index";

        DBWriter dbwRep(repDbSeq.c_str(), repDbSeqIdx.c_str(), static_cast<unsigned int>(par.threads), par.compressed,
                        reader.getDbtype());
        dbwRep.open();
        std::string seqsDbSeq = par.db3 + "_seq" + suffixes[prefix];
        std::string seqsDbSeqIdx = par.db3 + "_seq" + suffixes[prefix] + ".index";
        DBWriter dbwClu(seqsDbSeq.c_str(), seqsDbSeqIdx.c_str(), static_cast<unsigned int>(par.threads), par.compressed,
                        reader.getDbtype());
        dbwClu.open();
        Debug::Progress progress(clusterReader.getSize());
    #pragma omp parallel
        {
            unsigned int thread_idx = 0;
    #ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
    #endif
    #pragma omp for schedule(dynamic, 1)
            for (size_t id = 0; id < clusterReader.getSize(); id++) {
                progress.updateProgress();
                char *data = clusterReader.getData(id, thread_idx);
                size_t repKey = clusterReader.getDbKey(id);
                size_t repDataId = reader.getId(repKey);
                size_t repEntryLen = reader.getEntryLen(repDataId);
                dbwRep.writeData(reader.getData(repDataId, thread_idx), repEntryLen - 1, repKey, thread_idx);
                while (*data != '\0') {
                    // parse dbkey
                    size_t dbKey = Util::fast_atoi<unsigned int>(data);
                    if (dbKey == repKey) {
                        data = Util::skipLine(data);
                        continue;
                    }
                    size_t readerId = reader.getId(dbKey);
                    dbwClu.writeData(reader.getData(readerId, thread_idx),
                                     reader.getEntryLen(readerId) - 1, dbKey, thread_idx);
                    data = Util::skipLine(data);
                }
            }
        }
        dbwRep.close(true);
        dbwClu.close(true);
        reader.close();

        // merge index
        DBReader<unsigned int> dbrRep(repDbSeq.c_str(), repDbSeqIdx.c_str(), par.threads,
                                      DBReader<unsigned int>::USE_INDEX);
        dbrRep.open(DBReader<unsigned int>::NOSORT);
        DBReader<unsigned int> dbrSeq(seqsDbSeq.c_str(), seqsDbSeqIdx.c_str(), par.threads,
                                      DBReader<unsigned int>::USE_INDEX);
        dbrSeq.open(DBReader<unsigned int>::NOSORT);
        std::string seqsDbSeqIdxTmp = seqsDbSeqIdx + "_tmp";

        FILE *sIndex = FileUtil::openAndDelete(seqsDbSeqIdxTmp.c_str(), "w");
        std::vector<DBReader<unsigned int>::Index> allIndex(dbrSeq.getSize() + dbrRep.getSize());
        size_t dataSize = 0;
        for (size_t i = 0; i < dbrRep.getSize(); i++) {
            allIndex[i] = *dbrRep.getIndex(i);
            dataSize += allIndex[i].length;
        }
        for (size_t i = 0; i < dbrSeq.getSize(); i++) {
            DBReader<unsigned int>::Index *index = dbrSeq.getIndex(i);
            index->offset += dataSize;
            allIndex[dbrRep.getSize() + i] = *index;
        }
        SORT_PARALLEL(allIndex.begin(), allIndex.end(), DBReader<unsigned int>::Index::compareById);
        char buffer[1024];
        for (size_t i = 0; i < allIndex.size(); i++) {
            size_t len = DBWriter::indexToBuffer(buffer, allIndex[i].id, allIndex[i].offset, allIndex[i].length);
            size_t written = fwrite(buffer, sizeof(char), len, sIndex);
            if (written != len) {
                Debug(Debug::ERROR) << "Cannot write index file " << seqsDbSeqIdxTmp << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        if (fclose(sIndex) != 0) {
            Debug(Debug::ERROR) << "Cannot close index file " << seqsDbSeqIdxTmp << "\n";
            EXIT(EXIT_FAILURE);
        }
        FileUtil::move(seqsDbSeqIdxTmp.c_str(), seqsDbSeqIdx.c_str());
        FileUtil::symlinkAlias(repDbSeq, seqsDbSeq + ".0");
        FileUtil::move(seqsDbSeq.c_str(), (seqsDbSeq + ".1").c_str());
        dbrRep.close();
        dbrSeq.close();
    }
    clusterReader.close();
    DBReader<unsigned int>::copyDb(par.db2, par.db3 + "_clu");

    struct DBSuffix {
        DBFiles::Files flag;
        const char *suffix;
    };

    const DBSuffix suffices[] = {
            {DBFiles::LOOKUP,        ".lookup"},
            {DBFiles::SOURCE,        ".source"},
            {DBFiles::TAX_MAPPING,   "_mapping"},
            {DBFiles::TAX_NAMES,     "_names.dmp"},
            {DBFiles::TAX_NODES,     "_nodes.dmp"},
            {DBFiles::TAX_MERGED,    "_merged.dmp"},
            {DBFiles::TAX_MERGED,    "_taxonomy"},
    };

    for (size_t i = 0; i < ARRAY_SIZE(suffices); ++i) {
        std::string file = par.db1 + suffices[i].suffix;
        if (suffices[i].flag && FileUtil::fileExists(file.c_str())) {
            DBReader<unsigned int>::copyDb(file, par.db3 + suffices[i].suffix);
        }
    }
    for (size_t i = 0; i < ARRAY_SIZE(suffices); ++i) {
        std::string file = par.db3 + suffices[i].suffix;
        if (suffices[i].flag && FileUtil::fileExists(file.c_str())) {
            std::string fileToLinkTo = par.db3 + "_seq" + suffices[i].suffix;
            if (FileUtil::fileExists(fileToLinkTo.c_str())){
                DBReader<unsigned int>::removeDb(fileToLinkTo);
            }
            DBReader<unsigned int>::aliasDb(file, fileToLinkTo);
        }
    }
    return EXIT_SUCCESS;
}
