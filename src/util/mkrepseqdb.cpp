#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "FastSort.h"
#include "Parameters.h"

#ifdef OPENMP
#include <omp.h>
#endif

int mkrepseqdb(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::NOSORT);
    reader.readMmapedDataInMemory();

    DBReader<unsigned int> clusterReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    clusterReader.open(DBReader<unsigned int>::NOSORT);

    DBWriter dbwRep(par.db3.c_str(), par.db3Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed,  reader.getDbtype());
    dbwRep.open();
    std::string seqsDbSeq = par.db3 + "_seq";
    std::string seqsDbSeqIdx = par.db3 + "_seq.index";
    DBWriter dbwClu(seqsDbSeq.c_str(), seqsDbSeqIdx.c_str(), static_cast<unsigned int>(par.threads), par.compressed,  reader.getDbtype());
    dbwClu.open();
    Debug::Progress progress(clusterReader.getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::string resultBuffer;
        // write output file
#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < clusterReader.getSize(); id++) {
            progress.updateProgress();
            char *data = clusterReader.getData(id, thread_idx);
            size_t repKey = clusterReader.getDbKey(id);
            size_t repDataId = reader.getId(repKey);
            size_t repEntryLen = reader.getEntryLen(repDataId);
            dbwRep.writeData(reader.getData(repDataId, thread_idx), repEntryLen - 1, repKey, thread_idx);
            while(*data != '\0') {
                // parse dbkey
                size_t dbKey = Util::fast_atoi<unsigned int>(data);
                if(dbKey == repKey){
                    data = Util::skipLine(data);
                    continue;
                }
                size_t readerId =  reader.getId(dbKey);
                dbwClu.writeData(reader.getData(readerId, thread_idx),
                                 reader.getEntryLen(readerId) - 1, dbKey, thread_idx);
                data = Util::skipLine(data);
            }
            resultBuffer.clear();
        }
    }
    dbwRep.close(true);
    dbwClu.close(true);
    clusterReader.close();
    reader.close();

    // merge index
    DBReader<unsigned int> dbrRep(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX);
    dbrRep.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> dbrClu(seqsDbSeq.c_str(), seqsDbSeqIdx.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX);
    dbrClu.open(DBReader<unsigned int>::NOSORT);
    std::string seqsDbSeqIdxTmp = seqsDbSeqIdx + "_tmp";

    FILE *sIndex = FileUtil::openAndDelete(seqsDbSeqIdxTmp.c_str(), "w");
    std::vector<DBReader<unsigned int>::Index> allIndex(dbrClu.getSize() + dbrRep.getSize());
    size_t dataSize = 0;
    for(size_t i = 0; i < dbrRep.getSize(); i++) {
        allIndex[i] = *dbrRep.getIndex(i);
        dataSize += allIndex[i].length;
    }
    for(size_t i = 0; i < dbrClu.getSize(); i++) {
        DBReader<unsigned int>::Index * index = dbrClu.getIndex(i);
        index->offset += dataSize;
        allIndex[dbrRep.getSize() + i] = *index;
    }
    SORT_PARALLEL(allIndex.begin(), allIndex.end(), DBReader<unsigned int>::Index::compareById);
    char buffer[1024];
    for(size_t i = 0; i < allIndex.size(); i++){
        size_t len = DBWriter::indexToBuffer(buffer, allIndex[i].id, allIndex[i].offset, allIndex[i].length);
        size_t written = fwrite(buffer, sizeof(char), len, sIndex);
        if(written != len){
            Debug(Debug::ERROR) << "Cannot write index file " << seqsDbSeqIdxTmp << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    if (fclose(sIndex) != 0) {
        Debug(Debug::ERROR) << "Cannot close index file " << seqsDbSeqIdxTmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::move(seqsDbSeqIdxTmp.c_str(), seqsDbSeqIdx.c_str());
    FileUtil::symlinkAlias(par.db3, seqsDbSeq + ".0" );
    FileUtil::move(seqsDbSeq.c_str(), (seqsDbSeq + ".1").c_str());

    struct DBSuffix {
        DBFiles::Files flag;
        const char* suffix;
    };

    const DBSuffix suffices[] = {
            { DBFiles::HEADER,        "_h"                },
            { DBFiles::HEADER_INDEX,  "_h.index"          },
            { DBFiles::HEADER_DBTYPE, "_h.dbtype"         },
            { DBFiles::LOOKUP,        ".lookup"           },
            { DBFiles::SOURCE,        ".source"           },
            { DBFiles::TAX_MAPPING,   "_mapping"          },
            { DBFiles::TAX_NAMES,     "_names.dmp"        },
            { DBFiles::TAX_NODES,     "_nodes.dmp"        },
            { DBFiles::TAX_MERGED,    "_merged.dmp"       },
            { DBFiles::TAX_MERGED,    "_taxonomy"         },
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
            DBReader<unsigned int>::aliasDb(file, par.db3 + "_seq" + suffices[i].suffix);
        }
    }

    return EXIT_SUCCESS;
}
