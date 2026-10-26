#include "IndexReader.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Parameters.h"
#include "Util.h"

#include <unordered_set>

int mergedbs(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    if (par.filenames.size() <= 2) {
        Debug(Debug::ERROR) << "Need at least two databases for merging\n";
        EXIT(EXIT_FAILURE);
    }

    const std::vector<std::string> prefices = Util::split(par.mergePrefixes, ",");

    // When filtering, db1 is not only the key universe but also the per-key allow-list:
    // only merged lines whose target key (first column) is listed in the db1 entry are kept.
    const bool filterTarget = par.mergeFilterTarget;
    const int qDbrDataMode = filterTarget ? (DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA)
                                           : DBReader<DBKeyType>::USE_INDEX;
    const int preloadMode = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) ? IndexReader::PRELOAD_INDEX : 0;
    IndexReader qDbr(par.db1, 1, IndexReader::SEQUENCES, preloadMode, qDbrDataMode);

    // skip par.db{1,2}
    const size_t fileCount = par.filenames.size() - 2;
    DBReader<DBKeyType> **filesToMerge = new DBReader<DBKeyType>*[fileCount];
    for (size_t i = 0; i < fileCount; i++) {
        std::string indexName = par.filenames[i + 2] + ".index";
        filesToMerge[i] = new DBReader<DBKeyType>(par.filenames[i + 2].c_str(), indexName.c_str(), 1, DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_INDEX);
        filesToMerge[i]->open(DBReader<DBKeyType>::NOSORT);
    }

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), 1, par.compressed, filesToMerge[0]->getDbtype());
    writer.open();

    Debug(Debug::INFO) << "Merging the results to " << par.db2.c_str() << "\n";
    Debug::Progress progress(qDbr.sequenceReader->getSize());
    char dbKey[255];
    std::unordered_set<DBKeyType> allowedTargets;
    for (size_t id = 0; id < qDbr.sequenceReader->getSize(); id++) {
        progress.updateProgress();
        DBKeyType key = qDbr.sequenceReader->getDbKey(id);

        // build the per-key allow-list of target keys from the db1 entry
        if (filterTarget) {
            allowedTargets.clear();
            char *keyData = qDbr.sequenceReader->getData(id, 0);
            while (keyData != NULL && *keyData != '\0') {
                Util::parseKey(keyData, dbKey);
                allowedTargets.insert(Util::fast_atoi<DBKeyType>(dbKey));
                keyData = Util::skipLine(keyData);
            }
        }

        // get all data for the id from all files
        writer.writeStart(0);
        for (size_t i = 0; i < fileCount; i++) {
            size_t entryId = filesToMerge[i]->getId(key);
            if (entryId == DB_ENTRY_NOT_FOUND) {
                continue;
            }
            const char *data = filesToMerge[i]->getData(entryId, 0);
            if (data == NULL) {
                if (par.mergeStopEmpty == true) {
                    break;
                } else {
                    continue;
                }
            }
            if (i < prefices.size()) {
                writer.writeAdd(prefices[i].c_str(), prefices[i].size(), 0);
            }
            if (filterTarget == false) {
                writer.writeAdd(data, filesToMerge[i]->getEntryLen(entryId) - 1, 0);
            } else {
                // keep allow-listed targets; erase-on-write emits each target at most once
                // (dedups the self line and cross-file duplicates), so line count == cluster size.
                const char *line = data;
                while (*line != '\0') {
                    const char *next = Util::skipLine((char *) line);
                    Util::parseKey(line, dbKey);
                    std::unordered_set<DBKeyType>::iterator it = allowedTargets.find(Util::fast_atoi<DBKeyType>(dbKey));
                    if (it != allowedTargets.end()) {
                        writer.writeAdd(line, next - line, 0);
                        allowedTargets.erase(it);
                    }
                    line = next;
                }
            }
        }
        writer.writeEnd(key, 0);
    }
    writer.close();
    for (size_t i = 0; i < fileCount; i++) {
        filesToMerge[i]->close();
        delete filesToMerge[i];
    }
    delete[] filesToMerge;

    return EXIT_SUCCESS;
}
