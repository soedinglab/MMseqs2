// Computes either a PSSM or a MSA from clustering or alignment result
// For PSSMs: MMseqs just stores the position specific score in 1 byte

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <utility>

#include "Parameters.h"

#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

struct compareFirstEntry {
    bool
    operator()(const std::pair<std::string, unsigned int> &lhs, const std::pair<std::string, unsigned int> &rhs) const {
        return (lhs.first.compare(rhs.first) <= 0);
    }
};

struct compareKeyToFirstEntry {
    bool operator()(const std::pair<std::string, unsigned int> &lhs, const std::string &rhs) const {
        return  (lhs.first < rhs);
    }

    bool operator()(const std::string &lhs, const std::pair<std::string, unsigned int> &rhs) const {
        return  (lhs < rhs.first);
    }
};

int diffseqdbs(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 5);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    std::string headerDBold(par.db1);
    headerDBold.append("_h");
    std::string headerDBoldIndex(par.db1);
    headerDBoldIndex.append("_h.index");

    std::string headerDBnew(par.db2);
    headerDBnew.append("_h");
    std::string headerDBnewIndex(par.db2);
    headerDBnewIndex.append("_h.index");

    DBReader<unsigned int> oldReader(headerDBold.c_str(), headerDBoldIndex.c_str());
    oldReader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> newReader(headerDBnew.c_str(), headerDBnewIndex.c_str());
    newReader.open(DBReader<unsigned int>::NOSORT);

    std::ofstream removedSeqDBWriter, keptSeqDBWriter, newSeqDBWriter;
    removedSeqDBWriter.open(par.db3);
    keptSeqDBWriter.open(par.db4);
    newSeqDBWriter.open(par.db5);

    // Fill up the hash tables for the old and new DB
    size_t indexSizeOld = oldReader.getSize();
    // keys pairs are like : (headerID,key) where key is the ffindex key corresponding to the header
    std::pair<std::string, unsigned int> *keysOld
            = new std::pair<std::string, unsigned int>[indexSizeOld];

    #pragma omp parallel for schedule(static)
    for (size_t id = 0; id < indexSizeOld; ++id) {
        if (par.useSequenceId) {
            keysOld[id] = std::make_pair(
                    Util::parseFastaHeader(oldReader.getData(id)),
                    oldReader.getDbKey(id)
            );
        } else {
            keysOld[id] = std::make_pair(
                    oldReader.getData(id),
                    oldReader.getDbKey(id)
            );
        }
    }

    size_t indexSizeNew = newReader.getSize();
    std::pair<std::string, unsigned int> *keysNew
            = new std::pair<std::string, unsigned int>[indexSizeNew];

    #pragma omp parallel for schedule(static)
    for (size_t id = 0; id < indexSizeNew; ++id) {
        if (par.useSequenceId) {
            keysNew[id] = std::make_pair(
                    Util::parseFastaHeader(newReader.getData(id)),
                    newReader.getDbKey(id)
            );
        } else {
            keysNew[id] = std::make_pair(
                    newReader.getData(id),
                    newReader.getDbKey(id)
            );
        }
    }

    //sort by header for binary search
    std::stable_sort(keysNew, keysNew + indexSizeNew, compareFirstEntry());

    // default initialized with false
    bool* checkedNew = new bool[indexSizeNew]();
    // doesn't need to be initialized
    size_t *mappedIds = new size_t[indexSizeNew];

    bool* deletedIds = new bool[indexSizeOld]();

    #pragma omp parallel for schedule(static)
    for (size_t id = 0; id < indexSizeOld; ++id) {
        const std::string &keyToSearch = keysOld[id].first;
        std::pair<std::string, unsigned int> *mappedKey
                = std::lower_bound(keysNew, keysNew + indexSizeNew, keyToSearch, compareKeyToFirstEntry());

        if (mappedKey != (keysNew + indexSizeNew) && keyToSearch.compare(mappedKey->first) == 0) {
            // Found
            size_t indexInNewDB = (mappedKey - keysNew);
            checkedNew[indexInNewDB] = true;
            mappedIds[indexInNewDB] = id;
        } else {
            // Not found
            deletedIds[id] = true;
        }
    }

    for (size_t i = 0; i < indexSizeOld; ++i) {
        if(deletedIds[i]) {
            removedSeqDBWriter << keysOld[i].second << std::endl;
        }
    }
    removedSeqDBWriter.close();

    for (size_t id = 0; id < indexSizeNew; ++id) {
        if (checkedNew[id]) {
            keptSeqDBWriter << keysOld[mappedIds[id]].second << "\t" << keysNew[id].second << std::endl;
        } else {
            newSeqDBWriter << keysNew[id].second << std::endl;
        }
    }
    newSeqDBWriter.close();
    keptSeqDBWriter.close();

    delete[] deletedIds;
    delete[] mappedIds;
    delete[] checkedNew;
    delete[] keysNew;
    delete[] keysOld;

    newReader.close();
    oldReader.close();

    return EXIT_SUCCESS;
}

