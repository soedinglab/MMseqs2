#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "SubstitutionMatrix.h"
#include "tantan.h"
#include "Masker.h"

#ifdef OPENMP
#include <omp.h>
#endif

int makepaddedseqdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    const int mode = DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA;
    DBReader<unsigned int> dbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, mode);
    dbr.open(DBReader<unsigned int>::SORT_BY_LENGTH);

    DBReader<unsigned int> dbhr(par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads, mode);
    dbhr.open(DBReader<unsigned int>::NOSORT);

    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias);

    int dbType = DBReader<unsigned int>::setExtendedDbtype(dbr.getDbtype(), Parameters::DBTYPE_EXTENDED_GPU);
    DBWriter dbsw(par.db2.c_str(), par.db2Index.c_str(), par.threads, false, dbType);
    dbsw.open();
    DBWriter dbhw(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, false, Parameters::DBTYPE_GENERIC_DB);
    dbhw.open();

    // need to prune low scoring k-mers through masking

    Debug::Progress progress(dbr.getSize());
#pragma omp parallel
{
    unsigned int thread_idx = 0;
#ifdef OPENMP
    thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
    Masker masker(subMat);
    std::string result;
    result.reserve(par.maxSeqLen);

    const int ALIGN = 4;
    Sequence seq(dbr.getMaxSeqLen(), dbr.getDbtype(), &subMat,  0, false, false);

    size_t firstIt = SIZE_MAX;
    unsigned int seqKey = 0;

    size_t charSeqBufferSize = par.maxSeqLen + 1;
    unsigned char *charSequence = NULL;
    if (par.maskMode) {
        charSequence = (unsigned char*)malloc(charSeqBufferSize * sizeof(char));
    }

#pragma omp for schedule(static)
    for (size_t i = 0; i < dbr.getSize(); i++) {
        progress.updateProgress();

        if (firstIt == SIZE_MAX) {
            firstIt = i;
        }

        size_t id = dbr.getSize() - 1 - i;
        unsigned int key = dbr.getDbKey(id);
        char *data = dbr.getData(id, thread_idx);
        size_t seqLen = dbr.getSeqLen(id);
        seq.mapSequence(id, key, data, seqLen);

        if (charSequence != NULL) {
            if ((size_t)seq.L >= charSeqBufferSize) {
                charSeqBufferSize = seq.L * 1.5;
                charSequence = (unsigned char*)realloc(charSequence, charSeqBufferSize * sizeof(char));
            }
            memcpy(charSequence, seq.numSequence, seq.L);
            masker.maskSequence(seq, par.maskMode, par.maskProb, par.maskLowerCaseMode, par.maskNrepeats);
            for (int i = 0; i < seq.L; i++) {
                result.append(1, (seq.numSequence[i] == masker.maskLetterNum) ? charSequence[i] + 32 : charSequence[i]);
            }
        } else {
            for (int i = 0; i < seq.L; i++) {
                char aa = data[i];
                result.append(1, (islower(aa)) ? seq.numSequence[i] + 32 : seq.numSequence[i]);
            }
        }
        const size_t sequencepadding = (seq.L % ALIGN == 0) ? 0 : ALIGN - seq.L % ALIGN;
        result.append(sequencepadding, static_cast<char>(20));
        dbsw.writeData(result.c_str(), result.size(), key, thread_idx, false, false);

        // + 2 is needed for newline and null character
        size_t start = dbsw.getStart(thread_idx);
        if (start % 4 != 0) {
            Debug(Debug::ERROR) << "Misalligned entry\n";
            EXIT(EXIT_FAILURE);
        }
        dbsw.writeIndexEntry(firstIt + seqKey, start, seq.L + 2, thread_idx);

        unsigned int headerId = dbhr.getId(key);
        dbhw.writeData(dbhr.getData(headerId, thread_idx), dbhr.getEntryLen(headerId), firstIt + seqKey, thread_idx, false);

        seqKey++;
        result.clear();
    }
    if (charSequence != NULL) {
        free(charSequence);
    }
}
    dbsw.close(true, false);
    dbhw.close(true, false);
    dbhr.close();
    if (par.writeLookup == true) {
        DBReader<unsigned int> readerHeader(par.hdr2.c_str(), par.hdr2Index.c_str(), 1, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        readerHeader.open(DBReader<unsigned int>::NOSORT);
        // create lookup file
        std::string lookupFile = par.db2 + ".lookup";
        FILE* file = FileUtil::openAndDelete(lookupFile.c_str(), "w");
        std::string buffer;
        buffer.reserve(2048);
        DBReader<unsigned int>::LookupEntry entry;
        size_t totalSize = dbr.getSize();
        for (unsigned int id = 0; id < readerHeader.getSize(); id++) {
            char *header = readerHeader.getData(id, 0);
            entry.id = id;
            entry.entryName = Util::parseFastaHeader(header);
            entry.fileNumber = dbr.getDbKey(totalSize - 1 - id);
            readerHeader.lookupEntryToBuffer(buffer, entry);
            int written = fwrite(buffer.c_str(), sizeof(char), buffer.size(), file);
            if (written != (int)buffer.size()) {
                Debug(Debug::ERROR) << "Cannot write to lookup file " << lookupFile << "\n";
                EXIT(EXIT_FAILURE);
            }
            buffer.clear();
        }
        if (fclose(file) != 0) {
            Debug(Debug::ERROR) << "Cannot close file " << lookupFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        readerHeader.close();
    }
    dbr.close();
    return EXIT_SUCCESS;
}