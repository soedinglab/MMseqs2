#include <string>
#include <fstream>
#include <climits>
#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"
#include "tantan.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"


#ifdef OPENMP
#include <omp.h>

#endif

int masksequence(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::NOSORT);

    BaseMatrix *subMat;
    if (Parameters::isEqualDbtype(reader.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.c_str(), 1.0, 0.0);
    } else {
        // keep score bias at 0.0 (improved ROC)
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.c_str(), 2.0, 0.0);
    }
    size_t maxSeqLen = 0;
    for(size_t i = 0; i < reader.getSize(); i++){
        maxSeqLen = std::max(reader.getSeqLens(i), maxSeqLen);
    }
    // need to prune low scoring k-mers through masking
    ProbabilityMatrix probMatrix(*subMat);

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, reader.getDbtype());
    writer.open();
#pragma omp parallel
    {

        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        char *charSequence = new char[maxSeqLen];

#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < reader.getSize(); ++id) {
            char *seqData = reader.getData(id, thread_idx);
            unsigned int qKey = reader.getDbKey(id);
            unsigned int seqLen = 0;
            while (seqData[seqLen] != '\0') {
                charSequence[seqLen] = (char) subMat->aa2int[(int)seqData[seqLen]];
                seqLen++;
            }
            // s.print();
            tantan::maskSequences(charSequence,
                                  charSequence + seqLen,
                                  50 /*options.maxCycleLength*/,
                                  probMatrix.probMatrixPointers,
                                  0.005 /*options.repeatProb*/,
                                  0.05 /*options.repeatEndProb*/,
                                  0.9 /*options.repeatOffsetProbDecay*/,
                                  0, 0,
                                  0.5 /*options.minMaskProb*/,
                                  probMatrix.hardMaskTable);

            for (unsigned int pos = 0; pos < seqLen; pos++) {
                char aa = seqData[pos];
                charSequence[pos] = (charSequence[pos] == probMatrix.hardMaskTable[0]) ? tolower(aa) : toupper(aa);
            }
            writer.writeData(charSequence, seqLen, qKey, thread_idx);
        }
        delete [] charSequence;
    }
    writer.close(true);
    reader.close();

    FileUtil::symlinkAbs(par.hdr1, par.hdr2);
    FileUtil::symlinkAbs(par.hdr1Index, par.hdr2Index);
    DBWriter::writeDbtypeFile(par.hdr2.c_str(), Parameters::DBTYPE_GENERIC_DB, par.compressed);

    delete subMat;


    return EXIT_SUCCESS;
}
