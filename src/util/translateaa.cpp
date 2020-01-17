#include "SubstitutionMatrix.h"
#include "Sequence.h"
#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "TranslateNucl.h"

#ifdef OPENMP
#include <omp.h>
#endif

int translateaa(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_NUCLEOTIDES);
    writer.open();

    TranslateNucl translateNucl(static_cast<TranslateNucl::GenCode>(par.translationTable));
    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 2.0f, -0.0f);

    char lookupAA[21][3];
    const char nucLookup[4] = {'A', 'C', 'G', 'T'};
    char data[3];
    char writeAA[1];
    for (size_t i = 0; i < 20; i++) {
        bool found = false;
        for (size_t nuc1 = 0; nuc1 < 4 && found == false; nuc1++) {
            for (size_t nuc2 = 0; nuc2 < 4 && found == false; nuc2++) {
                for (size_t nuc3 = 0; nuc3 < 4 && found == false; nuc3++) {
                    data[0] = nucLookup[nuc1];
                    data[1] = nucLookup[nuc2];
                    data[2] = nucLookup[nuc3];
                    translateNucl.translate(writeAA, data, 3);
                    if (writeAA[0] == subMat.num2aa[i]) {
                        lookupAA[i][0] = data[0];
                        lookupAA[i][1] = data[1];
                        lookupAA[i][2] = data[2];
                        found = true;
                    }
                }
            }
        }
    }

    // add X aa
    lookupAA[20][0] = 'N';
    lookupAA[20][1] = 'N';
    lookupAA[20][2] = 'N';

#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        char *aa = new char[(par.maxSeqLen + 1) / 3 + 3 + 1];
        std::string nucSeq;
        nucSeq.reserve(10000);
        Sequence aaSequence(par.maxSeqLen + 1, Parameters::DBTYPE_AMINO_ACIDS, &subMat, 0, false, par.compBiasCorrection);

#pragma omp for schedule(dynamic, 5)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            unsigned int key = reader.getDbKey(i);
            char *data = reader.getData(i, thread_idx);
            aaSequence.mapSequence(0, key, data, reader.getSeqLen(i));

            // ignore null char at the end
            for (int pos = 0; pos < aaSequence.L; ++pos) {
                nucSeq.append(lookupAA[aaSequence.numSequence[pos]], 3);
            }

            nucSeq.append(1, '\n');
            writer.writeData(nucSeq.c_str(), nucSeq.size(), key, thread_idx);
            nucSeq.clear();
        }
        delete[] aa;
    }
    writer.close(true);
    reader.close();
    DBReader<unsigned int>::softlinkDb(par.db1, par.db2, DBFiles::SEQUENCE_ANCILLARY);

    return EXIT_SUCCESS;
}


