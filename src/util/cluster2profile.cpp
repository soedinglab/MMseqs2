// Computes either a PSSM or a MSA from clustering or alignment result
// For PSSMs: MMseqs just stores the position specific score in 1 byte
//

#include <string>
#include <vector>
#include <sstream>

#include "Matcher.h"
#include "SubstitutionMatrix.h"
#include "Parameters.h"
#include "Sequence.h"
#include "MultipleAlignment.h"
#include "PSSMCalculator.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Log.h"
#include "Util.h"
#include "Debug.h"

#ifdef OPENMP
#include <omp.h>
#endif

enum {
    MSA = 0,
    PSSM
};

int result2outputmode(int argn, const char **argv, int mode) {
    std::string usage("Converts a ffindex profile database to ffindex.\n");
    usage.append("USAGE: <queryDB> <targetDB> <resultDB> <outDB>\n");
    usage.append("\nDesigned and implemented by Martin Steinegger <martin.steinegger@mpibpc.mpg.de>\n");
    usage.append("\t & Milot Mirdita <milot@mirdita.de>");

    Parameters par;
    std::vector<MMseqsParameter> parameters;
    switch(mode) {
        case PSSM:
            parameters = par.result2profile;
            break;
        case MSA:
            parameters = par.result2msa;
            break;
        default:
            Debug(Debug::ERROR) << "Output mode not implemented!\n";
            EXIT(EXIT_FAILURE);

    }
    par.parseParameters(argn, argv, usage, parameters, 4);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<unsigned int> queryReader(par.db1.c_str(), par.db1Index.c_str());
    queryReader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int>* templateReader = &queryReader;

    bool sameDatabase = true;
    if (par.db1.compare(par.db2) != 0) {
        sameDatabase = false;
        templateReader = new DBReader<unsigned int>(DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str()));
        templateReader->open(DBReader<unsigned int>::NOSORT);
    }

    DBReader<unsigned int> clusterReader(par.db3.c_str(), par.db3Index.c_str());
    clusterReader.open(DBReader<unsigned int>::NOSORT);

    DBWriter::errorIfFileExist(par.db4.c_str());
    DBWriter::errorIfFileExist(par.db4Index.c_str());

    DBWriter writer(par.db4.c_str(), par.db4Index.c_str(), par.threads, DBWriter::BINARY_MODE);
    writer.open();

    // find longest sequence
    size_t maxSequenceLength = 0;
    for (size_t i = 0; i < queryReader.getSize(); i++) {
        maxSequenceLength = std::max((size_t) queryReader.getSeqLens()[i], maxSequenceLength);
    }
    for (size_t i = 0; i < templateReader->getSize(); i++) {
        maxSequenceLength = std::max((size_t) templateReader->getSeqLens()[i], maxSequenceLength);
    }

    //Find the max set size
    size_t maxSetSize = 0;
    for (size_t i = 0; i < clusterReader.getSize(); i++) {
        char *data = clusterReader.getData(i);
        size_t entry = 0;
        size_t position = 0;
        while (data[position] != '\0') {
            if (data[position] == '\n') {
                entry++;
            }
            position++;
        }
        maxSetSize = std::max(maxSetSize, entry);
    }
    maxSetSize += 1;

    bool noDeletion;
    switch (mode) {
        case MSA: noDeletion = true; break;
        case PSSM: noDeletion = !par.allowDeletion; break;
        default: EXIT(EXIT_FAILURE);
    }

    SubstitutionMatrix matrix(par.scoringMatrixFile.c_str(), 2.0f, -0.2f);
    Debug(Debug::INFO) << "Start computing " << (!mode ? "MSAs" : "profiles") << ".\n";
#pragma omp parallel
    {
        Matcher matcher(maxSequenceLength, &matrix, templateReader->getAminoAcidDBSize(), templateReader->getSize(),
                        par.compBiasCorrection);
        MultipleAlignment aligner(maxSequenceLength, maxSetSize, &matrix, &matcher);
        PSSMCalculator calculator(&matrix, maxSequenceLength);

        Sequence centerSequence(maxSequenceLength, matrix.aa2int, matrix.int2aa, Sequence::AMINO_ACIDS, 0, false);
        Sequence **sequences = new Sequence *[maxSetSize];
        for (size_t i = 0; i < maxSetSize; i++) {
            sequences[i] = new Sequence(maxSequenceLength, matrix.aa2int, matrix.int2aa, Sequence::AMINO_ACIDS, 0,
                                        false);
        }

#pragma omp for schedule(dynamic, 1000)
        for (size_t id = 0; id < clusterReader.getSize(); id++) {
            Log::printProgress(id);

            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif

            char *clusters = clusterReader.getData(id);
            unsigned int queryId = clusterReader.getDbKey(id);
            char *seqData = queryReader.getDataByDBKey(queryId);
            centerSequence.mapSequence(0, queryId, seqData);
            std::vector<Sequence *> seqSet;
            size_t position = 0;

            char dbKey[255 + 1];
            while (*clusters != '\0') {
                Util::parseKey(clusters, dbKey);
                unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
                if (key != queryId || !sameDatabase) {
                    char *dbSeqData = templateReader->getDataByDBKey(key);
                    sequences[position]->mapSequence(0, key, dbSeqData);
                    seqSet.push_back(sequences[position]);
                    position++;
                }
                clusters = Util::skipLine(clusters);
            }

            MultipleAlignment::MSAResult res = aligner.computeMSA(&centerSequence, seqSet, noDeletion);

            std::stringstream msa;
            std::string result;

            char *data;
            size_t dataSize;
            switch (mode) {
                case MSA:
                    for (size_t i = 0; i < res.setSize; i++) {
                        msa << ">" << sequences[i]->getDbKey() << "\n";
                        msa << std::string(res.msaSequence[i], 0, res.msaSequenceLength) << "\n";
                    }
                    result = msa.str();
                    data = (char *) result.c_str();
                    dataSize = result.length();
                    break;
                case PSSM:
                    data = (char *) calculator.computePSSMFromMSA(res.setSize, res.centerLength,
                                                                  res.msaSequence);
                    dataSize = res.centerLength * Sequence::PROFILE_AA_SIZE * sizeof(char);

                    for (size_t i = 0; i < dataSize; i++) {
                        // Avoid a null byte result
                        data[i] = data[i] ^ 0x80;
                    }
                    //pssm.printProfile(res.centerLength);
                    //calculator.printPSSM(res.centerLength);
                    break;
                default:
                    EXIT(EXIT_FAILURE);
            }

            writer.write(data, dataSize, SSTR(queryId).c_str(), thread_idx);
        }

        for (size_t i = 0; i < maxSetSize; i++) {
            delete sequences[i];
        }
        delete[] sequences;
    }

    // cleanup
    clusterReader.close();
    queryReader.close();
    if (!sameDatabase) {
        templateReader->close();
        delete templateReader;
    }
    writer.close();

    Debug(Debug::INFO) << "\nDone.\n";

    return EXIT_SUCCESS;
}

int result2profile(int argc, const char **argv) {
    return result2outputmode(argc, argv, PSSM);
}

int result2msa(int argc, const char **argv) {
    return result2outputmode(argc, argv, MSA);
}
