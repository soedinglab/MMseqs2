#include "Util.h"
#include "Parameters.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "StripedSmithWaterman.h"

#include <sys/time.h>

#ifdef OPENMP
#include <omp.h>
#endif

const char* binary_name = "test_ungappedcpuperf";
DEFAULT_PARAMETER_SINGLETON_INIT

#define AA_ALPHABET "ACDEFGHIKLMNPQRSTVWY"
#define ALPHABET_SIZE 20

void generateSequence(char *sequence, int length, unsigned int *seedp) {
    for (int i = 0; i < length; i++) {
        sequence[i] = AA_ALPHABET[rand_r(seedp) % ALPHABET_SIZE];
    }
    sequence[length] = '\0';
}

void generateNumSequence(unsigned char *sequence, int length, unsigned int *seedp) {
    for (int i = 0; i < length; i++) {
        sequence[i] = rand_r(seedp) % ALPHABET_SIZE;
    }
}

int main (int, const char**) {
    Parameters& par = Parameters::getInstance();
    par.initMatrices();

    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);
    int8_t* tinySubMat = new int8_t[subMat.alphabetSize * subMat.alphabetSize];
    for (int i = 0; i < subMat.alphabetSize; i++) {
        for (int j = 0; j < subMat.alphabetSize; j++) {
            tinySubMat[i*subMat.alphabetSize + j] = subMat.subMatrix[i][j];
        }
    }

    size_t targets = 5000000;

    std::vector<int> benchSizes = {
        32, 64, 96, 128, 160, 192, 224, 256, 288, 320, 352, 384, 416, 448, 480, 512, 576, 640, 704, 768, 832, 896, 960, 1024, 1152, 1280, 1408, 1536, 1664, 1792, 1920, 2048
    };

    for (auto seqLen : benchSizes) {
        char *seq = (char *)malloc((seqLen + 1) * sizeof(char));
        unsigned int qseed = 42;
        generateSequence(seq, seqLen, &qseed);

        size_t score = 0;
        double avgTime = 0;
        for (size_t rep = 0; rep < 5; rep++){
            double repTime = 0;
            size_t sanityCheck = 0;
#pragma omp parallel reduction(+:score,sanityCheck) reduction(max:repTime)
{
            unsigned int thread_idx = 0;
    #ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
    #endif
            size_t ignore, total;
            Util::decomposeDomain(targets, thread_idx, par.threads, &ignore, &total);
            sanityCheck += total;

            SmithWaterman aligner(seqLen, subMat.alphabetSize, false, 1.0, Parameters::DBTYPE_AMINO_ACIDS);
            Sequence qSeq(seqLen, Parameters::DBTYPE_AMINO_ACIDS, &subMat, 0, false, false);
            qSeq.mapSequence(0, 0, seq, seqLen);
            aligner.ssw_init(&qSeq, tinySubMat, &subMat);

            unsigned int tseed = 42 + thread_idx;
            unsigned char** targetSeqs = new unsigned char*[total];
            for (size_t i = 0; i < total; i++) {
                targetSeqs[i] = (unsigned char *)malloc(seqLen * sizeof(unsigned char));
                generateNumSequence(targetSeqs[i], seqLen, &tseed);
            }

            struct timeval start;
            struct timeval end;
            gettimeofday(&start, NULL);
            for (size_t i = 0; i < total; i++) {
                score = aligner.ungapped_alignment(targetSeqs[i], seqLen);
            }
            gettimeofday(&end, NULL);
            double diff = (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
            repTime = diff;

            for (size_t i = 0; i < total; i++) {
                free(targetSeqs[i]);
            }
            delete[] targetSeqs;
}
            Debug(Debug::INFO) << "total: " << sanityCheck << "\n";
            avgTime += repTime;
        }
        avgTime /= 5.0f;
        double cells = seqLen * seqLen * targets;
        Debug(Debug::INFO) << score << "\t" << seqLen << "\t" << (cells / (avgTime * 1000000000.0f)) << "\n";
        free(seq);
    }

    delete[] tinySubMat;
}

