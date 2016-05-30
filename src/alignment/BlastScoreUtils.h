//
// Created by mad on 10/4/15.
//

#ifndef MMSEQS_BLASTSCOREUTILS_H
#define MMSEQS_BLASTSCOREUTILS_H

#include <cmath>

class BlastScoreUtils {
public:
    struct BlastStat{
        const double lambda;
        const double K;
        const double H;
        const double alpha;
        const double beta;
        BlastStat(double lambda, double K, double H, double alpha, double beta)
                : lambda(lambda), K(K), H(H), alpha(alpha), beta(beta)
        {};
    };


#define BLOSUM45_VALUES_MAX 14
    static double blosum45_values[BLOSUM45_VALUES_MAX][8];
    static int blosum45_prefs[BLOSUM45_VALUES_MAX];


#define BLOSUM50_VALUES_MAX 16
    static double blosum50_values[BLOSUM50_VALUES_MAX][8];

    static int blosum50_prefs[BLOSUM50_VALUES_MAX];

#define BLOSUM62_VALUES_MAX 12
    static double blosum62_values[BLOSUM62_VALUES_MAX][8];
    static int blosum62_prefs[BLOSUM62_VALUES_MAX];


#define BLOSUM80_VALUES_MAX 10
    static double blosum80_values[BLOSUM80_VALUES_MAX][8];
    static int blosum80_prefs[BLOSUM80_VALUES_MAX];

#define BLOSUM90_VALUES_MAX 8
    static double blosum90_values[BLOSUM90_VALUES_MAX][8];

    static int blosum90_prefs[BLOSUM90_VALUES_MAX];

#define BLOSUM62_20_VALUES_MAX 65
    static double blosum62_20_values[BLOSUM62_20_VALUES_MAX][8];

    static int blosum62_20_prefs[BLOSUM62_20_VALUES_MAX];


/** Supported substitution and gap costs with corresponding quality values
 * for nucleotide sequence comparisons.
 * NB: the values 0 and 0 for the gap costs are treated as the defaults used for
 * the greedy gapped extension, i.e.
 * gap opening = 0,
 * gap extension = 1/2 match - mismatch.
 *
 * The fields are:
 *
 * 1. Gap opening cost,
 * 2. Gap extension cost,
 * 3. Lambda,
 * 4. K,
 * 5. H,
 * 6. Alpha,
 * 7. Beta,
 * 8. Theta
 */

/** Karlin-Altschul parameter values for substitution scores 1 and -5. */
    static const double blastn_values_1_5[][8];

/** Karlin-Altschul parameter values for substitution scores 1 and -4. */
    static const double blastn_values_1_4[][8];

/** Karlin-Altschul parameter values for substitution scores 2 and -7.
 * These parameters can only be applied to even scores. Any odd score must be
 * rounded down to the nearest even number before calculating the e-value.
 */
    static const double blastn_values_2_7[][8];
/** Karlin-Altschul parameter values for substitution scores 1 and -3. */
    static const double blastn_values_1_3[][8];

/** Karlin-Altschul parameter values for substitution scores 2 and -5.
 * These parameters can only be applied to even scores. Any odd score must be
 * rounded down to the nearest even number before calculating the e-value.
 */
    static const double blastn_values_2_5[][8];
/** Karlin-Altschul parameter values for substitution scores 1 and -2. */
    static const double blastn_values_1_2[][8];

/** Karlin-Altschul parameter values for substitution scores 2 and -3.
 * These parameters can only be applied to even scores. Any odd score must be
 * rounded down to the nearest even number before calculating the e-value.
 */
    static const double blastn_values_2_3[][8];
/** Karlin-Altschul parameter values for substitution scores 3 and -4. */
    static const double blastn_values_3_4[][8];
/** Karlin-Altschul parameter values for substitution scores 4 and -5. */
    static const double blastn_values_4_5[][8];
/** Karlin-Altschul parameter values for substitution scores 1 and -1. */
    static const double blastn_values_1_1[][8];
/** Karlin-Altschul parameter values for substitution scores 3 and -2. */
    static const double blastn_values_3_2[][8];
/** Karlin-Altschul parameter values for substitution scores 5 and -4. */
    static const double blastn_values_5_4[][8];

    static int BlastComputeLengthAdjustment(double K, double logK, double alpha_d_lambda, double beta,
                                            int query_length, size_t db_length, size_t db_num_seqs);

    static BlastStat getAltschulStatsForMatrix(std::string matrix, long gopen, long gextend);

    static double computeKmn(int qlen, double K, double lambda, double alpha, double beta, size_t dbLen,
                      size_t seqCnt);

    static inline double computeBitScore(double score, double lambdaLog2, double logKLog2) {
        return lambdaLog2 * score - logKLog2;
    }

    static inline double bitScoreToRawScore(double bitScore, double lambdaLog2, double logKLog2){
        return (bitScore + logKLog2) / lambdaLog2;
    }

    static inline double computeEvalue(double score, double Kmn, double lambda) {
        return Kmn * exp(-lambda * score);
    }
};


#endif //MMSEQS_BLASTSCOREUTILS_H
