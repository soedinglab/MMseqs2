//
// Created by mad on 10/4/15.
//

#include <algorithm>
#include <cmath>
#include <Debug.h>
#include <Util.h>
#include "BlastScoreUtils.h"
/* ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================*/


/**************************************************************************************

How the statistical parameters for the matrices are stored:
-----------------------------------------------------------
They parameters are stored in a two-dimensional array FloatHi (i.e.,
doubles, which has as it's first dimensions the number of different
gap existence and extension combinations and as it's second dimension 8.
The eight different columns specify:

1.) gap existence penalty (INT2_MAX denotes infinite).
2.) gap extension penalty (INT2_MAX denotes infinite).
3.) decline to align penalty (INT2_MAX denotes infinite).
4.) Lambda
5.) K
6.) H
7.) alpha
8.) beta

(Items 4-8 are explained in:
Altschul SF, Bundschuh R, Olsen R, Hwa T.
The estimation of statistical parameters for local alignment score distributions.
Nucleic Acids Res. 2001 Jan 15;29(2):351-61.).

Take BLOSUM45 (below) as an example.  Currently (5/17/02) there are
14 different allowed combinations (specified by "#define BLOSUM45_VALUES_MAX 14").
The first row in the array "blosum45_values" has INT2_MAX (i.e., INT2_MAX) for gap
existence, extension, and decline-to-align penalties.  For all practical purposes
this value is large enough to be infinite, so the alignments will be ungapped.
BLAST may also use this value (INT2_MAX) as a signal to skip gapping, so a
different value should not be used if the intent is to have gapless extensions.
The next row is for the gap existence penalty 13 and the extension penalty 3.
The decline-to-align penalty is only supported in a few cases, so it is normally
set to INT2_MAX.


How to add a new matrix to blastkar.c:
--------------------------------------
To add a new matrix to blastkar.c it is necessary to complete
four steps.  As an example consider adding the matrix
called TESTMATRIX

1.) add a define specifying how many different existence and extensions
penalties are allowed, so it would be necessary to add the line:

#define TESTMATRIX_VALUES_MAX 14

if 14 values were to be allowed.

2.) add a two-dimensional array to contain the statistical parameters:

static double testmatrix_values[TESTMATRIX_VALUES_MAX][8] ={ ...

3.) add a "prefs" array that should hint about the "optimal"
gap existence and extension penalties:

static int testmatrix_prefs[TESTMATRIX_VALUES_MAX] = {
BLAST_MATRIX_NOMINAL,
...
};

4.) Go to the function BlastLoadMatrixValues (in this file) and
add two lines before the return at the end of the function:

        matrix_info = MatrixInfoNew("TESTMATRIX", testmatrix_values, testmatrix_prefs, TESTMATRIX_VALUES_MAX);
        ValNodeAddPointer(&retval, 0, matrix_info);



***************************************************************************************/

#define BLAST_MATRIX_NOMINAL 0
#define BLAST_MATRIX_BEST 1
#define INT2_MAX 32767
double BlastScoreUtils::blosum45_values[BLOSUM45_VALUES_MAX][8] = {
        {(double) INT2_MAX, (double) INT2_MAX, (double) INT2_MAX, 0.2291, 0.0924, 0.2514, 0.9113, -5.7},
        {13, 3, (double) INT2_MAX, 0.207, 0.049, 0.14, 1.5, -22},
        {12, 3, (double) INT2_MAX, 0.199, 0.039, 0.11, 1.8, -34},
        {11, 3, (double) INT2_MAX, 0.190, 0.031, 0.095, 2.0, -38},
        {10, 3, (double) INT2_MAX, 0.179, 0.023, 0.075, 2.4, -51},
        {16, 2, (double) INT2_MAX, 0.210, 0.051, 0.14, 1.5, -24},
        {15, 2, (double) INT2_MAX, 0.203, 0.041, 0.12, 1.7, -31},
        {14, 2, (double) INT2_MAX, 0.195, 0.032, 0.10, 1.9, -36},
        {13, 2, (double) INT2_MAX, 0.185, 0.024, 0.084, 2.2, -45},
        {12, 2, (double) INT2_MAX, 0.171, 0.016, 0.061, 2.8, -65},
        {19, 1, (double) INT2_MAX, 0.205, 0.040, 0.11, 1.9, -43},
        {18, 1, (double) INT2_MAX, 0.198, 0.032, 0.10, 2.0, -43},
        {17, 1, (double) INT2_MAX, 0.189, 0.024, 0.079, 2.4, -57},
        {16, 1, (double) INT2_MAX, 0.176, 0.016, 0.063, 2.8, -67},
};

int BlastScoreUtils::blosum45_prefs[BLOSUM45_VALUES_MAX] = {
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_BEST,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL
};


#define BLOSUM50_VALUES_MAX 16
double BlastScoreUtils::blosum50_values[BLOSUM50_VALUES_MAX][8] = {
        {(double) INT2_MAX, (double) INT2_MAX, (double) INT2_MAX, 0.2318, 0.112, 0.3362, 0.6895, -4.0},
        {13, 3, (double) INT2_MAX, 0.212, 0.063, 0.19, 1.1, -16},
        {12, 3, (double) INT2_MAX, 0.206, 0.055, 0.17, 1.2, -18},
        {11, 3, (double) INT2_MAX, 0.197, 0.042, 0.14, 1.4, -25},
        {10, 3, (double) INT2_MAX, 0.186, 0.031, 0.11, 1.7, -34},
        {9, 3, (double) INT2_MAX, 0.172, 0.022, 0.082, 2.1, -48},
        {16, 2, (double) INT2_MAX, 0.215, 0.066, 0.20, 1.05, -15},
        {15, 2, (double) INT2_MAX, 0.210, 0.058, 0.17, 1.2, -20},
        {14, 2, (double) INT2_MAX, 0.202, 0.045, 0.14, 1.4, -27},
        {13, 2, (double) INT2_MAX, 0.193, 0.035, 0.12, 1.6, -32},
        {12, 2, (double) INT2_MAX, 0.181, 0.025, 0.095, 1.9, -41},
        {19, 1, (double) INT2_MAX, 0.212, 0.057, 0.18, 1.2, -21},
        {18, 1, (double) INT2_MAX, 0.207, 0.050, 0.15, 1.4, -28},
        {17, 1, (double) INT2_MAX, 0.198, 0.037, 0.12, 1.6, -33},
        {16, 1, (double) INT2_MAX, 0.186, 0.025, 0.10, 1.9, -42},
        {15, 1, (double) INT2_MAX, 0.171, 0.015, 0.063, 2.7, -76},
};

int BlastScoreUtils::blosum50_prefs[BLOSUM50_VALUES_MAX] = {
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_BEST,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL
};

double BlastScoreUtils::blosum62_values[BLOSUM62_VALUES_MAX][8] = {
        {(double) INT2_MAX, (double) INT2_MAX, (double) INT2_MAX, 0.3176, 0.134, 0.4012, 0.7916, -3.2},
        {11, 2, (double) INT2_MAX, 0.297, 0.082, 0.27, 1.1, -10},
        {10, 2, (double) INT2_MAX, 0.291, 0.075, 0.23, 1.3, -15},
        {9, 2, (double) INT2_MAX, 0.279, 0.058, 0.19, 1.5, -19},
        {8, 2, (double) INT2_MAX, 0.264, 0.045, 0.15, 1.8, -26},
        {7, 2, (double) INT2_MAX, 0.239, 0.027, 0.10, 2.5, -46},
        {6, 2, (double) INT2_MAX, 0.201, 0.012, 0.061, 3.3, -58},
        {13, 1, (double) INT2_MAX, 0.292, 0.071, 0.23, 1.2, -11},
        {12, 1, (double) INT2_MAX, 0.283, 0.059, 0.19, 1.5, -19},
        {11, 1, (double) INT2_MAX, 0.267, 0.041, 0.14, 1.9, -30},
        {10, 1, (double) INT2_MAX, 0.243, 0.024, 0.10, 2.5, -44},
        {9, 1, (double) INT2_MAX, 0.206, 0.010, 0.052, 4.0, -87},
};

int BlastScoreUtils::blosum62_prefs[BLOSUM62_VALUES_MAX] = {
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_BEST,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
};


#define BLOSUM80_VALUES_MAX 10
double BlastScoreUtils::blosum80_values[BLOSUM80_VALUES_MAX][8] = {
        {(double) INT2_MAX, (double) INT2_MAX, (double) INT2_MAX, 0.3430, 0.177, 0.6568, 0.5222, -1.6},
        {25, 2, (double) INT2_MAX, 0.342, 0.17, 0.66, 0.52, -1.6},
        {13, 2, (double) INT2_MAX, 0.336, 0.15, 0.57, 0.59, -3},
        {9, 2, (double) INT2_MAX, 0.319, 0.11, 0.42, 0.76, -6},
        {8, 2, (double) INT2_MAX, 0.308, 0.090, 0.35, 0.89, -9},
        {7, 2, (double) INT2_MAX, 0.293, 0.070, 0.27, 1.1, -14},
        {6, 2, (double) INT2_MAX, 0.268, 0.045, 0.19, 1.4, -19},
        {11, 1, (double) INT2_MAX, 0.314, 0.095, 0.35, 0.90, -9},
        {10, 1, (double) INT2_MAX, 0.299, 0.071, 0.27, 1.1, -14},
        {9, 1, (double) INT2_MAX, 0.279, 0.048, 0.20, 1.4, -19},
};

int BlastScoreUtils::blosum80_prefs[BLOSUM80_VALUES_MAX] = {
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_BEST,
        BLAST_MATRIX_NOMINAL
};

#define BLOSUM90_VALUES_MAX 8
double BlastScoreUtils::blosum90_values[BLOSUM90_VALUES_MAX][8] = {
        {(double) INT2_MAX, (double) INT2_MAX, (double) INT2_MAX, 0.3346, 0.190, 0.7547, 0.4434, -1.4},
        {9, 2, (double) INT2_MAX, 0.310, 0.12, 0.46, 0.67, -6},
        {8, 2, (double) INT2_MAX, 0.300, 0.099, 0.39, 0.76, -7},
        {7, 2, (double) INT2_MAX, 0.283, 0.072, 0.30, 0.93, -11},
        {6, 2, (double) INT2_MAX, 0.259, 0.048, 0.22, 1.2, -16},
        {11, 1, (double) INT2_MAX, 0.302, 0.093, 0.39, 0.78, -8},
        {10, 1, (double) INT2_MAX, 0.290, 0.075, 0.28, 1.04, -15},
        {9, 1, (double) INT2_MAX, 0.265, 0.044, 0.20, 1.3, -19},
};

int BlastScoreUtils::blosum90_prefs[BLOSUM90_VALUES_MAX] = {
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_BEST,
        BLAST_MATRIX_NOMINAL
};

#define BLOSUM62_20_VALUES_MAX 65
double BlastScoreUtils::blosum62_20_values[BLOSUM62_20_VALUES_MAX][8] = {
        {(double) INT2_MAX, (double) INT2_MAX, (double) INT2_MAX, 0.03391, 0.125, 0.4544, 0.07462, -3.2},
        {100, 12, (double) INT2_MAX, 0.0300, 0.056, 0.21, 0.14, -15},
        {95, 12, (double) INT2_MAX, 0.0291, 0.047, 0.18, 0.16, -20},
        {90, 12, (double) INT2_MAX, 0.0280, 0.038, 0.15, 0.19, -28},
        {85, 12, (double) INT2_MAX, 0.0267, 0.030, 0.13, 0.21, -31},
        {80, 12, (double) INT2_MAX, 0.0250, 0.021, 0.10, 0.25, -39},
        {105, 11, (double) INT2_MAX, 0.0301, 0.056, 0.22, 0.14, -16},
        {100, 11, (double) INT2_MAX, 0.0294, 0.049, 0.20, 0.15, -17},
        {95, 11, (double) INT2_MAX, 0.0285, 0.042, 0.16, 0.18, -25},
        {90, 11, (double) INT2_MAX, 0.0271, 0.031, 0.14, 0.20, -28},
        {85, 11, (double) INT2_MAX, 0.0256, 0.023, 0.10, 0.26, -46},
        {115, 10, (double) INT2_MAX, 0.0308, 0.062, 0.22, 0.14, -20},
        {110, 10, (double) INT2_MAX, 0.0302, 0.056, 0.19, 0.16, -26},
        {105, 10, (double) INT2_MAX, 0.0296, 0.050, 0.17, 0.17, -27},
        {100, 10, (double) INT2_MAX, 0.0286, 0.041, 0.15, 0.19, -32},
        {95, 10, (double) INT2_MAX, 0.0272, 0.030, 0.13, 0.21, -35},
        {90, 10, (double) INT2_MAX, 0.0257, 0.022, 0.11, 0.24, -40},
        {85, 10, (double) INT2_MAX, 0.0242, 0.017, 0.083, 0.29, -51},
        {115, 9, (double) INT2_MAX, 0.0306, 0.061, 0.24, 0.13, -14},
        {110, 9, (double) INT2_MAX, 0.0299, 0.053, 0.19, 0.16, -23},
        {105, 9, (double) INT2_MAX, 0.0289, 0.043, 0.17, 0.17, -23},
        {100, 9, (double) INT2_MAX, 0.0279, 0.036, 0.14, 0.20, -31},
        {95, 9, (double) INT2_MAX, 0.0266, 0.028, 0.12, 0.23, -37},
        {120, 8, (double) INT2_MAX, 0.0307, 0.062, 0.22, 0.14, -18},
        {115, 8, (double) INT2_MAX, 0.0300, 0.053, 0.20, 0.15, -19},
        {110, 8, (double) INT2_MAX, 0.0292, 0.046, 0.17, 0.17, -23},
        {105, 8, (double) INT2_MAX, 0.0280, 0.035, 0.14, 0.20, -31},
        {100, 8, (double) INT2_MAX, 0.0266, 0.026, 0.12, 0.23, -37},
        {125, 7, (double) INT2_MAX, 0.0306, 0.058, 0.22, 0.14, -18},
        {120, 7, (double) INT2_MAX, 0.0300, 0.052, 0.19, 0.16, -23},
        {115, 7, (double) INT2_MAX, 0.0292, 0.044, 0.17, 0.17, -24},
        {110, 7, (double) INT2_MAX, 0.0279, 0.032, 0.14, 0.20, -31},
        {105, 7, (double) INT2_MAX, 0.0267, 0.026, 0.11, 0.24, -41},
        {120,10,5, 0.0298, 0.049, 0.19, 0.16, -21},
        {115,10,5, 0.0290, 0.042, 0.16, 0.18, -25},
        {110,10,5, 0.0279, 0.033, 0.13, 0.21, -32},
        {105,10,5, 0.0264, 0.024, 0.10, 0.26, -46},
        {100,10,5, 0.0250, 0.018, 0.081, 0.31, -56},
        {125,10,4, 0.0301, 0.053, 0.18, 0.17, -25},
        {120,10,4, 0.0292, 0.043, 0.15, 0.20, -33},
        {115,10,4, 0.0282, 0.035, 0.13, 0.22, -36},
        {110,10,4, 0.0270, 0.027, 0.11, 0.25, -41},
        {105,10,4, 0.0254, 0.020, 0.079, 0.32, -60},
        {130,10,3, 0.0300, 0.051, 0.17, 0.18, -27},
        {125,10,3, 0.0290, 0.040, 0.13, 0.22, -38},
        {120,10,3, 0.0278, 0.030, 0.11, 0.25, -44},
        {115,10,3, 0.0267, 0.025, 0.092, 0.29, -52},
        {110,10,3, 0.0252, 0.018, 0.070, 0.36, -70},
        {135,10,2, 0.0292, 0.040, 0.13, 0.22, -35},
        {130,10,2, 0.0283, 0.034, 0.10, 0.28, -51},
        {125,10,2, 0.0269, 0.024, 0.077, 0.35, -71},
        {120,10,2, 0.0253, 0.017, 0.059, 0.43, -90},
        {115,10,2, 0.0234, 0.011, 0.043, 0.55, -121},
        {100,14,3, 0.0258, 0.023, 0.087, 0.33, -59},
        {105,13,3, 0.0263, 0.024, 0.085, 0.31, -57},
        {110,12,3, 0.0271, 0.028, 0.093, 0.29, -54},
        {115,11,3, 0.0275, 0.030, 0.10, 0.27, -49},
        {125,9,3, 0.0283, 0.034, 0.12, 0.23, -38},
        {130,8,3, 0.0287, 0.037, 0.12, 0.23, -40},
        {125,7,3, 0.0287, 0.036, 0.12, 0.24, -44},
        {140,6,3, 0.0285, 0.033, 0.12, 0.23, -40},
        {105,14,3, 0.0270, 0.028, 0.10, 0.27, -46},
        {110,13,3, 0.0279, 0.034, 0.10, 0.27, -50},
        {115,12,3, 0.0282, 0.035, 0.12, 0.24, -42},
        {120,11,3, 0.0286, 0.037, 0.12, 0.24, -44},
};

int BlastScoreUtils::blosum62_20_prefs[BLOSUM62_20_VALUES_MAX] = {
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_BEST,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL,
        BLAST_MATRIX_NOMINAL
};


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
const double BlastScoreUtils::blastn_values_1_5[][8] = {
        { 0, 0, 1.39, 0.747, 1.38, 1.00,  0, 100 },
        { 3, 3, 1.39, 0.747, 1.38, 1.00,  0, 100 }
};

/** Karlin-Altschul parameter values for substitution scores 1 and -4. */
const double BlastScoreUtils::blastn_values_1_4[][8] = {
        { 0, 0, 1.383, 0.738, 1.36, 1.02,  0, 100 },
        { 1, 2,  1.36,  0.67,  1.2,  1.1,  0,  98 },
        { 0, 2,  1.26,  0.43, 0.90,  1.4, -1,  91 },
        { 2, 1,  1.35,  0.61,  1.1,  1.2, -1,  98 },
        { 1, 1,  1.22,  0.35, 0.72,  1.7, -3,  88 }
};

/** Karlin-Altschul parameter values for substitution scores 2 and -7.
 * These parameters can only be applied to even scores. Any odd score must be
 * rounded down to the nearest even number before calculating the e-value.
 */
const double BlastScoreUtils::blastn_values_2_7[][8] = {
        { 0, 0,  0.69, 0.73, 1.34, 0.515,  0, 100 },
        { 2, 4,  0.68, 0.67,  1.2,  0.55,  0,  99 },
        { 0, 4,  0.63, 0.43, 0.90,   0.7, -1,  91 },
        { 4, 2, 0.675, 0.62,  1.1,   0.6, -1,  98 },
        { 2, 2,  0.61, 0.35, 0.72,   1.7, -3,  88 }
};

/** Karlin-Altschul parameter values for substitution scores 1 and -3. */
const double BlastScoreUtils::blastn_values_1_3[][8] = {
        { 0, 0, 1.374, 0.711, 1.31, 1.05,  0, 100 },
        { 2, 2,  1.37,  0.70,  1.2,  1.1,  0,  99 },
        { 1, 2,  1.35,  0.64,  1.1,  1.2, -1,  98 },
        { 0, 2,  1.25,  0.42, 0.83,  1.5, -2,  91 },
        { 2, 1,  1.34,  0.60,  1.1,  1.2, -1,  97 },
        { 1, 1,  1.21,  0.34, 0.71,  1.7, -2,  88 }
};

/** Karlin-Altschul parameter values for substitution scores 2 and -5.
 * These parameters can only be applied to even scores. Any odd score must be
 * rounded down to the nearest even number before calculating the e-value.
 */
const double BlastScoreUtils::blastn_values_2_5[][8] = {
        { 0, 0, 0.675, 0.65,  1.1,  0.6, -1, 99 },
        { 2, 4,  0.67, 0.59,  1.1,  0.6, -1, 98 },
        { 0, 4,  0.62, 0.39, 0.78,  0.8, -2, 91 },
        { 4, 2,  0.67, 0.61,  1.0, 0.65, -2, 98 },
        { 2, 2,  0.56, 0.32, 0.59, 0.95, -4, 82 }
};

/** Karlin-Altschul parameter values for substitution scores 1 and -2. */
const double BlastScoreUtils::blastn_values_1_2[][8] = {
        { 0, 0, 1.28, 0.46, 0.85, 1.5, -2, 96 },
        { 2, 2, 1.33, 0.62,  1.1, 1.2,  0, 99 },
        { 1, 2, 1.30, 0.52, 0.93, 1.4, -2, 97 },
        { 0, 2, 1.19, 0.34, 0.66, 1.8, -3, 89 },
        { 3, 1, 1.32, 0.57,  1.0, 1.3, -1, 99 },
        { 2, 1, 1.29, 0.49, 0.92, 1.4, -1, 96 },
        { 1, 1, 1.14, 0.26, 0.52, 2.2, -5, 85 }
};

/** Karlin-Altschul parameter values for substitution scores 2 and -3.
 * These parameters can only be applied to even scores. Any odd score must be
 * rounded down to the nearest even number before calculating the e-value.
 */
const double BlastScoreUtils::blastn_values_2_3[][8] = {
        { 0, 0,  0.55, 0.21, 0.46,  1.2, -5, 87 },
        { 4, 4,  0.63, 0.42, 0.84, 0.75, -2, 99 },
        { 2, 4, 0.615, 0.37, 0.72, 0.85, -3, 97 },
        { 0, 4,  0.55, 0.21, 0.46,  1.2, -5, 87 },
        { 3, 3, 0.615, 0.37, 0.68,  0.9, -3, 97 },
        { 6, 2,  0.63, 0.42, 0.84, 0.75, -2, 99 },
        { 5, 2, 0.625, 0.41, 0.78,  0.8, -2, 99 },
        { 4, 2,  0.61, 0.35, 0.68,  0.9, -3, 96 },
        { 2, 2, 0.515, 0.14, 0.33, 1.55, -9, 81 }
};

/** Karlin-Altschul parameter values for substitution scores 3 and -4. */
const double BlastScoreUtils::blastn_values_3_4[][8] = {
        { 6, 3, 0.389, 0.25, 0.56, 0.7, -5, 95},
        { 5, 3, 0.375, 0.21, 0.47, 0.8, -6, 92},
        { 4, 3, 0.351, 0.14, 0.35, 1.0, -9, 86},
        { 6, 2, 0.362, 0.16, 0.45, 0.8, -4, 88},
        { 5, 2, 0.330, 0.092, 0.28, 1.2, -13, 81},
        { 4, 2, 0.281, 0.046, 0.16, 1.8, -23, 69}
};

/** Karlin-Altschul parameter values for substitution scores 4 and -5. */
const double BlastScoreUtils::blastn_values_4_5[][8] = {
        { 0, 0, 0.22, 0.061, 0.22, 1.0, -15, 74 },
        { 6, 5, 0.28,  0.21, 0.47, 0.6 , -7, 93 },
        { 5, 5, 0.27,  0.17, 0.39, 0.7,  -9, 90 },
        { 4, 5, 0.25,  0.10, 0.31, 0.8, -10, 83 },
        { 3, 5, 0.23, 0.065, 0.25, 0.9, -11, 76 }
};

/** Karlin-Altschul parameter values for substitution scores 1 and -1. */
const double BlastScoreUtils::blastn_values_1_1[][8] = {
        { 3,  2, 1.09,  0.31, 0.55, 2.0,  -2, 99 },
        { 2,  2, 1.07,  0.27, 0.49, 2.2,  -3, 97 },
        { 1,  2, 1.02,  0.21, 0.36, 2.8,  -6, 92 },
        { 0,  2, 0.80, 0.064, 0.17, 4.8, -16, 72 },
        { 4,  1, 1.08,  0.28, 0.54, 2.0,  -2, 98 },
        { 3,  1, 1.06,  0.25, 0.46, 2.3,  -4, 96 },
        { 2,  1, 0.99,  0.17, 0.30, 3.3, -10, 90 }
};

/** Karlin-Altschul parameter values for substitution scores 3 and -2. */
const double BlastScoreUtils::blastn_values_3_2[][8] = {
        {  5,  5, 0.208, 0.030, 0.072, 2.9, -47, 77}
};

/** Karlin-Altschul parameter values for substitution scores 5 and -4. */
const double BlastScoreUtils::blastn_values_5_4[][8] = {
        { 10, 6, 0.163, 0.068, 0.16, 1.0, -19, 85 },
        {  8, 6, 0.146, 0.039, 0.11, 1.3, -29, 76 }
};

double BlastScoreUtils::computeKmn(int qlen, double K, double lambda, double alpha, double beta, size_t dbLen,
                                   size_t seqCnt) {



        double logK = log(K);

        size_t lenadj = BlastComputeLengthAdjustment(K, logK, alpha / lambda, beta, qlen, dbLen, seqCnt);
        //fprintf(stdout, "Params: lambda=%6.3g K=%6.3g alpha=%6.3g beta=%6.3g dbLen=%zu seqCnt=%d lenadj=%zu\n",
        //        lambda, K, alpha, beta, dbLen, seqCnt, lenadj);
        size_t m = std::max(1,qlen) - lenadj;
        size_t n = dbLen - seqCnt * lenadj;

        return K * (double)m * (double)n;
}

BlastScoreUtils::BlastStat BlastScoreUtils::getAltschulStatsForMatrix(std::string matrix, long gapOpen, long gapExtend) {
        double (*mat)[8];
        long val;

        if (matrix.compare("blosum45") == 0) {
                mat = blosum45_values;
                val = BLOSUM45_VALUES_MAX;
        }
        else if (matrix.compare("blosum50") == 0) {
                mat = blosum50_values;
                val = BLOSUM50_VALUES_MAX;
        }
        else if (matrix.compare("blosum62") == 0) {
                mat = blosum62_values;
                val = BLOSUM62_VALUES_MAX;
        }
        else if (matrix.compare("blosum62_20") == 0) {
                mat = blosum62_20_values;
                val = BLOSUM62_20_VALUES_MAX;
        }
        else if (matrix.compare("blosum80") == 0) {
                mat = blosum80_values;
                val = BLOSUM80_VALUES_MAX;
        }
        else if (matrix.compare("blosum90") == 0) {
                mat = blosum90_values;
                val = BLOSUM90_VALUES_MAX;
        }
        else {
                Debug(Debug::ERROR) << "Could not find statistics for " << matrix << " gapopen " << gapOpen << " and gapextend " << gapExtend << "\n";
                EXIT(EXIT_FAILURE);
        }
        for (long i = 0; i < val; i++) {
                if ((fabs(mat[i][0] - ((double) gapOpen)) < 0.1) &&
                    (fabs(mat[i][1] - ((double) gapExtend)) < 0.1)) {
                        return BlastStat(mat[i][3], mat[i][4], mat[i][5], mat[i][6], mat[i][7]);
                }
        }
        Debug(Debug::ERROR) << "Could not find statistics for " << matrix << " gapopen " << gapOpen << " and gapextend " << gapExtend << "\n";
        EXIT(EXIT_FAILURE);

        return BlastStat(0.0,0.0,0.0,0.0,0.0);
}


/**
 * Computes the adjustment to the lengths of the query and database sequences
 * that is used to compensate for edge effects when computing evalues.
 *
 * The length adjustment is an integer-valued approximation to the fixed
 * point of the function
 *
 *    f(ell) = beta +
 *               (alpha/lambda) * (log K + log((m - ell)*(n - N ell)))
 *
 * where m is the query length n is the length of the database and N is the
 * number of sequences in the database. The values beta, alpha, lambda and
 * K are statistical, Karlin-Altschul parameters.
 *
 * The value of the length adjustment computed by this routine, A,
 * will always be an integer smaller than the fixed point of
 * f(ell). Usually, it will be the largest such integer.  However, the
 * computed length adjustment, A, will also be so small that
 *
 *    K * (m - A) * (n - N * A) > min(m,n).
 *
 * Moreover, an iterative method is used to compute A, and under
 * unusual circumstances the iterative method may not converge.
 *
 * @param K      the statistical parameter K
 * @param logK   the natural logarithm of K
 * @param alpha_d_lambda    the ratio of the statistical parameters
 *                          alpha and lambda (for ungapped alignments, the
 *                          value 1/H should be used)
 * @param beta              the statistical parameter beta (for ungapped
 *                          alignments, beta == 0)
 * @param query_length      the length of the query sequence
 * @param db_length         the length of the database
 * @param db_num_seq        the number of sequences in the database
 * @param length_adjustment the computed value of the length adjustment [out]
 *
 * @return   0 if length_adjustment is known to be the largest integer less
 *           than the fixed point of f(ell); 1 otherwise.
 */
int
BlastScoreUtils::BlastComputeLengthAdjustment(double K, double logK, double alpha_d_lambda, double beta,
                                              int query_length, size_t db_length, size_t db_num_seqs)
{
        int i;                     /* iteration index */
        int length_adjustment = 0;
        const int maxits = 100;     /* maximum allowed iterations */
#ifdef ORIGINAL_NCBI_CODE
        double m = query_length, n = db_length, N = db_num_seqs;
#else
        double m = query_length, n = (double) db_length, N = db_num_seqs;
#endif

        double ell;            /* A float value of the length adjustment */
        double ss;             /* effective size of the search space */
        double ell_min = 0, ell_max;   /* At each iteration i,
                                         * ell_min <= ell <= ell_max. */
        bool converged    = false;       /* True if the iteration converged */
        double ell_next = 0;   /* Value the variable ell takes at iteration
                                 * i + 1 */
        /* Choose ell_max to be the largest nonnegative value that satisfies
         *
         *    K * (m - ell) * (n - N * ell) > max(m,n)
         *
         * Use quadratic formula: 2 c /( - b + sqrt( b*b - 4 * a * c )) */
        { /* scope of a, mb, and c, the coefficients in the quadratic formula
       * (the variable mb is -b) */
                double a  = N;
                double mb = m * N + n;
                double c  = n * m - std::max(m, n) / K;

                if(c < 0) {
                        length_adjustment = 0;
                        return 1;
                } else {
                        ell_max = 2 * c / (mb + std::sqrt(mb * mb - 4 * a * c));
                }
        } /* end scope of a, mb and c */

        for(i = 1; i <= maxits; i++) {      /* for all iteration indices */
                double ell_bar;    /* proposed next value of ell */
                ell      = ell_next;
                ss       = (m - ell) * (n - N * ell);
                ell_bar  = alpha_d_lambda * (logK + log(ss)) + beta;
                if(ell_bar >= ell) { /* ell is no bigger than the true fixed point */
                        ell_min = ell;
                        if(ell_bar - ell_min <= 1.0) {
                                converged = true;
                                break;
                        }
#ifdef ORIGINAL_NCBI_CODE
                        if(ell_min == ell_max) { /* There are no more points to check */
#else
                        if(ell_min >= ell_max) { /* There are no more points to check */
#endif
                                break;
                        }
                } else { /* else ell is greater than the true fixed point */
                        ell_max = ell;
                }
                if(ell_min <= ell_bar && ell_bar <= ell_max) {
                        /* ell_bar is in range. Accept it */
                        ell_next = ell_bar;
                } else { /* else ell_bar is not in range. Reject it */
                        ell_next = (i == 1) ? ell_max : (ell_min + ell_max) / 2;
                }
        } /* end for all iteration indices */
        if(converged) { /* the iteration converged */
                /* If ell_fixed is the (unknown) true fixed point, then we
                 * wish to set (*length_adjustment) to floor(ell_fixed).  We
                 * assume that floor(ell_min) = floor(ell_fixed) */
                length_adjustment = (int) ell_min;
                /* But verify that ceil(ell_min) != floor(ell_fixed) */
                ell = ceil(ell_min);
                if( ell <= ell_max ) {
                        ss = (m - ell) * (n - N * ell);
                        if(alpha_d_lambda * (logK + log(ss)) + beta >= ell) {
                                /* ceil(ell_min) == floor(ell_fixed) */
                                length_adjustment = (int) ell;
                        }
                }
        } else { /* else the iteration did not converge. */
                /* Use the best value seen so far */
                length_adjustment = (int) ell_min;
        }
        if(converged == 0){
//                Debug(Debug::WARNING) << "BlastComputeLengthAdjustment did not converge\n";
                ;
        }

        return length_adjustment;
}
