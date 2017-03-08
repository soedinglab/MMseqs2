// Copyright 2010 Martin C. Frith
// tantan is distributed under the GNU General Public License, either
//        version 3 of the License, or (at your option) any later version.  For
//        details, see COPYING.txt.
//
// If you use tantan in your research, please cite:
// "A new repeat-masking method enables specific detection of homologous
// sequences", MC Frith, Nucleic Acids Research 2011 39(4):e23.
//
// tantan's website is: http://www.cbrc.jp/tantan/
//
// If you have any questions, comments, or problems concerning tantan,
// please email: tantan (ATmark) cbrc (dot) jp.
// These are routines for masking simple regions (low-complexity and
// short-period tandem repeats) in biological sequences.  To
// understand them in detail, see the published article (in
// preparation).

// Typically, you would just use the maskSequences routine.  The other
// routines are more specialized.  The inputs to maskSequences are as
// follows.

// seqBeg: pointer to the start of the sequence.
// seqEnd: pointer to one-past-the-end of the sequence.
// maxRepeatOffset: the maximum tandem-repeat period-size to consider.

// likelihoodRatioMatrix: a matrix of the form Qxy / (Px * Py), where
// Qxy is the probability of seeing letters x and y in equivalent
// positions of a tandem repeat, and Px and Py are the background
// probabilities of the letters.  This matrix is related to a scoring
// matrix (e.g. blosum62) by:
// likelihoodRatioMatrix[x][y] = exp(lambda * scoringMatrix[x][y]).

// The uchars in the sequence will be used as indexes into the matrix
// (e.g. likelihoodRatioMatrix[x][y]).  So typically the uchars are
// small integers (e.g. 0, 1, 2, 3).  The matrix needs to have entries
// for any uchars that can occur.

// repeatProb: the probability of a repetitive segment starting per position.
// repeatEndProb: the probability of a repetitive segment ending per position.

// repeatOffsetProbDecay: the probability of a period-(i) repeat
// divided by the probability of a period-(i-1) repeat.

// firstGapProb: the probability of initiating an insertion or
// deletion in a repetitive region.

// otherGapProb: the probability of extending an insertion or deletion
// by one more letter.

// minMaskProb: mask letters whose posterior probability of being
// repetitive is >= this.

// maskTable: how to do the masking.  Letter x will be changed to
// maskTable[x].  So maskTable needs to have entries for any uchar
// that can occur.

// Typical usage:
// tantan::maskSequences(seqBeg, seqEnd, 100, likelihoodRatioMatrix,
// 0.005, 0.05, 0.9, 0, 0, 0.5, maskTable)

#ifndef TANTAN_HH
#define TANTAN_HH

namespace tantan {

typedef const double *const_double_ptr;

int maskSequences(char *seqBeg,
                   char *seqEnd,
                   int maxRepeatOffset,
                   const const_double_ptr *likelihoodRatioMatrix,
                   double repeatProb,
                   double repeatEndProb,
                   double repeatOffsetProbDecay,
                   double firstGapProb,
                   double otherGapProb,
                   double minMaskProb,
                   const char *maskTable);

// The following routine gets the posterior probability that each
// letter is repetitive.  It stores the results in "probabilities",
// which must point to enough pre-allocated space to fit the results.

void getProbabilities(const char *seqBeg,
                      const char *seqEnd,
                      int maxRepeatOffset,
                      const const_double_ptr *likelihoodRatioMatrix,
                      double repeatProb,
                      double repeatEndProb,
                      double repeatOffsetProbDecay,
                      double firstGapProb,
                      double otherGapProb,
                      float *probabilities);

// The following routine masks each letter whose corresponding entry
// in "probabilities" is >= minMaskProb.

int maskProbableLetters(char *seqBeg,
                         char *seqEnd,
                         const float *probabilities,
                         double minMaskProb,
                         const char *maskTable);

}

#endif
