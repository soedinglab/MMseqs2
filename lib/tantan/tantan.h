// Author: Martin C. Frith 2010
// SPDX-License-Identifier: MPL-2.0

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

typedef unsigned char uchar;
typedef const double *const_double_ptr;

void maskSequences(uchar *seqBeg,
                   uchar *seqEnd,
                   int maxRepeatOffset,
                   const const_double_ptr *likelihoodRatioMatrix,
                   double repeatProb,
                   double repeatEndProb,
                   double repeatOffsetProbDecay,
                   double firstGapProb,
                   double otherGapProb,
                   double minMaskProb,
                   const uchar *maskTable);

// The following routine gets the posterior probability that each
// letter is repetitive.  It stores the results in "probabilities",
// which must point to enough pre-allocated space to fit the results.

void getProbabilities(const uchar *seqBeg,
                      const uchar *seqEnd,
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

void maskProbableLetters(uchar *seqBeg,
                         uchar *seqEnd,
                         const float *probabilities,
                         double minMaskProb,
                         const uchar *maskTable);

// The following routine counts the expected number of transitions
// from the background (non-repeat) state to other states.  It adds
// the results to "transitionCounts", which must point to
// pre-initialized space for (maxRepeatOffset+1) items.  The
// background->background transition count is stored in
// transitionCounts[0].  The background->(period-i repeat) transition
// count is stored in transitionCounts[i].

// (In this routine, the HMM begin and end states are counted as
// background states.  Thus, begin->X is added to background->X, and
// X->end is added to X->background.)

void countTransitions(const uchar *seqBeg,
                      const uchar *seqEnd,
                      int maxRepeatOffset,
                      const const_double_ptr *likelihoodRatioMatrix,
                      double repeatProb,
                      double repeatEndProb,
                      double repeatOffsetProbDecay,
                      double firstGapProb,
                      double otherGapProb,
                      double *transitionCounts);

}

#endif
