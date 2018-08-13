// Copyright 2010 Martin C. Frith

#include "tantan.h"

#include <algorithm>  // fill, max
#include <cassert>
#include <cmath>  // pow, abs
#include <iostream>  // cerr
#include <numeric>  // accumulate
#include <vector>

#define BEG(v) ((v).empty() ? 0 : &(v).front())
#define END(v) ((v).empty() ? 0 : &(v).back() + 1)

namespace tantan {

    void multiplyAll(std::vector<double> &v, double factor) {
        for (std::vector<double>::iterator i = v.begin(); i < v.end(); ++i)
            *i *= factor;
    }

    double firstRepeatOffsetProb(double probMult, int maxRepeatOffset) {
        if (probMult < 1 || probMult > 1)
            return (1 - probMult) / (1 - std::pow(probMult, maxRepeatOffset));
        else
            return 1.0 / maxRepeatOffset;
    }

    void checkForwardAndBackwardTotals(double fTot, double bTot) {
        double x = std::abs(fTot);
        double y = std::abs(bTot);

        // ??? Is 1e6 suitable here ???
        if (std::abs(fTot - bTot) > std::max(x, y) / 1e6)
            std::cerr << "tantan: warning: possible numeric inaccuracy\n"
                      << "tantan:          forward algorithm total: " << fTot << "\n"
                      << "tantan:          backward algorithm total: " << bTot << "\n";
    }

    struct Tantan {
        enum { scaleStepSize = 16 };

        const char *seqBeg;  // start of the sequence
        const char *seqEnd;  // end of the sequence
        const char *seqPtr;  // current position in the sequence

        int maxRepeatOffset;

        const const_double_ptr *likelihoodRatioMatrix;

        double b2b;  // transition probability from background to background
        double f2b;  // transition probability from foreground to background
        double g2g;  // transition probability from gap/indel to gap/indel
        //double f2g;  // transition probability from foreground to gap/indel
        //double g2f;  // transition probability from gap/indel to foreground
        double oneGapProb;  // f2g * g2f
        double endGapProb;  // f2g * 1
        double f2f0;  // foreground to foreground, if there are 0 indel transitions
        double f2f1;  // foreground to foreground, if there is 1 indel transition
        double f2f2;  // foreground to foreground, if there are 2 indel transitions
        double b2fDecay;
        double b2fGrowth;
        double b2fFirst;  // background state to first foreground state
        double b2fLast;  // background state to last foreground state

        double backgroundProb;
        std::vector<double> foregroundProbs;
        std::vector<double> insertionProbs;

        std::vector<double> scaleFactors;

        Tantan(const char *seqBeg,
               const char *seqEnd,
               int maxRepeatOffset,
               const const_double_ptr *likelihoodRatioMatrix,
               double repeatProb,
               double repeatEndProb,
               double repeatOffsetProbDecay,
               double firstGapProb,
               double otherGapProb) {
            assert(maxRepeatOffset > 0);
            assert(repeatProb >= 0 && repeatProb < 1);
            // (if repeatProb==1, then any sequence is impossible)
            assert(repeatEndProb >= 0 && repeatEndProb <= 1);
            assert(repeatOffsetProbDecay > 0 && repeatOffsetProbDecay <= 1);
            assert(otherGapProb >= 0 && otherGapProb <= 1);
            assert(firstGapProb >= 0);
            assert(repeatEndProb + firstGapProb * 2 <= 1);

            this->seqBeg = seqBeg;
            this->seqEnd = seqEnd;
            this->seqPtr = seqBeg;
            this->maxRepeatOffset = maxRepeatOffset;
            this->likelihoodRatioMatrix = likelihoodRatioMatrix;

            b2b = 1 - repeatProb;
            f2b = repeatEndProb;
            g2g = otherGapProb;
            //f2g = firstGapProb;
            //g2f = 1 - otherGapProb;
            oneGapProb = firstGapProb * (1 - otherGapProb);
            endGapProb = firstGapProb * 1;
            f2f0 = 1 - repeatEndProb;
            f2f1 = 1 - repeatEndProb - firstGapProb;
            f2f2 = 1 - repeatEndProb - firstGapProb * 2;

            b2fDecay = repeatOffsetProbDecay;
            b2fGrowth = 1 / repeatOffsetProbDecay;

            b2fFirst = repeatProb * firstRepeatOffsetProb(b2fDecay, maxRepeatOffset);
            b2fLast = repeatProb * firstRepeatOffsetProb(b2fGrowth, maxRepeatOffset);

            foregroundProbs.resize(maxRepeatOffset);
            insertionProbs.resize(maxRepeatOffset - 1);

            scaleFactors.resize((seqEnd - seqBeg) / scaleStepSize);
        }

        void initializeForwardAlgorithm() {
            backgroundProb = 1.0;
            std::fill(foregroundProbs.begin(), foregroundProbs.end(), 0.0);
            std::fill(insertionProbs.begin(), insertionProbs.end(), 0.0);
        }

        double forwardTotal() {
            double fromForeground = std::accumulate(foregroundProbs.begin(),
                                                    foregroundProbs.end(), 0.0);
            fromForeground *= f2b;
            double total = backgroundProb * b2b + fromForeground;
            assert(total > 0);
            return total;
        }

        void initializeBackwardAlgorithm() {
            backgroundProb = b2b;
            std::fill(foregroundProbs.begin(), foregroundProbs.end(), f2b);
            std::fill(insertionProbs.begin(), insertionProbs.end(), 0.0);
        }

        double backwardTotal() {
            assert(backgroundProb > 0);
            return backgroundProb;
        }

        void calcForwardTransitionProbsWithGaps() {
            double fromBackground = backgroundProb * b2fLast;
            double *foregroundPtr = &foregroundProbs.back();
            double f = *foregroundPtr;
            double fromForeground = f;

            if (insertionProbs.empty()) {
                *foregroundPtr = fromBackground + f * f2f0;
            } else {
                double *insertionPtr = &insertionProbs.back();
                double i = *insertionPtr;
                *foregroundPtr = fromBackground + f * f2f1 + i * endGapProb;
                double d = f;
                --foregroundPtr;
                fromBackground *= b2fGrowth;

                while (foregroundPtr > &foregroundProbs.front()) {
                    f = *foregroundPtr;
                    fromForeground += f;
                    i = *(insertionPtr - 1);
                    *foregroundPtr = fromBackground + f * f2f2 + (i + d) * oneGapProb;
                    *insertionPtr = f + i * g2g;
                    d = f + d * g2g;
                    --foregroundPtr;
                    --insertionPtr;
                    fromBackground *= b2fGrowth;
                }

                f = *foregroundPtr;
                fromForeground += f;
                *foregroundPtr = fromBackground + f * f2f1 + d * endGapProb;
                *insertionPtr = f;
            }

            fromForeground *= f2b;
            backgroundProb = backgroundProb * b2b + fromForeground;
        }

        void calcBackwardTransitionProbsWithGaps() {
            double toBackground = f2b * backgroundProb;
            double *foregroundPtr = &foregroundProbs.front();
            double f = *foregroundPtr;
            double toForeground = f;

            if (insertionProbs.empty()) {
                *foregroundPtr = toBackground + f2f0 * f;
            } else {
                double *insertionPtr = &insertionProbs.front();
                double i = *insertionPtr;
                *foregroundPtr = toBackground + f2f1 * f + i;
                double d = endGapProb * f;
                ++foregroundPtr;
                toForeground *= b2fGrowth;

                while (foregroundPtr < &foregroundProbs.back()) {
                    f = *foregroundPtr;
                    toForeground += f;
                    i = *(insertionPtr + 1);
                    *foregroundPtr = toBackground + f2f2 * f + (i + d);
                    double oneGapProb_f = oneGapProb * f;
                    *insertionPtr = oneGapProb_f + g2g * i;
                    d = oneGapProb_f + g2g * d;
                    ++foregroundPtr;
                    ++insertionPtr;
                    toForeground *= b2fGrowth;
                }

                f = *foregroundPtr;
                toForeground += f;
                *foregroundPtr = toBackground + f2f1 * f + d;
                *insertionPtr = endGapProb * f;
            }

            toForeground *= b2fLast;
            backgroundProb = b2b * backgroundProb + toForeground;
        }

        void calcForwardTransitionProbs() {
            if (endGapProb > 0) return calcForwardTransitionProbsWithGaps();

            double fromBackground = backgroundProb * b2fLast;
            double fromForeground = 0;
            double *foregroundPtr = END(foregroundProbs);
            double *foregroundBeg = BEG(foregroundProbs);

            while (foregroundPtr > foregroundBeg) {
                --foregroundPtr;
                double f = *foregroundPtr;
                fromForeground += f;
                *foregroundPtr = fromBackground + f * f2f0;
                fromBackground *= b2fGrowth;
            }

            fromForeground *= f2b;
            backgroundProb = backgroundProb * b2b + fromForeground;
        }

        void calcBackwardTransitionProbs() {
            if (endGapProb > 0) return calcBackwardTransitionProbsWithGaps();

            double toBackground = f2b * backgroundProb;
            double toForeground = 0;
            double *foregroundPtr = BEG(foregroundProbs);
            double *foregroundEnd = END(foregroundProbs);

            while (foregroundPtr < foregroundEnd) {
                toForeground *= b2fGrowth;
                double f = *foregroundPtr;
                toForeground += f;
                *foregroundPtr = toBackground + f2f0 * f;
                ++foregroundPtr;
            }

            toForeground *= b2fLast;
            backgroundProb = b2b * backgroundProb + toForeground;
        }

        void calcEmissionProbs() {
            const double *lrRow = likelihoodRatioMatrix[static_cast<int>(*seqPtr)];

            bool isNearSeqBeg = (seqPtr - seqBeg < maxRepeatOffset);
            const char *seqStop = isNearSeqBeg ? seqBeg : seqPtr - maxRepeatOffset;

            double *foregroundPtr = BEG(foregroundProbs);
            const char *offsetPtr = seqPtr;

            while (offsetPtr > seqStop) {
                --offsetPtr;
                *foregroundPtr *= lrRow[static_cast<int>(*offsetPtr)];
                ++foregroundPtr;
            }

            while (foregroundPtr < END(foregroundProbs)) {
                *foregroundPtr *= 0;
                ++foregroundPtr;
            }
        }

        void rescale(double scale) {
            backgroundProb *= scale;
            multiplyAll(foregroundProbs, scale);
            multiplyAll(insertionProbs, scale);
        }

        void rescaleForward() {
            if ((seqPtr - seqBeg) % scaleStepSize == scaleStepSize - 1) {
                assert(backgroundProb > 0);
                double scale = 1 / backgroundProb;
                scaleFactors[(seqPtr - seqBeg) / scaleStepSize] = scale;
                rescale(scale);
            }
        }

        void rescaleBackward() {
            if ((seqPtr - seqBeg) % scaleStepSize == scaleStepSize - 1) {
                double scale = scaleFactors[(seqPtr - seqBeg) / scaleStepSize];
                rescale(scale);
            }
        }

        void calcRepeatProbs(float *letterProbs) {
            initializeForwardAlgorithm();

            while (seqPtr < seqEnd) {
                calcForwardTransitionProbs();
                calcEmissionProbs();
                rescaleForward();
                *letterProbs = static_cast<float>(backgroundProb);
                ++letterProbs;
                ++seqPtr;
            }

            double z = forwardTotal();

            initializeBackwardAlgorithm();

            while (seqPtr > seqBeg) {
                --seqPtr;
                --letterProbs;
                double nonRepeatProb = *letterProbs * backgroundProb / z;
                // Convert nonRepeatProb to a float, so that it is more likely
                // to be exactly 1 when it should be, e.g. for the 1st letter of
                // a sequence:
                *letterProbs = 1 - static_cast<float>(nonRepeatProb);
                rescaleBackward();
                calcEmissionProbs();
                calcBackwardTransitionProbs();
            }

            double z2 = backwardTotal();
            checkForwardAndBackwardTotals(z, z2);
        }

    };

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
                       const char *maskTable) {
        std::vector<float> p(seqEnd - seqBeg);
        float *probabilities = BEG(p);

        getProbabilities(seqBeg, seqEnd, maxRepeatOffset,
                         likelihoodRatioMatrix, repeatProb, repeatEndProb,
                         repeatOffsetProbDecay, firstGapProb, otherGapProb,
                         probabilities);

        return maskProbableLetters(seqBeg, seqEnd, probabilities, minMaskProb, maskTable);
    }

    void getProbabilities(const char *seqBeg,
                          const char *seqEnd,
                          int maxRepeatOffset,
                          const const_double_ptr *likelihoodRatioMatrix,
                          double repeatProb,
                          double repeatEndProb,
                          double repeatOffsetProbDecay,
                          double firstGapProb,
                          double otherGapProb,
                          float *probabilities) {
        Tantan tantan(seqBeg, seqEnd, maxRepeatOffset, likelihoodRatioMatrix,
                      repeatProb, repeatEndProb, repeatOffsetProbDecay,
                      firstGapProb, otherGapProb);
        tantan.calcRepeatProbs(probabilities);
    }

    int maskProbableLetters(char *seqBeg,
                             char *seqEnd,
                             const float *probabilities,
                             double minMaskProb,
                             const char *maskTable) {
        int masked = 0;
        while (seqBeg < seqEnd) {
            if (*probabilities >= minMaskProb){
                masked++;
                *seqBeg = maskTable[static_cast<int>(*seqBeg)];
            }
            ++probabilities;
            ++seqBeg;
        }
        return masked;
    }

}
