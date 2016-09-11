#ifndef MATH_H
#define MATH_H

#include <cmath>
#include <climits>
#include <cfloat>
#include <vector>
#include <limits>

class MathUtil {
public:
    static inline int ipow(int base, int exponent) {
        int res = 1;
        for (int i = 0; i < exponent; i++)
            res = res * base;
        return res;
    }

    static inline short sadd16_signed(short x, short y) {
        unsigned short ux = x;
        unsigned short uy = y;
        unsigned short res = ux + uy;

        /* Calculate overflowed result. (Don't change the sign bit of ux) */
        ux = (ux >> 15) + SHRT_MAX;

        /* Force compiler to use cmovns instruction */
        if ((short) ((ux ^ uy) | ~(uy ^ res)) >= 0) {
            res = ux;
        }

        return res;
    }

    static inline short ssub16_signed(short x, short y) {
        unsigned short ux = x;
        unsigned short uy = y;
        unsigned short res = ux - uy;

        ux = (ux >> 15) + SHRT_MAX;

        /* Force compiler to use cmovns instruction */
        if ((short) ((ux ^ uy) & (ux ^ res)) < 0) {
            res = ux;
        }

        return res;
    }

    // doesnt produce NANs careful
    static inline float flog2(float x) {
        if (x <= 0)
            return -128;
        int *px = (int *) (&x);        // store address of float as pointer to long int
        float e = (float) (((*px & 0x7F800000) >> 23) -
                           0x7f); // shift right by 23 bits and subtract 127 = 0x7f => exponent
        *px = ((*px & 0x007FFFFF) | 0x3f800000);  // set exponent to 127 (i.e., 0)
        x -= 1.0;         // and calculate x-1.0
        x *= (1.441740
              + x * (-0.7077702 +
                     x * (0.4123442 + x * (-0.1903190 + x * 0.0440047)))); // 5'th order polynomial approx. of log(1+x)
        return x + e;
    }

    static inline double fpow2(float x) {
        if (x >= FLT_MAX_EXP) return FLT_MAX;
        if (x <= FLT_MIN_EXP) return 0.0f;

        int *px = (int *) (&x);        // store address of float as pointer to long int
        float tx = (x - 0.5f) + (3 << 22); // temporary value for truncation: x-0.5 is added to a large integer (3<<22),
        // 3<<22 = (1.1bin)*2^23 = (1.1bin)*2^(150-127),
        // which, in internal bits, is written 0x4b400000 (since 10010110bin = 150)
        int lx = *((int *) &tx) - 0x4b400000;   // integer value of x
        float dx = x - (float) (lx);             // float remainder of x
//   x = 1.0f + dx*(0.69606564f           // cubic apporoximation of 2^x for x in the range [0, 1]
//            + dx*(0.22449433f           // Gives relative deviation < 1.5E-4
//            + dx*(0.07944023f)));       // Speed: 1.9E-8s
        x = 1.0f + dx * (0.693019f // polynomial approximation of 2^x for x in the range [0, 1]
                         + dx * (0.241404f             // Gives relative deviation < 4.6E-6
                                 + dx * (0.0520749f            // Speed: 2.1E-8s
                                         + dx * 0.0134929f)));
//   x = 1.0f + dx*(0.693153f             // polynomial apporoximation of 2^x for x in the range [0, 1]
//            + dx*(0.240153f             // Gives relative deviation < 2.3E-7
//            + dx*(0.0558282f            // Speed: 2.3E-8s
//            + dx*(0.00898898f
//            + dx* 0.00187682f ))));
        *px += (lx << 23);                      // add integer power of 2 to exponent
        return x;
    }

    static inline unsigned int concatenate(unsigned int x, unsigned int y) {
        unsigned int pow = 10;
        while (y >= pow)
            pow *= 10;
        return x * pow + y;
    }

    static inline double log2(double x) {
        return log10(x) / 0.301029996;
    }

    static inline unsigned short sadd16(const unsigned short a, const unsigned short b) {
        return (a > 0xFFFF - b) ? 0xFFFF : a + b;
    }

    // Compute the sum of bits of one or two integers
    static inline int popCount(int i) {
        i = i - ((i >> 1) & 0x55555555);
        i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
        return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
    }

    static inline float getCoverage(size_t start, size_t end, size_t length) {
        return static_cast<float>(end - start + 1) / static_cast<float>(length);
    }

// compute look up table based on stirling approximation
    static void computeFactorial(double *output, const size_t range) {
        output[0] = log(1.0);
        for(size_t score = 1; score < range; score++){
            const double scoreDbl = static_cast<double>(score);

            const double S_fact = std::min(std::numeric_limits<double>::max(),
                                           sqrt(2 * M_PI * scoreDbl) * pow(scoreDbl / exp(1), scoreDbl) * exp(1 / (12  * scoreDbl)));
            output[score] = log(S_fact);
        }
    }
};

#endif
