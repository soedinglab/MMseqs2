#ifndef MATH_H
#define MATH_H

#include <cmath>
#include <cstddef>
#include <climits>
#include <cfloat>
#include <vector>
#include <limits>

// Not provided by Cygwin
#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif

#if defined(__has_attribute)
#  define HAS_ATTRIBUTE(x) __has_attribute(x)
#else
#  define HAS_ATTRIBUTE(x) (0)
#endif


#ifndef MAY_ALIAS
#if HAS_ATTRIBUTE(__may_alias__)
#  define MAY_ALIAS(x) x __attribute__((__may_alias__))
#else
#  define MAY_ALIAS(x) x
#endif
#endif

class MathUtil {
public:

    template<typename T>
    static T ceilIntDivision(T x, T y){
        return 1 + ((x - 1) / y);
    }

    static bool isWithinModRang(unsigned int start, unsigned int end, unsigned int n, unsigned int mod){
        unsigned int startRange = start % mod;  // % = mod
        unsigned int endRange = end % mod;
        unsigned int number = n % mod;

        return (number - startRange) % mod <=  (endRange - startRange) % mod;

    }

    template<typename T>
    static inline T ipow(int base, int exponent) {
        T res = 1;
        for (int i = 0; i < exponent; i++)
            res = res * base;
        return res;
    }

    static inline unsigned int log10base(unsigned int v){
        return (v >= 1000000000u) ? 1000000000u : (v >= 100000000u) ? 100000000u :
                                        (v >= 10000000u) ? 10000000u : (v >= 1000000u) ? 1000000u:
                                                               (v >= 100000u) ?  100000u: (v >= 10000u) ? 10000u:
                                                                                    (v >= 1000u) ? 1000u: (v >= 100u) ? 100u: (v >= 10u) ? 10u: 1u;
    }

    static inline unsigned int log10(unsigned int v){
        return (v >= 1000000000u) ? 9 : (v >= 100000000u) ? 8 :
                                        (v >= 10000000u) ? 7 : (v >= 1000000u) ? 6 :
                                                               (v >= 100000u) ? 5 : (v >= 10000u) ? 4 :
                                                                                    (v >= 1000u) ? 3 : (v >= 100u) ? 2 : (v >= 10u) ? 1u : 0u;
    }

    static bool AreSame(float a, float b)
    {
        return fabs(a - b) < std::numeric_limits<float>::epsilon();
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
        MAY_ALIAS(int) *px = (int *) (&x);        // store address of float as pointer to long int
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

        MAY_ALIAS(int) *px = (int *) (&x);        // store address of float as pointer to long int
        float tx = (x - 0.5f) + (3 << 22); // temporary value for truncation: x-0.5 is added to a large integer (3<<22),
        // 3<<22 = (1.1bin)*2^23 = (1.1bin)*2^(150-127),
        // which, in internal bits, is written 0x4b400000 (since 10010110bin = 150)
        MAY_ALIAS(int) *ix = (int *) (&tx);   // integer value of x
        int lx = *ix - 0x4b400000;
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

    static inline float getCoverage(size_t start, size_t end, size_t length) {
        return static_cast<float>(end - start + 1) / static_cast<float>(length);
    }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wconstant-conversion"
#pragma GCC diagnostic ignored "-Woverflow"
    /** A single gain expressed as minifloat */
#define EXPONENT_BITS   3
#define MANTISSA_BITS   5
#define EXPONENT_MAX    ((1 << EXPONENT_BITS) - 1)
#define EXCESS          ((1 << EXPONENT_BITS) - 2)
#define MANTISSA_MAX    ((1 << MANTISSA_BITS) - 1)
#define HIDDEN_BIT      (1 << MANTISSA_BITS)
#define ONE_FLOAT       ((float) (1 << (MANTISSA_BITS + 1)))
#define MINIFLOAT_MAX   ((EXPONENT_MAX << MANTISSA_BITS) | MANTISSA_MAX)
#define MINIFLOAT_MIN   1
#if EXPONENT_BITS + MANTISSA_BITS != 8
#error EXPONENT_BITS and MANTISSA_BITS must sum to 16
#endif

    static char convertFloatToChar(float v)
    {
        if (std::isnan(v) || v <= 0.0f) {
            return 0;
        }
        if (v >= 2.0f) {
            return MINIFLOAT_MAX;
        }
        int exp;
        float r = frexpf(v, &exp);
        if ((exp += EXCESS) > EXPONENT_MAX) {
            return MINIFLOAT_MAX;
        }
        if (-exp >= MANTISSA_BITS) {
            return 0;
        }
        int mantissa = (int) (r * ONE_FLOAT);
        return exp > 0 ? (exp << MANTISSA_BITS) | (mantissa & ~HIDDEN_BIT) :
               (mantissa >> (1 - exp)) & MANTISSA_MAX;
    }

    static float convertCharToFloat(char a)
    {
        int mantissa = a & MANTISSA_MAX;
        int exponent = (a >> MANTISSA_BITS) & EXPONENT_MAX;
        return ldexpf((exponent > 0 ? HIDDEN_BIT | mantissa : mantissa << 1) / ONE_FLOAT,
                      exponent - EXCESS);
    }
#pragma GCC diagnostic pop


    static char convertNeffToChar(const float neff) {
        float retVal = std::min(255.0f, 1.0f+64.0f*flog2(neff) );
        return std::max(static_cast<unsigned char>(1), static_cast<unsigned char>(retVal + 0.5) );
    }

    static float convertNeffToFloat(unsigned char neffToScale) {
        float retNeff = fpow2((static_cast<float>(neffToScale)-1.0f)/64.0f);;
        return retNeff;
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


    // Normalize a float array such that it sums to one
    // If it sums to 0 then assign def_array elements to array (optional)
    static inline float NormalizeTo1(float* array, int length, const double* def_array = NULL) {
        float sum = 0.0f;
        for (int k = 0; k < length; k++){
            sum += array[k];
        }
        if (sum != 0.0f) {
            float fac = 1.0 / sum;
            for (int i = 0; i < length; i++) {
                array[i] *= fac;
            }
        } else if (def_array) {
            for (int i = 0; i < length; i++) {
                array[i] = def_array[i];
            }
        }
        return sum;
    }


};

#endif
