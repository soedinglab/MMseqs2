#ifndef UTIL_H
#define UTIL_H

#include <cstddef>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <limits.h>

#include <iostream>
#include <float.h>

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

#ifdef OPENMP
#include <omp.h>
#endif

#include "MMseqsMPI.h"

#ifndef EXIT
#define EXIT(exitCode) exit(exitCode)
#endif

#define CHECK_BIT(var,pos) ((var) & (1<<(pos)))
#define SSTR( x ) dynamic_cast< std::ostringstream& >( \
( std::ostringstream().flush() << std::dec << x ) ).str()

#define ARRAY_SIZE(a) (sizeof(a) / sizeof(a[0]))

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif



class Util {
public:
    static void decompose_domain(size_t domain_size, size_t world_rank,
            size_t world_size, size_t *subdomain_start,
            size_t *subdomain_size);
    static void rankedDescSort20(short *val, unsigned int *index);
    static void decomposeDomainByAminoaAcid(size_t aaSize, unsigned int *seqSizes, size_t count,
            size_t worldRank, size_t worldSize, size_t *start, size_t *end);

    static size_t count_lines(const char * file, size_t endPos );

        static inline int fast_atoi( const char * str )
    {
        int val = 0;
        while (*str >= '0' && *str <= '9') {
            val = val*10 + (*str++ - '0');
        }
        return val;
    }
    
    // this is needed because with GCC4.7 omp_get_num_threads() returns just 1.
    static int omp_thread_count() {
        int n = 0;
        #pragma omp parallel reduction(+:n)
            n += 1;
        return n;
    }

    static int ipow (int base, int exponent){
        int res = 1;
        for (int i = 0; i < exponent; i++)
            res = res*base;
        return res;
    }

    static bool startWith(std::string prefix, std::string str){
        return (!str.compare(0, prefix.size(), prefix));
    }
    
    static std::vector<std::string> split(std::string str,std::string sep){
        char buffer[1024];
        snprintf(buffer, 1024, "%s", str.c_str());
        char* cstr = (char *) &buffer;
        char* current;
        std::vector<std::string> arr;
        current=strtok(cstr,sep.c_str());
        while(current!=NULL){
            arr.push_back(current);
            current=strtok(NULL,sep.c_str());
        }
        return arr;
    }

    
    
    static inline char * skipLine(char * data){
        while( *data !='\n' ) { data++; }
        return (data+1);
    }
    
    static inline size_t skipWhitespace(char * data){
        size_t counter = 0;
        while( (data[counter] == ' ' || data[counter] == '\t') == true ) {
            counter++;
        }
        return counter;
    }
    
    static inline size_t skipNoneWhitespace(char * data){
        //A value different from zero (i.e., true) if indeed c is a white-space character. Zero (i.e., false) otherwise.
        size_t counter = 0;
        while(( data[counter] == ' '  || data[counter] == '\t'
             || data[counter] == '\n' || data[counter] == '\0' ) == false ) {
            counter++;
        }
        return counter;
    }
    
    
    static ffindex_index_t* openIndex(const char* indexFileName){
        // count the number of entries in the clustering
        char line [1000];
        int cnt = 0;
        std::ifstream index_file(indexFileName);
        if (index_file.is_open()) {
            while ( index_file.getline (line, 1000) ){
                cnt++;
            }
            index_file.close();
        }
        else{
            std::cerr << "Could not open ffindex index file " << indexFileName << "\n";
            EXIT(EXIT_FAILURE);
        }
        // open clustering ffindex
        FILE* indexFile = fopen(indexFileName, "r");
        if( indexFile == NULL) { fferror_print(__FILE__, __LINE__, "DBReader", indexFileName);  EXIT(EXIT_FAILURE); }
        
        ffindex_index_t* index = ffindex_index_parse(indexFile, cnt);
        return index;
    }
    
    
    static inline size_t getWordsOfLine(char * data, char ** words, size_t maxElement ){
        size_t elementCounter = 0;
        while(*data !=  '\n' && *data != '\0'){
            data += skipWhitespace(data);
            words[elementCounter] = data;
            elementCounter++;
            if(elementCounter >= maxElement)
                return elementCounter;
            data += skipNoneWhitespace(data);
        }
        return elementCounter;
    }
    
    
    static inline unsigned short sadd16(const unsigned short  a, const unsigned short  b)
    { return (a > 0xFFFF - b) ? 0xFFFF : a + b; }
    
    static inline short sadd16_signed(short x, short y)
    {
        unsigned short ux = x;
        unsigned short uy = y;
        unsigned short res = ux + uy;
        
        /* Calculate overflowed result. (Don't change the sign bit of ux) */
        ux = (ux >> 15) + SHRT_MAX;
        
        /* Force compiler to use cmovns instruction */
        if ((short) ((ux ^ uy) | ~(uy ^ res)) >= 0)
        {
            res = ux;
        }
        
        return res;
    }
    
    static inline short ssub16_signed (short x, short y)
    {
        unsigned short ux = x;
        unsigned short uy = y;
        unsigned short res = ux - uy;
        
        ux = (ux >> 15) + SHRT_MAX;
        
        /* Force compiler to use cmovns instruction */
        if ((short)((ux ^ uy) & (ux ^ res)) < 0)
        {
            res = ux;
        }
        
        return res;
    }

    static std::string parseFastaHeader(std::string header);

    static inline char toUpper(char character){
        character += ('a' <= character && character <= 'z') ? ('A' - 'a') : 0;
        return character;
    }

    static inline void parseKey(char *data, char * key) {
        char * startPosOfKey = data;
        char * endPosOfId    = data + Util::skipNoneWhitespace(data);
        ptrdiff_t keySize =  (endPosOfId - startPosOfKey);
        strncpy(key, data, keySize);
        key[keySize] = '\0';
    }

    static FILE* openFileOrDie(const char * fileName, const char * mode);

	static inline void parseByColumnNumber(char *data, char * key, int position) {
        char * startPosOfKey = data;
        for (int i = 0; i < position; ++i) {
            startPosOfKey = startPosOfKey + Util::skipNoneWhitespace(startPosOfKey);
            startPosOfKey = startPosOfKey + Util::skipWhitespace(startPosOfKey);

        }

        char * endPosOfId    = startPosOfKey + Util::skipNoneWhitespace(startPosOfKey);
        ptrdiff_t keySize =  (endPosOfId - startPosOfKey);
        strncpy(key, startPosOfKey, keySize);
        key[keySize] = '\0';
    }

    static inline float
    flog2(float x)
    {
        if (x <= 0)
            return -128;
        int *px = (int*) (&x);        // store address of float as pointer to long int
        float e = (float) (((*px & 0x7F800000) >> 23) - 0x7f); // shift right by 23 bits and subtract 127 = 0x7f => exponent
        *px = ((*px & 0x007FFFFF) | 0x3f800000);  // set exponent to 127 (i.e., 0)
        x -= 1.0;         // and calculate x-1.0
        x *= (1.441740
                + x * (-0.7077702 + x * (0.4123442 + x * (-0.1903190 + x * 0.0440047)))); // 5'th order polynomial approx. of log(1+x)
        return x + e;
    }

    static inline double fpow2(float x) {
        if (x>=FLT_MAX_EXP) return FLT_MAX;
        if (x<=FLT_MIN_EXP) return 0.0f;

        int *px = (int*) (&x);        // store address of float as pointer to long int
        float tx = (x - 0.5f) + (3 << 22); // temporary value for truncation: x-0.5 is added to a large integer (3<<22),
        // 3<<22 = (1.1bin)*2^23 = (1.1bin)*2^(150-127),
        // which, in internal bits, is written 0x4b400000 (since 10010110bin = 150)
        int lx = *((int*) &tx) - 0x4b400000;   // integer value of x
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
};
#endif
