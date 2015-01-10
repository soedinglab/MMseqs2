#ifndef UTIL_H
#define UTIL_H
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <limits.h>



#include <stdlib.h>
#include <iostream>
#include <IOKit/IODataQueueClient.h>

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef OPENMP
#include <omp.h>
#endif

#ifdef HAVE_MPI
#define EXIT(exitCode) MPI_Finalize(); exit(exitCode)
#else
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

static const float _2p23 = 8388608.0f;
typedef struct {
    unsigned int  precision_m;
    unsigned int* pTable_m;
} PowFast;




static void powFastSetTable(unsigned int* const pTable, const unsigned int  precision )
{
    /* step along table elements and x-axis positions */
    float zeroToOne = 1.0f / ((float)(1 << precision) * 2.0f);
    int   i;
    for( i = 0;  i < (1 << precision);  ++i )
    {
        /* make y-axis value for table element */
        const float f = ((float)pow( 2.0f, zeroToOne ) - 1.0f) * _2p23;
        pTable[i] = (unsigned int)( f < _2p23 ? f : (_2p23 - 1.0f) );
        
        zeroToOne += 1.0f / (float)(1 << precision);
    }
}


static void* powFastCreate( unsigned int precision )
{
    PowFast* pPowFast = 0;
    
    precision = (precision <= 18u) ? precision : 18u;
    
    pPowFast = (PowFast*)malloc( sizeof(PowFast) + ((1 << precision) *
                                                    sizeof(*(pPowFast->pTable_m))) );
    if( pPowFast )
    {
        pPowFast->precision_m = precision;
        pPowFast->pTable_m   = (unsigned int*)((char*)pPowFast + sizeof(PowFast));
        
        powFastSetTable( pPowFast->pTable_m, pPowFast->precision_m );
    }
    
    return pPowFast;
}

inline static float powFastLookup ( const float val, const float ilog2,
  unsigned int* const pTable, const unsigned int precision )
{
    /* build float bits */
    const int i = (int)( (val * (_2p23 * ilog2)) + (127.0f * _2p23) );
    
    /* replace mantissa with lookup */
    const int it = (i & 0xFF800000) | pTable[(i & 0x7FFFFF) >> (23 - precision)];
    
    /* convert bits to float */
    union { int i; float f; } pun;
    return pun.i = it,  pun.f;
}


static const void* powFastAdj = powFastCreate( 11 );

static inline float powFast2 ( float f)
{ 
    const PowFast* ppf = (const PowFast*)powFastAdj;
    
    return powFastLookup( f, 1.0f, ppf->pTable_m, ppf->precision_m );
}







class Util {
public:
    static void decompose_domain(size_t domain_size, size_t world_rank,
            size_t world_size, size_t *subdomain_start,
            size_t *subdomain_size);
    static void rankedDescSort20(short *val, unsigned int *index);
    static void decomposeDomainByAminoaAcid(size_t aaSize, unsigned short *seqSizes, size_t count,
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
        while(( data[counter] == ' ' || data[counter] == '\t' || data[counter] == '\n') == false ) {
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
};
#endif
