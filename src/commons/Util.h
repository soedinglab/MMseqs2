#ifndef UTIL_H
#define UTIL_H
#include <math.h>
#include <stdlib.h>
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#ifdef HAVE_MPI
#define EXIT(exitCode) MPI_Finalize(); exit(exitCode)
#else
#define EXIT(exitCode) exit(exitCode)
#endif

#define CHECK_BIT(var,pos) ((var) & (1<<(pos)))
#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
( std::ostringstream() << std::dec << x ) ).str()



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

#include <stdlib.h>
#include <iostream>

class Util {
public:
	static void * mem_align(size_t bound, size_t size);
    static void decompose_domain(int domain_size, int world_rank,
                                 int world_size, int* subdomain_start,
                                 int* subdomain_size);
    static void rankedDescSort20(short * val, int * index);

    
    static inline int fast_atoi( const char * str )
    {
        int val = 0;
        while (*str >= '0' && *str <= '9') {
            val = val*10 + (*str++ - '0');
        }
        return val;
    }
    

    
 
    
    static inline char * skipLine(char * data){
        while( *data !='\n' ) { data++; }
        return (data+1);
    }
    
    static inline char * skipWhitespace(char * data){
        while( isspace(*data) == true ) {
            data++;
        }
        return data;
    }
    
    static inline char * skipNoneWhitespace(char * data){
        //A value different from zero (i.e., true) if indeed c is a white-space character. Zero (i.e., false) otherwise.
        while( isspace(*data) == false ) {
            data++;
        }
        return data;
    }
    
    
    static inline size_t getWordsOfLine(char * data, char ** words, size_t maxElement ){
        size_t elementCounter = 0;
        while(*data !=  '\n' && *data != '\0'){
            data = skipWhitespace(data);
            words[elementCounter++] = data;
            if(elementCounter >= maxElement)
                break;
            data = skipNoneWhitespace(data);
        }
        return elementCounter;
    }


};
#endif
