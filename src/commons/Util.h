#ifndef UTIL_H
#define UTIL_H

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
        while( isspace(*data) ) { data++; }
        return data;
    }
    
    static inline char * skipNoneWhitespace(char * data){
        while( !isspace(*data) ) { data++; }
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
