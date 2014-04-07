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
class Util {
public:
	static void * mem_align(size_t bound, size_t size);
    static void decompose_domain(int domain_size, int world_rank,
                                 int world_size, int* subdomain_start,
                                 int* subdomain_size);    
    static char * skipLine(char * data);
    static char * skipWhitespace(char * data);
    static char * skipNoneWhitespace(char * data);
    static size_t getWordsOfLine(char * data, char ** words, size_t maxElement );
    

};
#endif
