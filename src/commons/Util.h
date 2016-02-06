#ifndef UTIL_H
#define UTIL_H

#include <cstddef>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <map>

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

#include "MMseqsMPI.h"

#ifndef EXIT
#define EXIT(exitCode) exit(exitCode)
#endif

#define CHECK_BIT(var,pos) ((var) & (1<<(pos)))
#if __cplusplus <= 199711L
#define SSTR( x ) \
dynamic_cast< std::ostringstream& >( \
( std::ostringstream().flush() << std::dec << x ) ).str()
#else
#define SSTR( x ) std::to_string(x)
#endif

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
    static void decomposeDomain(size_t domain_size, size_t world_rank,
                                size_t world_size, size_t *subdomain_start,
                                size_t *subdomain_size);
    static void rankedDescSort20(short *val, unsigned int *index);
    static void decomposeDomainByAminoaAcid(size_t aaSize, unsigned int *seqSizes, size_t count,
            size_t worldRank, size_t worldSize, size_t *start, size_t *end);

    static size_t countLines(const char *data, size_t length);

    static inline int fast_atoi( const char * str )
    {
        int val = 0;
        while (*str >= '0' && *str <= '9') {
            val = val*10 + (*str++ - '0');
        }
        return val;
    }

    static bool startWith(std::string prefix, std::string str){
        return (!str.compare(0, prefix.size(), prefix));
    }

    static std::vector<std::string> split(std::string str, std::string sep) {
		char buffer[1024];
		snprintf(buffer, 1024, "%s", str.c_str());
		char *cstr = (char *)&buffer;
		char *current;
		char *rest;
		std::vector<std::string> arr;
		current = strtok_r(cstr, sep.c_str(), &rest);
		while (current != NULL) {
			arr.push_back(current);
			current = strtok_r(NULL, sep.c_str(), &rest);
		}
		return arr;
    }

    static inline char * skipLine(char * data){
        while( *data !='\n' ) { data++; }
        return (data+1);
    }

    static inline bool getLine(const char* data, size_t dataLength, char* buffer, size_t bufferLength) {
        size_t keySize = 0;
        while (((data[keySize] != '\n') && (data[keySize] != '\0')) && keySize < dataLength) {
            keySize++;
        }
        size_t maxLength = std::min(keySize + 1, bufferLength);
        strncpy(buffer, data, maxLength);
        buffer[maxLength - 1] = '\0';

        return bufferLength > dataLength;
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

    static std::pair<std::string, std::string> createTmpFileNames(std::string db, std::string dbindex, int numb){
        std::string splitSuffix = std::string("_tmp_") + SSTR(numb);
        std::string dataFile  = db + splitSuffix;
        std::string indexFile = dbindex + splitSuffix;
        return std::make_pair(dataFile, indexFile);
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

    // Compute the sum of bits of one or two integers
    inline static int NumberOfSetBits(int i)
    {
        i = i - ((i >> 1) & 0x55555555);
        i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
        return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
    }

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

    static std::string base_name(std::string const & path, std::string const & delims)
    {
        return path.substr(path.find_last_of(delims) + 1);
    }

    static std::string remove_extension(std::string const & filename)
    {
        typename std::string::size_type const p(filename.find_last_of('.'));
        return p > 0 && p != std::string::npos ? filename.substr(0, p) : filename;
    }

    static std::map<std::string, size_t> readMapping(const char *fastaFile);


};
#endif
