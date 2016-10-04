#ifndef UTIL_H
#define UTIL_H

#include <cstddef>
#include <cstring>
#include <vector>
#include <sstream>
#include <map>

#include "MMseqsMPI.h"
#include "Sequence.h"
#include "BaseMatrix.h"

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

#ifndef __has_builtin
#define __has_builtin(x) 0
#endif

#if defined(__GNUC__) || __has_builtin(__builtin_expect)
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
    static void decomposeDomainByAminoAcid(size_t aaSize, unsigned int *seqSizes, size_t count,
                                           size_t worldRank, size_t worldSize, size_t *start, size_t *end);
    static void decomposeDomainSizet(size_t aaSize, size_t *seqSizes, size_t count,
                                     size_t worldRank, size_t worldSize, size_t *start, size_t *size);
    static size_t getTotalSystemMemory();
    static size_t get_phys_pages();
    static size_t countLines(const char *data, size_t length);

    static inline int fast_atoi( const char * str )
    {
        int val = 0;
        while (*str >= '0' && *str <= '9') {
            val = val*10 + (*str++ - '0');
        }
        return val;
    }

    static bool startWith(const std::string &prefix, const std::string &str, const size_t offset = 0){
        if (str.length() < prefix.length()) {
            return false;
        }

        return (!str.compare(offset, prefix.length(), prefix));
    }

    static bool endsWith(const std::string &suffix, const std::string &str){
        if (str.length() < suffix.length()) {
            return false;
        }

        return (!str.compare(str.length() - suffix.length(), suffix.length(), suffix));
    }

    static std::vector<std::string> split(const std::string &str, const std::string &sep);

    static inline char * skipLine(char * data){
        while( *data !='\n' ) { data++; }
        return (data+1);
    }

    static inline char * seekToNextEntry(char * data){
        while( *data !='\0' ) { data++; }
        return (data+1);
    }

    static bool getLine(const char* data, size_t dataLength, char* buffer, size_t bufferLength);

    static inline size_t skipWhitespace(char * data){
        size_t counter = 0;
        while( (data[counter] == ' ' || data[counter] == '\t') == true ) {
            counter++;
        }
        return counter;
    }

    // Return the index i such that data[i] <- '\0' makes a string terminated
    // by a non Whitespace character
    static inline size_t getLastNonWhitespace(char * data, size_t len){
        size_t counter = len;

        if (counter && data[counter] == '\0')
            counter--;

        while( (data[counter] == ' ' || data[counter] == '\t')) {
            if(!counter)
                return 0;
            counter--;
        }

        return counter + 1;
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

    static std::pair<std::string, std::string> createTmpFileNames(const std::string &db,
                                                                  const std::string &dbindex, int count){
        std::string suffix = std::string("_tmp_") + SSTR(count);
        std::string data  = db + suffix;
        std::string index = dbindex + suffix;
        return std::make_pair(data, index);
    }

    static std::pair<std::string, std::string> databaseNames(const std::string &basename) {
        std::string index = basename;
        index.append(".index");
        return std::make_pair(basename, index);
    };

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
        if(elementCounter < maxElement)
            words[elementCounter] = data;

        return elementCounter;
    }


    static std::pair<ssize_t,ssize_t> getFastaHeaderPosition(const std::string& header);
    static std::string parseFastaHeader(const std::string& header);

    static inline char toUpper(char character){
        character += ('a' <= character && character <= 'z') ? ('A' - 'a') : 0;
        return character;
    }

    static void parseKey(char *data, char * key);

    static void parseByColumnNumber(char *data, char * key, int position);

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

    static std::map<unsigned int, std::string> readLookup(const std::string& fastaFile);

    static void checkAllocation(void *pointer, std::string message);

    template <typename Iterator, typename Container>
    static bool isLastIterator(Iterator iterator, const Container& container) {
        return (iterator != container.end()) && (++iterator == container.end());
    }

    static std::string csvEscape(const std::string &s) {
        size_t n = s.size(), wp = 0;
        std::vector<char> result(n * 2);
        for (size_t i = 0; i < n; i++) {
            if (s[i] == '\n' || s[i] == '\t') {
                result[wp++] = '\\';
                result[wp++] = 'n';
                continue;
            }
            result[wp++] = s[i];
        }
        return std::string(&result[0], &result[wp]);
    }

    static size_t maskLowComplexity(BaseMatrix * mat, Sequence *sequence, int seqLen, int windowSize, int maxAAinWindow, int alphabetSize, int maskValue);

    static void filterRepeates(int *seq, int seqLen, char *mask, int p, int W, int MM);

    static void filterByBiasCorrection(Sequence *s, int seqLen, BaseMatrix *m, char *mask, int scoreThr);

    static std::string removeAfterFirstSpace(std::string in) {
        in.erase(in.find_first_of(" "));
        return in;
    }

    static size_t overlappingKmers(int seqLen, unsigned int kmerSize) {
        int kmersPerSize = seqLen - static_cast<int>(kmerSize);
        return  (kmersPerSize >= 0) ? kmersPerSize + 1 :  0;
    }

    static int omp_thread_count();
};
#endif
