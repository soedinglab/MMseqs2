#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <cstddef>
#include <cstring>
#include <vector>
#include <limits>
#include <map>
#include <cstdint>
#include "MMseqsMPI.h"

#ifndef EXIT
#define EXIT(exitCode) do { int __status = (exitCode); std::cerr.flush(); std::cout.flush(); exit(__status); } while(0)
#endif

#define BIT_SET(a,b) ((a) | (1ULL<<(b)))
#define BIT_CLEAR(a,b) ((a) & ~(1ULL<<(b)))
#define BIT_CHECK(a,b) (!!((a) & (1ULL<<(b))))

template<typename T>
struct assert_false : std::false_type
{ };

template<typename T>
std::string SSTR(T) {
    static_assert(assert_false<T>::value , "Not implemented for requested type");
    return "";
}

template<> std::string SSTR(const char*);
template<> std::string SSTR(char*);
template<> std::string SSTR(bool);
template<> std::string SSTR(const char[]);
template<> std::string SSTR(const std::string&);
template<> std::string SSTR(std::string);
template<> std::string SSTR(char);
template<> std::string SSTR(short);
template<> std::string SSTR(unsigned short);
template<> std::string SSTR(int);
template<> std::string SSTR(unsigned int);
template<> std::string SSTR(long);
template<> std::string SSTR(unsigned long);
template<> std::string SSTR(long long);
template<> std::string SSTR(unsigned long long);
template<> std::string SSTR(double);
template<> std::string SSTR(float);

#define ARRAY_SIZE(a) (sizeof(a) / sizeof(a[0]))

#if defined(__has_builtin)
#  define HAS_BUILTIN(x) __has_builtin(x)
#else
#  define HAS_BUILTIN(x) (0)
#endif

#if HAS_BUILTIN(__builtin_expect)
#  define LIKELY(x) __builtin_expect((x),1)
#  define UNLIKELY(x) __builtin_expect((x),0)
#else
#  define LIKELY(x) (x)
#  define UNLIKELY(x) (x)
#endif

#if defined(__has_attribute)
#  define HAS_ATTRIBUTE(x) __has_attribute(x)
#else
#  define HAS_ATTRIBUTE(x) (0)
#endif

#if HAS_ATTRIBUTE(__unused__)
#  define MAYBE_UNUSED(x) x __attribute__((__unused__))
#else
#  define MAYBE_UNUSED(x) x
#endif

class Util {
public:
    static void decomposeDomain(size_t domain_size, size_t world_rank,
                                size_t world_size, size_t *subdomain_start,
                                size_t *subdomain_size);

    static void rankedDescSort20(short *val, unsigned int *index);

    static size_t getTotalSystemMemory();
    static size_t getPageSize();
    static size_t getTotalMemoryPages();
    static uint64_t getL2CacheSize();

    static char touchMemory(const char* memory, size_t size);

    static size_t countLines(const char *data, size_t length);

    static size_t ompCountLines(const char *data, size_t length, unsigned int threads);

    static int readMapping(std::string mappingFilename, std::vector<std::pair<unsigned int, unsigned int> > & mapping);

    template<typename T>
    static inline T fast_atoi( const char * str )
    {
        T val = 0;
        int sign=1;
        if(std::numeric_limits<T>::is_signed == true){
            if(*str == '-') {
                sign = -1;
                str++;
            }
        }
        while (*str >= '0' && *str <= '9') {
            val = val*10 + (*str++ - '0');
        }
        return sign * val;
    }

    static char* fastSeqIdToBuffer(float seqId, char* buffer);

    static bool isNumber(const std::string& s)
    {
        std::string::const_iterator it = s.begin();
        while (it != s.end() && std::isdigit(*it)) ++it;
        return !s.empty() && it == s.end();
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

    static std::string& implode(const std::vector<std::string> &splits, char sep, std::string& str) {
        for (std::vector<std::string>::const_iterator it = splits.begin(); it != splits.end(); ++it) {
            str += (*it);
            if (it + 1 != splits.end()) {
                str += sep;
            }
        }
        return str;
    }

    static std::string implode(const std::vector<std::string> &splits, char sep) {
        std::string str;
        implode(splits, sep, str);
        return str;
    }

    static inline char * skipLine(char * data){
        while( *data !='\n' ) { data++; }
        return (data+1);
    }

    static inline char * seekToNextEntry(char * data){
        while( *data !='\0' ) { data++; }
        return (data+1);
    }

    static bool getLine(const char* data, size_t dataLength, char* buffer, size_t bufferLength);

    static inline size_t skipWhitespace(const char* data) {
        size_t counter = 0;
        while ((data[counter] == ' ' || data[counter] == '\t') == true) {
            counter++;
        }
        return counter;
    }

    static inline size_t skipTab(const char* data) {
        size_t counter = 0;
        while (data[counter] == '\t') {
            counter++;
        }
        return counter;
    }

    // Return the index i such that data[i] <- '\0' makes a string terminated
    // by a non Whitespace character
    static inline size_t getLastNonWhitespace(char * data, size_t len){
        size_t counter = len;

        if (!counter)
            return 0;

        if (counter && data[counter] == '\0')
            counter--;

        while( (data[counter] == ' ' || data[counter] == '\t')) {
            if(!counter)
                return 0;
            counter--;
        }

        return counter + 1;
    }

    static inline size_t skipNoneWhitespace(const char * data){
        //A value different from zero (i.e., true) if indeed c is a white-space character. Zero (i.e., false) otherwise.
        size_t counter = 0;
        while((data[counter] == ' '  || data[counter] == '\t' || data[counter] == '\n' || data[counter] == '\0') == false) {
            counter++;
        }
        return counter;
    }

    static inline size_t skipNonTab(const char * data) {
        size_t counter = 0;
        while (data[counter] != '\t' && data[counter] != '\n' && data[counter] != '\0') {
            counter++;
        }
        return counter;
    }

    static std::pair<std::string, std::string> createTmpFileNames(const std::string &db,
                                                                  const std::string &dbindex, int count){
        // TODO take only db and not db and dbindex
        // TODO check naming of old database paths
        std::string suffix = std::string("_tmp_") + SSTR(count);
        std::string data = db + suffix;
        
        std::string index = "";
        if (dbindex.compare(db + ".index") == 0) {
            index.append(db + suffix + ".index");
        } else {
            //Debug(Debug::WARNING) << "Database index name is deprecated (" << dbindex << ")\n";
            index.append(dbindex + suffix);
        }

        return std::make_pair(data, index);
    }

    static std::pair<std::string, std::string> databaseNames(const std::string &basename) {
        std::string index = basename;
        index.append(".index");
        return std::make_pair(basename, index);
    };

    static inline size_t getWordsOfLine(const char* data, const char** words, size_t maxElement ){
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

    static inline size_t getFieldsOfLine(const char* data, const char** words, size_t maxElement) {
        size_t elementCounter = 0;
        while (*data != '\n' && *data != '\0') {
            data += skipTab(data);
            words[elementCounter] = data;
            elementCounter++;
            if (elementCounter >= maxElement)
                return elementCounter;
            data += skipNonTab(data);
        }
        if (elementCounter < maxElement)
            words[elementCounter] = data;

        return elementCounter;
    }


    static std::pair<ssize_t,ssize_t> getFastaHeaderPosition(const std::string & header);
    static std::string parseFastaHeader(const char * header);

    static inline char toUpper(char character){
        character += ('a' <= character && character <= 'z') ? ('A' - 'a') : 0;
        return character;
    }

    static void parseKey(const char *data, char * key);

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

    static std::string removeAfterFirstSpace(std::string in) {
        in.erase(in.find_first_of(" "));
        return in;
    }

    static size_t overlappingKmers(int seqLen, unsigned int kmerSize) {
        int kmersPerSize = seqLen - static_cast<int>(kmerSize);
        return  (kmersPerSize >= 0) ? kmersPerSize + 1 :  0;
    }


    template <typename T>
    static size_t hash(T * x, size_t length){
        const size_t INITIAL_VALUE = 0;
        const size_t A = 31;
        size_t h = INITIAL_VALUE;
        for (size_t i = 0; i < length; ++i){
            h = ((h*A) + x[i]);
        }
        return h;
    }

    static int omp_thread_count();

    static std::string removeWhiteSpace(std::string in);

    static std::map<unsigned int, std::string> readLookup(const std::string& lookupFile,
                                                          const bool removeSplit = false);

    static bool canBeCovered(const float covThr, const int covMode, float queryLength, float targetLength);

    static bool hasCoverage(float covThr, int covMode, float queryCov, float targetCov);

    static int swapCoverageMode(int covMode);

    static float computeSeqId(int seqIdMode, int aaIds, int qLen, int tLen, int alnLen);

    static uint64_t revComplement(const uint64_t kmer, const int k);

    static bool hasAlignmentLength(int alnLenThr, int alnLen) {
        return alnLen >= alnLenThr;
    }

    static size_t computeMemory(size_t limit);
};
#endif
