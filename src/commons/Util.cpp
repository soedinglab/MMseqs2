#include "Util.h"
#include "Debug.h"
#include "kseq.h"
#include "FileUtil.h"
#include "BaseMatrix.h"
#include "SubstitutionMatrix.h"
#include "Sequence.h"
#include "Parameters.h"
#include <sys/time.h>
#include <sys/resource.h>

#include <unistd.h>
#ifdef __APPLE__
#include <sys/param.h>
#include <sys/sysctl.h>
#endif

#include <fstream>
#include <algorithm>

KSEQ_INIT(int, read)

size_t Util::countLines(const char *data, size_t length) {
    size_t newlines = 0;
    for (size_t i = 0; i < length; i++ ) {
        if (data[i] == '\n' ) {
            newlines++;
        }
    }
    return newlines;
}

void Util::decomposeDomain(size_t domain_size, size_t world_rank,
                           size_t world_size, size_t *subdomain_start,
                           size_t *subdomain_size) {
    if (world_size > domain_size) {
        // Don't worry about this special case. Assume the domain size
        // is greater than the world size.
        Debug(Debug::ERROR) << "World Size: " << world_size << " aaSize: " << domain_size << "\n";
        EXIT(EXIT_FAILURE);
    }
    *subdomain_start = domain_size / world_size * world_rank;
    *subdomain_size = domain_size / world_size;
    if (world_rank == world_size - 1) {
        // Give remainder to last process
        *subdomain_size += domain_size % world_size;
    }
}

std::map<std::string, size_t> Util::readMapping(const char *fastaFile) {
    std::map<std::string, size_t> map;
    FILE *fasta_file = FileUtil::openFileOrDie(fastaFile, "r", true);
    kseq_t *seq = kseq_init(fileno(fasta_file));
    size_t i = 0;
    while (kseq_read(seq) >= 0) {
        std::string key = Util::parseFastaHeader(seq->name.s);
        if (map.find(key) == map.end()) {
            map[key] = i;
            i++;
        } else {
            Debug(Debug::ERROR) << "Duplicated key " << key << " in function readMapping.\n";
            EXIT(EXIT_FAILURE);
        }
    }
    kseq_destroy(seq);
    fclose(fasta_file);
    return map;
}



template <typename T>
void Util::decomposeDomainByAminoAcid(size_t aaSize, T seqSizes, size_t count,
                                      size_t worldRank, size_t worldSize, size_t *start, size_t *size){
    if (worldSize > aaSize) {
        // Assume the domain size is greater than the world size.
        Debug(Debug::ERROR) << "World Size: " << worldSize << " aaSize: " << aaSize << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (worldSize == 1) {
        *start = 0;
        *size = count;
        return;
    }

    size_t aaPerSplit =  aaSize / worldSize;
    size_t currentRank = 0;
    size_t currentSize = 0;
    *start = 0;
    *size = 0;
    for(size_t i = 0; i < count; i++ ){
        if(currentSize > aaPerSplit){
            currentSize = 0;
            currentRank++;
            if (currentRank > worldRank) {
                *size = (i) - *start;
                break;
            } else if (currentRank == worldRank) {
                *start = i;
                *size = count - *start;
                if (worldRank == worldSize - 1) {
                    break;
                }
            }
        }
        currentSize += seqSizes[i];
    }
}

template
void Util::decomposeDomainByAminoAcid<unsigned int*>(size_t aaSize, unsigned int *seqSizes,
                                                     size_t count, size_t worldRank, size_t worldSize,
                                                     size_t *start, size_t *size);

template
void Util::decomposeDomainByAminoAcid<size_t*>(size_t aaSize, size_t *seqSizes,
                                               size_t count, size_t worldRank, size_t worldSize,
                                               size_t *start, size_t *size);




// http://jgamble.ripco.net/cgi-bin/nw.cgi?inputs=20&algorithm=batcher&output=svg
// sorting networks
void Util::rankedDescSort32(short *val, unsigned int *index){
#define SWAP(x,y){\
if( val[x] < val[y] ){   \
short tmp1 = val[x];    \
val[x] = val[y];     \
val[y] = tmp1;        \
unsigned int tmp2 = index[x];      \
index[x] = index[y]; \
index[y] = tmp2;      \
} \
}
    SWAP(0,16);SWAP(1,17);SWAP(2,18);SWAP(3,19);SWAP(4,20);SWAP(5,21);SWAP(6,22);SWAP(7,23);SWAP(8,24);SWAP(9,25);SWAP(10,26);SWAP(11,27);
    SWAP(12,28);SWAP(13,29);SWAP(14,30);SWAP(15,31);SWAP(0,8);SWAP(1,9);SWAP(2,10);SWAP(3,11);SWAP(4,12);SWAP(5,13);SWAP(6,14);SWAP(7,15);
    SWAP(16,24);SWAP(17,25);SWAP(18,26);SWAP(19,27);SWAP(20,28);SWAP(21,29);SWAP(22,30);SWAP(23,31);SWAP(8,16);SWAP(9,17);SWAP(10,18);
    SWAP(11,19);SWAP(12,20);SWAP(13,21);SWAP(14,22);SWAP(15,23);SWAP(0,4);SWAP(1,5);SWAP(2,6);SWAP(3,7);SWAP(24,28);SWAP(25,29);SWAP(26,30);
    SWAP(27,31);SWAP(8,12);SWAP(9,13);SWAP(10,14);SWAP(11,15);SWAP(16,20);SWAP(17,21);SWAP(18,22);SWAP(19,23);SWAP(0,2);SWAP(1,3);SWAP(28,30);
    SWAP(29,31);SWAP(4,16);SWAP(5,17);SWAP(6,18);SWAP(7,19);SWAP(12,24);SWAP(13,25);SWAP(14,26);SWAP(15,27);SWAP(0,1);SWAP(30,31);
    SWAP(4,8);SWAP(5,9);SWAP(6,10);SWAP(7,11);SWAP(12,16);SWAP(13,17);SWAP(14,18);SWAP(15,19);SWAP(20,24);SWAP(21,25);SWAP(22,26);SWAP(23,27);
    SWAP(4,6);SWAP(5,7);SWAP(8,10);SWAP(9,11);SWAP(12,14);SWAP(13,15);SWAP(16,18);SWAP(17,19);SWAP(20,22);SWAP(21,23);SWAP(24,26);SWAP(25,27);
    SWAP(2,16);SWAP(3,17);SWAP(6,20);SWAP(7,21);SWAP(10,24);SWAP(11,25);SWAP(14,28);SWAP(15,29);
    SWAP(2,8);SWAP(3,9);SWAP(6,12);SWAP(7,13);SWAP(10,16);SWAP(11,17);SWAP(14,20);SWAP(15,21);SWAP(18,24);SWAP(19,25);SWAP(22,28);SWAP(23,29);
    SWAP(2,4);SWAP(3,5);SWAP(6,8);SWAP(7,9);SWAP(10,12);SWAP(11,13);SWAP(14,16);SWAP(15,17);SWAP(18,20);SWAP(19,21);SWAP(22,24);SWAP(23,25);
    SWAP(26,28);SWAP(27,29);SWAP(2,3);SWAP(4,5);SWAP(6,7);SWAP(8,9);SWAP(10,11);SWAP(12,13);SWAP(14,15);SWAP(16,17);SWAP(18,19);SWAP(20,21);
    SWAP(22,23);SWAP(24,25);SWAP(26,27);SWAP(28,29);SWAP(1,16);SWAP(3,18);SWAP(5,20);SWAP(7,22);SWAP(9,24);SWAP(11,26);SWAP(13,28);SWAP(15,30);
    SWAP(1,8);SWAP(3,10);SWAP(5,12);SWAP(7,14);SWAP(9,16);SWAP(11,18);SWAP(13,20);SWAP(15,22);SWAP(17,24);SWAP(19,26);SWAP(21,28);SWAP(23,30);
    SWAP(1,4);SWAP(3,6);SWAP(5,8);SWAP(7,10);SWAP(9,12);SWAP(11,14);SWAP(13,16);SWAP(15,18);SWAP(17,20);SWAP(19,22);SWAP(21,24);SWAP(23,26);
    SWAP(25,28);SWAP(27,30);SWAP(1,2);SWAP(3,4);SWAP(5,6);SWAP(7,8);SWAP(9,10);SWAP(11,12);SWAP(13,14);SWAP(15,16);SWAP(17,18);SWAP(19,20);
    SWAP(21,22);SWAP(23,24);SWAP(25,26);SWAP(27,28);SWAP(29,30);
#undef SWAP
}


// http://jgamble.ripco.net/cgi-bin/nw.cgi?inputs=20&algorithm=batcher&output=svg
// sorting networks
void Util::rankedDescSort20(short *val, unsigned int *index){
#define SWAP(x,y){\
if( val[x] < val[y] ){   \
short tmp1 = val[x];    \
val[x] = val[y];     \
val[y] = tmp1;        \
unsigned int tmp2 = index[x];      \
index[x] = index[y]; \
index[y] = tmp2;      \
} \
}
    SWAP(0,16);SWAP(1,17);SWAP(2,18);SWAP(3,19);SWAP(4,12);SWAP(5,13);SWAP(6,14);SWAP(7,15);
    SWAP(0,8);SWAP(1,9);SWAP(2,10);SWAP(3,11);
    SWAP(8,16);SWAP(9,17);SWAP(10,18);SWAP(11,19);SWAP(0,4);SWAP(1,5);SWAP(2,6);SWAP(3,7);
    SWAP(8,12);SWAP(9,13);SWAP(10,14);SWAP(11,15);SWAP(4,16);SWAP(5,17);SWAP(6,18);SWAP(7,19);SWAP(0,2);SWAP(1,3);
    SWAP(4,8);SWAP(5,9);SWAP(6,10);SWAP(7,11);SWAP(12,16);SWAP(13,17);SWAP(14,18);SWAP(15,19);SWAP(0,1);
    SWAP(4,6);SWAP(5,7);SWAP(8,10);SWAP(9,11);SWAP(12,14);SWAP(13,15);SWAP(16,18);SWAP(17,19);
    SWAP(2,16);SWAP(3,17);SWAP(6,12);SWAP(7,13);SWAP(18,19);
    SWAP(2,8);SWAP(3,9);SWAP(10,16);SWAP(11,17);
    SWAP(2,4);SWAP(3,5);SWAP(6,8);SWAP(7,9);SWAP(10,12);SWAP(11,13);SWAP(14,16);SWAP(15,17);
    SWAP(2,3);SWAP(4,5);SWAP(6,7);SWAP(8,9);SWAP(10,11);SWAP(12,13);SWAP(14,15);SWAP(16,17);
    SWAP(1,16);SWAP(3,18);SWAP(5,12);SWAP(7,14);
    SWAP(1,8);SWAP(3,10);SWAP(9,16);SWAP(11,18);
    SWAP(1,4);SWAP(3,6);SWAP(5,8);SWAP(7,10);SWAP(9,12);SWAP(11,14);SWAP(13,16);SWAP(15,18);
    SWAP(1,2);SWAP(3,4);SWAP(5,6);SWAP(7,8);SWAP(9,10);SWAP(11,12);SWAP(13,14);SWAP(15,16);SWAP(17,18);
#undef SWAP
}

// find start and end position of an identifier in a FASTA header
std::pair<ssize_t,ssize_t> Util::getFastaHeaderPosition(const std::string& header) {
    const std::pair<size_t, size_t> errorPosition = std::make_pair(-1, -1);
    if (header.length() == 0)
        return errorPosition;

    size_t offset = 0;
    if (Util::startWith("consensus_", header)) {
        offset = 10;
    }

    struct Databases {
        std::string prefix;
        unsigned int length;
        unsigned int verticalBarPos;
    };

    const struct Databases databases[] = {
            { "uc",   2, 0}, // Uniclust
            { "cl|",   3, 1},
            { "sp|",   3, 1}, // Swiss prot
            { "tr|",   3, 1}, // trembl
            { "gb|",   3, 1}, // GenBank
            { "ref|",  4, 1}, // NCBI Reference Sequence
            { "pdb|",  4, 1}, // Brookhaven Protein Data Bank
            { "bbs|",  4, 1}, // GenInfo Backbone Id
            { "lcl|",  4, 1}, // Local Sequence identifier
            { "pir||", 5, 1}, // NBRF PIR
            { "prf||", 5, 1}, // Protein Research Foundation
            { "gnl|",  4, 2}, // General database identifier
            { "pat|",  4, 2}, // Patents
            { "gi|",   3, 3}  // NCBI GI
    };
    const unsigned int database_count = 14;

    for (size_t i = 0; i < database_count; ++i) {
        if (Util::startWith(databases[i].prefix, header, offset)) {
            size_t start = offset + databases[i].length;
            if (databases[i].verticalBarPos > 1) {
                for (size_t j = 0; j < databases[i].verticalBarPos - 1; ++j) {
                    size_t end = header.find_first_of('|', start);
                    if (end != std::string::npos) {
                        start = end + 1;
                    } else {
                        return errorPosition;
                    }
                }
            }

            size_t end = header.find_first_of('|', start);
            if (end != std::string::npos) {
                return std::make_pair(start, end);
            } else {
                end = header.find_first_of(" \n", start);
                if (end != std::string::npos) {
                    return std::make_pair(start, end);
                } else {
                    // return until the end of the line
                    return std::make_pair(start, header.length());
                }
            }
        }
    }

    // if we can not find one of the existing database ids,
    // we use the first part of the string or the whole string
    size_t end = header.find_first_of(" \n", offset);
    if (end != std::string::npos) {
        return std::make_pair(offset, end);
    } else {
        // return until the end of the line
        return std::make_pair(offset, header.length());
    }
}


std::string Util::parseFastaHeader(const std::string& header) {
    std::pair<ssize_t, ssize_t> pos = Util::getFastaHeaderPosition(header);
    if(pos.first == -1 && pos.second == -1)
        return "";

    return header.substr(pos.first, pos.second - pos.first);
}

void Util::parseByColumnNumber(char *data, char *key, int position) {
    char *startPosOfKey = data;
    for (int i = 0; i < position; ++i) {
        startPosOfKey = startPosOfKey + Util::skipNoneWhitespace(startPosOfKey);
        startPosOfKey = startPosOfKey + Util::skipWhitespace(startPosOfKey);
    }
    char *endPosOfId = startPosOfKey + Util::skipNoneWhitespace(startPosOfKey);
    ptrdiff_t keySize = (endPosOfId - startPosOfKey);
    strncpy(key, startPosOfKey, keySize);
    key[keySize] = '\0';
}

void Util::parseKey(char *data, char *key) {
    char *startPosOfKey = data;
    char *endPosOfId = data + Util::skipNoneWhitespace(data);
    ptrdiff_t keySize = (endPosOfId - startPosOfKey);
    strncpy(key, data, keySize);
    key[keySize] = '\0';
}

std::vector<std::string> Util::split(const std::string &str, const std::string &sep) {
    std::vector<std::string> arr;

    char *cstr = strdup(str.c_str());
    const char* csep = sep.c_str();
    char *rest;
    char *current = strtok_r(cstr, csep, &rest);
    while (current != NULL) {
        arr.emplace_back(current);
        current = strtok_r(NULL, csep, &rest);
    }
    free(cstr);

    return arr;
}

bool Util::getLine(const char* data, size_t dataLength, char* buffer, size_t bufferLength) {
    size_t keySize = 0;
    while (((data[keySize] != '\n') && (data[keySize] != '\0')) && keySize < dataLength) {
        keySize++;
    }
    size_t maxLength = std::min(keySize + 1, bufferLength);
    strncpy(buffer, data, maxLength);
    buffer[maxLength - 1] = '\0';

    bool didCutoff = (keySize + 1) > bufferLength;
    return didCutoff == false;
}

void Util::checkAllocation(void *pointer, std::string message) {
    if(pointer == NULL){
        Debug(Debug::ERROR) << message << "\n";
        EXIT(EXIT_FAILURE);
    }
}

size_t Util::getPageSize() {
    return sysconf(_SC_PAGE_SIZE);
}

size_t Util::getTotalMemoryPages() {
#if __APPLE__
    size_t mem;
    size_t len = sizeof(mem);
    sysctlbyname("hw.memsize", &mem, &len, NULL, 0);
    static size_t phys_pages = mem / Util::getPageSize();
#else
    static size_t phys_pages = sysconf(_SC_PHYS_PAGES);
#endif
    return phys_pages;
}

size_t Util::getTotalSystemMemory()
{
    // check for real physical memory
    long pages = getTotalMemoryPages();
    long page_size = getPageSize();
    uint64_t sysMemory = pages * page_size;
    // check for ulimit
//    struct rlimit limit;
//    getrlimit(RLIMIT_MEMLOCK, &limit);
    return sysMemory;
}

char Util::touchMemory(char *memory, size_t size) {
    size_t pageSize = getPageSize();
    char bytes = 0;
    for(size_t i = 0; i < size; i+=pageSize){
        bytes ^= memory[i];
    }

    return bytes;
}


/* Scaled log-likelihood ratios for coiled-coil heptat repeat */
const float     ccoilmat[23][7] =
        {
                {249, 310, 74, 797, -713, 235, -102},           // 0 A
                {-85, -6214, -954, -820, -1980, -839, -2538},   // 1 C
                {-2688, 743, 498, -1703, -409, 458, 337},       // 2 D
                {-1269, 1209, 1097, -236, 1582, 1006, 1338},    // 3 E
                {-713, -2590, -939, -447, -2079, -2513, -3270}, // 4 F
                {-2476, -1537, -839, -2198, -1877, -1002, -2079}, // 5 G
                {-527, -436, -537, -171, -1180, -492, -926},      // 6 H
                {878, -1343, -1064, -71, -911, -820, -1241},      // 7 I
                {209, 785, 597, -492, 739, 522, 706},             // 8 K
                {1097, -1313, -1002, 1348, -673, -665, -576},     // 9 L
                {770, -502, -816, 365, -499, -783, -562},         // 10 M
                {207, 520, 768, -1624, 502, 887, 725},            // 11 N
                {-5521, -2225, -4017, -5115, -4605, -5521, -4961},// 12 P
                {-1167, 828, 845, -209, 953, 767, 949},           // 13 Q
                {13, 389, 571, -2171, 511, 696, 611},             // 14 R
                {-1102, -283, -72, -858, -309, -221, -657},       // 15 S
                {-1624, -610, -435, -385, -99, -441, -213},       // 16 T
                {421, -736, -1049, -119, -1251, -1049, -1016},    // 17 V
                {-2718, -2748, -2733, -291, -5115, -2162, -4268}, // 18 W
                {276, -2748, -2513, 422, -1589, -2137, -2343},    // 19 Y
                {0, 0, 0, 0, 0, 0, 0}                             // 20 X
        };

/* Sample Means for 100-residue windows */
const float     aamean100[20] =
        {
                7.44,5.08,4.69,5.36,1.54,3.93,6.24,6.34,2.24,6.09,
                9.72, 6.00,2.39,4.30,4.69,7.23,5.61,1.25,3.31,6.53
        };
/* Sample Standard Deviations for 100-residue windows */
const float     aasd100[20] =
        {
                4.02,3.05,2.87,2.71,1.88,2.63,3.46,3.45,1.79,3.19,
                3.77,3.64,1.71,2.62,3.00,3.63,2.83,1.32,2.18,2.92
        };


std::map<unsigned int, std::string> Util::readLookup(const std::string& file, const bool removeSplit) {
    std::map<unsigned int, std::string> mapping;
    if (file.length() > 0) {
        std::ifstream mappingStream(file);
        if (mappingStream.fail()) {
            Debug(Debug::ERROR) << "File " << file << " not found!\n";
            EXIT(EXIT_FAILURE);
        }

        std::string line;
        while (std::getline(mappingStream, line)) {
            std::vector<std::string> split = Util::split(line, "\t");
            unsigned int id = strtoul(split[0].c_str(), NULL, 10);

            std::string& name = split[1];

            size_t pos;
            if (removeSplit && (pos = name.find_last_of('_')) != std::string::npos) {
                name = name.substr(0, pos);
            }

            mapping.emplace(id, name);
        }
    }

    return mapping;
}

std::map<std::string, unsigned int> Util::readLookupReverse(const std::string& file, const bool removeSplit) {
    std::map<std::string, unsigned int> mapping;
    if (file.length() > 0) {
        std::ifstream mappingStream(file);
        if (mappingStream.fail()) {
            Debug(Debug::ERROR) << "File " << file << " not found!\n";
            EXIT(EXIT_FAILURE);
        }

        std::string line;
        while (std::getline(mappingStream, line)) {
            std::vector<std::string> split = Util::split(line, "\t");
            unsigned int id = strtoul(split[0].c_str(), NULL, 10);
            std::string& name = split[1];

            size_t pos;
            if (removeSplit && (pos = name.find_last_of('_')) != std::string::npos) {
                name = name.substr(0, pos);
            }

            mapping.emplace(name, id);
        }
    }

    return mapping;
}

int Util::omp_thread_count() {
    int n = 0;
#pragma omp parallel reduction(+:n)
    n += 1;
    return n;
}

std::string Util::removeWhiteSpace(std::string in) {
    in.erase(std::remove_if(in.begin(), in.end(), isspace), in.end());
    return in;
}

bool Util::canBeCovered(const float covThr, const int covMode, float queryLength, float targetLength) {
    switch(covMode){
        case Parameters::COV_MODE_BIDIRECTIONAL:
            return ((queryLength / targetLength >= covThr) || (targetLength / queryLength >= covThr));
        case Parameters::COV_MODE_QUERY:
            return ((targetLength / queryLength) >= covThr);
        default:
            return true;
    }
}

bool Util::hasCoverage(float covThr, int covMode, float queryCov, float targetCov){
    switch(covMode){
        case Parameters::COV_MODE_BIDIRECTIONAL:
            return ((queryCov >= covThr) && (targetCov >= covThr));
        case Parameters::COV_MODE_QUERY:
            return (queryCov >= covThr);
        case Parameters::COV_MODE_TARGET:
            return (targetCov >= covThr);
        default:
            return true;
    }
}
