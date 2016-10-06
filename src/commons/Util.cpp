#include "Util.h"
#include "Debug.h"
#include "kseq.h"
#include "FileUtil.h"
#include "BaseMatrix.h"
#include "SubstitutionMatrix.h"
#include "Sequence.h"

#include <unistd.h>
#ifdef __APPLE__
#include <sys/param.h>
#include <sys/sysctl.h>
#endif

#include <fstream>

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
    FILE * fasta_file = FileUtil::openFileOrDie(fastaFile, "r", true);
    kseq_t *seq = kseq_init(fileno(fasta_file));
    size_t i = 0;
    while (kseq_read(seq) >= 0) {
        std::string key = Util::parseFastaHeader(seq->name.s);
        if(map.find(key) == map.end()){
            map[key] = i;
            i++;
        }else{
            Debug(Debug::ERROR) << "Duplicated key "<< key <<" in function readMapping.\n";
            EXIT(EXIT_FAILURE);
        }
    }
    kseq_destroy(seq);
    return map;
}



void Util::decomposeDomainSizet(size_t aaSize, size_t *seqSizes, size_t count,
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
    for(size_t i = 0; i < count; i++ ){
        if(currentSize > aaPerSplit){
            currentSize = 0;
            currentRank++;
            if(currentRank > worldRank){
                *size = (i) - *start ;
                break;
            }else if (currentRank == worldRank){
                *start = i;
                if(worldRank == worldSize-1){
                    *size = count - *start;
                    break;
                }
            }
        }
        currentSize += seqSizes[i];
    }
}



void Util::decomposeDomainByAminoAcid(size_t aaSize, unsigned int *seqSizes, size_t count,
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
    for(size_t i = 0; i < count; i++ ){
        if(currentSize > aaPerSplit){
            currentSize = 0;
            currentRank++;
            if(currentRank > worldRank){
                *size = (i) - *start ;
                break;
            }else if (currentRank == worldRank){
                *start = i;
                if(worldRank == worldSize-1){
                    *size = count - *start;
                    break;
                }
            }
        }
        currentSize += seqSizes[i];
    }
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

size_t Util::get_phys_pages () {
#if __APPLE__
    uint64_t mem;
    size_t len = sizeof(mem);
    sysctlbyname("hw.memsize", &mem, &len, NULL, 0);
    static unsigned phys_pages = mem/sysconf(_SC_PAGE_SIZE);
#else
    static unsigned phys_pages = sysconf(_SC_PHYS_PAGES);
#endif
    return phys_pages;
}

size_t Util::getTotalSystemMemory()
{
    long pages = get_phys_pages();
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
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

void Util::filterRepeates(int *seq, int seqLen, char *mask, int p, int W, int MM){
    char id[1000];
    for (int j=0; j<W; j++){
        id[j]=0;
    }
    int sum = 0;         // number of identities
    int j = 0;
    for (int i = p; i < seqLen; i++)
    {
        sum -= id[j];
        id[j] = ( seq[i - p] == seq[i] ? 1 : 0 );
        sum += id[j];
        if (sum >= W - MM)
            for (int k = std::max(0,i - W - p + 1); k <= i; k++){
                mask[k] = 1;
            }
        if (++j >= W) {
            j=0;
        }
    }
}
size_t Util::maskLowComplexity(BaseMatrix * m, Sequence *s,
                               int seqLen,
                               int windowSize,
                               int maxAAinWindow,
                               int alphabetSize, int maskValue) {
    size_t aafreq[21];

    char * mask = new char[seqLen];
    memset(mask, 0, seqLen * sizeof(char));

    // Filter runs of 4 identical residues
    filterRepeates(s->int_sequence, seqLen, mask, 1, 4, 0);
//
//    // Filter runs of 4 doublets with a maximum of one mismatch
    filterRepeates(s->int_sequence, seqLen, mask, 2, 8, 1);
//
//    // Filter runs of 4 triplets with a maximum of two mismatches
    filterRepeates(s->int_sequence, seqLen, mask, 3, 9, 2);
    filterByBiasCorrection(s, seqLen, m, mask, 70);

    // filter low complex
    for (int i = 0; i < seqLen - windowSize; i++)
    {
        for (int j = 0; j < alphabetSize; j++){
            aafreq[j] = 0;
        }
        for (int j = 0; j < windowSize; j++){
            aafreq[s->int_sequence[i + j]]++;
        }
        int n = 0;
        for (int j = 0; j < alphabetSize; j++){
            if (aafreq[j]){
                n++; // count # amino acids
            }
        }
        if (n <= maxAAinWindow)
        {
            for (int j = 0; j < windowSize; j++)
                mask[i + j] = 1;
        }
    }

    // filter coil
    // look at 3 possible coils -> 21 pos (3 * 7)
    for (int i = 0; i < seqLen - 21; i++)
    {
        int tot = 0;
        for (int l = 0; l < 21; l++)
            tot += ccoilmat[s->int_sequence[i + l]][l % 7];
        if (tot > 10000)
        {
            for (int l = 0; l < 21; l++)
                mask[i + l] = 1;
        }
    }
    size_t maskedResidues = 0;
    for (int i = 0; i < seqLen; i++){
        maskedResidues += (mask[i]);
        s->int_sequence[i] = (mask[i]) ? maskValue : s->int_sequence[i];
    }
    delete [] mask;
    return maskedResidues;
}

std::map<unsigned int, std::string> Util::readLookup(const std::string& file) {
    std::map<unsigned int, std::string> mapping;
    if (file.length() > 0) {
        std::fstream mappingStream(file);
        if (mappingStream.fail()) {
            Debug(Debug::ERROR) << "File " << file << " not found!\n";
            EXIT(EXIT_FAILURE);
        }

        std::string line;
        while (std::getline(mappingStream, line)) {
            std::vector<std::string> split = Util::split(line, "\t");
            unsigned int id = strtoul(split[0].c_str(), NULL, 10);
            mapping.emplace(id, split[1]);
        }
    }

    return mapping;
}

void Util::filterByBiasCorrection(Sequence *s, int seqLen, BaseMatrix *m, char *mask, int scoreThr) {
    char * compositionBias = new char[seqLen];
    float * tmp_composition_bias = new float[seqLen];
    SubstitutionMatrix::calcLocalAaBiasCorrection(m, s->int_sequence, seqLen, tmp_composition_bias);
    //const int8_t seed_6_spaced[] = {1, 1, 0, 1, 0, 1, 0, 0, 1, 1}; // better than 11101101
    const int8_t pos[] = {0, 1, 3, 5, 8, 9}; // better than 11101101

    for(int i = 0; i < seqLen; i++){
        compositionBias[i] = (int8_t) (tmp_composition_bias[i] < 0.0)
                             ? tmp_composition_bias[i] - 0.5: tmp_composition_bias[i] + 0.5;
    }
    for(int i = 0; i < seqLen - 10; i++){
        int kmerScore = 0;
        for (unsigned int kmerPos = 0; kmerPos < 6; kmerPos++) {
            unsigned int aa_pos = i + pos[kmerPos];
            unsigned int aa = s->int_sequence[aa_pos];
            kmerScore += m->subMatrix[aa][aa] + compositionBias[aa_pos];
        }
        if(kmerScore <= scoreThr) {
            for (unsigned int kmerPos = 0; kmerPos < 6; kmerPos++) {
                unsigned int aa_pos = i +  pos[kmerPos];
                mask[aa_pos] = 1;
            }
        }
    }
    delete [] compositionBias;
    delete [] tmp_composition_bias;
}

int Util::omp_thread_count() {
        int n = 0;
#pragma omp parallel reduction(+:n)
        n += 1;
        return n;
}







