#include "Util.h"
#include "Debug.h"
#include "kseq.h"
#include <unistd.h>

#include  <stdio.h>
KSEQ_INIT(int, read)

#include <sys/stat.h>

size_t Util::count_lines(const char * file, size_t endPos ) {
    size_t newlines = 0;
    for ( size_t i = 0; i < endPos; i++ ) {
        if ( file[i] == '\n' ) {
            newlines++;
        }
    }
    return newlines;
}

void Util::decompose_domain(size_t domain_size, size_t world_rank,
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
    FILE * fasta_file = fopen(fastaFile, "r");
    if(fasta_file == NULL) { perror(fastaFile);  }
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
    return map;
}

void Util::decomposeDomainByAminoaAcid(size_t aaSize, unsigned int *seqSizes, size_t count,
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

std::string Util::parseFastaHeader(std::string header){
    if(header.length() == 0)
        return "";

    std::vector<std::string> arr = Util::split(header,"|");
    if(arr.size() > 1) {
        if (Util::startWith("cl|",   header) ||
                Util::startWith("sp|",   header) ||
                Util::startWith("tr|",   header) ||
                Util::startWith("ref|",  header) ||
                Util::startWith("pdb|",  header) ||
                Util::startWith("bbs|",  header) ||
                Util::startWith("lcl|",  header) ||
                Util::startWith("pir||", header) ||
                Util::startWith("prf||", header)) {
            return arr[1];
        }
        else if (Util::startWith("gnl|", header) || Util::startWith("pat|", header))
            return arr[2];
        else if (Util::startWith("gi|", header))
            return arr[3];

    }
    arr = Util::split(header," ");
    return arr[0];
}


FILE* Util::openFileOrDie(const char * fileName, const char * mode) {
    struct stat st;
    FILE* file;
    if(stat(fileName, &st) == 0) { errno = EEXIST; perror(fileName); EXIT(EXIT_FAILURE); }
    file = fopen(fileName, mode);
    if(file == NULL) { perror(fileName); EXIT(EXIT_FAILURE); }
    return file;
}



