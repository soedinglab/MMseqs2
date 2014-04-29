#include "Util.h"



void * Util::mem_align(size_t boundary, size_t size) {
  void *pointer;
  if (posix_memalign(&pointer,boundary,size) != 0) {
	std::cerr<<"Error: Could not allocate memory by memalign. Please report this bug to developers\n";
	exit(3);
   }
   return pointer;
}

void Util::decompose_domain(int domain_size, int world_rank,
                      int world_size, int* subdomain_start,
                      int* subdomain_size) {
    if (world_size > domain_size) {
        // Don't worry about this special case. Assume the domain size
        // is greater than the world size.
        EXIT(1);
    }
    *subdomain_start = domain_size / world_size * world_rank;
    *subdomain_size = domain_size / world_size;
    if (world_rank == world_size - 1) {
        // Give remainder to last process
        *subdomain_size += domain_size % world_size;
    }
}

void Util::decomposeDomainByAminoaAcid(int aaSize, unsigned short * seqSizes, size_t count,
                            int worldRank, int worldSize, int *start, int *size){
    if (worldSize > aaSize) {
        // Assume the domain size is greater than the world size.
        EXIT(1);
    }
    if (worldSize == 1) {
        *start = 0;
        *size = count;
        return;
    }
    
    int aaPerSplitt =  aaSize / worldSize;
    int currentRank = 0;
    int currentSize = 0;
    *start = 0;
    for(int i = 0; i < count; i++ ){
        if(currentSize > aaPerSplitt){
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
    std::cout << currentSize << std::endl;
}


// http://jgamble.ripco.net/cgi-bin/nw.cgi?inputs=20&algorithm=bosenelson&output=svg
// // sorting networks
void Util::rankedDescSort20(short * val, int * index){
#define SWAP(x,y){\
    if( val[x] < val[y] ){   \
        int tmp = val[x];    \
        val[x] = val[y];     \
        val[y] = tmp;        \
        tmp = index[x];      \
        index[x] = index[y]; \
        index[y] = tmp;      \
    } \
}
    SWAP(0,1);SWAP(3,4);SWAP(5,6);SWAP(8,9);SWAP(10,11);SWAP(13,14);SWAP(15,16);SWAP(18,19);
    SWAP(2,4);SWAP(7,9);SWAP(12,14);SWAP(17,19);
    SWAP(2,3);SWAP(1,4);SWAP(7,8);SWAP(6,9);SWAP(12,13);SWAP(11,14);SWAP(17,18);SWAP(16,19);
    SWAP(0,3);SWAP(5,8);SWAP(4,9);SWAP(10,13);SWAP(15,18);SWAP(14,19);
    SWAP(0,2);SWAP(1,3);SWAP(5,7);SWAP(6,8);SWAP(10,12);SWAP(11,13);SWAP(15,17);SWAP(16,18);SWAP(9,19);
    SWAP(1,2);SWAP(6,7);SWAP(0,5);SWAP(3,8);SWAP(11,12);SWAP(16,17);SWAP(10,15);SWAP(13,18);
    SWAP(1,6);SWAP(2,7);SWAP(4,8);SWAP(11,16);SWAP(12,17);SWAP(14,18);SWAP(0,10);
    SWAP(1,5);SWAP(3,7);SWAP(11,15);SWAP(13,17);SWAP(8,18);
    SWAP(4,7);SWAP(2,5);SWAP(3,6);SWAP(14,17);SWAP(12,15);SWAP(13,16);SWAP(1,11);SWAP(9,18);
    SWAP(4,6);SWAP(3,5);SWAP(14,16);SWAP(13,15);SWAP(1,10);SWAP(2,12);SWAP(7,17);
    SWAP(4,5);SWAP(14,15);SWAP(3,13);SWAP(2,10);SWAP(6,16);SWAP(8,17);
    SWAP(4,14);SWAP(3,12);SWAP(5,15);SWAP(9,17);SWAP(8,16);
    SWAP(4,13);SWAP(3,11);SWAP(6,15);SWAP(9,16);
    SWAP(4,12);SWAP(3,10);SWAP(7,15);
    SWAP(4,11);SWAP(8,15);SWAP(7,12);
    SWAP(4,10);SWAP(9,15);SWAP(6,11);SWAP(8,13);
    SWAP(5,10);SWAP(9,14);SWAP(8,12);
    SWAP(6,10);SWAP(9,13);SWAP(8,11);
    SWAP(9,12);SWAP(7,10);
    SWAP(9,11);SWAP(8,10);
    SWAP(9,10);
#undef SWAP
#undef min
#undef max
#undef min2
#undef max2
}

