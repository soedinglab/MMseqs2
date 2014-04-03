#include "Util.h"
#include <iostream>



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

