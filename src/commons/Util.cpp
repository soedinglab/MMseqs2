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
