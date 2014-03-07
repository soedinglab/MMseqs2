#ifndef UTIL_H
#define UTIL_H
#include <stdlib.h>
class Util {
public:
	static void * mem_align(size_t bound, size_t size);
};
#endif
