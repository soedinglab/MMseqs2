#ifndef UTIL_H
#define UTIL_H

#define CHECK_BIT(var,pos) ((var) & (1<<(pos)))

#include <stdlib.h>
class Util {
public:
	static void * mem_align(size_t bound, size_t size);
};
#endif
