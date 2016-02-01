#include "Log.h"
#include "Debug.h"

#include <cstdio>

void Log::printProgress(size_t id){
    if (id % 1000000 == 0 && id > 0){
        Debug(Debug::INFO) << "\t" << (id/1000000) << " Mio. sequences processed\n";
        fflush(stdout);
    }
    else if (id % 10000 == 0 && id > 0) {
        Debug(Debug::INFO) << ".";
        fflush(stdout);
    }
}
