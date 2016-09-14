#include "Debug.h"


int Debug::debugLevel = Debug::INFO;

Debug::Debug( int level)
{
    this->level = level;
}

void Debug::setDebugLevel (int i) {
    debugLevel = i;
}

void Debug::printProgress(size_t id){
    if (id % 1000000 == 0 && id > 0){
        Debug(INFO) << "\t" << (id / 1000000) << " Mio. sequences processed\n";
        fflush(stdout);
    }
    else if (id % 10000 == 0 && id > 0) {
        Debug(INFO) << ".";
        fflush(stdout);
    }
}