#include "Debug.h"


int Debug::debugLevel = Debug::INFO;

Debug::Debug( int level)
{
    this->level = level;
}

void Debug::setDebugLevel (int i) {
    debugLevel = i;
}
