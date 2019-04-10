#ifndef DEBUG_H
#define DEBUG_H

#include <iostream>
#include "Util.h"

#include <cstddef>

class Debug
{

public:
    static const int NOTHING = 0;
    static const int ERROR = 1;
    static const int WARNING = 2;
    static const int INFO = 3;

    static int debugLevel;
    Debug( int level ) : level(level) {};

    ~Debug(){
        if (level <= ERROR && level <= debugLevel){
            std::cout << std::flush;
            std::cerr << buffer;
            std::cerr << std::flush;
        }
        else if(level > ERROR && level <= debugLevel)
        {
            std::cout << buffer;
            std::cout << std::flush;
        }
    }


    template<typename T>
    Debug& operator<<( T t)
    {
        buffer.append(SSTR(t));
        return *this;
    }

    Debug& operator<<(double val){
        char str[64];
        snprintf(str, sizeof(str), "%f", val);
        buffer.append(str);
        return *this;
    }

    Debug& operator<<(float val){
        char str[64];
        snprintf(str, sizeof(str), "%f", val);
        buffer.append(str);
        return *this;
    }
    static void setDebugLevel(int i);

    static void printProgress(size_t id);

private:
    const int level;
    std::string buffer;

};



#endif
