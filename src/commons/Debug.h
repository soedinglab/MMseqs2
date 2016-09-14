#ifndef DEBUG_H
#define DEBUG_H

#include <iostream>
#include <cstddef>

class Debug
{

public:
    static const int NOTHING = 0;
    static const int ERROR = 1;
    static const int WARNING = 2;
    static const int INFO = 3;

    static int debugLevel;

    explicit Debug( int level );

    template<typename T>
    Debug& operator<<( T t)
    {
        if (level <= ERROR && level <= debugLevel){
            std::cerr << t;
            std::cerr << std::flush;
            return *this;
        }
        else if(level > ERROR && level <= debugLevel)
        {
            std::cout << t;
            std::cout << std::flush;
            return *this;
        }
        else{
            return *this;
        }
    }
    static void setDebugLevel(int i);

    static void printProgress(size_t id);

private:
    int level;
};


#endif
