#ifndef DEBUG_H
#define DEBUG_H

#include <iostream>

class Debug
{

    public:
        static const int NOTHING = 0;
        static const int ERROR = 1;
        static const int WARNING = 2;
        static const int INFO = 3;

        static int debugLevel;

        Debug( int level ); 

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
    private:
        int level;
};


#endif
