#ifndef DEBUG_H
#define DEBUG_H

#include <iostream>
#include "Util.h"
#include <ostream>
#include <unistd.h>
#include <cstddef>

class Debug
{

public:
    static const int NOTHING = 0;
    static const int ERROR = 1;
    static const int WARNING = 2;
    static const int INFO = 3;

    enum Color {
        FG_RED      = 31,
        FG_GREEN    = 32,
        FG_YELLOW    = 33,
        FG_BLUE     = 34,
        FG_DEFAULT  = 39,
        BG_RED      = 41,
        BG_GREEN    = 42,
        BG_BLUE     = 44,
        BG_DEFAULT  = 49
    };

    static int debugLevel;


    Debug( int level ) : level(level), errIsTty(isatty(fileno(stderr))), outIsTty(isatty(fileno(stdout))) {};

    ~Debug(){
        if (level <= ERROR && level <= debugLevel){
            std::cout << std::flush;
            if(errIsTty){
                std::cerr << "\033[" << Color::FG_RED << "m" << buffer;
            }else{
                std::cerr << buffer;
            }
            std::cerr << std::flush;
        } else if(level == WARNING && level <= debugLevel){
            if(outIsTty){
                std::cout << "\033[" << Color::FG_YELLOW << "m" << buffer;
            }else{
                std::cout << buffer;
            }
            std::cout << std::flush;
        } else if(level > WARNING && level <= debugLevel) {
            std::cout << buffer;
//            std::cout << std::flush;
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

    class Progress{
        private:
            size_t currentPos;
            size_t totalEntries;
            float prevPrintedProgress;
            const static int BARWIDTH = 70;
        public:
            Progress(size_t totalEntries) :  currentPos(0), totalEntries(totalEntries){
                prevPrintedProgress = 0.0;
            }

            void updateProgress(){
                size_t id = __sync_add_and_fetch(&currentPos, 1);
                float progress = (static_cast<float>(id) / static_cast<float>(totalEntries-1));
                if(progress-prevPrintedProgress > 0.001 || id == (totalEntries - 1)  ){

                    Debug(Debug::INFO)  << "[";
                    int pos = BARWIDTH * progress;
                    for (int i = 0; i < BARWIDTH; ++i) {
                        if (i < pos) Debug(Debug::INFO) << "=";
                        else if (i == pos) Debug(Debug::INFO)  << ">";
                        else Debug(Debug::INFO)  << " ";
                    }
                    char buffer[32];
                    int n = sprintf(buffer, "%.2f", progress * 100.0f);
                    std::string progressPercent(buffer, n);
                    Debug(Debug::INFO)  << "] " << progressPercent << " %\r";
                    std::cout.flush();
                    prevPrintedProgress=progress;
                }
                if(id == (totalEntries - 1) ){
                    Debug(Debug::INFO)  << "\n";
                }
            }
    };
private:
    const int level;
    std::string buffer;
    const int errIsTty;
    const int outIsTty;

};



#endif
