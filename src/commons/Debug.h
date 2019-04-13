#ifndef DEBUG_H
#define DEBUG_H

#include <iostream>
#include "Util.h"
#include "Timer.h"
#include <unistd.h>
#include <stdlib.h>
#include <cstddef>
#include <sys/stat.h>

class TtyCheck {
public:
    TtyCheck()  {
        tty = false;

        bool stdoutIsTty = isatty(fileno(stdout));
        bool stderrtIsTty = isatty(fileno(stderr));
        struct stat stats;
        fstat(fileno(stdin), &stats);
        bool isChr = S_ISCHR (stats.st_mode) == true; // is terminal
        bool isFifo = S_ISFIFO(stats.st_mode) == false; // is no pipe
        bool isReg = S_ISREG(stats.st_mode) == false;
        if (isFifo && stdoutIsTty && stderrtIsTty && isReg && isChr) {
            tty = true;
        }

        char* ttyEnv = getenv("TTY");
        if(ttyEnv != NULL && strcasecmp(ttyEnv, "1") == 0) {
            tty = true;
        }
        if(ttyEnv != NULL && strcasecmp(ttyEnv, "0") == 0) {
            tty = false;
        }

    };
    bool tty;
};

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


    explicit Debug( int level ) : level(level) {
        static TtyCheck check;
        interactive = check.tty;
    };

    ~Debug(){
        if (level <= ERROR && level <= debugLevel){
            std::cout << std::flush;
            if(interactive){
                std::cerr << "\033[" << Color::FG_RED << "m" << buffer << "\033[" << Color::FG_DEFAULT << "m";;
            }else{
                std::cerr << buffer;
            }
            std::cerr << std::flush;
        } else if(level == WARNING && level <= debugLevel){
            if(interactive){
                std::cout << "\033[" << Color::FG_YELLOW << "m" << buffer << "\033[" << Color::FG_DEFAULT << "m";
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
        size_t prevPrintedId;
        size_t totalEntries;
        bool interactive;
        Timer timer;

        const static int BARWIDTH = 70;

    public:
        Progress(size_t totalEntries)
                :  currentPos(0), prevPrintedId(0), totalEntries(totalEntries){
            static TtyCheck check;
            interactive = check.tty;
        }

        Progress() : currentPos(0),  prevPrintedId(0), totalEntries(SIZE_MAX){
            static TtyCheck check;
            interactive = check.tty;
        }

        void updateProgress(){
            size_t id = __sync_fetch_and_add(&currentPos, 1);
            // if no active terminal exists write dots
            if(interactive == false){
                if(totalEntries==SIZE_MAX) {
                    if(id==0) {
                        Debug(INFO) << '[';
                    }
                    if (id % 1000000 == 0 && id > 0){
                        Debug(INFO) << "\t" << (id / 1000000) << " Mio. sequences processed\n";
                        fflush(stdout);
                    }
                    else if (id % 10000 == 0 && id > 0) {
                        Debug(INFO) << "=";
                        fflush(stdout);
                    }
                }else{
                    if(id==0) {
                        Debug(INFO) << '[';
                    }
                    float progress = (static_cast<float>(id) / static_cast<float>(totalEntries-1));
                    if(id > 0) {
                        float prevPrintedProgress = (static_cast<float>(id-1) / static_cast<float>(totalEntries-1));
                        int prevPos = BARWIDTH * prevPrintedProgress;
                        int pos     = BARWIDTH * progress;
                        for (int write = prevPos; write < pos; write++) {
                            Debug(INFO) << '=';
                            fflush(stdout);
                        }
                    }

                    if(id == (totalEntries - 1) ){
                        Debug(INFO) << "] ";
                        Debug(INFO) << timer.lapProgress();
                        Debug(INFO) << "\n";
                    }
                }
            }else{
                if(totalEntries==SIZE_MAX){
                    if(id > prevPrintedId + 100) {
                        std::string line;
                        line.push_back('[');
                        line.append(SSTR(id+1));
                        line.append("] ");
                        line.append(timer.lapProgress());
                        line.push_back('\r');
                        Debug(Debug::INFO) << line;
                        fflush(stdout);
                        prevPrintedId = id;
                    }
                }else{
                    float progress = (static_cast<float>(id) / static_cast<float>(totalEntries-1));
                    float prevPrintedProgress = (static_cast<float>(prevPrintedId) / static_cast<float>(totalEntries-1));
                    if(progress-prevPrintedProgress > 0.001 || id == (totalEntries - 1)  || id == 0 ){
                        std::string line;
                        line.push_back('[');
                        int pos = BARWIDTH * progress;
                        for (int i = 0; i < BARWIDTH; ++i) {
                            if (i < pos) {
                                line.push_back('=');
                            }else if (i == pos) {
                                line.push_back('>');
                            } else {
                                line.push_back(' ');
                            }
                        }
                        char buffer[32];
                        int n = sprintf(buffer, "%.2f", progress * 100.0f);
                        line.append("] ");
                        line.append(buffer, n);
                        line.append(" % ");
                        line.append(timer.lapProgress());
                        //printf("%zu\t%zu\t%f\n", id, totalEntries, progress);
                        line.push_back((id == (totalEntries - 1) ) ? '\n' : '\r' );
                        Debug(Debug::INFO) << line;
                        fflush(stdout);
                        prevPrintedId=id;
                    }
                }
            }
        }
    };




private:
    const int level;
    std::string buffer;
    bool interactive;
};



#endif
