#ifndef MMSEQS_TIMER_H
#define MMSEQS_TIMER_H

#include <sys/time.h>
#include <string>
#include <sstream>

class Timer {
public:
    Timer() {
        reset();
    };

    std::string lapProgress() {
        struct timeval end;
        gettimeofday(&end, NULL);
        std::ostringstream ss;
        double timediff = (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
        time_t sec = (time_t)timediff;
        time_t msec = (time_t)((timediff - sec) * 1e3);

        if(sec >= 3600){
            ss << (sec / 3600) << "h ";
        }
        if(sec >= 60){
            ss << (sec % 3600 / 60) << "m ";
        }
        ss << (sec % 60) << "s " << msec << "ms";
        return ss.str();
    }

    std::string lap() {
        struct timeval end;
        gettimeofday(&end, NULL);
        std::ostringstream ss;
        double timediff = (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
        time_t sec = (time_t)timediff;
        time_t msec = (time_t)((timediff - sec) * 1e3);
        ss << (sec / 3600) << "h " << (sec % 3600 / 60) << "m " << (sec % 60) << "s " << msec << "ms";
        return ss.str();
    }

    double getTimediff(){
        struct timeval end;
        gettimeofday(&end, NULL);
        double timediff = (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
        return timediff;
    }

    void reset() {
        gettimeofday(&start, NULL);
    }
private:
    struct timeval start;
};

#endif
