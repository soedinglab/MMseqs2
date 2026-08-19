#ifndef MMSEQS_TIMER_H
#define MMSEQS_TIMER_H

#include <sys/time.h>
#include <string>

class Timer {
public:
    Timer() {
        reset();
    };

    std::string lapProgress();

    std::string lap();

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
