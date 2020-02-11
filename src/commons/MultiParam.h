#ifndef MULTIPARAM_H
#define MULTIPARAM_H

/*
 * MultiParam: class to store sequence type specific parameter values
 * written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
 */

#include <string>
#include <cstring>
#include <limits.h>

template <class T>
class MultiParam {

public:
    T aminoacids;
    T nucleotides;

    MultiParam(T aminoacids, T nucleotides);
    MultiParam(const char* parametercstring);
    static std::string format(const MultiParam<T> multiparam);
    MultiParam& operator=(T value);
};


#endif