#ifndef MULTIPARAM_H
#define MULTIPARAM_H

/*
 * MultiParam: class to store sequence type specific parameter values
 * written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
 */

#include <string>
#include <cstring>
#include <limits.h>
#include <stdlib.h>

template <class T>
class MultiParam {

public:
    T aminoacids;
    T nucleotides;
    MultiParam(){};
    MultiParam(T aminoacids, T nucleotides);
    MultiParam(const char* parametercstring);
    ~MultiParam(){};

    static std::string format(const MultiParam<T> &multiparam);

    MultiParam<T>& operator=(T value) {
        nucleotides = value;
        aminoacids = value;
        return *this;
    }
};

template <>
class MultiParam<char*> {

public:
    char* aminoacids;
    char* nucleotides;

    MultiParam(const char* aminoacids, const char* nucleotides);
    MultiParam(const char* parametercstring);
    MultiParam(const MultiParam<char*> & par) {
        nucleotides = strdup(par.nucleotides);
        aminoacids = strdup(par.aminoacids);
    };

    explicit MultiParam<char*>(const std::string& parameterstring) : MultiParam<char*>(parameterstring.c_str()) {}
    ~MultiParam();

    static std::string format(const MultiParam<char*> &multiparam);
    //MultiParam<char>& operator=(char* value);



    //ScoreMatrixFile(const ScoreMatrixFile& copy) : ScoreMatrixFile(copy.aminoacids, copy.nucleotides) {}

    MultiParam<char*>& operator=(const MultiParam<char*>& other) {
        if (nucleotides != NULL) {
            free(nucleotides);
        }
        if (aminoacids != NULL) {
            free(aminoacids);
        }
        nucleotides = strdup(other.nucleotides);
        aminoacids = strdup(other.aminoacids);
        return *this;
    }


    bool operator==(const char* other) const;
    bool operator==(const std::string& other) const;
    bool operator==(const MultiParam<char*>& other) const;

    bool operator!=(const char* other) const {
        return !(operator==(other));
    }
    bool operator!=(const std::string& other) const {
        return !(operator==(other));
    }
    bool operator!=(const MultiParam<char*>& other) const {
        return !(operator==(other));
    }


};


#endif
