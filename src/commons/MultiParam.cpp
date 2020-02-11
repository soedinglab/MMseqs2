
#include "MultiParam.h"
#include <stdio.h>
#include <stdlib.h>
#include <cfloat>

template <typename T>
MultiParam<T>::MultiParam(T aminoacids, T nucleotides) {
    this->nucleotides = nucleotides;
    this->aminoacids = aminoacids;
}

template <typename T>
std::string MultiParam<T>::format(const MultiParam<T> multiparam) {
    if (multiparam.nucleotides == multiparam.aminoacids) {
        return std::to_string(multiparam.nucleotides);
    } else {
        return std::string("nucl:") + std::to_string(multiparam.nucleotides) + ",aa:" + std::to_string(multiparam.aminoacids);
    }
}


template <typename T>
MultiParam<T>& MultiParam<T>::operator=(T value) {
    nucleotides = value;
    aminoacids = value;
    return *this;
}

template<>
MultiParam<int>::MultiParam(const char* parametercstring) {
    if (strchr(parametercstring, ',') != NULL) {
        if (sscanf(parametercstring, "aa:%d,nucl:%d", &aminoacids, &nucleotides) != 2 &&
            sscanf(parametercstring, "nucl:%d,aa:%d", &nucleotides, &aminoacids) != 2) {
            nucleotides = INT_MAX;
            aminoacids = INT_MAX;
        }
    } else {

        if (sscanf(parametercstring, "%d", &aminoacids) != 1) {
            nucleotides = INT_MAX;
            aminoacids = INT_MAX;
        }
        else
            nucleotides = aminoacids;
    }
}

template<>
MultiParam<float>::MultiParam(const char* parametercstring) {
    if (strchr(parametercstring, ',') != NULL) {
        if (sscanf(parametercstring, "aa:%f,nucl:%f", &aminoacids, &nucleotides) != 2 &&
            sscanf(parametercstring, "nucl:%f,aa:%f", &nucleotides, &aminoacids) != 2) {
            nucleotides = FLT_MAX;
            aminoacids = FLT_MAX;
        }
    } else {

        if (sscanf(parametercstring, "%f", &aminoacids) != 1) {
            nucleotides = FLT_MAX;
            aminoacids = FLT_MAX;
        }
        else
            nucleotides = aminoacids;
    }
}

template class MultiParam<int>;
template class MultiParam<float>;