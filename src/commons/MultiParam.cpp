#include "MultiParam.h"
#include <stdio.h>
#include <cfloat>

#include "Util.h"

template <typename T>
MultiParam<T>::MultiParam(T aminoacids, T nucleotides) {
    this->nucleotides = nucleotides;
    this->aminoacids = aminoacids;
}

template <typename T>
std::string MultiParam<T>::format(const MultiParam<T> &multiparam) {
    if (multiparam.nucleotides == multiparam.aminoacids) {
        return SSTR(multiparam.nucleotides);
    } else {
        return std::string("nucl:") + SSTR(multiparam.nucleotides) + ",aa:" + SSTR(multiparam.aminoacids);
    }
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

/* explicit implementation for MultiParam<char*> */

MultiParam<char*>::MultiParam(const char*  aminoacids, const char*  nucleotides) {
    this->nucleotides = strdup(nucleotides);
    this->aminoacids = strdup(aminoacids);
}

MultiParam<char*>::MultiParam(const char* filename) {
    if (strchr(filename, ',') != NULL) {
        size_t len = strlen(filename);
        aminoacids = (char*) malloc(len * sizeof(char));
        nucleotides = (char*) malloc(len * sizeof(char));
        if (sscanf(filename, "aa:%[^,],nucl:%s", aminoacids, nucleotides) != 2 && sscanf(filename, "nucl:%[^,],aa:%s", nucleotides, aminoacids) != 2) {
            free((char*)nucleotides);
            free((char*)aminoacids);
            nucleotides = strdup("INVALID");
            aminoacids = strdup("INVALID");
        }
    } else {
        nucleotides = strdup(filename);
        aminoacids = strdup(filename);
    }
}

MultiParam<char*>::~MultiParam() {
    free(nucleotides);
    free(aminoacids);
}

std::string MultiParam<char*>::format(const MultiParam<char*> &file) {
    /*if (strncmp(file.nucleotides, file.aminoacids, strlen(file.aminoacids)) == 0) {
        return file.nucleotides;
    } else {*/
        return std::string("nucl:") + file.nucleotides + ",aa:" + file.aminoacids;
    //}
}


bool MultiParam<char*>::operator==(const char* other) const {
    return strncmp(other, nucleotides, strlen(nucleotides)) == 0 || strncmp(other, aminoacids, strlen(aminoacids)) == 0;
}

bool MultiParam<char*>::operator==(const std::string& other) const {
    return strncmp(other.c_str(), nucleotides, strlen(nucleotides)) == 0 || strncmp(other.c_str(), aminoacids, strlen(aminoacids)) == 0;
}

bool MultiParam<char*>::operator==(const MultiParam<char*>& other) const {
    return strncmp(other.nucleotides, nucleotides, strlen(nucleotides)) == 0 && strncmp(other.aminoacids, aminoacids, strlen(aminoacids)) == 0;
}
