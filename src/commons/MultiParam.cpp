#include "MultiParam.h"

template<typename T>
MultiParam<T>::MultiParam(const char *parametercstring) {
    values.first = T::max;
    values.second = T::max;
    if (strchr(parametercstring, ',') == NULL) {
        if (assign(parametercstring, values.second) == false) {
            values.first = values.second = T::max;
        } else {
            values.first = values.second;
        }
        return;
    }
    std::vector<std::string> params = Util::split(parametercstring, ",");
    if (params.size() != 2) {
        return;
    }
    for (size_t i = 0; i < params.size(); i++) {
        if (params[i].rfind(T::constFirst + ":", 0) == 0) {
            std::vector<std::string> pair = Util::split(params[i], ":");
            if (pair.size() != 2 || assign(pair[1], values.first) == false) {
                values.first = T::max;
            }
        }
        if (params[i].rfind(T::constSecond + ":", 0) == 0) {
            std::vector<std::string> pair = Util::split(params[i], ":");
            if (pair.size() != 2 || assign(pair[1], values.second) == false) {
                values.second = T::max;
            }
        }
    }
}

template<typename T>
std::string MultiParam<T>::format(const MultiParam<T> &file) {
    return T::constFirst + ":" + SSTR(file.values.first) + "," + T::constSecond + ":" + SSTR(file.values.second);
}

template<> const int NuclAA<int>::max(INT_MAX);
template<> const std::string NuclAA<int>::constFirst = "aa";
template<> const std::string NuclAA<int>::constSecond = "nucl";
template class MultiParam<NuclAA<int>>;

template<> const float NuclAA<float>::max(FLT_MAX);
template<> const std::string NuclAA<float>::constFirst = "aa";
template<> const std::string NuclAA<float>::constSecond = "nucl";
template class MultiParam<NuclAA<float>>;

template<> const std::string NuclAA<std::string>::max("INVALID");
template<> const std::string NuclAA<std::string>::constFirst = "aa";
template<> const std::string NuclAA<std::string>::constSecond = "nucl";
template class MultiParam<NuclAA<std::string>>;

template<> const int SeqProf<int>::max(INT_MAX);
template<> const std::string SeqProf<int>::constFirst = "seq";
template<> const std::string SeqProf<int>::constSecond = "prof";
template class MultiParam<SeqProf<int>>;

const float PseudoCounts::max(FLT_MAX);
const std::string PseudoCounts::constFirst = "substitution";
const std::string PseudoCounts::constSecond = "context";
template class MultiParam<PseudoCounts>;
