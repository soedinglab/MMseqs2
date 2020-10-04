#include "MultiParam.h"
template <> const int NuclAA<int>::max(INT_MAX);
template <> const std::string NuclAA<int>::constFirst = "nucl";
template <> const std::string NuclAA<int>::parseStr = "%d";
template <> const std::string NuclAA<int>::constSecond = "aa";

template <> const float NuclAA<float>::max(FLT_MAX);
template <> const std::string NuclAA<float>::constFirst = "nucl";
template <> const std::string NuclAA<float>::parseStr = "%f";
template <> const std::string NuclAA<float>::constSecond = "aa";

template <> const std::string NuclAA<std::string>::max("INVALID");
template <> const std::string NuclAA<std::string>::constFirst = "nucl";
template <> const std::string NuclAA<std::string>::parseStr = "%s";
template <> const std::string NuclAA<std::string>::constSecond = "aa";

const float PseudoCounts::max(FLT_MAX);
const std::string PseudoCounts::constFirst = "substitution";
const std::string PseudoCounts::parseStr = "%f";
const std::string PseudoCounts::constSecond = "context";



template <typename T>
MultiParam<T>::MultiParam(const char* parametercstring){
    if (strchr(parametercstring, ',') != NULL) {
        std::vector<std::string> params = Util::split(parametercstring,",");
        if(params.size() == 2){
            for(size_t i = 0; i < params.size();i++) {
                if (params[i].rfind(T::constFirst +":", 0) == 0) {
                    std::vector<std::string> pair = Util::split(params[i],":");
                    values.first = T::max;
                    if(pair.size() == 2) {
                        if(assign(pair[1], values.first) != 1){
                            values.first = T::max;
                        }
                    }
                }
                if (params[i].rfind(T::constSecond +":", 0) == 0) {
                    std::vector<std::string> pair = Util::split(params[i],":");
                    values.second = T::max;
                    if(pair.size() == 2) {
                        if(assign(pair[1], values.second)!=1){
                            values.second = T::max;
                        }
                    }
                }
            }
        } else {
            values.first = T::max;
            values.second = T::max;
        }
    } else {
        if (sscanf(parametercstring, T::parseStr.c_str(), &values.second) != 1) {
            values.first = T::max;
            values.second = T::max;
        }
        else
            values.first = values.second;
    }
}

template class MultiParam<NuclAA<int>>;
template class MultiParam<NuclAA<float>>;
template class MultiParam<NuclAA<std::string>>;
template class MultiParam<PseudoCounts>;