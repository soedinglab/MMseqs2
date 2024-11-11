#ifndef TARGET_SUBJECT_IDS_CUH
#define TARGET_SUBJECT_IDS_CUH

#include "config.hpp"

#include <algorithm>
#include <string>
#include <fstream>
#include <vector>

namespace cudasw4{

struct TargetSubjectIds{
    std::vector<ReferenceIdT> subjectIds;

    TargetSubjectIds() = default;
    TargetSubjectIds(const std::string& filename){
        std::ifstream is(filename);
        if(!is){
            throw std::runtime_error("File " + filename + " could not be opened");
        }
        std::string line;
        while(std::getline(is, line)){
            try{
                std::int64_t val = std::stoull(line);
                if(0 <= val && val < std::int64_t(MaxSequencesInDB::value())){
                    subjectIds.push_back(val);
                }else{
                    std::cerr << "Invalid reference id '" << line << "'. Skipping line...\n";
                }
            }catch(std::invalid_argument&){
                std::cerr << "Could not convert '" << line << "' to number. Skipping line...\n";
            }
        }

        std::sort(subjectIds.begin(), subjectIds.end());
    }

    TargetSubjectIds(std::vector<ReferenceIdT> ids) : subjectIds(std::move(ids)){
        std::sort(subjectIds.begin(), subjectIds.end());
    }

    void removeOutOfBoundsTargets(size_t databaseSize){
        subjectIds.erase(
            std::remove_if(subjectIds.begin(), subjectIds.end(), [&](size_t id){return id >= databaseSize; }),
            subjectIds.end()
        );
    }

    size_t size() const{
        return subjectIds.size();
    }

    auto begin() const{
        return subjectIds.begin();
    }

    auto end() const{
        return subjectIds.end();
    }
};

}







#endif