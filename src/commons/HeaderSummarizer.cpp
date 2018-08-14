#include "HeaderSummarizer.h"
#include "Util.h"
#include "Debug.h"
#include "PatternCompiler.h"

#include <set>
#include <algorithm>

struct UniprotHeader {
    std::string dbType;
    std::string identifier;
    std::string proteinName;
    std::string organismName;
    unsigned int existence;
    unsigned int priority;

    UniprotHeader(const std::string& dbType,
                  const std::string& identifier,
                  const std::string& proteinName,
                  const std::string& organismName,
                  unsigned int existence
    ) : dbType(dbType),
        identifier(identifier),
        proteinName(proteinName),
        organismName(organismName),
        existence(existence) {
        updatePriority();
    };

    PatternCompiler& isUninformative() {
        static PatternCompiler uninformative("hypothetical|unknown|putative|predicted|unnamed|probable|partial|possible|uncharacterized|fragment");
        return uninformative;
    }

    void updatePriority() {
        priority = 0;

        if(isUninformative().isMatch(identifier.c_str()))
            return;

        if(dbType == "sp") {
            priority = 4;
        } else if (dbType == "tr") {
            priority = 1;
        }

        priority += (std::min(existence, 5u) - 5u);
    }

    friend bool operator<(const UniprotHeader& h1, const UniprotHeader& h2)
    {
        return h1.priority < h2.priority;
    }
};

struct MetaclustHeader {
    std::string dbType;
    std::string identifier;
    int priority;

    MetaclustHeader(const std::string& dbType,
                  const std::string& identifier
    ) : dbType(dbType),
        identifier(identifier)
    {
        updatePriority();
    };

    void updatePriority() {
        priority = 0;

        if(dbType == "UPI") {
            priority = 4;
        } else {
            priority = 1;
        }
    }

    friend bool operator<(const MetaclustHeader& h1, const MetaclustHeader& h2)
    {
        return h1.priority < h2.priority;
    }
};


std::string UniprotHeaderSummarizer::summarize(const std::vector<std::string>& headers) {
    std::vector<UniprotHeader> headerQueue;

    std::string representingIdentifier;

    unsigned int clusterMembers = static_cast<unsigned int>(headers.size());

    for(std::vector<std::string>::const_iterator it = headers.begin(); it != headers.end(); ++it) {
        const std::string& header = *it;

        size_t start = 0;
        size_t end = header.find('|');
        if(end == std::string::npos)
            continue;
        std::string dbType = header.substr(start, end);

        start = end + 1;
        end = header.find('|', start);
        if(end == std::string::npos)
            continue;
        std::string identifier = header.substr(start, end - start);
        if(it == headers.begin()) {
            representingIdentifier = identifier;
        }

        start = header.find(' ', end);
        if(start == std::string::npos)
            continue;
        start++;

        end = header.find(" OS=", start);
        if(end == std::string::npos)
            continue;
        std::string proteinName = header.substr(start, end - start);

        start = header.find('=', end);
        if(start == std::string::npos)
            continue;
        start++;

        end = header.find(" GN=", start);
        if(end == std::string::npos) {
            end = header.find(" PE=", start);
            if(end == std::string::npos) {
                continue;
            }
        }
        std::string organismName = header.substr(start, end - start);

        start = header.find("PE=", end);
        if(start == std::string::npos)
            continue;
        start += 3;
        end = header.find(" SV=", start);
        if(end == std::string::npos)
            continue;

        std::string existenceString = header.substr(start, end - start);
        unsigned int existence = static_cast<unsigned int>(strtoul(existenceString.c_str(), NULL, 10));

        headerQueue.emplace_back(dbType, identifier, proteinName, organismName, existence);
    }

    std::make_heap(headerQueue.begin(), headerQueue.end());

    const unsigned int maxDescriptions = 5;

    std::ostringstream summarizedHeader;
    summarizedHeader << "Representative=" << representingIdentifier.c_str();
    summarizedHeader << " n=" << clusterMembers;

    std::set<std::string> usedDescriptions;
    summarizedHeader << " Descriptions=[";
    unsigned int descriptionCount = 0;
    for (std::vector<UniprotHeader>::const_iterator it = headerQueue.begin(); it != headerQueue.end(); ++it) {
        if (descriptionCount > maxDescriptions)
            break;

        const UniprotHeader& header = *it;

        if (usedDescriptions.find(header.proteinName) != usedDescriptions.end()) {
            continue;
        }

        summarizedHeader << header.proteinName;

        usedDescriptions.emplace(header.proteinName);

        descriptionCount++;

        if(Util::isLastIterator(it, headerQueue) == false && descriptionCount <= maxDescriptions) {
            summarizedHeader << "|";
        }
    }
    summarizedHeader << "]";

    summarizedHeader<< " Members=";
    for (std::vector<UniprotHeader>::const_iterator it = headerQueue.begin(); it != headerQueue.end(); ++it) {
        const UniprotHeader& header = *it;
        summarizedHeader << header.identifier;
        if(Util::isLastIterator(it, headerQueue) == false) {
            summarizedHeader << ",";
        }
    }
    summarizedHeader << "\n";

    return summarizedHeader.str();
}



std::string MetaclustHeaderSummarizer::summarize(const std::vector<std::string>& headers) {
    std::vector<MetaclustHeader> headerQueue;

    std::string representingIdentifier;

    unsigned int clusterMembers = static_cast<unsigned int>(headers.size());

    for(std::vector<std::string>::const_iterator it = headers.begin(); it != headers.end(); ++it) {
        const std::string& header = *it;

        size_t start = 0;
        //UPI0008DB4360
        size_t end = header.find("UPI");
        std::string dbType = "lessImportant";
        if(end != std::string::npos){
            dbType = "UPI";
        }

        end = header.find(' ', start);
        if(end == std::string::npos)
            continue;
        std::string identifier = header.substr(0, end);
        if(it == headers.begin()) {
            representingIdentifier = identifier;
        }
        headerQueue.emplace_back(dbType, identifier);
    }

    std::make_heap(headerQueue.begin(), headerQueue.end());

    std::ostringstream summarizedHeader;
    summarizedHeader << "Representative=" << representingIdentifier.c_str();
    summarizedHeader << " n=" << clusterMembers;

    summarizedHeader<< " Members=";
    // +1 to skip representative
    for (std::vector<MetaclustHeader>::const_iterator it = headerQueue.begin(); it != headerQueue.end(); ++it) {
        const MetaclustHeader& header = *it;
        if(header.identifier.compare(representingIdentifier) == 0){
            continue;
        }
        summarizedHeader << header.identifier;
        if(Util::isLastIterator(it, headerQueue) == false) {
            summarizedHeader << ",";
        }
    }
    std::string header = summarizedHeader.str();
    if (header[header.size() -1] == ','){
        header[header.size() -1] = '\n';
    }else{
        header.push_back('\n');
    }

    return header;
}
