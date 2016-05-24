#ifndef MMSEQS_HEADERSUMMARIZER_H
#define MMSEQS_HEADERSUMMARIZER_H

#include <string>
#include <vector>

class HeaderSummarizer {
public:
    virtual std::string summarize(const std::vector<std::string>& headers, const std::string& prefix) = 0;
};

class UniprotHeaderSummarizer : public HeaderSummarizer {
public:
    std::string summarize(const std::vector<std::string>& headers, const std::string& prefix);
};


#endif //MMSEQS_HEADERSUMMARIZER_H
