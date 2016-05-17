#ifndef MMSEQS_PATTERNCOMPILER_H
#define MMSEQS_PATTERNCOMPILER_H
#include <regex.h>

class PatternCompiler {
public:
    PatternCompiler(const char* pattern);
    ~PatternCompiler();

    bool isMatch(const char* string);
private:
    regex_t regex;
};


#endif //MMSEQS_PATTERNCOMPILER_H
