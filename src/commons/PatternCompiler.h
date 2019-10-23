#ifndef MMSEQS_PATTERNCOMPILER_H
#define MMSEQS_PATTERNCOMPILER_H
#include <regex.h>
#include "PatternCompiler.h"
#include "Debug.h"
#include "Util.h"


class PatternCompiler {
public:
    PatternCompiler(const char* pattern)  {
        if (regcomp(&regex, pattern, REG_EXTENDED | REG_NEWLINE) != 0 ){
            Debug(Debug::ERROR) << "Error in regex " << pattern << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

    ~PatternCompiler() {
        regfree(&regex);
    }

    bool isMatch(const char *target) {
        return regexec(&regex, target, 0, NULL, 0) == 0;
    }

    bool isMatch(const char *target, size_t nmatch, regmatch_t *pmatch) {
        return regexec(&regex, target, nmatch, pmatch, 0) == 0;
    }


private:
    regex_t regex;
};


#endif //MMSEQS_PATTERNCOMPILER_H
