#include "PatternCompiler.h"
#include "Debug.h"
#include "Util.h"

PatternCompiler::PatternCompiler(const char* pattern)  {
    if (regcomp(&regex, pattern, REG_EXTENDED | REG_NEWLINE) != 0 ){
        Debug(Debug::ERROR) << "Error in regex " << pattern << "\n";
        EXIT(EXIT_FAILURE);
    }
}

PatternCompiler::~PatternCompiler() {
    regfree(&regex);
}

bool PatternCompiler::isMatch(const char *target) {
    return regexec(&regex, target, 0, NULL, 0) == 0;
}
