#ifndef MMSEQS_COMMANDCALLER_H
#define MMSEQS_COMMANDCALLER_H

#include <cstddef>

class CommandCaller {
public:
    CommandCaller();

    void addVariable(const char* key, const char* value);

    int callProgram(const char* program, size_t argc, const char **argv);

    // Does not return on success
    void execProgram(const char* program, size_t argc, const char **argv);
};

#endif //MMSEQS_COMMANDCALLER_H
