#ifndef MMSEQS_COMMANDCALLER_H
#define MMSEQS_COMMANDCALLER_H

#include <string>

class CommandCaller {
public:
    void addVariable(std::string key, std::string value);

    int callProgram(std::string program, size_t argc, const char **argv);

    // Does not return on success
    void execProgram(std::string program, size_t argc, const char **argv);

    int spawnProgram(std::string program, size_t argc, const char **argv);
};

#endif //MMSEQS_COMMANDCALLER_H
