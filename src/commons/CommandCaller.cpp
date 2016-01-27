#include "CommandCaller.h"
#include "Util.h"

#include <cstdlib>
#include <unistd.h>
#include <sstream>


void CommandCaller::addVariable(std::string key, std::string value) {
    setenv(key.c_str(), value.c_str(), true);
}

int CommandCaller::callProgram(std::string program, size_t argc, const char **argv) {
    std::stringstream argString(program);
    for (size_t i = 0; i < argc; i++) {
        argString << " " << argv[i];
    }

    if (std::system(argString.str().c_str()) != EXIT_SUCCESS) {
        EXIT(EXIT_FAILURE);
    }

    return 0;
}

void CommandCaller::execProgram(std::string program, size_t argc, const char **argv) {
    char *name = strdup(program.c_str());

    // hack: our argv string does not contain a program name anymore, readd it
    char *pArgv[argc + 2];
    pArgv[0] = name;
    for (size_t i = 0; i < argc; ++i) {
        pArgv[i + 1] = (char *) argv[i];
    }
    pArgv[argc + 1] = NULL;

    int res = execvp(program.c_str(), pArgv);

    // should not be reached in the normal case
    free(name);
    EXIT(res);
}
