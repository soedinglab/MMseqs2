#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"

#define _GNU_SOURCE

#include <sched.h>
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

bool CommandCaller::resetAffinity() {
#ifdef __linux__
    cpu_set_t *mask;
    size_t size;
    // 256 cores should be more than enough
    int nrcpus = 256;
    int i;

    // allow running on all processors
    mask = CPU_ALLOC(nrcpus);
    for (i = 0; i < nrcpus; i++)
        CPU_SET_S(i, size, mask);
    if (sched_setaffinity(0, size, mask) == -1) {
        return false;
    }
    CPU_FREE(mask);
#endif

    return true;
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

    if(!resetAffinity()) {
        Debug(Debug::WARNING) << "Could not reset CPU affinity. Performance of child programs might be degraded!\n";
    }

    int res = execvp(program.c_str(), pArgv);

    // should not be reached in the normal case
    free(name);
    EXIT(res);
}
