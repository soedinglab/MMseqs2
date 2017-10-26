#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"

#include <strings.h>
#include <cstdlib>
#include <unistd.h>
#include <sstream>

#ifdef OPENMP
#include <omp.h>
#endif


CommandCaller::CommandCaller() {
#ifdef OPENMP
#if _OPENMP >= 201307
    if(omp_get_proc_bind() != omp_proc_bind_false){
#else
    char* procBind = getenv("OMP_PROC_BIND");
    if(procBind != NULL && strcasecmp(procBind, "false") != 0  && strcasecmp(procBind, "0") != 0) {
#endif
        Debug(Debug::WARNING) << "Error: Calling program has OMP_PROC_BIND set in its environment. Please unset OMP_PROC_BIND.\n";
        EXIT(EXIT_FAILURE);
    }
#endif
}

void CommandCaller::addVariable(const char* key, const char* value) {
    setenv(key, value, true);
}

int CommandCaller::callProgram(const char* program, size_t argc, const char **argv) {
    std::stringstream argStream(program);
    for (size_t i = 0; i < argc; i++) {
        argStream << " " << argv[i];
    }

    std::string argString = argStream.str();
    if (std::system(argString.c_str()) != EXIT_SUCCESS) {
        EXIT(EXIT_FAILURE);
    }

    return 0;
}

void CommandCaller::execProgram(const char* program, const std::vector<std::string> &argv) {
    // hack: our argv string does not contain a program name anymore, readd it
    const char **pArgv = new const char*[argv.size() + 2];
    pArgv[0] = program;
    for (size_t i = 0; i < argv.size(); ++i) {
        pArgv[i + 1] = argv[i].c_str();
    }
    pArgv[argv.size() + 1] = NULL;

    int res = execvp(program, (char * const *) pArgv);

    if (res == -1) {
        Debug(Debug::ERROR) << "Failed to execute " << program << " with error " << errno << ".\n";
    }

    // should not be reached in the normal case
    delete[] pArgv;
    EXIT(EXIT_FAILURE);
}
