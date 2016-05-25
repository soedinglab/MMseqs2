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
        Debug(Debug::WARNING) << "Warning: Calling program has OMP_PROC_BIND set in its environment. Performance may be degraded!\n";
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

void CommandCaller::execProgram(const char* program, size_t argc, const char **argv) {
    // hack: our argv string does not contain a program name anymore, readd it
    const char **pArgv = new const char*[argc + 2];
    pArgv[0] = program;
    for (size_t i = 0; i < argc; ++i) {
        pArgv[i + 1] = argv[i];
    }
    pArgv[argc + 1] = NULL;

    int res = execvp(program, (char * const *) pArgv);

    // should not be reached in the normal case
    delete[] pArgv;
    EXIT(res);
}
