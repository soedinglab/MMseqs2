//
// Created by mad on 5/9/15.
//

#include "CommandCaller.h"

#include "iostream"
#include "Util.h"
#include <cstdlib>
#include <sstream>


CommandCaller::CommandCaller() {

}

void CommandCaller::addVariable(std::string key, std::string value) {
    setenv(key.c_str(), value.c_str(), true);
}

int CommandCaller::callProgram(std::string program,const char **argv, size_t argc) {
    std::stringstream argString;
    argString << program << " ";
    for(size_t i = 0; i < argc; i++){
        argString << argv[i] << " ";
    };
    if(std::system(argString.str().c_str())){
        EXIT(EXIT_FAILURE);
    }
    return 0;
}
