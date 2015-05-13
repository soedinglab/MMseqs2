//
// Created by mad on 5/9/15.
//

#ifndef MMSEQS_COMMANDCALLER_H
#define MMSEQS_COMMANDCALLER_H

#include <string>
#include <map>

class CommandCaller {
public:
    CommandCaller();
    void addVariable(std::string key, std::string value);
    int callProgram(std::string program, const char **argv, size_t argc);

private:
    std::map<std::string, std::string> variables;
};


#endif //MMSEQS_COMMANDCALLER_H
