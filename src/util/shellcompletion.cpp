#include "Debug.h"
#include "Parameters.h"
#include "Util.h"
#include "Command.h"

#include <vector>

extern std::vector<struct Command> commands;

int shellcompletion(int argc, const char** argv, const Command&) {
    // mmseqs programs
    if(argc == 0) {
        for (size_t i = 0; i < commands.size(); i++) {
            struct Command &p = commands[i];
            if(p.mode == COMMAND_HIDDEN)
                continue;
            Debug(Debug::INFO) << p.cmd << " ";
        }
        Debug(Debug::INFO) << "\n";
    }

    // mmseqs parameters for given program
    if(argc == 1) {
        for (size_t i = 0; i < commands.size(); i++) {
            struct Command &p = commands[i];
            if(strcmp(p.cmd, argv[0]) != 0) {
                continue;
            }
            if(p.params == NULL) {
                continue;
            }
            for(std::vector<MMseqsParameter *>::const_iterator it = p.params->begin(); it != p.params->end(); ++it) {
                Debug(Debug::INFO) << (*it)->name << " ";
            }
            Debug(Debug::INFO) << "\n";
            break;
        }
        Debug(Debug::INFO) << "\n";
    }

    EXIT(EXIT_SUCCESS);
}
