#include "Debug.h"
#include "Util.h"
#include "Command.h"
#include "DistanceCalculator.h"
#include "FileUtil.h"
#include "Timer.h"

#include <iomanip>

extern const char *binary_name;
extern const char *tool_name;
extern const char *tool_introduction;
extern const char *main_author;
extern const char *version;
extern const char *show_extended_help;
extern const char *show_bash_info;
extern bool hide_base_commands;

extern std::vector<std::vector<Command>*> commands;
extern std::vector<Categories> categories;
extern void (*validatorUpdate)(void);
extern void (*initCommands)(void);

const Command* getCommandByName(const char *s) {
    // allow base commands also to be called with a prefix, e.g. "mmseqs base:createdb"
    // this allows inheriting programs to find shadowed base modules
    const char *prefix = "base:";
    const char *check = strncmp(s, prefix, strlen(prefix)) == 0 ? s + strlen(prefix) : s;

    for (size_t i = commands.size() - 1; true; i--) {
        for (size_t j = 0; j < commands[i]->size(); j++) {
            const Command &p = (*commands[i])[j];
            if (!strcmp(s, p.cmd))
                return &p;
            if (i == 0 && !strcmp(check, p.cmd))
                return &p;
        }
        if (i == 0) {
            break;
        }
    }
    return NULL;
}

int runCommand(const Command *p, int argc, const char **argv) {
    Timer timer;
    int status = p->commandFunction(argc, argv, *p);
    Debug(Debug::INFO) << "Time for processing: " << timer.lap() << "\n";
    return status;
}

void printUsage(bool showExtended) {
    std::stringstream usage;

    usage << tool_introduction << "\n\n";
    usage << tool_name << " Version: " << version << "\n";
    usage << "© " << main_author << "\n\n";
    usage << "usage: " << binary_name << " <command> [<args>]" << "\n";

    std::vector<int> showCategoryHeader(categories.size(), 0);
    for (size_t i = 0; i < categories.size(); ++i) {
        size_t start = 0;
        if (hide_base_commands) {
            start = commands.size() - 1;
        }
        for (size_t j = commands.size() - 1; j >= start; j--) {
            for (size_t k = 0; k < commands[j]->size(); k++) {
                const Command &p = (*commands[j])[k];
                if (p.mode & categories[i].mode) {
                    showCategoryHeader[i] = 1;
                    break;
                }
            }
            if (j == 0) {
                break;
            }
        }
    }


    for (size_t i = 0; i < categories.size(); ++i) {
        if (showExtended == false
            && (categories[i].mode & COMMAND_MAIN) == 0
            && (categories[i].mode & COMMAND_EASY) == 0
            && (categories[i].mode & COMMAND_DATABASE_CREATION) == 0
            && (categories[i].mode & COMMAND_FORMAT_CONVERSION) == 0
            ) {
            continue;
        }

        if (showCategoryHeader[i] == 0) {
            continue;
        }

        usage << "\n" << std::setw(20) << categories[i].title << "\n";
        size_t start = 0;
        if (hide_base_commands) {
            start = commands.size() - 1;
        }
        for (size_t j = commands.size() - 1; j >= start; j--) {
            for (size_t k = 0; k < commands[j]->size(); k++) {
                const Command &p = (*commands[j])[k];
                if (showExtended == false && (p.mode & COMMAND_EXPERT) != 0) {
                    continue;
                }
                if (p.mode & categories[i].mode) {
                    usage << std::left << std::setw(20) << "  " + std::string(p.cmd) << "\t" << p.description << "\n";
                }
            }
            if (j == 0) {
                break;
            }
        }
    }

    if (show_extended_help != NULL) {
        if (showExtended == false) {
            usage << "\nAn extended list of all modules can be obtained by calling '" << binary_name << " -h'.\n";
        }
    }
    if (show_bash_info != NULL) {
        usage  << "\nBash completion for modules and parameters can be installed by adding \"source MMSEQS_HOME/util/bash-completion.sh\" to your \"$HOME/.bash_profile\".\nInclude the location of the " << tool_name << " binary in your \"$PATH\" environment variable.";
    }
    Debug(Debug::INFO) << usage.str() << "\n";
}

int shellcompletion(int argc, const char **argv) {
    // mmseqs programs
    if (argc == 0) {
        size_t start = 0;
        if (hide_base_commands) {
            start = commands.size() - 1;
        }
        for (size_t i = commands.size() - 1; i >= start; i--) {
            for (size_t j = 0; j < commands[i]->size(); j++) {
                const Command &p = (*commands[i])[j];
                if (p.mode & COMMAND_HIDDEN)
                    continue;
                Debug(Debug::INFO) << p.cmd << " ";
            }
            if (i == 0) {
                break;
            }
        }
        Debug(Debug::INFO) << "\n";
    }

    // mmseqs parameters for given program
    if (argc == 1) {
        size_t start = 0;
        if (hide_base_commands) {
            start = commands.size() - 1;
        }
        for (size_t i = commands.size() - 1; i >= start; i--) {
            for (size_t j = 0; j < commands[i]->size(); j++) {
                const Command &p = (*commands[i])[j];
                if (strcmp(p.cmd, argv[0]) != 0) {
                    continue;
                }
                if (p.params == NULL) {
                    continue;
                }
                for (std::vector<MMseqsParameter *>::const_iterator it = p.params->begin(); it != p.params->end(); ++it) {
                    Debug(Debug::INFO) << (*it)->name << " ";
                }
                Debug(Debug::INFO) << "\n";
                break;
            }
            if (i == 0) {
                break;
            }
        }
        Debug(Debug::INFO) << "\n";
    }
    return EXIT_SUCCESS;
}

int main(int argc, const char **argv) {
    if (initCommands != NULL) {
        initCommands();
    }

    if (argc < 2) {
        printUsage(false);
        return EXIT_SUCCESS;
    }

    if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0){
        printUsage(true);
        return EXIT_SUCCESS;
    }

    if(validatorUpdate != NULL){
        (*validatorUpdate)();
    }
    FileUtil::fixRlimitNoFile();

    setenv("MMSEQS", argv[0], true);
    const Command *c = NULL;
    if (strncmp(argv[1], "shellcompletion", strlen("shellcompletion")) == 0) {
        return shellcompletion(argc - 2, argv + 2);
    } else if ((c = getCommandByName(argv[1])) != NULL) {
        EXIT(runCommand(c, argc - 2, argv + 2));
    } else {
        printUsage(true);
        Debug(Debug::INFO) << "\nInvalid Command: " << argv[1] << "\n";

        // Suggest some command that the user might have meant
        size_t indexI = SIZE_MAX;
        size_t indexJ = SIZE_MAX;
        int maxDistance = 0;
        size_t start = 0;
        if (hide_base_commands) {
            start = commands.size() - 1;
        }
        for (size_t i = commands.size() - 1; i >= start; i--) {
            for (size_t j = 0; j < commands[i]->size(); j++) {
                const Command &p = (*commands[i])[j];
                if (p.mode & COMMAND_HIDDEN) {
                    continue;
                }

                int distance = DistanceCalculator::localLevenshteinDistance(argv[1], p.cmd);
                if (distance > maxDistance) {
                    maxDistance = distance;
                    indexI = i;
                    indexJ = j;
                }
            }
            if (i == 0) {
                break;
            }
        }

        if (indexI != SIZE_MAX && indexJ != SIZE_MAX) {
            Debug(Debug::WARNING) << "Did you mean \"" << argv[0] << " " << (*commands[indexI])[indexJ].cmd << "\"?\n";
        }

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
