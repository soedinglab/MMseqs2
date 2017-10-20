#include "Debug.h"
#include "Util.h"
#include "Command.h"
#include "DistanceCalculator.h"

#include <CpuInfo.h>
#include <iomanip>

extern const char* tool_name;
extern const char* tool_introduction;
extern const char* main_author;
extern const char* version;
extern std::vector<struct Command> commands;
extern std::vector<Categories> categories;

void checkCpu() {
    CpuInfo info;
    if(info.HW_x64 == false) {
        Debug(Debug::ERROR) << "64 bit system is required to run MMseqs.\n";
        EXIT(EXIT_FAILURE);
    }
#ifdef SEE
    if(info.HW_SSE41 == false) {
        Debug(Debug::ERROR) << "SSE4.1 is required to run MMseqs.\n";
        EXIT(EXIT_FAILURE);
    }
#endif
#ifdef AVX2
    if(info.HW_AVX2 == false){
        Debug(Debug::ERROR) << "Your machine does not support AVX2.\n";
        if(info.HW_SSE41 == true) {
            Debug(Debug::ERROR) << "Please compile with SSE4.1 cmake -DHAVE_SSE4_1=1 \n";
        }else{
            Debug(Debug::ERROR) << "SSE 4.1 is the minimum requirement to run MMseqs.\n";
        }
        EXIT(EXIT_FAILURE);
    }
#endif
}

int getCommandIndex(const char *s) {
    for (size_t i = 0; i < commands.size(); i++) {
        struct Command &p = commands[i];
        if (!strcmp(s, p.cmd))
            return i;
    }
    return -1;
}

int runCommand(const Command &p, int argc, const char **argv) {
    return p.commandFunction(argc, argv, p);
}

void printUsage() {
    std::stringstream usage;
    usage << tool_introduction << "\n\n";
    usage << tool_name << " Version: " << version << "\n";
    usage << "© " << main_author << "\n";

    for(size_t i = 0; i < categories.size(); ++i) {
        usage << "\n" << std::setw(20) << categories[i].title << "\n";
        for (size_t j = 0; j < commands.size(); j++) {
            struct Command &p = commands[j];
            if (p.mode == categories[i].mode) {
                usage << std::left << std::setw(20) << "  " + std::string(p.cmd) << "\t" << p.shortDescription << "\n";
            }
        }
    }

    usage << "\nBash completion for tools and parameters can be installed by adding \"source path/to/mmseqs/util/bash-completion.sh\" to your \"$HOME/.bash_profile\".\n"
            "Include the location of the " << tool_name << " binary is in your \"$PATH\" environment variable.";

    Debug(Debug::INFO) << usage.str() << "\n";
}

int main(int argc, const char **argv) {
    checkCpu();
    if (argc < 2) {
        printUsage();
        EXIT(EXIT_FAILURE);
    }
    setenv("MMSEQS", argv[0], true);
    int i;
    if ((i = getCommandIndex(argv[1])) != -1) {
        const struct Command &p = commands[i];
        EXIT(runCommand(p, argc - 2, argv + 2));
    } else {
        printUsage();
        Debug(Debug::ERROR) << "\nInvalid Command: " << argv[1] << "\n";

        // Suggest some command that the user might have meant
        size_t index = SIZE_MAX;
        size_t minDistance = SIZE_MAX;
        for (size_t i = 0; i < commands.size(); ++i) {
            struct Command &p = commands[i];
            if(p.mode == COMMAND_HIDDEN)
                continue;

            size_t distance = DistanceCalculator::levenshteinDistance(argv[1], p.cmd);
            if(distance < minDistance) {
                minDistance = distance;
                index = i;
            }
        }

        if(index != SIZE_MAX) {
            Debug(Debug::ERROR) << "Did you mean \"" <<  argv[0] << " " << commands[index].cmd << "\"?\n";
        }

        EXIT(EXIT_FAILURE);
    }

    return 0;
}
