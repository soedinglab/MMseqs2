#include "Command.h"
#include "Debug.h"

extern const char* version;

int versionstring(int argc, const char **argv, const Command& command) {
    Debug(Debug::INFO) << version << "\n";
}
