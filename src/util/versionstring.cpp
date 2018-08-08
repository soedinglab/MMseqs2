#include "Command.h"
#include "Debug.h"
#include "Util.h"

extern const char* version;

int versionstring(int, const char**, const Command&) {
    Debug(Debug::INFO) << version << "\n";
    EXIT(EXIT_SUCCESS);
}
