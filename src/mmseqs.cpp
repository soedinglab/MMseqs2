#include "Command.h"
#include "DownloadDatabase.h"
#include "Prefiltering.h"
#include "Parameters.h"

const char* binary_name = "mmseqs";
const char* tool_name = "MMseqs2";
const char* tool_introduction = "MMseqs2 (Many against Many sequence searching) is an open-source software suite for very fast, \nparallelized protein sequence searches and clustering of huge protein sequence data sets.\n\nPlease cite: M. Steinegger and J. Soding. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, doi:10.1038/nbt.3988 (2017).";
const char* main_author = "Martin Steinegger (martin.steinegger@snu.ac.kr)";
const char* show_extended_help = "1";
const char* show_bash_info = "1";
extern const char* MMSEQS_CURRENT_INDEX_VERSION;
const char* index_version_compatible = MMSEQS_CURRENT_INDEX_VERSION;
bool hide_base_commands = false;
void (*validatorUpdate)(void) = 0;

extern std::vector<Command> baseCommands;
void init() {
    registerCommands(&baseCommands);
}
void (*initCommands)(void) = init;

DEFAULT_PARAMETER_SINGLETON_INIT

std::vector<DatabaseDownload> externalDownloads = {};
std::vector<KmerThreshold> externalThreshold = {};

bool hide_base_downloads = false;
