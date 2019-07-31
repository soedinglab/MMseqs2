#include "Command.h"

const char* binary_name = "mmseqs";
const char* tool_name = "MMseqs2";
const char* tool_introduction = "MMseqs2 (Many against Many sequence searching) is an open-source software suite for very fast, \nparallelized protein sequence searches and clustering of huge protein sequence data sets.\n\nPlease cite: M. Steinegger and J. Soding. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, doi:10.1038/nbt.3988 (2017).";
const char* main_author = "Martin Steinegger (martin.steinegger@mpibpc.mpg.de)";
const char* show_extended_help = "1";
const char* show_bash_info = "1";
bool hide_base_commands = false;
std::vector<Command> commands = {};
