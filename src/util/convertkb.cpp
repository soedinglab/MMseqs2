#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "FileUtil.h"
#include "Log.h"
#include "Debug.h"

int convertkb(int argn, const char **argv)
{
    std::string usage;
    usage.append("Turns an UniprotKB file into a ffindex database.\n");
    usage.append("USAGE: <uniprotKB> <outDB>\n");
    usage.append("\nDesigned and implemented by Milot Mirdita <milot@mirdita.de>.\n");

    Parameters par;
    par.parseParameters(argn, argv, usage, par.onlyverbosity, 2);

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str());
    writer.open();

    FILE* kbFile = FileUtil::openFileOrDie(par.db1.c_str(), "r", true);

    bool isInEntry = false;

    std::string identifier;
    std::stringstream s;
    char * line = NULL;
    size_t n = 0;
    ssize_t read = 0;
    while ((read = getline(&line, &n, kbFile)) != -1) {
        if(read < 2) {
            Debug(Debug::WARNING) << "Invalid line" << "\n";
            continue;
        }

        if(strncmp("ID", line, 2) == 0) {
            isInEntry = true;
            std::vector<std::string> words = Util::split(line, " ");
            if(words.size() < 2) {
                Debug(Debug::WARNING) << "Invalid entry" << "\n";
                isInEntry = false;
                continue;
            }
            identifier = words[1];
        }

        if(isInEntry) {
            s << line;
        }

        if(strncmp("//", line, 2) == 0) {
            isInEntry = false;
            std::string entry = s.str();
            s.str("");
            s.clear();
            writer.write(entry.c_str(), entry.length(), identifier.c_str());
        }
    }

    fclose(kbFile);
    writer.close();
    return EXIT_SUCCESS;
}
