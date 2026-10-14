#include "MMseqsNode.h"
#include "Debug.h"
#include "FileUtil.h"
#include "Parameters.h"
#include "Util.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <unistd.h>

unsigned int MMseqsNode::rank = 0;
unsigned int MMseqsNode::count = 1;

void MMseqsNode::init() {
    const char *procId = getenv("SLURM_PROCID");
    const char *taskCount = getenv("SLURM_STEP_NUM_TASKS");
    if (procId != NULL && taskCount != NULL) {
        rank = strtoul(procId, NULL, 10);
        count = strtoul(taskCount, NULL, 10);
    }
    if (count == 0 || rank >= count) {
        Debug(Debug::ERROR) << "Invalid SLURM_PROCID " << rank << " for " << count << " tasks\n";
        EXIT(EXIT_FAILURE);
    }
    if (isMaster() == false) {
        Parameters &par = Parameters::getInstance();
        par.verbosity = Debug::ERROR;
        Debug::setDebugLevel(Debug::ERROR);
    }
}

void MMseqsNode::publishFile(const std::string &file, const std::string &content) {
    const std::string tmpFile = file + ".tmp";
    FILE *handle = FileUtil::openAndDelete(tmpFile.c_str(), "w");
    if (fwrite(content.c_str(), sizeof(char), content.size(), handle) != content.size()) {
        Debug(Debug::ERROR) << "Cannot write " << tmpFile << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fflush(handle) != 0 || fsync(fileno(handle)) != 0) {
        Debug(Debug::ERROR) << "Cannot flush " << tmpFile << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fclose(handle) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << tmpFile << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (rename(tmpFile.c_str(), file.c_str()) != 0) {
        Debug(Debug::ERROR) << "Cannot rename " << tmpFile << " to " << file << "\n";
        EXIT(EXIT_FAILURE);
    }
}

std::string MMseqsNode::waitForFile(const std::string &file) {
    while (FileUtil::fileExists(file.c_str()) == false) {
        sleep(1);
    }
    std::ifstream handle(file.c_str(), std::ios::binary);
    if (handle.is_open() == false) {
        Debug(Debug::ERROR) << "Cannot read " << file << "\n";
        EXIT(EXIT_FAILURE);
    }
    std::stringstream content;
    content << handle.rdbuf();
    return content.str();
}

std::string MMseqsNode::gatherFiles(const std::string &prefix) {
    std::string content;
    for (unsigned int node = 0; node < count; node++) {
        const std::string file = prefix + SSTR(node);
        content += waitForFile(file);
        FileUtil::remove(file.c_str());
    }
    return content;
}
