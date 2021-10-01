#include "Util.h"
#include "Parameters.h"
#include "Debug.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include "DownloadDatabase.h"
#include <cassert>

extern std::vector<DatabaseDownload> downloads;

const int PAD_LEFT = 0;
const int PAD_RIGHT = 1;
void appendPadded(std::string& dst, const std::string& value, size_t n, int direction = PAD_LEFT, char padding = ' ') {
    if (n < value.size()) {
        dst.append(value);
        return;
    }
    if (direction == PAD_RIGHT) {
        dst.append(n - value.size(), padding);
    }
    dst.append(value);
    if (direction == PAD_LEFT) {
        dst.append(n - value.size(), padding);
    }
}

std::string listDatabases(const Command &command, bool detailed) {
    size_t nameWidth = 4, urlWidth = 3, dbTypeWidth = 4;
    for (size_t i = 0; i < downloads.size(); ++i) {
        nameWidth = std::max(nameWidth, strlen(downloads[i].name));
        urlWidth = std::max(urlWidth, strlen(downloads[i].url));
        dbTypeWidth = std::max(dbTypeWidth, strlen(Parameters::getDbTypeName(downloads[i].dbType)));
    }

    std::string description;
    description.reserve(1024);
    if (detailed) {
        description += " By ";
        description += command.author;
        description += "\n";
    }

    description += "\n  ";
    appendPadded(description, "Name", nameWidth);
    description.append(1, '\t');
    appendPadded(description, "Type", dbTypeWidth);
    description.append(1, '\t');
    appendPadded(description, "Taxonomy", 8);
    description.append(1, '\t');
    appendPadded(description, "Url", urlWidth);
    description.append(1, '\n');

    for (size_t i = 0; i < downloads.size(); ++i) {
        description.append("- ");
        appendPadded(description, downloads[i].name, nameWidth);
        description.append(1, '\t');
        appendPadded(description, Parameters::getDbTypeName(downloads[i].dbType), dbTypeWidth);
        description.append(1, '\t');
        appendPadded(description, (downloads[i].hasTaxonomy ? "yes" : "-"), 8, PAD_RIGHT);
        description.append(1, '\t');
        // last field in line should not be padded
        //appendPadded(description, downloads[i].url, urlWidth);
        description.append(downloads[i].url);
        description.append(1, '\n');
        if (detailed) {
            if (strlen(downloads[i].description) > 0) {
                description.append(2, ' ');
                description.append(downloads[i].description);
                description.append(1, '\n');
            }
            if (strlen(downloads[i].citation) > 0) {
                description.append("  Cite: ");
                description.append(downloads[i].citation);
                description.append(1, '\n');
            }
            description.append(1, '\n');
        }
    }

    return description;
}

int databases(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    std::string description = listDatabases(command, par.help);
    if (par.filenames.size() == 0 || par.help) {
        par.printUsageMessage(command, par.help ? MMseqsParameter::COMMAND_EXPERT : 0, description.c_str());
        EXIT(EXIT_SUCCESS);
    }

    ssize_t downloadIdx = -1;
    for (size_t i = 0; i < downloads.size(); ++i) {
        if (par.db1 == std::string(downloads[i].name)) {
            downloadIdx = i;
            break;
        }
    }
    if (downloadIdx == -1) {
        par.printUsageMessage(command, par.help ? MMseqsParameter::COMMAND_EXPERT : 0, description.c_str());
        Debug(Debug::ERROR) << "Selected database " << par.db1 << " was not found\n";
        EXIT(EXIT_FAILURE);
    }
    par.printParameters(command.cmd, argc, argv, par.databases);
    std::string tmpDir = par.db3;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.databases));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    for (size_t i = 0; i < downloads[downloadIdx].environment.size(); ++i) {
        cmd.addVariable(downloads[downloadIdx].environment[i].key, downloads[downloadIdx].environment[i].value);
    }
    cmd.addVariable("TAXONOMY", downloads[downloadIdx].hasTaxonomy ? "TRUE" : NULL);
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("VERB_PAR", par.createParameterString(par.onlyverbosity).c_str());
    cmd.addVariable("COMP_PAR", par.createParameterString(par.verbandcompression).c_str());
    // aria2c gives an (undocumented error with more than 16 connections)
    cmd.addVariable("ARIA_NUM_CONN", SSTR(std::min(16, par.threads)).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("THREADS_COMP_PAR", par.createParameterString(par.threadsandcompression).c_str());
    std::string program = tmpDir + "/download.sh";
    FileUtil::writeFile(program, downloads[downloadIdx].script, downloads[downloadIdx].scriptLength);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    EXIT(EXIT_FAILURE);
}
