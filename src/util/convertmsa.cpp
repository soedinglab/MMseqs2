#include "Parameters.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "GzReader.h"

#include <algorithm>
#include <map>


int convertmsa(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    GzReader in(par.db1);
    if (in.fail()) {
        Debug(Debug::ERROR) << "File " << par.db1 << " not found!\n";
        return EXIT_FAILURE;
    }

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), 1, par.compressed, Parameters::DBTYPE_MSA_DB);
    writer.open();

    std::string line;
    unsigned int i = 0;
    bool inEntry = false;
    const char *buffer[255];
    std::vector<std::string> seqOrder;
    std::map<std::string, std::string> sequences;
    std::string identifier;
    std::string result;
    result.reserve(10 * 1024 * 1024);

    Debug::Progress progress;
    while (in.getline(line)) {
        size_t lineLength = line.length();
        if (lineLength < 1) {
            continue;
        }

        if (inEntry == false && line == "# STOCKHOLM 1.0") {
            progress.updateProgress();
            inEntry = true;
            continue;
        }

        if (inEntry == true  && line == "//") {
            inEntry = false;
            result.clear();
            size_t j = 0;
            for (std::vector<std::string>::const_iterator it = seqOrder.begin(); it != seqOrder.end(); ++it) {
                const std::string &accession = *it;
                const std::string &sequence = sequences[*it];
                result.append(">");
                if (j == 0 && identifier.length() > 0) {
                    result.append(identifier);
                    result.append(" ");
                }
                result.append(accession);
                result.append("\n");
                result.append(sequence);
                result.append("\n");
                j++;
            }

            writer.writeData(result.c_str(), result.length(), i++, 0);

            seqOrder.clear();
            sequences.clear();
            identifier = "";
            continue;
        }

        if (inEntry == false) {
            continue;
        }

        size_t columns = Util::getWordsOfLine(const_cast<char*>(line.c_str()), buffer, 255);
        if (line[0] == '#') {
            if (Util::startWith("#=GF", line)) {
                if (columns < 3) {
                    Debug(Debug::ERROR) << "Invalid annotation!\n";
                    inEntry = false;
                    continue;
                }

                if (par.identifierField == 1 && strncmp("AC", buffer[1], 2) == 0) {
                    const char *id = buffer[2];
                    size_t length = Util::skipNoneWhitespace(buffer[2]);
                    identifier = std::string(id, length);
                } else if (par.identifierField == 0 && strncmp("ID", buffer[1], 2) == 0) {
                    const char *id = buffer[2];
                    size_t length = Util::skipNoneWhitespace(buffer[2]);
                    identifier = std::string(id, length);
                }
            }
        } else {
            if (columns < 2) {
                Debug(Debug::ERROR) << "Invalid sequence!\n";
                inEntry = false;
                continue;
            }
            char accession[255];
            Util::parseKey(buffer[0], accession);

            const char* seq = buffer[1];
            size_t length = Util::skipNoneWhitespace(buffer[1]);

            std::map<std::string, std::string>::iterator it = sequences.find(accession);
            if (it == sequences.end()) {
                std::string sequence(seq, length);
                sequence.reserve(1024);
                std::replace(sequence.begin(), sequence.end(), '.', '-');
                sequences.emplace(accession, sequence);
                seqOrder.emplace_back(accession);
            } else {
                (*it).second.append(std::string(seq, length));
            }
        }
    }
    writer.close();

    return EXIT_SUCCESS;
}
