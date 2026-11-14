#ifndef EM_TARGET_GROUPS_H
#define EM_TARGET_GROUPS_H

// --target-groups for `mmseqs em --ref-group 1`: which targets of a target DB belong to one genome.
#include "Debug.h"
#include "FileUtil.h"
#include "KSeqWrapper.h"
#include "Util.h"

#include <dirent.h>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace EMTargetGroups {

const unsigned int NONE = std::numeric_limits<unsigned int>::max();

// calls f(key, name, file number) for every line of a DB's .lookup; false if there is no .lookup
template <typename F>
bool forEachLookup(const std::string &db, F f) {
    std::ifstream in(db + ".lookup");
    if (!in) {
        return false;
    }
    std::string line;
    while (std::getline(in, line)) {
        const size_t a = line.find('\t');
        if (a == std::string::npos) {
            continue;
        }
        const size_t b = line.find('\t', a + 1);
        const unsigned int file = b == std::string::npos ? 0u : static_cast<unsigned int>(std::strtoul(line.c_str() + b + 1, NULL, 10));
        f(static_cast<unsigned int>(std::strtoul(line.c_str(), NULL, 10)),
          line.substr(a + 1, b == std::string::npos ? std::string::npos : b - a - 1), file);
    }
    return true;
}

// a file name without compression and FASTA extensions; isFasta tells whether it had a FASTA extension
inline std::string stripFastaName(const std::string &path, bool *isFasta) {
    std::string b = FileUtil::baseName(path);
    const char *comp[] = {".gz", ".bz2", ".zst", ".xz"};
    for (size_t i = 0; i < sizeof(comp) / sizeof(comp[0]); ++i) {
        const std::string e(comp[i]);
        if (b.size() > e.size() && b.compare(b.size() - e.size(), e.size(), e) == 0) {
            b.erase(b.size() - e.size());
            break;
        }
    }
    const char *ext[] = {".fasta", ".fna", ".faa", ".ffn", ".frn", ".fas", ".fa", ".seq", ".fsa"};
    *isFasta = false;
    for (size_t i = 0; i < sizeof(ext) / sizeof(ext[0]); ++i) {
        const std::string e(ext[i]);
        if (b.size() > e.size() && b.compare(b.size() - e.size(), e.size(), e) == 0) {
            b.erase(b.size() - e.size());
            *isFasta = true;
            break;
        }
    }
    return b;
}

// the FASTA files a --target-groups value names: a directory (its FASTA files) or a text file listing FASTA paths;
// empty if it is neither
inline std::vector<std::string> genomeFastaFiles(const std::string &spec) {
    std::vector<std::string> files;
    bool isFasta;
    if (FileUtil::directoryExists(spec.c_str())) {
        DIR *dir = opendir(spec.c_str());
        if (dir != NULL) {
            struct dirent *e;
            while ((e = readdir(dir)) != NULL) {
                const std::string name(e->d_name);
                if (name.empty() || name[0] == '.') {
                    continue;
                }
                const std::string path = spec + "/" + name;
                stripFastaName(path, &isFasta);
                if (!FileUtil::directoryExists(path.c_str()) && isFasta) {
                    files.push_back(path);
                }
            }
            closedir(dir);
        }
        std::sort(files.begin(), files.end());
        return files;
    }
    std::ifstream in(spec);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line.find('\t') != std::string::npos || !FileUtil::fileExists(line.c_str())) {
            return std::vector<std::string>();          // a TSV, not a list of files
        }
        files.push_back(line);
    }
    return files;
}

struct Groups {
    std::vector<unsigned int> groupOfKey;     // target key -> genome, NONE = not assigned (its own group)
    std::vector<std::string> names;           // genome names
    size_t mapped = 0;                        // targets assigned to a genome
    std::string how;                          // where the genomes came from, for the log
};

// genomes of the targets of targetDb (keys 0..maxKey) as --target-groups spec defines them
inline Groups load(const std::string &targetDb, const std::string &spec, unsigned int maxKey) {
    Groups out;
    out.groupOfKey.assign(static_cast<size_t>(maxKey) + 1, NONE);
    std::unordered_map<std::string, unsigned int> index;
    auto genome = [&](const std::string &name) -> unsigned int {
        std::unordered_map<std::string, unsigned int>::const_iterator it = index.find(name);
        if (it != index.end()) {
            return it->second;
        }
        const unsigned int id = static_cast<unsigned int>(out.names.size());
        index.emplace(name, id);
        out.names.push_back(name);
        return id;
    };
    auto assign = [&](unsigned int key, unsigned int g) {
        if (key <= maxKey) {
            out.groupOfKey[key] = g;
            ++out.mapped;
        }
    };
    bool hasLookup = true;
    bool isFasta;
    if (spec == "source") {
        // createdb of several files: .lookup column 3 = file number, .source = file number -> file name
        std::unordered_map<unsigned int, std::string> genomeOfFileNumber;
        std::ifstream src(targetDb + ".source");
        if (!src) {
            Debug(Debug::ERROR) << targetDb << ".source is missing: --target-groups source needs a target DB made by createdb.\n";
            EXIT(EXIT_FAILURE);
        }
        std::string line;
        while (std::getline(src, line)) {
            const size_t t = line.find('\t');
            if (t != std::string::npos) {
                genomeOfFileNumber[static_cast<unsigned int>(std::strtoul(line.c_str(), NULL, 10))] = stripFastaName(line.substr(t + 1), &isFasta);
            }
        }
        if (genomeOfFileNumber.size() < 2) {
            Debug(Debug::WARNING) << targetDb << " was made from " << genomeOfFileNumber.size()
                                  << " input file(s): with --target-groups source every target of a file is one genome; "
                                  << "give createdb one FASTA per genome, or set --target-groups.\n";
        }
        hasLookup = forEachLookup(targetDb, [&](unsigned int key, const std::string &, unsigned int file) {
            std::unordered_map<unsigned int, std::string>::const_iterator it = genomeOfFileNumber.find(file);
            if (it != genomeOfFileNumber.end()) {
                assign(key, genome(it->second));
            }
        });
        out.how = "createdb input file (" + SSTR(genomeOfFileNumber.size()) + " files)";
    } else if (spec == "taxid") {
        std::unordered_map<unsigned int, std::string> taxOfKey;
        std::ifstream map(targetDb + "_mapping");
        if (!map) {
            Debug(Debug::ERROR) << targetDb << "_mapping is missing: --target-groups taxid needs a target DB annotated by createtaxdb.\n";
            EXIT(EXIT_FAILURE);
        }
        std::string line;
        while (std::getline(map, line)) {
            const size_t t = line.find('\t');
            if (t != std::string::npos) {
                taxOfKey[static_cast<unsigned int>(std::strtoul(line.c_str(), NULL, 10))] = line.substr(t + 1);
            }
        }
        hasLookup = forEachLookup(targetDb, [&](unsigned int key, const std::string &, unsigned int) {
            std::unordered_map<unsigned int, std::string>::const_iterator it = taxOfKey.find(key);
            if (it != taxOfKey.end()) {
                assign(key, genome(it->second));
            }
        });
        out.how = "taxonomy (" + targetDb + "_mapping)";
    } else if (spec.compare(0, 6, "field:") == 0) {
        // field:<separator>:<n>  genome = the first n fields of the identifier (a trailing :<m> is accepted and ignored)
        const size_t c1 = spec.find(':', 6);
        const std::string sep = spec.substr(6, c1 == std::string::npos ? std::string::npos : c1 - 6);
        const size_t c2 = c1 == std::string::npos ? std::string::npos : spec.find(':', c1 + 1);
        const int nG = c1 == std::string::npos ? 1 : std::atoi(spec.substr(c1 + 1, c2 == std::string::npos ? std::string::npos : c2 - c1 - 1).c_str());
        if (sep.empty() || nG < 1) {
            Debug(Debug::ERROR) << "--target-groups " << spec << ": expected field:<separator>:<n>, e.g. field:|:1\n";
            EXIT(EXIT_FAILURE);
        }
        hasLookup = forEachLookup(targetDb, [&](unsigned int key, const std::string &name, unsigned int) {
            size_t pos = 0;
            for (int i = 0; i < nG; ++i) {
                pos = name.find(sep, pos);
                if (pos == std::string::npos) {
                    if (i + 1 < nG) {
                        return;                         // fewer fields than n: not assigned
                    }
                    break;
                }
                if (i + 1 < nG) {
                    pos += sep.size();
                }
            }
            assign(key, genome(pos == std::string::npos ? name : name.substr(0, pos)));
        });
        out.how = "the first " + SSTR(nG) + " '" + sep + "'-separated field(s) of the target identifier";
    } else {
        std::vector<std::string> fasta = genomeFastaFiles(spec);
        std::unordered_map<std::string, unsigned int> byName;
        if (!fasta.empty()) {
            // one FASTA file per genome: every sequence identifier (as createdb reads it) -> the file's genome
            size_t duplicated = 0;
            for (size_t i = 0; i < fasta.size(); ++i) {
                const unsigned int g = genome(stripFastaName(fasta[i], &isFasta));
                KSeqWrapper *kseq = KSeqFactory(fasta[i].c_str());
                while (kseq->ReadEntry()) {
                    const std::string id = Util::parseFastaHeader(kseq->entry.name.s);
                    if (!id.empty() && !byName.emplace(id, g).second) {
                        ++duplicated;
                    }
                }
                delete kseq;
            }
            if (duplicated > 0) {
                Debug(Debug::WARNING) << duplicated << " sequence identifiers occur in more than one genome FASTA file; the first file wins. "
                                      << "Build the target DB from these files with createdb and use --target-groups source instead.\n";
            }
            out.how = SSTR(fasta.size()) + " genome FASTA files (" + spec + ")";
        } else if (FileUtil::fileExists(spec.c_str())) {
            std::ifstream in(spec);
            std::string line;
            while (std::getline(in, line)) {
                if (line.empty() || line[0] == '#') {
                    continue;
                }
                const size_t a = line.find('\t');
                if (a == std::string::npos) {
                    continue;
                }
                const size_t b = line.find('\t', a + 1);
                byName[line.substr(0, a)] = genome(line.substr(a + 1, b == std::string::npos ? std::string::npos : b - a - 1));
            }
            out.how = "the table " + spec;
        } else {
            Debug(Debug::ERROR) << "--target-groups " << spec << ": not 'source', 'taxid', 'field:...', a directory of genome "
                                << "FASTA files, a file listing them, or a target->genome table.\n";
            EXIT(EXIT_FAILURE);
        }
        hasLookup = forEachLookup(targetDb, [&](unsigned int key, const std::string &name, unsigned int) {
            std::unordered_map<std::string, unsigned int>::const_iterator it = byName.find(name);
            if (it != byName.end()) {
                assign(key, it->second);
            }
        });
    }
    if (!hasLookup) {
        Debug(Debug::ERROR) << targetDb << ".lookup is missing: --target-groups needs the target identifiers.\n";
        EXIT(EXIT_FAILURE);
    }
    if (out.mapped == 0) {
        Debug(Debug::WARNING) << "--target-groups " << spec << " assigned no target of " << targetDb
                              << ": every target is its own group.\n";
    }
    return out;
}

}  // namespace EMTargetGroups

#endif
