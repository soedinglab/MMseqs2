#ifndef MMSEQS_NODEPLACEMENT_H
#define MMSEQS_NODEPLACEMENT_H

#include "Parameters.h"
#include "Debug.h"
#include "Util.h"

#include <climits>
#include <cstddef>
#include <cctype>
#include <cstring>
#include <string>
#include <functional>
#include <unistd.h>
#include <vector>

class Lin8DbIndex;

struct NodePlacement {
    unsigned int index;
    unsigned int count;

    // with one node there is no node to name, so the sentence opens with its verb instead
    std::string says(const char *verb) const;
    // one node's share of anything is all of it
    std::string share(size_t mine, size_t total) const;

    static NodePlacement resolve(const Parameters &par);
};

inline std::string NodePlacement::says(const char *verb) const {
    if (count == 1) {
        std::string opening(verb);
        opening[0] = static_cast<char>(toupper(opening[0]));
        return opening + " ";
    }
    return "Node " + SSTR(index) + " of " + SSTR(count) + " " + verb + " ";
}

inline std::string NodePlacement::share(size_t mine, size_t total) const {
    return mine == total ? SSTR(total) : SSTR(mine) + " of " + SSTR(total);
}

inline NodePlacement NodePlacement::resolve(const Parameters &par) {
    NodePlacement placement;
    placement.index = 0;
    placement.count = 1;
    // the short host name, so a fully qualified name still matches a plain name in --node-list
    char host[HOST_NAME_MAX + 1];
    memset(host, 0, sizeof(host));
    if (gethostname(host, HOST_NAME_MAX) != 0) {
        Debug(Debug::ERROR) << "Cannot read the host name to place this node\n";
        EXIT(EXIT_FAILURE);
    }
    char *dot = strchr(host, '.');
    if (dot != NULL) {
        *dot = '\0';
    }
    if (par.linclusterdbNodeList.empty() == false) {
        const std::vector<std::string> names = Util::split(par.linclusterdbNodeList, ",");
        if (names.empty()) {
            Debug(Debug::ERROR) << "--node-list " << par.linclusterdbNodeList << " is empty\n";
            EXIT(EXIT_FAILURE);
        }
        placement.count = static_cast<unsigned int>(names.size());
        unsigned int matched = UINT_MAX;
        for (size_t i = 0; i < names.size(); i++) {
            const std::string listed = names[i].substr(0, names[i].find('.'));
            const bool mine = listed == host;
            Debug(Debug::INFO) << "Node " << i << " is " << names[i] << (mine ? "  <- this host\n" : "\n");
            if (mine) {
                matched = static_cast<unsigned int>(i);
            }
        }
        if (par.linclusterdbNodeId >= 0) {
            placement.index = static_cast<unsigned int>(par.linclusterdbNodeId);
            if (matched != UINT_MAX && matched != placement.index) {
                Debug(Debug::ERROR) << "Host " << host << " is node " << matched
                                    << " in --node-list, but --node-id says " << placement.index
                                    << ". Give this host the node id that matches its place in the list\n";
                EXIT(EXIT_FAILURE);
            }
        } else if (matched != UINT_MAX) {
            placement.index = matched;
        } else {
            Debug(Debug::ERROR) << "Host " << host << " is not in --node-list "
                                << par.linclusterdbNodeList << ", pass --node-id instead\n";
            EXIT(EXIT_FAILURE);
        }
    } else {
        if (par.linclusterdbNodeCount > 0) {
            placement.count = static_cast<unsigned int>(par.linclusterdbNodeCount);
        }
        if (par.linclusterdbNodeId >= 0) {
            placement.index = static_cast<unsigned int>(par.linclusterdbNodeId);
        }
    }
    if (placement.index >= placement.count) {
        Debug(Debug::ERROR) << "Node id " << placement.index << " is outside the node count "
                            << placement.count << "\n";
        EXIT(EXIT_FAILURE);
    }
    Debug(Debug::INFO) << "This is host " << host;
    if (placement.count > 1) {
        Debug(Debug::INFO) << ", node " << placement.index << " of " << placement.count;
    }
    Debug(Debug::INFO) << "\n";
    return placement;
}

std::vector<size_t> nodeFileSlots(const Lin8DbIndex &index, const NodePlacement &node,
                                  const std::function<uint64_t(uint32_t)> &costOfLength = NULL);

#endif
