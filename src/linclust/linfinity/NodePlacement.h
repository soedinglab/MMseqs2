#ifndef MMSEQS_NODEPLACEMENT_H
#define MMSEQS_NODEPLACEMENT_H

#include "Parameters.h"
#include "Debug.h"
#include "Util.h"

#include <climits>
#include <cstddef>
#include <cstring>
#include <string>
#include <functional>
#include <unistd.h>
#include <vector>

class SequenceLocator;

struct NodePlacement {
    unsigned int index;
    unsigned int count;

    static NodePlacement resolve(const Parameters &par);
};

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
    Debug(Debug::INFO) << "This is host " << host << ", node " << placement.index
                       << " of " << placement.count << "\n";
    return placement;
}

std::vector<size_t> nodeFileSlots(const SequenceLocator &runs, const NodePlacement &node,
                                  const std::function<uint64_t(uint32_t)> &costOfLength = NULL);

#endif
