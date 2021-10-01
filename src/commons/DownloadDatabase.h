//
// Created by Martin Steinegger on 10/1/21.
//

#ifndef MMSEQS_DOWNLOADDATABASE_H
#define MMSEQS_DOWNLOADDATABASE_H
#include <vector>
#include <string>


struct EnvironmentEntry {
    const char* key;
    const char* value;
};

struct DatabaseDownload {
    const char *name;
    const char *description;
    const char *citation;
    const char *url;
    bool hasTaxonomy;
    int dbType;
    const unsigned char *script;
    size_t scriptLength;
    std::vector<EnvironmentEntry> environment;
};


#endif //MMSEQS_DOWNLOADDATABASE_H
