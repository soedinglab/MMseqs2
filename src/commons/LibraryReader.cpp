//
// Created by mad on 11/29/16.
//
#include "LibraryReader.h"
#include <cstring>
#include "Debug.h"

// Returns pointer to first non-white-space character in str OR to NULL if none
// found
const char* LibraryReader::strscn(const char* str) {
    if (!str) return NULL;
    const char* ptr = str;
    while (*ptr != '\0' && isspace(*ptr)) ++ptr;
    return (*ptr == '\0') ? NULL : ptr;
}

// Returns true iff next non-blank line in file 'fp' contains string 'id'.
bool LibraryReader::StreamStartsWith(std::stringstream &in, const char* id) {
    std::string line;
    std::string idStr(id);
    if(in.good()){
        std::getline(in, line);
    }

    return (!line.compare(0, idStr.size(), idStr));
}


// Parse serialization record and return integer value following label 'str' in
// line read from file pointer 'fp'.
int LibraryReader::ReadInt(const char* line,
                           const char* label,
                           const char* errmsg = NULL) {
    int rv = 0;
    if (strstr(line, label)) {
        const char* ptr = line + strlen(label);
        rv = atoi(ptr);
    } else if (errmsg) {
        Debug(Debug::WARNING) << errmsg;
    }
    return rv;
}

// Parse serialization record and return double value following label 'label' in
// line read from file pointer 'fp'.
double LibraryReader::ReadDouble(const char* line,
                                 const char* label,
                                 const char* errmsg = NULL) {
    double rv = 0;
    if (strstr(line, label)) {
        rv = atof(line + strlen(label));
    } else if (errmsg) {
        Debug(Debug::WARNING) << errmsg;
    }
    return rv;
}


// Parse serialization record and return string following label 'label' in
// line read from file pointer 'fp'.
std::string  LibraryReader::ReadString(const char* line,
                                       const char* label,
                                       const char* errmsg = NULL) {
    std::string rv;
    if (strstr(line, label)) {
        const char* ptr = strscn(line + strlen(label));
        rv = ptr;
    } else if (errmsg) {
        Debug(Debug::WARNING) << errmsg;
    }
    return rv;
}

std::string LibraryReader::getline(std::stringstream &in) {
    std::string ret;
    if(in.good()){
        std::getline(in, ret);
    }
    return ret;
}



// Parse serialization record and return bool value following label 'str' in
// line read from file pointer 'fp'.
bool LibraryReader::ReadBool(const char* line,
                             const char* label,
                             const char* errmsg = NULL) {
    bool rv = false;
    if (strstr(line, label)) {
        const char* ptr = line + strlen(label);
        if (strchr(ptr, 'T') != NULL || strchr(ptr, '1') != NULL)
            rv = true;
        else if (strchr(ptr, 'F') != NULL || strchr(ptr, '0') != NULL)
            rv = false;
        else if (errmsg)
            Debug(Debug::WARNING) << errmsg;
    } else if (errmsg) {
        Debug(Debug::WARNING) << errmsg;
    }
    return rv;
}


std::vector<std::string> LibraryReader::tokenize(const char* str, char sep = ' ') {

    std::stringstream tok;
    std::vector<std::string> tokens;
    size_t pos=0;

    while (str[pos] != '\0')
    {
        if (str[pos] == sep)
        {
            tokens.push_back(tok.str());
            tok.str("");
        }

        // go to the next non-empty field
        while (str[pos] != '\0' && str[pos] == sep)
            pos++;

        tok << str[pos];
        pos++;
    }
    tokens.push_back(tok.str());
    return tokens;
}

