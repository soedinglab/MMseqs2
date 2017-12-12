//
// Created by mad on 11/29/16.
//

#ifndef MMSEQS_LIBRARYREADER_H
#define MMSEQS_LIBRARYREADER_H


#include <sstream>
#include <vector>

class LibraryReader {
public:
    bool StreamStartsWith(std::stringstream &in, const char* id);
    int ReadInt(const char* line,const char* label,const char* errmsg);
    double ReadDouble(const char* line,const char* label,const char* errmsg);
    std::string ReadString(const char* line,const char* label,const char* errmsg);
    bool ReadBool(const char* line,const char* label,const char* errmsg);
    const char* strscn(const char* str) ;
    static std::vector<std::string> tokenize(const char* str, char sep);
    std::string getline(std::stringstream &in);
};


#endif //MMSEQS_LIBRARYREADER_H
