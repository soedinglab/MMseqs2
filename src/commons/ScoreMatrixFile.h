#ifndef SCOREMATRIXFILE_H
#define SCOREMATRIXFILE_H

#include <string>

class ScoreMatrixFile {
public:
    explicit ScoreMatrixFile(const char* filename);
    explicit ScoreMatrixFile(const std::string& filename) : ScoreMatrixFile(filename.c_str()) {}
    ScoreMatrixFile(const char* aminoacids, const char* nucleotides);
    ScoreMatrixFile(const ScoreMatrixFile& copy) : ScoreMatrixFile(copy.aminoacids, copy.nucleotides) {}
    ScoreMatrixFile& operator=(const ScoreMatrixFile& other);
    ~ScoreMatrixFile();

    bool operator==(const char* other) const;
    bool operator==(const std::string& other) const;
    bool operator==(const ScoreMatrixFile& other) const;

    bool operator!=(const char* other) const {
        return !(operator==(other));
    }
    bool operator!=(const std::string& other) const {
        return !(operator==(other));
    }
    bool operator!=(const ScoreMatrixFile& other) const {
        return !(operator==(other));
    }

    static std::string format(const ScoreMatrixFile &file);

    char* nucleotides;
    char* aminoacids;
};

#endif
