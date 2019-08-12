#include "ScoreMatrixFile.h"
#include "Debug.h"

#include <cstring>
#include <cstdlib>

ScoreMatrixFile::ScoreMatrixFile(const char* filename) {
    if (strchr(filename, ',') != NULL) {
        size_t len = strlen(filename);
        aminoacids = (char*) malloc(len * sizeof(char));
        nucleotides = (char*) malloc(len * sizeof(char));
        if (sscanf(filename, "aa:%[^,],nucl:%s", aminoacids, nucleotides) != 2 && sscanf(filename, "nucl:%[^,],aa:%s", nucleotides, aminoacids) != 2) {
            free(nucleotides);
            free(aminoacids);
            nucleotides = strdup("INVALID");
            aminoacids = strdup("INVALID");
        }
    } else {
        nucleotides = strdup(filename);
        aminoacids = strdup(filename);
    }
}

ScoreMatrixFile::ScoreMatrixFile(const char* aminoacids, const char* nucleotides) {
    this->nucleotides = strdup(nucleotides);
    this->aminoacids = strdup(aminoacids);
}

ScoreMatrixFile& ScoreMatrixFile::operator=(const ScoreMatrixFile& other) {
    if (nucleotides != NULL) {
        free(nucleotides);
    }
    if (aminoacids != NULL) {
        free(aminoacids);
    }
    nucleotides = strdup(other.nucleotides);
    aminoacids = strdup(other.aminoacids);
    return *this;
}

ScoreMatrixFile::~ScoreMatrixFile() {
    free(nucleotides);
    free(aminoacids);
}

bool ScoreMatrixFile::operator==(const char* other) const {
    return strncmp(other, nucleotides, strlen(nucleotides)) == 0 || strncmp(other, aminoacids, strlen(aminoacids)) == 0;
}

bool ScoreMatrixFile::operator==(const std::string& other) const {
    return strncmp(other.c_str(), nucleotides, strlen(nucleotides)) == 0 || strncmp(other.c_str(), aminoacids, strlen(aminoacids)) == 0;
}

bool ScoreMatrixFile::operator==(const ScoreMatrixFile& other) const {
    return strncmp(other.nucleotides, nucleotides, strlen(nucleotides)) == 0 && strncmp(other.aminoacids, aminoacids, strlen(aminoacids)) == 0;
}

std::string ScoreMatrixFile::format(const ScoreMatrixFile &file) {
    if (strncmp(file.nucleotides, file.aminoacids, strlen(file.aminoacids)) == 0) {
        return file.nucleotides;
    } else {
        return std::string("nucl:") + file.nucleotides + ",aa:" + file.aminoacids;
    }
}
