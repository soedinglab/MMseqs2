#include "CandidateEdge.h"

#include "Debug.h"
#include "FileUtil.h"
#include "Util.h"

#include <cerrno>
#include <cstring>

std::string EdgeWriter::partitionPath(const std::string &dir, unsigned int partition) {
    return dir + "/p" + SSTR(partition) + ".edges";
}

EdgeWriter::EdgeWriter(const std::string &path, size_t bufferRecords)
    : path(path), file(NULL), bufferRecords(bufferRecords), edgeCount(0), closed(false) {
    buffer.reserve(bufferRecords);
}

EdgeWriter::~EdgeWriter() {
    close();
}

void EdgeWriter::flush() {
    if (buffer.empty()) {
        return;
    }
    if (file == NULL) {
        file = fopen(path.c_str(), "wb");
        if (file == NULL) {
            Debug(Debug::ERROR) << "Cannot open edge file " << path << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    if (fwrite(buffer.data(), sizeof(CandidateEdge), buffer.size(), file) != buffer.size()) {
        Debug(Debug::ERROR) << "Cannot write " << buffer.size() << " edges to " << path << ": "
                            << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    buffer.clear();
}

void EdgeWriter::append(const CandidateEdge &edge) {
    buffer.push_back(edge);
    if (buffer.size() >= bufferRecords) {
        flush();
    }
    edgeCount++;
}

void EdgeWriter::close() {
    // The destructor calls this too. Without the guard a second call would fall
    // into the empty-file branch below and truncate what the first call wrote.
    if (closed) {
        return;
    }
    closed = true;
    flush();
    if (file != NULL) {
        // A partition that produced no edge still gets an empty file, so a reader
        // can tell "this partition was reduced and had nothing" from "this
        // partition was never reduced".
        if (fclose(file) != 0) {
            Debug(Debug::ERROR) << "Cannot close edge file " << path << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        file = NULL;
    } else {
        FILE *empty = fopen(path.c_str(), "wb");
        if (empty == NULL) {
            Debug(Debug::ERROR) << "Cannot create edge file " << path << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(empty);
    }
}
