#ifndef PHP_A3M_READER_H
#define PHP_A3M_READER_H

#include <cstddef>
#include <vector>
#include <string>

class A3mReader {
public:
    A3mReader(std::string a3m);

    std::string getFasta();

private:
    void addSequence(const std::string& sequence);

    bool columnHasInsertion(size_t col);

    std::vector<std::string> headers;
    std::vector<std::vector<char>> entries;
    size_t length;
};

#endif
