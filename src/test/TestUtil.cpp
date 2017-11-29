#include <string>
#include <iostream>

#include "Util.h"

const char* binary_name = "test_util";

int main (int argc, const char **argv) {
    unsigned int sizes[5] = {1882, 150, 630, 929, 167};

    for (size_t i = 0; i < 5; ++i) {
        size_t start = 0;
        size_t end = i;

        Util::decomposeDomainByAminoAcid(3758, sizes, 5, i, 5, &start, &end);
        std::cout << start << " " << end << std::endl;
    }

    for (size_t i = 0; i < 5; ++i) {
        size_t start = 0;
        size_t end = i;

        Util::decomposeDomain(5, i, 5, &start, &end);
        std::cout << start << " " << end << std::endl;
    }
}

