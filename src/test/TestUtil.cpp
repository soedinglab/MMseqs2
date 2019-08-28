#include <string>
#include <iostream>
#include <cassert>

#include "Debug.h"
#include "Util.h"

const char* binary_name = "test_util";

int main (int, const char**) {
    assert(SSTR((unsigned int)1) == "1");
    assert(SSTR((int)1) == "1");
    assert(SSTR((unsigned long long int)1) == "1");
    assert(SSTR((long long int)1) == "1");
    assert(SSTR((unsigned int)1) == "1");
    assert(SSTR((unsigned int)1) != "2");
    assert(SSTR('c') == "c");
    assert(SSTR(0.00314) == "3.140E-03");
    assert(SSTR((double)0.00314) == "3.140E-03");
    assert(SSTR("TEST") == "TEST");

//    unsigned int sizes[5] = {1882, 150, 630, 929, 167};

//    for (size_t i = 0; i < 5; ++i) {
//        size_t start = 0;
//        size_t end = i;
//
//        Util::decomposeDomainByAminoAcid(3758, sizes, 5, i, 5, &start, &end);
//        std::cout << start << " " << end << std::endl;
//    }

    for (size_t i = 0; i < 5; ++i) {
        size_t start = 0;
        size_t end = i;

        Util::decomposeDomain(5, i, 5, &start, &end);
        std::cout << start << " " << end << std::endl;
    }
}

