#include <string>
#include <iostream>
#include <cassert>

#include "Debug.h"
#include "Util.h"
#include "SubstitutionMatrix.h"

const char *binary_name = "test_lambdacalculator";

int main(int argc, const char ** argv) {
    if (argc < 2) {
        return EXIT_FAILURE;
    }

    SubstitutionMatrix mat(argv[1], 2.0f, 0.0f);
}