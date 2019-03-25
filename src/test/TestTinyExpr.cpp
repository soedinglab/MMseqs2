#include <iostream>
#include <cassert>

#include "ExpressionParser.h"

const char* binary_name = "test_tinyexpr";

int main (int, const char**) {
    ExpressionParser expression("sqrt($11^2+$2^2)");
    if (expression.isOk()) {
        std::vector<int> indices = expression.findBindableIndices();
        assert(indices[0] == 11);
        assert(indices[1] == 2);
        expression.bind(11, 3);
        expression.bind(2, 4);
        double result = expression.evaluate();
        std::cout << result << std::endl;
        assert(result == 5);
    } else {
        std::cerr << "Failed to parse expression" << std::endl;
        assert(false);
    }
}

