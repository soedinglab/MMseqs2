#ifndef LAMBDA_CALC_H
#define LAMBDA_CALC_H
#include <vector>

double calculate_lambda(
    const double** mat_b,
    const int alpha_size,
    std::vector<double>& p,
    std::vector<double>& q
);

#endif