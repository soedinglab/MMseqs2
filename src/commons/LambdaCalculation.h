#ifndef LAMBDA_CALC_H
#define LAMBDA_CALC_H
#include <vector>

double calculate_lambda(
    const double** raw_mat_b,
    const int alpha_size,
    std::vector<double>& p,
    std::vector<double>& q
);

double calculate_lambda_direct(
    const double** raw_mat_b,
    const int alpha_size,
    std::vector<double>& p,
    std::vector<double>& q
);

#endif