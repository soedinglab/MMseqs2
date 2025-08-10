#ifndef LAMBDA_CALC_H
#define LAMBDA_CALC_H
#include <vector>

double calculate_lambda(
    const double** raw_mat_b,
    const int alpha_size,
    std::vector<double>& final_p,
    std::vector<double>& final_q,
    double eps = 1e-12
);

#endif