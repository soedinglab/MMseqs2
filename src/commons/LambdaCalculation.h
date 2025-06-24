#ifndef LAMBDA_CALC_H
#define LAMBDA_CALC_H
#include <vector>

// Enum to specify matrix units
enum class MatrixUnits {
    HALF_BIT,    // Typical for BLOSUM matrices
    BIT,         // Full bit units
    NATURAL_LOG  // Natural log units (required for lambda calculation)
};

// Forward declaration
struct Matrix;

double calculate_lambda(
    const double** raw_mat_b,
    const int alpha_size,
    std::vector<double>& final_p,
    std::vector<double>& final_q,
    double eps = 1e-12
);

#endif