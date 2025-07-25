#include "LambdaCalculation.h"
#include "Debug.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <numeric>
#include <vector>
#include <random>

#define LambdaCalculation_DEBUG 0

struct LambdaCalculationConstants {
    static constexpr double MATRIX_SINGULARITY_THRESHOLD = 1e-15;
    static constexpr double CONVERGENCE_TOLERANCE = 1e-10;
    static constexpr double NEWTON_CONVERGENCE_TOLERANCE = 1e-10;
    static constexpr double LAMBDA_LOWER_BOUND = 1e-6;
    static constexpr double BISECTION_CONVERGENCE_TOLERANCE = 1e-12;
    static constexpr double POWER_ITERATION_TOLERANCE = 1e-12;
    static constexpr double DAMPING_FACTOR = 0.8;
    static constexpr double LAMBDA_UPPER_BOUND_SAFETY = 0.95;
    static constexpr double PROBABILITY_VALIDATION_TOLERANCE = 1e-10;
    static constexpr double RESTRICTION_VALIDATION_TOLERANCE = 1e-6;
    static constexpr double ZERO_ROW_COL_THRESHOLD = 1e-10;

    static constexpr int MAX_NEWTON_ITERATIONS = 50;
    static constexpr int MAX_BISECTION_ITERATIONS = 100;
    static constexpr int MAX_POWER_ITERATIONS = 1000;
    static constexpr int MAX_BRACKET_ATTEMPTS = 50;

    static constexpr double BRACKET_INCREMENT_FACTOR = 1.5;
    static constexpr double BRACKET_START_LAMBDA = 0.01;
};

struct Probs {
    int alphabet_size;
    std::vector<double> values;

    // Constructor: uniform_init = true for uniform, random_init = true for random valid probabilities
    Probs(const int alphabet_size, bool uniform_init = false, bool random_init = false)
        : alphabet_size(alphabet_size), values(alphabet_size, 0.0) {
        if (uniform_init) {
            double uniform_prob = 1.0 / static_cast<double>(alphabet_size);
            std::fill(values.begin(), values.end(), uniform_prob);
        } else if (random_init) {
            // Generate random valid probability vector (sum to 1, all >= 0)
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(0.0, 1.0);

            double sum = 0.0;
            for (int i = 0; i < alphabet_size; ++i) {
                values[i] = dis(gen);
                sum += values[i];
            }
            if (sum > 0.0) {
                for (int i = 0; i < alphabet_size; ++i) {
                    values[i] /= sum;
                }
            } else {
                // fallback to uniform if all randoms are zero (very unlikely)
                double uniform_prob = 1.0 / static_cast<double>(alphabet_size);
                std::fill(values.begin(), values.end(), uniform_prob);
            }
        }
    }

    // Copy constructor from another Probs
    Probs(const Probs& other)
        : alphabet_size(other.alphabet_size), values(other.values) {}

    // Assignment from another Probs
    Probs& operator=(const Probs& other) {
        if (this != &other) {
            alphabet_size = other.alphabet_size;
            values = other.values;
        }
        return *this;
    }

    // Initialize from another Probs (explicit function)
    void initialize_from(const Probs& other) {
        alphabet_size = other.alphabet_size;
        values = other.values;
    }

    double &operator[](int i) { return values[i]; }
    const double &operator[](int i) const { return values[i]; }
};

struct Matrix {
    int row_dim;
    int col_dim;
    std::vector<double> values;

    Matrix(const int row_dim, const int col_dim)
        : row_dim(row_dim), col_dim(col_dim), values(row_dim * col_dim, 0.0) {}

    double &at(int i, int j) {
        return values[i * col_dim + j];
    }
    const double &at(int i, int j) const {
        return values[i * col_dim + j];
    }

    void copy_from(const double **mat) {
        for (int i = 0; i < row_dim; ++i) {
            for (int j = 0; j < col_dim; ++j) {
                at(i, j) = mat[i][j];
            }
        }
    }

    // Fills the given matrix A with: A[i][j] = p[i] * q[j] * exp(lambda * S[i][j])
    void fill_A_matrix(Matrix& A, const Probs& p, const Probs& q, double lambda) const {
        for (int i = 0; i < row_dim; ++i) {
            for (int j = 0; j < col_dim; ++j) {
                A.at(i, j) = p[i] * q[j] * std::exp(lambda * this->at(i, j));
            }
        }
    }

    bool validate() const {
        if (row_dim <= 0 || col_dim <= 0) {
            return false;
        }
        if (values.size() != static_cast<size_t>(row_dim * col_dim)) {
            return false;
        }
        return std::none_of(values.begin(), values.end(), [](double v) {
            return std::isnan(v) || std::isinf(v);
        });
    }

    // Fill this matrix with the joint probability distribution p[i] * q[j]
    // Assumes this matrix is already sized to (p.alphabet_size, q.alphabet_size)
    void fill_joint_from_probs(const Probs &p, const Probs &q) {
        for (int i = 0; i < row_dim; ++i) {
            for (int j = 0; j < col_dim; ++j) {
                this->at(i, j) = p[i] * q[j];
            }
        }
    }

    // Function to compute the initial score distribution Ps
    // Given a score matrix S (this), and probability vectors p and q,
    // returns a map from score value s to probability Ps
    std::map<int, double>
    get_score_distribution(const std::vector<double> &p, const std::vector<double> &q) const {
        std::map<int, double> score_dist;
        // Find min and max score values (assume integer scores)
        int min_score = static_cast<int>(std::floor(values[0]));
        int max_score = static_cast<int>(std::ceil(values[0]));
        for (int i = 0; i < row_dim; ++i) {
            for (int j = 0; j < col_dim; ++j) {
                int s = static_cast<int>(std::round(at(i, j)));
                if (s < min_score) {
                    min_score = s;
                }
                if (s > max_score) {
                    max_score = s;
                }
            }
        }
        // For each (i, j), accumulate p[i] * q[j] for the corresponding score s
        double total = 0.0;
        for (int i = 0; i < row_dim; ++i) {
            for (int j = 0; j < col_dim; ++j) {
                int s = static_cast<int>(std::round(at(i, j)));
                double prob = p[i] * q[j];
                score_dist[s] += prob;
                total += prob;
            }
        }
        // Normalize so that the sum of probabilities is 1.0
        if (total > 0.0) {
            for (auto &kv : score_dist) {
                kv.second /= total;
            }
        }
        return score_dist;
    }

    double calculate_determinant(void) {
        // Only for square matrices
        if (row_dim != col_dim) {
            return 0.0;
        }
        int n = row_dim;
        // Make a copy of the matrix values for row operations
        std::vector<double> mat(values);
        double det = 1.0;
        int sign = 1;

        for (int i = 0; i < n; ++i) {
            // Find pivot
            int pivot = i;
            for (int j = i + 1; j < n; ++j) {
                if (std::abs(mat[j * n + i]) > std::abs(mat[pivot * n + i])) {
                    pivot = j;
                }
            }
            if (std::abs(mat[pivot * n + i]) < LambdaCalculationConstants::MATRIX_SINGULARITY_THRESHOLD) {
                return 0.0;
            }
            if (pivot != i) {
                // Swap rows
                for (int k = 0; k < n; ++k) {
                    std::swap(mat[i * n + k], mat[pivot * n + k]);
                }
                sign *= -1;
            }
            det *= mat[i * n + i];
            // Eliminate below
            for (int j = i + 1; j < n; ++j) {
                double factor = mat[j * n + i] / mat[i * n + i];
                for (int k = i; k < n; ++k) {
                    mat[j * n + k] -= factor * mat[i * n + k];
                }
            }
        }
        return det * sign;
    }

    // LU factorization with partial pivoting
    // Returns false if matrix is singular
    bool lu_factorize(std::vector<int>& pivot_indices, std::vector<double>& lu_matrix) const {
        if (row_dim != col_dim) {
            return false;	// Only square matrices
        }

        int n = row_dim;
        lu_matrix = values;	// Copy matrix values
        pivot_indices.resize(n);

        // Initialize pivot indices
        for (int i = 0; i < n; ++i) {
            pivot_indices[i] = i;
        }

        for (int k = 0; k < n - 1; ++k) {
            // Find pivot
            int pivot = k;
            double max_val = std::abs(lu_matrix[k * n + k]);

            for (int i = k + 1; i < n; ++i) {
                double val = std::abs(lu_matrix[i * n + k]);
                if (val > max_val) {
                    max_val = val;
                    pivot = i;
                }
            }

            // Check for singularity
            if (max_val < LambdaCalculationConstants::MATRIX_SINGULARITY_THRESHOLD) {
                return false;
            }

            // Swap rows if needed
            if (pivot != k) {
                std::swap(pivot_indices[k], pivot_indices[pivot]);
                for (int j = 0; j < n; ++j) {
                    std::swap(lu_matrix[k * n + j], lu_matrix[pivot * n + j]);
                }
            }

            // Eliminate below pivot
            for (int i = k + 1; i < n; ++i) {
                double factor = lu_matrix[i * n + k] / lu_matrix[k * n + k];
                lu_matrix[i * n + k] = factor;	// Store multiplier in L part

                for (int j = k + 1; j < n; ++j) {
                    lu_matrix[i * n + j] -= factor * lu_matrix[k * n + j];
                }
            }
        }

        return true;
    }

    // Forward substitution: solve L * y = b where L is lower triangular
    void forward_substitution(const std::vector<double>& lu_matrix,
                             const std::vector<int>& pivot_indices,
                             const std::vector<double>& b,
                             std::vector<double>& y) const {
        int n = row_dim;
        y.resize(n);

        // Apply permutation to b
        for (int i = 0; i < n; ++i) {
            y[i] = b[pivot_indices[i]];
        }

        // Forward substitution
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < i; ++j) {
                y[i] -= lu_matrix[i * n + j] * y[j];
            }
        }
    }

    // Back substitution: solve U * x = y where U is upper triangular
    void back_substitution(const std::vector<double>& lu_matrix,
                        const std::vector<double>& y,
                        std::vector<double>& x) const {
        int n = row_dim;
        x.resize(n);

        // Back substitution
        for (int i = n - 1; i >= 0; --i) {
            x[i] = y[i];
            for (int j = i + 1; j < n; ++j) {
                x[i] -= lu_matrix[i * n + j] * x[j];
            }
            x[i] /= lu_matrix[i * n + i];
        }
    }

    // Solve A * x = b using LU factorization
    bool solve_linear_system(const std::vector<double>& b, std::vector<double>& x) const {
        std::vector<int> pivot_indices;
        std::vector<double> lu_matrix;

        if (!lu_factorize(pivot_indices, lu_matrix)) {
            return false;	// Singular matrix
        }

        std::vector<double> y;
        forward_substitution(lu_matrix, pivot_indices, b, y);
        back_substitution(lu_matrix, y, x);

        return true;
    }
};

struct ScoreGroup {
    double score;
    double prob_sum;

    ScoreGroup(double s, double p) : score(s), prob_sum(p) {}
};

class LambdaDebug {
private:
    // Helper to print variadic arguments recursively (C++11/14 compatible)
    template<typename T>
    static void print_args(T&& arg) {
        if (enabled) {
            Debug(Debug::ERROR) << arg;
        }
    }

    template<typename T, typename... Args>
    static void print_args(T&& first, Args&&... rest) {
        if (enabled) {
            Debug(Debug::ERROR) << first;
            print_args(rest...);
        }
    }

public:
    template<typename... Args>
    static void log(Args&&... args) {
        if (enabled) {
            Debug(Debug::ERROR) << "Lambda: ";
            print_args(args...);
            Debug(Debug::ERROR) << "\n";
        }
    }

    static void log_matrix_info(const Matrix& matrix, const std::string& name) {
        if (enabled) {
            log(name, " matrix (first 5x5):");
            for (int i = 0; i < std::min(5, matrix.row_dim); ++i) {
                for (int j = 0; j < std::min(5, matrix.col_dim); ++j) {
                    Debug(Debug::ERROR) << matrix.at(i, j) << "\t";
                }
                Debug(Debug::ERROR) << "\n";
            }
        }
    }

    static void log_probabilities(const Probs& probs, const std::string& name, int count = 5) {
        if (enabled) {
            Debug(Debug::ERROR) << "Lambda: " << name << " (first " << count << "): ";
            for (int i = 0; i < std::min(count, probs.alphabet_size); ++i) {
                Debug(Debug::ERROR) << probs[i] << " ";
            }
            Debug(Debug::ERROR) << "\n";
        }
    }

    static void log_probability_changes(const std::vector<double>& old_probs, const Probs& new_probs,
                                         const std::string& name, int count = 5) {
        if (enabled) {
            Debug(Debug::ERROR) << "Lambda: " << name << " changes (first " << count << "): ";
            for (int i = 0; i < std::min(count, new_probs.alphabet_size); ++i) {
                double change = new_probs[i] - old_probs[i];
                Debug(Debug::ERROR) << name << "[" << i << "]:" << old_probs[i]
                                    << "->" << new_probs[i] << "(Δ=" << change << ") ";
            }
            Debug(Debug::ERROR) << "\n";
        }
    }

    static void log_joint_sample(const Matrix& joint_probs, int size = 3) {
        if (enabled) {
            log("Sample joint probabilities (first ", size, "x", size, "):");
            for (int i = 0; i < std::min(size, joint_probs.row_dim); ++i) {
                for (int j = 0; j < std::min(size, joint_probs.col_dim); ++j) {
                    Debug(Debug::ERROR) << joint_probs.at(i, j) << "\t";
                }
                Debug(Debug::ERROR) << "\n";
            }
        }
    }

    static bool validate_probability_distribution(const Probs& probs, const std::string& name) {
        double sum = std::accumulate(probs.values.begin(), probs.values.end(), 0.0);
        bool valid = true;

        if (enabled) {
            // Check sum is approximately 1.0
            if (std::abs(sum - 1.0) > LambdaCalculationConstants::PROBABILITY_VALIDATION_TOLERANCE) {
                log(name, " sum=", sum, " (should be 1.0)");
                valid = false;
            }

            // Check all values are non-negative and finite
            for (int i = 0; i < probs.alphabet_size; ++i) {
                if (probs[i] < 0.0) {
                    log(name, "[", i, "]=", probs[i], " (should be >= 0)");
                    valid = false;
                }
                if (std::isnan(probs[i]) || std::isinf(probs[i])) {
                    log(name, "[", i, "] is NaN or Inf");
                    valid = false;
                }
            }

            if (valid) {
                log(name, " is valid");
            }
        }

        return valid;
    }

    static const bool enabled = LambdaCalculation_DEBUG;
};

class ScoreMatrixValidator {
public:
    static bool validate_score_matrix(const Matrix& matrix) {
        LambdaDebug::log("Validating score matrix");

        // Check basic matrix properties
        if (!matrix.validate()) {
            LambdaDebug::log("Score matrix failed basic validation");
            return false;
        }

        // Check if matrix is square
        if (matrix.row_dim != matrix.col_dim) {
            LambdaDebug::log("Score matrix is not square: ", matrix.row_dim, "x", matrix.col_dim);
            return false;
        }

        // Check if matrix is symmetric (substitution matrices should be symmetric)
        for (int i = 0; i < matrix.row_dim; ++i) {
            for (int j = 0; j < matrix.col_dim; ++j) {
                if (std::abs(matrix.at(i, j) - matrix.at(j, i)) > LambdaCalculationConstants::CONVERGENCE_TOLERANCE) {
                    LambdaDebug::log("Score matrix is not symmetric at (", i, ",", j, "): ",
                                     matrix.at(i, j), " vs ", matrix.at(j, i));
                    return false;
                }
            }
        }

        // Check if diagonal elements are positive (self-alignment should be favorable)
        for (int i = 0; i < matrix.row_dim; ++i) {
            if (matrix.at(i, i) <= 0.0) {
                LambdaDebug::log("Diagonal element at (", i, ",", i, ") is not positive: ", matrix.at(i, i));
                return false;
            }
        }

        // Check if the matrix has reasonable score range
        double min_score = matrix.values[0];
        double max_score = matrix.values[0];
        for (double score : matrix.values) {
            min_score = std::min(min_score, score);
            max_score = std::max(max_score, score);
        }

        if (max_score - min_score <= 0.0) {
            LambdaDebug::log("Score matrix has no variation: min=", min_score, ", max=", max_score);
            return false;
        }

        LambdaDebug::log("Score matrix validation passed: range=[", min_score, ",", max_score, "]");
        return true;
    }
};

class MatrixAnalyzer {
public:
    // Group scores by their values to optimize restriction value calculations
    static std::vector<ScoreGroup> group_scores_by_value(const Matrix &score_matrix,
                                                        const Probs &p, const Probs &q) {
        std::map<double, double> score_to_prob;
        const double tolerance = LambdaCalculationConstants::CONVERGENCE_TOLERANCE;

        LambdaDebug::log("Grouping scores by value for optimization");

        for (int i = 0; i < score_matrix.row_dim; ++i) {
            for (int j = 0; j < score_matrix.col_dim; ++j) {
                double score_val = score_matrix.at(i, j);

                // Find existing score within tolerance or create new one
                bool found = false;
                for (auto &pair : score_to_prob) {
                    if (std::abs(pair.first - score_val) < tolerance) {
                        pair.second += p[i] * q[j];
                        found = true;
                        break;
                    }
                }

                if (!found) {
                    score_to_prob[score_val] += p[i] * q[j];
                }
            }
        }

        std::vector<ScoreGroup> score_groups;
        score_groups.reserve(score_to_prob.size());
        for (const auto &pair : score_to_prob) {
            score_groups.emplace_back(pair.first, pair.second);
        }

        LambdaDebug::log("Created ", score_groups.size(), " score groups from ",
                        score_matrix.row_dim * score_matrix.col_dim, " matrix entries");

        return score_groups;
    }

    // Optimized restriction value calculation using score groups
    static double restriction_value_grouped(const std::vector<ScoreGroup> &score_groups, double lambda) {
        double sum = 0.0;
        for (const auto &group : score_groups) {
            sum += group.prob_sum * std::exp(lambda * group.score);
        }
        return sum - 1.0;
    }

    // Optimized restriction value derivative calculation using score groups
    static double restriction_value_derivative_grouped(const std::vector<ScoreGroup> &score_groups, double lambda) {
        double sum = 0.0;
        for (const auto &group : score_groups) {
            sum += group.prob_sum * group.score * std::exp(lambda * group.score);
        }
        return sum;
    }

    static bool matrix_solvable_for_lambda(const Matrix &mat, double &lambda_upper_bound) {
        LambdaDebug::log("Checking matrix solvability");

        // Check for zero rows/columns
        std::vector<double> row_sums(mat.row_dim, 0.0);
        std::vector<double> col_sums(mat.col_dim, 0.0);

        for (int i = 0; i < mat.row_dim; ++i) {
            for (int j = 0; j < mat.col_dim; ++j) {
                row_sums[i] += mat.at(i, j);
                col_sums[j] += mat.at(i, j);
            }
        }

        if (LambdaDebug::enabled) {
            // Build row sums string for debugging
            std::string row_sums_str;
            for (int i = 0; i < std::min(5, mat.row_dim); ++i) {
                row_sums_str += std::to_string(row_sums[i]) + " ";
            }
            LambdaDebug::log("Row sums (first 5): ", row_sums_str);

            // Build col sums string for debugging
            std::string col_sums_str;
            for (int i = 0; i < std::min(5, mat.col_dim); ++i) {
                col_sums_str += std::to_string(col_sums[i]) + " ";
            }
            LambdaDebug::log("Col sums (first 5): ", col_sums_str);
        }

        if (std::any_of(row_sums.begin(), row_sums.end(),
                        [](double sum) { return std::abs(sum) < LambdaCalculationConstants::ZERO_ROW_COL_THRESHOLD; }) ||
            std::any_of(col_sums.begin(), col_sums.end(),
                        [](double sum) { return std::abs(sum) < LambdaCalculationConstants::ZERO_ROW_COL_THRESHOLD; })) {
            LambdaDebug::log("Matrix has zero row/column sums");
            return false;
        }

        // Check that each row and column has both positive and negative values
        std::vector<double> row_maxes(mat.row_dim, 0.0);
        std::vector<double> col_maxes(mat.col_dim, 0.0);

        for (int i = 0; i < mat.row_dim; ++i) {
            bool row_pos = false, row_neg = false;
            for (int j = 0; j < mat.col_dim; ++j) {
                double val = mat.at(i, j);
                if (val > 0.0) {
                    row_pos = true;
                    row_maxes[i] = std::max(row_maxes[i], val);
                    col_maxes[j] = std::max(col_maxes[j], val);
                } else if (val < 0.0) {
                    row_neg = true;
                }
            }
            if (!(row_pos && row_neg)) {
                LambdaDebug::log("Row ", i, " lacks both positive and negative values");
                return false;
            }
        }

        // Check columns
        for (int j = 0; j < mat.col_dim; ++j) {
            bool col_pos = false, col_neg = false;
            for (int i = 0; i < mat.row_dim; ++i) {
                double val = mat.at(i, j);
                if (val > 0.0)
                col_pos = true;
                else if (val < 0.0)
                col_neg = true;
            }
            if (!(col_pos && col_neg)) {
                LambdaDebug::log("Column ", j, " lacks both positive and negative values");
                return false;
            }
        }

        lambda_upper_bound =
            std::min(*std::min_element(row_maxes.begin(), row_maxes.end()),
                    *std::min_element(col_maxes.begin(), col_maxes.end()));

        LambdaDebug::log("Computed lambda_upper_bound=", lambda_upper_bound);
        LambdaDebug::log("Matrix is solvable");

        return true;
    }

    // Function to compute restriction value: sum_{i,j} p[i] * q[j] * exp(lambda * S[i,j]) - 1
    static double restriction_value_matrix(const Matrix &score_matrix, const Probs &p,
                                    const Probs &q, double lambda) {
        double sum = 0.0;
        for (int i = 0; i < score_matrix.row_dim; ++i) {
            for (int j = 0; j < score_matrix.col_dim; ++j) {
                sum += p[i] * q[j] * std::exp(lambda * score_matrix.at(i, j));
            }
        }
        return sum - 1.0;
    }

    // Function to compute derivative of restriction value
    static double restriction_value_derivative_matrix(const Matrix &score_matrix,
                                                 const Probs &p, const Probs &q,
                                                 double lambda) {
        double sum = 0.0;
        for (int i = 0; i < score_matrix.row_dim; ++i) {
            for (int j = 0; j < score_matrix.col_dim; ++j) {
                sum += p[i] * q[j] * score_matrix.at(i, j) * std::exp(lambda * score_matrix.at(i, j));
            }
        }
        return sum;
    }
};

bool bracket_lambda(Matrix score_matrix, const Probs p, const Probs q,
                    double &lambda_low, double &lambda_high,
                    double lambda_upper_bound) {
    // For a proper scoring matrix, f(0) = 0 and f'(0) < 0
    // So we need to find a small positive lambda where f(lambda) > 0

    // Use score grouping optimization for faster restriction value computation
    std::vector<ScoreGroup> score_groups = MatrixAnalyzer::group_scores_by_value(score_matrix, p, q);

    double f_zero = MatrixAnalyzer::restriction_value_grouped(score_groups, 0.0);
    LambdaDebug::log("f(0) = ", f_zero);

    // Check that we're close to f(0) = 0 (within tolerance)
    if (std::abs(f_zero) > LambdaCalculationConstants::PROBABILITY_VALIDATION_TOLERANCE) {
        LambdaDebug::log("f(0) = ", f_zero, " is not close to 0, substitution matrix may not be properly normalized");
        return false;
    }

    // Check the derivative at 0 (should be negative for a proper scoring matrix)
    double f_prime_zero = MatrixAnalyzer::restriction_value_derivative_grouped(score_groups, 0.0);
    LambdaDebug::log("f'(0) = ", f_prime_zero);

    if (f_prime_zero >= 0) {
        LambdaDebug::log("f'(0) = ", f_prime_zero, " >= 0, expected negative derivative");
        return false;
    }

    // Start with lambda_low = 0 and find lambda_high where f(lambda_high) > 0
    lambda_low = 0.0;
    lambda_high = LambdaCalculationConstants::BRACKET_START_LAMBDA; // Start with a small positive value

    const double increment_factor = LambdaCalculationConstants::BRACKET_INCREMENT_FACTOR;
    const int max_attempts = LambdaCalculationConstants::MAX_BRACKET_ATTEMPTS;

    for (int i = 0; i < max_attempts; ++i) {
        double f_high = MatrixAnalyzer::restriction_value_grouped(score_groups, lambda_high);

        LambdaDebug::log("Trying lambda_high = ", lambda_high, ", f(lambda_high) = ", f_high);

        if (f_high > 0.0) {
            LambdaDebug::log("Successfully bracketed lambda: [", lambda_low, ", ", lambda_high, "]");
            return true;
        }

        lambda_high *= increment_factor;

        // Don't exceed the theoretical upper bound
        if (lambda_high > lambda_upper_bound * LambdaCalculationConstants::LAMBDA_UPPER_BOUND_SAFETY) {
            lambda_high = lambda_upper_bound * LambdaCalculationConstants::LAMBDA_UPPER_BOUND_SAFETY;
            double f_high_final = MatrixAnalyzer::restriction_value_grouped(score_groups, lambda_high);

            LambdaDebug::log("At upper bound lambda_high = ", lambda_high, ", f(lambda_high) = ", f_high_final);

            if (f_high_final > 0.0) {
                LambdaDebug::log("Successfully bracketed lambda at upper bound: [", lambda_low, ", ", lambda_high, "]");
                return true;
            } else {
                LambdaDebug::log("Failed to bracket lambda even at upper bound");
                return false;
            }
        }
    }

    return false;
}

double find_lambda_bisection(const Matrix score_matrix, const Probs p,
                             const Probs q, double lambda_low,
                             double lambda_high, double epsilon = LambdaCalculationConstants::BISECTION_CONVERGENCE_TOLERANCE,
                             int max_iters = LambdaCalculationConstants::MAX_BISECTION_ITERATIONS) {
    LambdaDebug::log("Starting bisection: [", lambda_low, ", ", lambda_high, "]");

    // Use score grouping optimization for faster restriction value computation
    std::vector<ScoreGroup> score_groups = MatrixAnalyzer::group_scores_by_value(score_matrix, p, q);

    for (int iter = 0; iter < max_iters; ++iter) {
        if (lambda_high - lambda_low < epsilon) {
            LambdaDebug::log("Converged: interval width = ", (lambda_high - lambda_low));
            break;
        }

        double mid = 0.5 * (lambda_low + lambda_high);
        double f_mid = MatrixAnalyzer::restriction_value_grouped(score_groups, mid);

        LambdaDebug::log("iter ", iter + 1, ": mid=", mid, " f(mid)=", f_mid);

        // Update bounds
        if (f_mid > 0.0) {
            lambda_high = mid;
        } else {
            lambda_low = mid;
        }

        if (std::abs(f_mid) < epsilon) {
            LambdaDebug::log("Converged: |f(mid)| = ", std::abs(f_mid));
            return mid;
        }
    }

    return 0.5 * (lambda_low + lambda_high);
}

// Newton-Raphson method for finding lambda
// Uses both function value and derivative for faster convergence
double find_lambda_newton_raphson(const Matrix score_matrix, const Probs p,
                                    const Probs q, double lambda_low,
                                    double lambda_high, double epsilon = LambdaCalculationConstants::NEWTON_CONVERGENCE_TOLERANCE,
                                    int max_iters = LambdaCalculationConstants::MAX_BISECTION_ITERATIONS) {
    LambdaDebug::log("Starting Newton-Raphson: initial bracket [", lambda_low, ", ", lambda_high, "]");

    // Use score grouping optimization for faster restriction value computation
    std::vector<ScoreGroup> score_groups = MatrixAnalyzer::group_scores_by_value(score_matrix, p, q);

    // Start with midpoint of the bracket
    double lambda = 0.5 * (lambda_low + lambda_high);

    for (int iter = 0; iter < max_iters; ++iter) {
        double f_val = MatrixAnalyzer::restriction_value_grouped(score_groups, lambda);
        double f_prime = MatrixAnalyzer::restriction_value_derivative_grouped(score_groups, lambda);

        LambdaDebug::log("Newton-Raphson iter ", iter + 1, ": lambda=", lambda, " f(lambda)=", f_val, " f'(lambda)=", f_prime);

        // Check for convergence
        if (std::abs(f_val) < epsilon) {
            LambdaDebug::log("Newton-Raphson converged: |f(lambda)| = ", std::abs(f_val));
            return lambda;
        }

        // Check for zero derivative (would cause division by zero)
        if (std::abs(f_prime) < LambdaCalculationConstants::MATRIX_SINGULARITY_THRESHOLD) {
            LambdaDebug::log("Newton-Raphson: derivative too small, falling back to bisection");
            return find_lambda_bisection(score_matrix, p, q, lambda_low, lambda_high, epsilon, max_iters);
        }

        // Newton-Raphson step: lambda_new = lambda_old - f(lambda_old) / f'(lambda_old)
        double newton_step = -f_val / f_prime;
        double lambda_new = lambda + newton_step;

        LambdaDebug::log("Newton-Raphson: full step would be ", newton_step, " -> lambda=", lambda_new);

        // Use adaptive damping to stay within bounds
        double damping_factor = 1.0;
        if (lambda_new < lambda_low || lambda_new > lambda_high) {
            // Calculate damping factor to stay within bounds
            double max_allowed_step = 0.0;
            if (newton_step > 0) {
                max_allowed_step = (lambda_high - lambda) * LambdaCalculationConstants::DAMPING_FACTOR;	// Stay within bounds
            } else {
                max_allowed_step = (lambda_low - lambda) * LambdaCalculationConstants::DAMPING_FACTOR;	 // Stay within bounds
            }

            if (std::abs(newton_step) > std::abs(max_allowed_step)) {
                damping_factor = max_allowed_step / newton_step;
            }

            LambdaDebug::log("Newton-Raphson: damping factor=", damping_factor);
        }

        lambda_new = lambda + damping_factor * newton_step;
        LambdaDebug::log("Newton-Raphson: damped step=", (damping_factor * newton_step), " -> lambda=", lambda_new);

        // Check for convergence in lambda
        if (std::abs(lambda_new - lambda) < epsilon) {
            LambdaDebug::log("Newton-Raphson converged: lambda change = ", std::abs(lambda_new - lambda));
            return lambda_new;
        }

        // Update bracket bounds to maintain the bracket property
        double f_new = MatrixAnalyzer::restriction_value_grouped(score_groups, lambda_new);
        if (f_new > 0.0) {
            lambda_high = lambda_new;
        } else {
            lambda_low = lambda_new;
        }

        lambda = lambda_new;
    }

    LambdaDebug::log("Newton-Raphson: Max iterations reached, final lambda=", lambda);

    return lambda;
}

class LambdaFinder {
public:
    static bool bracket_lambda(Matrix score_matrix, const Probs p, const Probs q,
                                double &lambda_low, double &lambda_high,
                                double lambda_upper_bound) {
        // For a proper scoring matrix, f(0) = 0 and f'(0) < 0
        // So we need to find a small positive lambda where f(lambda) > 0

        // Use score grouping optimization for faster restriction value computation
        std::vector<ScoreGroup> score_groups = MatrixAnalyzer::group_scores_by_value(score_matrix, p, q);

        double f_zero = MatrixAnalyzer::restriction_value_grouped(score_groups, 0.0);

        LambdaDebug::log("f(0) = ", f_zero);

        // Check that we're close to f(0) = 0 (within tolerance)
        if (std::abs(f_zero) > LambdaCalculationConstants::PROBABILITY_VALIDATION_TOLERANCE) {
            LambdaDebug::log("f(0) = ", f_zero, " is not close to 0, substitution matrix may not be properly normalized");
            return false;
        }

        // Check the derivative at 0 (should be negative for a proper scoring matrix)
        double f_prime_zero = MatrixAnalyzer::restriction_value_derivative_grouped(score_groups, 0.0);
        LambdaDebug::log("f'(0) = ", f_prime_zero);

        if (f_prime_zero >= 0) {
            LambdaDebug::log("f'(0) = ", f_prime_zero, " >= 0, expected negative derivative");
            return false;
        }

        // Start with lambda_low = 0 and find lambda_high where f(lambda_high) > 0
        lambda_low = 0.0;
        lambda_high = LambdaCalculationConstants::BRACKET_START_LAMBDA; // Start with a small positive value

        const double increment_factor = LambdaCalculationConstants::BRACKET_INCREMENT_FACTOR;
        const int max_attempts = LambdaCalculationConstants::MAX_BRACKET_ATTEMPTS;

        for (int i = 0; i < max_attempts; ++i) {
            double f_high = MatrixAnalyzer::restriction_value_grouped(score_groups, lambda_high);

            LambdaDebug::log("Trying lambda_high = ", lambda_high, ", f(lambda_high) = ", f_high);

            if (f_high > 0.0) {
                LambdaDebug::log("Successfully bracketed lambda: [", lambda_low, ", ", lambda_high, "]");
                return true;
            }

            lambda_high *= increment_factor;

            // Don't exceed the theoretical upper bound
            if (lambda_high > lambda_upper_bound * LambdaCalculationConstants::LAMBDA_UPPER_BOUND_SAFETY) {
                lambda_high = lambda_upper_bound * LambdaCalculationConstants::LAMBDA_UPPER_BOUND_SAFETY;
                double f_high_final = MatrixAnalyzer::restriction_value_grouped(score_groups, lambda_high);

                LambdaDebug::log("At upper bound lambda_high = ", lambda_high, ", f(lambda_high) = ", f_high_final);

                if (f_high_final > 0.0) {
                    LambdaDebug::log("Successfully bracketed lambda at upper bound: [", lambda_low, ", ", lambda_high, "]");
                    return true;
                } else {
                    LambdaDebug::log("Failed to bracket lambda even at upper bound");
                    return false;
                }
            }
        }

        return false;
    }
};

class EigenvalueSolver {
public:
    // Solve for the dominant eigenvector of A where A[i,j] = exp(lambda * S[i,j])
    // Using the method suggested by collaborators: solve A * x = e_j for each unit vector e_j
    // Then combine the results to get the eigenvector
    static std::vector<double> solve_eigenvector_lu_decomposition(const Matrix& score_matrix, double lambda) {
        int n = score_matrix.row_dim;
        std::vector<double> empty_result;

        if (n <= 0 || lambda <= 0.0) {
            LambdaDebug::log("Invalid parameters: n=", n, ", lambda=", lambda);
            return empty_result;
        }

        LambdaDebug::log("Solving eigenvector using LU decomposition with lambda=", lambda);

        // Build the matrix A where A[i,j] = exp(lambda * S[i,j])
        Matrix A(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                A.at(i, j) = std::exp(lambda * score_matrix.at(i, j));
            }
        }

        // Initialize all background frequencies to zero
        std::vector<double> eigenvector(n, 0.0);

        // For each residue type j, set col[j] = 1.0 and solve A * x = col
        for (int j = 0; j < n; ++j) {
            std::vector<double> col(n, 0.0);
            col[j] = 1.0;	// Unit vector with 1.0 at position j

            std::vector<double> solution;
            if (A.solve_linear_system(col, solution)) {
                // Add this solution to our eigenvector estimate
                for (int i = 0; i < n; ++i) {
                    eigenvector[i] += solution[i];
                }
                LambdaDebug::log("Successfully solved for residue ", j);
            } else {
                LambdaDebug::log("Failed to solve for residue ", j);
                return empty_result;
            }
        }

        // Check if we got a valid result
        double sum_p = std::accumulate(eigenvector.begin(), eigenvector.end(), 0.0);
        if (sum_p <= 0.0) {
            LambdaDebug::log("Invalid eigenvector sum: ", sum_p);
            return empty_result;
        }

        // Normalize the eigenvector to be a probability vector
        for (int i = 0; i < n; ++i) {
            eigenvector[i] /= sum_p;
        }

        // Verify all probabilities are non-negative
        for (int i = 0; i < n; ++i) {
            if (eigenvector[i] < 0.0) {
                LambdaDebug::log("Negative probability at position ", i, ": ", eigenvector[i]);
                return empty_result;
            }
        }

        LambdaDebug::log("LU decomposition eigenvector solution completed successfully");
        return eigenvector;
    }
};

// ===== SIMPLE PROBABILITY SOLVER =====

class ProbabilitySolver {
public:
    // Algebraic method to solve for probabilities given scoring matrix and lambda
    // Returns a probability vector p where p = q and sum_{i,j} p[i] * p[j] * exp(lambda * S[i,j]) = 1
    // Uses eigenvalue decomposition for efficient algebraic solution
    static std::vector<double> solve_probabilities_from_lambda(
        const Matrix& score_matrix, double lambda,
        double tolerance = LambdaCalculationConstants::POWER_ITERATION_TOLERANCE,
        int max_iterations = LambdaCalculationConstants::MAX_POWER_ITERATIONS) {
        int n = score_matrix.row_dim;
        std::vector<double> empty_result;

        if (n <= 0 || lambda <= 0.0) {
            LambdaDebug::log("Invalid parameters: n=", n, ", lambda=", lambda);
            return empty_result;
        }

        LambdaDebug::log("Solving for probabilities algebraically with lambda=", lambda);

        // Build the exponential matrix A[i,j] = exp(lambda * S[i,j])
        Matrix A(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                A.at(i, j) = std::exp(lambda * score_matrix.at(i, j));
            }
        }

        // The constraint p = q with sum_{i,j} p[i] * p[j] * A[i,j] = 1 leads to:
        // p[i] = sum_j p[j] * A[i,j] / Z, where Z = sum_{i,j} p[i] * p[j] * A[i,j]
        // This is equivalent to finding the dominant eigenvector of A, then scaling

        // Use power iteration to find the dominant eigenvector
        std::vector<double> p(n, 1.0 / n);	// Start with uniform distribution

        for (int iter = 0; iter < max_iterations; ++iter) {
            // Compute A * p
            std::vector<double> Ap(n, 0.0);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    Ap[i] += A.at(i, j) * p[j];
                }
            }

            // Normalize to get unit vector (for numerical stability)
            double norm = std::sqrt(std::inner_product(Ap.begin(), Ap.end(), Ap.begin(), 0.0));
            if (norm <= 0.0) {
                LambdaDebug::log("Zero norm in power iteration");
                return empty_result;
            }

            for (int i = 0; i < n; ++i) {
                Ap[i] /= norm;
            }

            // Check convergence
            double max_change = 0.0;
            for (int i = 0; i < n; ++i) {
                double change = std::abs(Ap[i] - p[i]);
                max_change = std::max(max_change, change);
            }

            p = Ap;

            if (max_change < tolerance) {
                LambdaDebug::log("Power iteration converged after ", iter + 1, " iterations");
                break;
            }
        }

        // Now we have the dominant eigenvector, but we need to scale it to satisfy the constraint
        // The constraint is: sum_{i,j} p[i] * p[j] * A[i,j] = 1
        // First, normalize p to be a probability vector
        double sum_p = std::accumulate(p.begin(), p.end(), 0.0);
        if (sum_p <= 0.0) {
            LambdaDebug::log("Invalid probability sum: ", sum_p);
            return empty_result;
        }

        for (int i = 0; i < n; ++i) {
            p[i] /= sum_p;
        }

        // Compute the constraint value with the normalized probabilities
        double constraint_val = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                constraint_val += p[i] * p[j] * A.at(i, j);
            }
        }

        LambdaDebug::log("Constraint value before scaling: ", constraint_val);

        // Scale probabilities to satisfy the constraint exactly
        // If sum_{i,j} p[i] * p[j] * A[i,j] = C, then we need to scale by 1/sqrt(C)
        if (constraint_val <= 0.0) {
            LambdaDebug::log("Invalid constraint value: ", constraint_val);
            return empty_result;
        }

        double scale_factor = 1.0 / std::sqrt(constraint_val);
        for (int i = 0; i < n; ++i) {
            p[i] *= scale_factor;
        }

        // Verify the scaling worked
        double final_sum = std::accumulate(p.begin(), p.end(), 0.0);
        double final_constraint = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                final_constraint += p[i] * p[j] * A.at(i, j);
            }
        }

        LambdaDebug::log("Final probability sum: ", final_sum);
        LambdaDebug::log("Final constraint value: ", final_constraint);

        // Check if we have a valid probability vector
        if (final_sum <= 0.0) {
            LambdaDebug::log("Invalid final probability sum");
            return empty_result;
        }

        // Renormalize to ensure sum = 1
        for (int i = 0; i < n; ++i) {
            p[i] /= final_sum;
        }

        LambdaDebug::log("Algebraic solution completed successfully");
        return p;
    }
};

double calculate_lambda(const double **raw_mat_b, const int alpha_size,
                        std::vector<double> &final_p,
                        std::vector<double> &final_q,
                        double eps) {
    if (alpha_size <= 0) {
        LambdaDebug::log("Invalid alphabet size: ", alpha_size);
        return -1.0;
    }

    LambdaDebug::log("Starting lambda calculation with eigenvalue-solved background probabilities");
    LambdaDebug::log("alphabet_size=", alpha_size, ", eps=", eps);

    Matrix score_matrix(alpha_size, alpha_size);
    score_matrix.copy_from(raw_mat_b);

    // Validate and debug the score matrix
    if (!ScoreMatrixValidator::validate_score_matrix(score_matrix)) {
        LambdaDebug::log("Score matrix validation failed!");
        return -1.0;
    }

    LambdaDebug::log_matrix_info(score_matrix, "Score");

    // Check average score (should be negative for a proper scoring matrix)
    double avg_score = 0.0;
    for (int i = 0; i < alpha_size; ++i) {
        for (int j = 0; j < alpha_size; ++j) {
            avg_score += score_matrix.at(i, j);
        }
    }
    avg_score /= (alpha_size * alpha_size);
    LambdaDebug::log("Average score: ", avg_score);

    // Check if matrix is solvable
    double lambda_upper_bound = 0.0;
    if (!MatrixAnalyzer::matrix_solvable_for_lambda(score_matrix, lambda_upper_bound)) {
        LambdaDebug::log("Matrix is not solvable for lambda");
        return -1.0;
    }
    if (lambda_upper_bound <= 0.0) {
        LambdaDebug::log("Invalid lambda upper bound: ", lambda_upper_bound);
        return -1.0;
    }

    LambdaDebug::log("Lambda upper bound: ", lambda_upper_bound);

    // Initialize background probabilities (will be solved via eigenvalue decomposition)
    Probs p(alpha_size, true);	// start with uniform as initial guess
    Probs q(p);

    LambdaDebug::log("Initial background probabilities (uniform):");
    LambdaDebug::log_probabilities(p, "p (initial)");
    LambdaDebug::log_probabilities(q, "q (initial)");

    // Validate initial probabilities
    if (!LambdaDebug::validate_probability_distribution(p, "background p") ||
        !LambdaDebug::validate_probability_distribution(q, "background q")) {
        LambdaDebug::log("Invalid background probabilities");
        return -1.0;
    }

    // First, bracket lambda
    double lambda_low = LambdaCalculationConstants::LAMBDA_LOWER_BOUND;
    double lambda_high = LambdaCalculationConstants::BRACKET_START_LAMBDA;
    if (!LambdaFinder::bracket_lambda(score_matrix, p, q, lambda_low, lambda_high, lambda_upper_bound)) {
        LambdaDebug::log("Failed to bracket lambda");
        return -1.0;
    }

    LambdaDebug::log("Lambda bracket: [", lambda_low, ", ", lambda_high, "]");
    // Use Newton-Raphson method with LU decomposition to find the root
    // where restriction_value = 0
    LambdaDebug::log("Starting Newton-Raphson search for lambda with LU decomposition");

    double lambda_current = 0.5 * (lambda_low + lambda_high);	// Start in the middle
    double lambda = lambda_current;	// Initialize lambda for later use
    double epsilon = LambdaCalculationConstants::NEWTON_CONVERGENCE_TOLERANCE;
    int max_iterations = LambdaCalculationConstants::MAX_NEWTON_ITERATIONS;

    for (int iter = 0; iter < max_iterations; ++iter) {
        // Solve for probabilities at current lambda
        std::vector<double> solved_probs = EigenvalueSolver::solve_eigenvector_lu_decomposition(score_matrix, lambda_current);

        // If LU decomposition fails, try power iteration as fallback
        if (solved_probs.empty()) {
            LambdaDebug::log("LU decomposition failed at lambda=", lambda_current, ", trying power iteration");
            solved_probs = ProbabilitySolver::solve_probabilities_from_lambda(score_matrix, lambda_current);

            if (solved_probs.empty()) {
                LambdaDebug::log("Both methods failed at lambda=", lambda_current);
                return -1.0;
            }
        }

        // Convert to Probs objects
        Probs p_current(alpha_size);
        Probs q_current(alpha_size);
        for (int j = 0; j < alpha_size; ++j) {
            p_current[j] = solved_probs[j];
            q_current[j] = solved_probs[j]; // p = q for symmetric case
        }

        // Use score grouping optimization for faster restriction value computation
        std::vector<ScoreGroup> score_groups = MatrixAnalyzer::group_scores_by_value(score_matrix, p_current, q_current);

        // Compute function value and derivative using optimized grouped calculations
        double f_val = MatrixAnalyzer::restriction_value_grouped(score_groups, lambda_current);
        double f_prime = MatrixAnalyzer::restriction_value_derivative_grouped(score_groups, lambda_current);

        LambdaDebug::log("Newton-Raphson iter ", iter + 1, ": lambda=", lambda_current, " f(lambda)=", f_val, " f'(lambda)=", f_prime);

        // Check for convergence
        if (std::abs(f_val) < epsilon) {
            LambdaDebug::log("Newton-Raphson converged with lambda=", lambda_current, ", |f(lambda)|=", std::abs(f_val));
            p = p_current;
            q = q_current;
            lambda = lambda_current;
            break;
        }

        // Check for zero derivative
        if (std::abs(f_prime) < LambdaCalculationConstants::MATRIX_SINGULARITY_THRESHOLD) {
            LambdaDebug::log("Newton-Raphson: derivative too small, stopping");
            break;
        }

        // Newton-Raphson step with damping to stay within bounds
        double newton_step = -f_val / f_prime;
        double damping_factor = 1.0;
        double lambda_new = lambda_current + newton_step;

        // Apply damping if we go outside bounds
        if (lambda_new < lambda_low || lambda_new > lambda_high) {
            if (newton_step > 0) {
                damping_factor = std::min(1.0, (lambda_high - lambda_current) * LambdaCalculationConstants::DAMPING_FACTOR / newton_step);
            } else {
                damping_factor = std::min(1.0, (lambda_low - lambda_current) * LambdaCalculationConstants::DAMPING_FACTOR / newton_step);
            }
            lambda_new = lambda_current + damping_factor * newton_step;
        }
        LambdaDebug::log("Newton-Raphson step: ", newton_step, " damped: ", damping_factor * newton_step, " -> lambda=", lambda_new);

        // Check for convergence in lambda
        if (std::abs(lambda_new - lambda_current) < epsilon) {
            LambdaDebug::log("Newton-Raphson converged in lambda change: ", std::abs(lambda_new - lambda_current));
            // Use the current probabilities since we converged
            p = p_current;
            q = q_current;
            lambda = lambda_new;
            break;
        }
        lambda_current = lambda_new;

        // Final iteration check
        if (iter == max_iterations - 1) {
            LambdaDebug::log("Newton-Raphson: Maximum iterations reached");
            p = p_current;
            q = q_current;
            lambda = lambda_current;
        }
    }

    LambdaDebug::log_probabilities(p, "Updated p");
    LambdaDebug::log_probabilities(q, "Updated q");

    LambdaDebug::log("Final lambda=", lambda);

    // Validate final result using optimized calculation
    std::vector<ScoreGroup> final_score_groups = MatrixAnalyzer::group_scores_by_value(score_matrix, p, q);
    double final_restriction = MatrixAnalyzer::restriction_value_grouped(final_score_groups, lambda);
    LambdaDebug::log("Final restriction value: ", final_restriction);

    if (std::abs(final_restriction) > LambdaCalculationConstants::RESTRICTION_VALIDATION_TOLERANCE) {
        LambdaDebug::log("Failed to satisfy restriction equation, |f(lambda)| = ", std::abs(final_restriction));
        return -1.0;
    }

    if (lambda <= 0.0) {
        LambdaDebug::log("Lambda is non-positive: ", lambda);
        return -1.0;
    }

    // Return the solved background probabilities
    final_p = p.values;
    final_q = q.values;

    LambdaDebug::log("Lambda calculation successful");
    LambdaDebug::log("Final lambda=", lambda);
    LambdaDebug::log("Background probabilities solved via eigenvalue decomposition");

    return lambda;
}