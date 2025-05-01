#include "LambdaCalculation.h"
#include "Debug.h"
#include <algorithm>
#include <cmath>
#include <vector>

#define LambdaCalculation_DEBUG false

// A struct to represent a matrix using a 1-D vector,
// which provides a nice interface and better memory access pattern
struct Matrix {
    int row_dim;
    int col_dim;
    std::vector<double> values;
    
    Matrix(const int row_dim, const int col_dim) 
        : row_dim(row_dim), col_dim(col_dim), values(row_dim * col_dim, 0.0) {}
    
    double& at(int i, int j) {
        return values[i * col_dim + j];
    }
    
    const double& at(int i, int j) const {
        return values[i * col_dim + j];
    }

    // Copy from a raw pointer matrix
    void copy_from(const double** mat) {
        for (int i = 0; i < row_dim; ++i) {
            for (int j = 0; j < col_dim; ++j) {
                at(i, j) = mat[i][j];
            }
        }
    }

    // Verify all values are within reasonable numeric limits of 0 <= abs(val) <= 10^8
    bool is_well_conditioned() const {
        double max_val = 0.0;
        double min_val = std::numeric_limits<double>::max();
        for (const double& val : values) {
            if (std::abs(val) > 0.0) {
                max_val = std::max(max_val, std::abs(val));
                min_val = std::min(min_val, std::abs(val));
            }
        }
        return (min_val > 0.0) && (max_val / min_val < 1e8);
    }

    // Ensure the matrix has row_dim and col_dim > 0,
    // that its values fit in a size_t, and that there are no NaNs / infs
    bool validate() const {
        if (row_dim <= 0 || col_dim <= 0) {
            return false;
        }
        if (values.size() != static_cast<size_t>(row_dim * col_dim)) {
            return false;
        }
        
        // Check for NaN or Inf values
        return std::none_of(values.begin(), values.end(), 
                           [](double v) { return std::isnan(v) || std::isinf(v); });
    }
};

// Helper struct to compute Karlin-Altschul's Lambda,
// which is a single scalar value that is used downstream to
// normalize things like the e-value and bit score.
// 
// Lambda is first guessed, then honed with bisection to a smaller range,
// then pinned using Newton's method. Newton's method is quadratic, which is why
// we first use bisection.
struct Lambda {
    double lower_bound = 0.0;
    double upper_bound = -1.0;
    double value = -1.0;
    double epsilon = 1e-3;
    std::vector<double> convergence_history;

    void bisection_search(const Matrix& score_matrix, Matrix& A_matrix);
    void newton_refine(const Matrix& score_matrix, Matrix& A_matrix);

    void record_iteration(double det) {
        convergence_history.push_back(std::abs(det));
    }
    
    bool is_converging() const {
        if (convergence_history.size() < 3) {
            return true;
        }
        return convergence_history.back() < convergence_history[convergence_history.size()-2];
    }
};

bool matrix_solvable_for_lambda(const Matrix& mat_b, Lambda& lambda) {
    #ifdef LambdaCalculation_DEBUG
    // Print matrix values
    Debug(Debug::ERROR) << "Matrix values:" << "\n";
    for (int i = 0; i < mat_b.row_dim; ++i) {
        for (int j = 0; j < mat_b.col_dim; ++j) {
            Debug(Debug::ERROR) << mat_b.at(i, j) << "\t";
        }
        Debug(Debug::ERROR) << "\n";
    }
    Debug(Debug::ERROR) << "\n";
    #endif
    
    // Add check for zero rows/columns
    std::vector<double> row_sums(mat_b.row_dim, 0.0);
    std::vector<double> col_sums(mat_b.col_dim, 0.0);
    
    for (int i = 0; i < mat_b.row_dim; ++i) {
        for (int j = 0; j < mat_b.col_dim; ++j) {
            row_sums[i] += mat_b.at(i,j);
            col_sums[j] += mat_b.at(i,j);
        }
    }
    
    // Check for all-zero rows/columns
    if (std::any_of(row_sums.begin(), row_sums.end(), 
                    [](double sum) { return std::abs(sum) < 1e-10; }) ||
        std::any_of(col_sums.begin(), col_sums.end(), 
                    [](double sum) { return std::abs(sum) < 1e-10; })) {
        Debug(Debug::ERROR) << "Failed: Found zero row/column sum\n";
        return false;
    }
    
    std::vector<bool> col_pos(mat_b.col_dim, false);
    std::vector<bool> col_neg(mat_b.col_dim, false);
    std::vector<double> row_max_pos(mat_b.row_dim, -1.0);
    std::vector<double> col_max_pos(mat_b.col_dim, -1.0);
    
    // First pass: check row conditions and find row maximums
    for (int i = 0; i < mat_b.row_dim; ++i) {
        bool row_pos = false;
        bool row_neg = false;
        for (int j = 0; j < mat_b.col_dim; ++j) {
            double val = mat_b.at(i, j);
            row_pos = row_pos | (val > 0.0);
            row_neg = row_neg | (val < 0.0);
            col_pos[j] = col_pos[j] || (val > 0.0);
            col_neg[j] = col_neg[j] || (val < 0.0);
            
            // Track maximum positive values
            if (val > 0.0) {
                row_max_pos[i] = (row_max_pos[i] < 0.0) ? val : std::min(val, row_max_pos[i]);
                col_max_pos[j] = (col_max_pos[j] < 0.0) ? val : std::min(val, col_max_pos[j]);
            }
        }
        if (!(row_pos && row_neg)) {
            Debug(Debug::ERROR) << "Failed: Row " << i << " does not have both positive and negative values\n";
            return false;
        }
    }
    
    // Check column conditions
    if (!(std::all_of(col_pos.begin(), col_pos.end(), [](bool v) { return v; }) &&
          std::all_of(col_neg.begin(), col_neg.end(), [](bool v) { return v; }))) {
        Debug(Debug::ERROR) << "Failed: Not all columns have both positive and negative values\n";
        return false;
    }
    
    // Find the largest positive value as the upper bound
    lambda.upper_bound = -1.0;
    for (int i = 0; i < mat_b.row_dim; ++i) {
        for (int j = 0; j < mat_b.col_dim; ++j) {
            double val = mat_b.at(i, j);
            if (val > 0.0) {
                lambda.upper_bound = (lambda.upper_bound < 0.0) ? val : std::max(lambda.upper_bound, val);
            }
        }
    }
    
    #ifdef LambdaCalculation_DEBUG
    Debug(Debug::ERROR) << "Matrix is solvable. Upper bound: " << lambda.upper_bound << "\n";
    #endif
    return true;
}

void exponentiate_matrix(const Matrix& score_matrix, const double lambda, Matrix& result) {
    for (int i = 0; i < score_matrix.row_dim; ++i) {
        for (int j = 0; j < score_matrix.col_dim; ++j) {
            result.at(i, j) = std::exp(lambda * score_matrix.at(i, j));
        }
    }
}

double compute_determinant(const Matrix& mat) {
    if (mat.row_dim != mat.col_dim) {
        return 0.0;
    }

    // Create a copy for LU decomposition
    Matrix work_matrix(mat.row_dim, mat.col_dim);
    std::copy(mat.values.begin(), mat.values.end(), work_matrix.values.begin());

    double determinant = 1.0;
    
    // Add scaling factor for better numerical stability
    double scale = 0.0;
    for (const double& val : mat.values) {
        scale = std::max(scale, std::abs(val));
    }
    if (scale > 0.0) {
        for (double& val : work_matrix.values) {
            val /= scale;
        }
        // Adjust determinant for scaling
        determinant = std::pow(scale, mat.row_dim);
    }

    // Perform LU decomposition with partial pivoting
    std::vector<int> pivot_indices(mat.row_dim);
    for (int i = 0; i < mat.row_dim; i++) {
        pivot_indices[i] = i;
    }

    for (int i = 0; i < mat.row_dim; i++) {
        // Find pivot
        double pivot = work_matrix.at(i, i);
        int pivot_row = i;
        
        for (int j = i + 1; j < mat.row_dim; j++) {
            if (std::abs(work_matrix.at(j, i)) > std::abs(pivot)) {
                pivot = work_matrix.at(j, i);
                pivot_row = j;
            }
        }

        if (std::abs(pivot) < 1e-10) {
            return 0.0;
        }

        // Swap rows if necessary
        if (pivot_row != i) {
            for (int j = 0; j < mat.row_dim; j++) {
                std::swap(work_matrix.at(i, j), work_matrix.at(pivot_row, j));
            }
            std::swap(pivot_indices[i], pivot_indices[pivot_row]);
            determinant *= -1.0;
        }

        determinant *= work_matrix.at(i, i);

        // Eliminate below
        for (int j = i + 1; j < mat.row_dim; j++) {
            double factor = work_matrix.at(j, i) / work_matrix.at(i, i);
            work_matrix.at(j, i) = factor;  // Store the multiplier
            for (int k = i + 1; k < mat.row_dim; k++) {
                work_matrix.at(j, k) -= factor * work_matrix.at(i, k);
            }
        }
    }

    return determinant;
}

void compute_joint_probabilities(
    const Matrix& score_matrix, 
    double lambda,
    Matrix& joint_prob_matrix) {
    // Protect from overflow by scaling the matrix to log space
    const double max_exp = 700.0;  // log(DBL_MAX) ≈ 709
    double max_score = 0.0;
    for (const double& val : score_matrix.values) {
        max_score = std::max(max_score, std::abs(val));
    }
    
    if (lambda * max_score > max_exp) {
        // Scale lambda to prevent overflow
        lambda = max_exp / max_score; 
    }
    
    // Compute in log space first
    double max_val = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < score_matrix.row_dim; ++i) {
        for (int j = 0; j < score_matrix.col_dim; ++j) {
            double log_val = lambda * score_matrix.at(i,j);
            max_val = std::max(max_val, log_val);
            joint_prob_matrix.at(i,j) = log_val;
        }
    }
    
    // Shift and exponentiate
    double Z = 0.0;
    for (int i = 0; i < score_matrix.row_dim; ++i) {
        for (int j = 0; j < score_matrix.col_dim; ++j) {
            joint_prob_matrix.at(i,j) = std::exp(joint_prob_matrix.at(i,j) - max_val);
            Z += joint_prob_matrix.at(i,j);
        }
    }
    
    // Normalize
    for (double& val : joint_prob_matrix.values) {
        val /= Z;
    }
}

void compute_background_probabilities(
    const Matrix& joint_prob_matrix,
    std::vector<double>& p,
    std::vector<double>& q) {
    p.resize(joint_prob_matrix.row_dim, 0.0);
    q.resize(joint_prob_matrix.col_dim, 0.0);
    
    for (int i = 0; i < joint_prob_matrix.row_dim; ++i) {
        for (int j = 0; j < joint_prob_matrix.col_dim; ++j) {
            p[i] += joint_prob_matrix.at(i, j);
            q[j] += joint_prob_matrix.at(i, j);
        }
    }
}

double calculate_lambda(
    const double** raw_mat_b,
    const int alpha_size,
    std::vector<double>& p,
    std::vector<double>& q) {
    if (alpha_size <= 0) {
        return -1.0;
    }
    
    // error: cannot use 'try' with exceptions disabled
    // try {
        // Convert input matrix to our Matrix class
        Matrix score_matrix(alpha_size, alpha_size);
        score_matrix.copy_from(raw_mat_b);
        
        Lambda lambda_calc;
        
        // Check solvability and get upper bound
        if (!matrix_solvable_for_lambda(score_matrix, lambda_calc)) {
            return -1.0;
        }
        
        // Initialize lambda to midpoint
        lambda_calc.value = lambda_calc.upper_bound / 2.0;
        
        // Create matrices for the scaled determinant and the joint probabilities
        Matrix A_matrix(alpha_size, alpha_size);
        Matrix joint_probs(alpha_size, alpha_size);
        
        // Find lambda using bisection search
        lambda_calc.bisection_search(score_matrix, A_matrix);
        
        // Refine using Newton's method
        lambda_calc.newton_refine(score_matrix, A_matrix);
        
        // Compute final probabilities
        compute_joint_probabilities(score_matrix, lambda_calc.value, joint_probs);
        compute_background_probabilities(joint_probs, p, q);
        
        return lambda_calc.value;
    // } catch (const std::exception& e) {
    //     // Handle numerical errors gracefully
    //     return -1.0;
    // }
}

void Lambda::bisection_search(const Matrix& score_matrix, Matrix& A_matrix) {
    // Check both bounds to ensure proper bracketing
    exponentiate_matrix(score_matrix, lower_bound, A_matrix);
    double det_low = compute_determinant(A_matrix);
    
    exponentiate_matrix(score_matrix, upper_bound, A_matrix);
    double det_high = compute_determinant(A_matrix);
    
    // Verify opposite signs
    if (det_low * det_high >= 0.0) {
        // No root bracketed, try to adjust bounds
        if (std::abs(det_low) < std::abs(det_high)) {
            upper_bound = (upper_bound + lower_bound) / 2.0;
        } else {
            lower_bound = (upper_bound + lower_bound) / 2.0;
        }
        
        // Recompute determinants with adjusted bounds
        exponentiate_matrix(score_matrix, lower_bound, A_matrix);
        det_low = compute_determinant(A_matrix);
        
        exponentiate_matrix(score_matrix, upper_bound, A_matrix);
        det_high = compute_determinant(A_matrix);
        
        // If still no sign change, try a more aggressive approach
        if (det_low * det_high >= 0.0) {
            // Try a wider range
            upper_bound *= 2.0;
            exponentiate_matrix(score_matrix, upper_bound, A_matrix);
            det_high = compute_determinant(A_matrix);
            
            if (det_low * det_high >= 0.0) {
                // If still no sign change, try a different approach
                // Use a fixed upper bound based on the maximum score
                double max_score = 0.0;
                for (const double& val : score_matrix.values) {
                    max_score = std::max(max_score, std::abs(val));
                }
                upper_bound = 1.0 / max_score;  // A reasonable starting point
                
                exponentiate_matrix(score_matrix, upper_bound, A_matrix);
                det_high = compute_determinant(A_matrix);
            }
        }
    }
    
    value = (lower_bound + upper_bound) / 2.0;
    exponentiate_matrix(score_matrix, value, A_matrix);
    double prev_det = compute_determinant(A_matrix);
    
    // We'll use a slightly larger epsilon for bisection
    const double bisection_epsilon = epsilon * 10.0;
    
    // Add early convergence check
    if (std::abs(prev_det) < bisection_epsilon) {
        return;  // Already close enough to zero
    }
    
    int max_iterations = 100;  // Prevent infinite loops
    int iteration = 0;
    
    while ((upper_bound - lower_bound) > bisection_epsilon && iteration < max_iterations) {
        // Compute A(λ) = exp(λsij)
        exponentiate_matrix(score_matrix, value, A_matrix);
        
        // Compute determinant
        double det = compute_determinant(A_matrix);
        
        // Record for convergence tracking
        record_iteration(det);
        
        // Update bounds based on determinant sign change
        if (det * prev_det < 0.0) {
            // Sign change detected - we've bracketed the root
            upper_bound = value;
        } else {
            lower_bound = value;
        }
        
        prev_det = det;
        value = (lower_bound + upper_bound) / 2.0;
        
        // Check for convergence
        if (std::abs(det) < bisection_epsilon) {
            break;
        }
        
        iteration++;
    }
}

void Lambda::newton_refine(const Matrix& score_matrix, Matrix& A_matrix) {
    const int max_iterations = 30;  // Increased from 20
    const double newton_epsilon = epsilon;  // Tighter tolerance for final refinement
    Matrix A_prime(score_matrix.row_dim, score_matrix.col_dim);  // For derivative
    
    // Record initial value for convergence tracking
    double initial_value = value;
    
    for (int iter = 0; iter < max_iterations; ++iter) {
        // Compute A(λ) and A'(λ)
        for (int i = 0; i < score_matrix.row_dim; ++i) {
            for (int j = 0; j < score_matrix.col_dim; ++j) {
                double exp_term = std::exp(value * score_matrix.at(i, j));
                A_matrix.at(i, j) = exp_term;
                A_prime.at(i, j) = score_matrix.at(i, j) * exp_term;
            }
        }
        
        // Compute determinant
        double det = compute_determinant(A_matrix);
        
        // Record for convergence tracking
        record_iteration(det);
        
        // Check if we've converged
        if (std::abs(det) < newton_epsilon) {
            break;
        }
        
        // Compute derivative analytically:
        // f'(λ) = sum(sij * exp(λsij)) * det(A(λ))
        double trace_term = 0.0;
        for (int i = 0; i < score_matrix.row_dim; ++i) {
            for (int j = 0; j < score_matrix.col_dim; ++j) {
                trace_term += score_matrix.at(i,j) * A_matrix.at(i,j);
            }
        }
        double derivative = trace_term * det;
        
        // Check for near-zero derivative
        if (std::abs(derivative) < 1e-10) {
            // Derivative too small, use a small step in the right direction
            derivative = (det > 0) ? 1e-10 : -1e-10;
        }
        
        // Newton step
        double delta = det / derivative;
        
        // Limit step size to prevent overshooting
        double max_step = (upper_bound - lower_bound) * 0.1;
        if (std::abs(delta) > max_step) {
            delta = (delta > 0 ? max_step : -max_step);
        }
        
        // Update lambda
        double new_value = value - delta;
        
        // Keep lambda within bounds
        new_value = std::max(lower_bound, std::min(upper_bound, new_value));
        
        // Check for convergence
        if (std::abs(new_value - value) < newton_epsilon) {
            value = new_value;
            break;
        }
        
        // Check if we're oscillating
        if (iter > 5 && !is_converging()) {
            // If not converging, try a different approach
            // Use a fixed step size in the direction of the root
            double step = (det > 0) ? -0.01 : 0.01;
            new_value = value + step;
            
            // Keep within bounds
            new_value = std::max(lower_bound, std::min(upper_bound, new_value));
        }
        
        value = new_value;
    }
    
    // If we've moved too far from the initial value, revert
    if (std::abs(value - initial_value) > 0.5) {
        value = initial_value;
    }
}
