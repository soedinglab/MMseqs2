#include "LambdaCalculation.h"
#include "Debug.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <numeric>
#include <set> // TODO: might try using an unordered set. We don't need repeated values, we do one insertion per iteration, and we do one full search per iteration.

// #define LambdaCalculation_DEBUG 1

struct Probs {
    int alphabet_size;
    std::vector<double> values;
    Probs(const int alphabet_size, bool uniform_init = false, bool random_init = false) : alphabet_size(alphabet_size), values(alphabet_size, 0.0){
        if (uniform_init){
            for (auto& val : values){
                val = 1.0 / static_cast<double>(alphabet_size);
            }
        }
        else if (random_init){
            double sum = 0.0;
            for (auto& val : values){
                val = static_cast<double>(std::rand() % 100);
                sum += static_cast<double>(val);
            }
            for (auto& val : values){
                val = val / sum;
            }
            
        }
    }
    double& operator[](int i) {
        return values[i];
    }

    const double& operator[](int i) const {
        return values[i];
    }
};



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

    void set_uniform() {
        int divisor = this->col_dim * this->row_dim;
        double val = 1.0 / static_cast<double>(divisor);
        for (int i = 0; i < this->row_dim; ++i){
            for (int j = 0; j < this->col_dim; ++j){
                this->at(i, j) = val;
            }
        }
    }

    void print_to_err(){
        for (int i = 0; i < this->row_dim; ++i){
            for (int j = 0; j < this->col_dim; ++j){
                Debug(Debug::ERROR) << this->at(i, j);
                if (j < this->col_dim -1){
                    Debug(Debug::ERROR) << "\t";
                }
            }
            Debug(Debug::ERROR) << "\n";
        }
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
    double lower_bound = 1e-10;
    double upper_bound = -1.0;
    double value = -1.0;
    double epsilon = 1e-121;
    int min_iters = 1;
    int max_iters = 100;
    bool converged = false;
    std::set<double> convergence_history;

    void bisection_search(const Matrix& score_matrix, const Probs p, const Probs q);
    void newton_refine(const Matrix& score_matrix,
                                const Probs p,
                                const Probs q);
};

inline void swap(double& a, double& b){
    double tmp = a;
    a = b;
    b = tmp;
}

inline double restriction_value(const Probs p,
                            const Probs q,
                            const double lambda,
                            const Matrix score_matrix){
    double sum = 0.0;
    for (int i = 0; i < score_matrix.row_dim; ++i){
        for (int j = 0; j < score_matrix.col_dim; ++j){
            // TODO: precompute p * q, though one extra mult probably isn't killing us.
            sum += p[i] * q[j] * exp(lambda * score_matrix.at(i, j));
        }
    }
    return sum - 1.0;
}

inline double restriction_value_first_derivative(const Probs p,
                                            const Probs q,
                                            const double lambda,
                                            const Matrix score_matrix){
    double sum = 0.0;
    for (int i = 0; i < score_matrix.row_dim; ++i){
        for (int j = 0; j < score_matrix.col_dim; ++j){
            // TODO: precompute p * q
            sum += p[i] * q[j] * score_matrix.at(i,j) * exp(lambda * score_matrix.at(i, j));
        }
    }
    return sum - 1.0;
}

inline bool sign_change(const double a, const double b){
    /**
     * Checks if a and b have different sign by multiplying
     * them and checking if the product is less than 0.
     */
    return a * b < 0.0;
}

template <typename T>
inline int signer(T val) {
    return (T(0) < val) - (val < T(0));
}

bool matrix_solvable_for_lambda(const Matrix& mat_b, double& lambda_upper_bound) {
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

    // We must find the smallest maximum positive score in a row and column
    std::vector<double> row_maxes (mat_b.row_dim, 0);
    std::vector<double> col_maxes (mat_b.col_dim, 0);
    
    // Check row/col conditions and find row maximums
    for (int i = 0; i < mat_b.row_dim; ++i) {
        bool row_pos = false;
        bool row_neg = false;
        for (int j = 0; j < mat_b.col_dim; ++j) {
            double val = mat_b.at(i, j);
            row_pos = row_pos | (val > 0.0);
            row_neg = row_neg | (val < 0.0);
            col_pos[j] = col_pos[j] || (val > 0.0);
            col_neg[j] = col_neg[j] || (val < 0.0);
            if (val > 0.0){
                row_maxes[i] = std::max(row_maxes[i], val);
                col_maxes[j] = std::max(col_maxes[j], val);
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
    // Take the minimum of the maxes, across both rows and columns, as the lambda upper bound.
    lambda_upper_bound = std::min(
        *(std::min_element(row_maxes.begin(), row_maxes.end())),
        *(std::min_element(col_maxes.begin(), col_maxes.end()))
    );

    
    if (count_positive > 0) {
        // Use the average positive score to estimate a reasonable upper bound
        // double avg_positive = sum_positive / count_positive;
        // For BLOSUM62-like matrices, lambda tends to be around 0.3-0.4
        // So we'll set the upper bound to 1.0 to ensure we include the expected range
        // lambda.upper_bound = 1.0;
    }
    
    #ifdef LambdaCalculation_DEBUG
    Debug(Debug::ERROR) << "Matrix is solvable. Upper bound: " << lambda_upper_bound << "\n";
    #endif
    return true;
}

void exponentiate_matrix(const Matrix& score_matrix, const double lambda, Matrix& result) {
    #ifdef LambdaCalculation_DEBUG
    Debug(Debug::ERROR) << "Computing with lambda = " << lambda << "\n";
    #endif
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
    
    // Base cases for 1x1 and 2x2 matrices
    if (mat.row_dim == 1) {
        return mat.at(0, 0);
    }
    if (mat.row_dim == 2) {
        return mat.at(0, 0) * mat.at(1, 1) - mat.at(0, 1) * mat.at(1, 0);
    }

    // Create working copy of matrix
    Matrix lu(mat.row_dim, mat.col_dim);
    std::copy(mat.values.begin(), mat.values.end(), lu.values.begin());

    // Track row permutations
    std::vector<int> pivot_indices(mat.row_dim);
    for (int i = 0; i < mat.row_dim; i++) {
        pivot_indices[i] = i;
    }

    double det = 1.0;
    int sign = 1;

    // LU decomposition with partial pivoting
    for (int k = 0; k < mat.row_dim - 1; k++) {
        // Find pivot
        double max_val = std::abs(lu.at(k, k));
        int pivot_row = k;
        
        for (int i = k + 1; i < mat.row_dim; i++) {
            double val = std::abs(lu.at(i, k));
            if (val > max_val) {
                max_val = val;
                pivot_row = i;
            }
        }

        // Check for singularity
        if (max_val < 1e-10) {
            return 0.0;
        }

        // Swap rows if necessary
        if (pivot_row != k) {
            for (int j = 0; j < mat.col_dim; j++) {
                std::swap(lu.at(k, j), lu.at(pivot_row, j));
            }
            std::swap(pivot_indices[k], pivot_indices[pivot_row]);
            sign *= -1;
        }

        // Eliminate below
        for (int i = k + 1; i < mat.row_dim; i++) {
            double factor = lu.at(i, k) / lu.at(k, k);
            lu.at(i, k) = factor;  // Store the multiplier
            
            for (int j = k + 1; j < mat.col_dim; j++) {
                lu.at(i, j) -= factor * lu.at(k, j);
            }
        }
    }

    // Compute determinant from diagonal elements
    for (int i = 0; i < mat.row_dim; i++) {
        det *= lu.at(i, i);
    }

    return det * sign;
}

void compute_joint_probabilities(
    const Matrix& score_matrix, 
    double lambda,
    Matrix& joint_prob_matrix,
    bool use_log_space = false) {
    // Protect from overflow by scaling the matrix
    const double max_exp = 700.0;  // log(DBL_MAX) ≈ 709
    double max_score = 0.0;
    for (const double& val : score_matrix.values) {
        max_score = std::max(max_score, std::abs(val));
    }
    
    if (lambda * max_score > max_exp) {
        // Scale lambda to prevent overflow
        lambda = max_exp / max_score; 
    }
    
    if (use_log_space) {
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
    } else {
        // Compute exponential values directly
        double Z = 0.0;
        for (int i = 0; i < score_matrix.row_dim; ++i) {
            for (int j = 0; j < score_matrix.col_dim; ++j) {
                joint_prob_matrix.at(i,j) = std::exp(lambda * score_matrix.at(i,j));
                Z += joint_prob_matrix.at(i,j);
            }
        }
        
        // Normalize
        for (double& val : joint_prob_matrix.values) {
            val /= Z;
        }
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
        
        
        // Check solvability and get upper bound
        double lambda_upper_bound = 0.0;
        if (!matrix_solvable_for_lambda(score_matrix, lambda_upper_bound)) {
            return -1.0;
        }
        
        // Initialize lambda to midpoint
        Lambda lambda_calc;
        lambda_calc.upper_bound = lambda_upper_bound;
        lambda_calc.value = lambda_calc.upper_bound / 2.0;
        

        // Initialize probabilities
        Probs local_p(alpha_size, false, true);
        Probs local_q(alpha_size, false, true);

        // Find lambda using bisection search and refine using Newton's method
        lambda_calc.bisection_search(score_matrix, local_p, local_q);
        lambda_calc.newton_refine(score_matrix, local_p, local_q);
        Matrix joint_probs(alpha_size, alpha_size);
        
        // Compute final probabilities
        compute_joint_probabilities(score_matrix, lambda_calc.value, joint_probs);
        compute_background_probabilities(joint_probs, p, q);
        
        // Normalize probabilities
        double sum_p = 0.0, sum_q = 0.0;
        for (int i = 0; i < alpha_size; ++i) {
            sum_p += new_p[i];
            sum_q += new_q[i];
        }
        
        for (int i = 0; i < alpha_size; ++i) {
            p[i] = new_p[i] / sum_p;
            q[i] = new_q[i] / sum_q;
        }
        
        iteration++;
    }
    
    return lambda;
}

double calculate_lambda_direct(
    const double** raw_mat_b,
    const int alpha_size,
    std::vector<double>& p,
    std::vector<double>& q) {
    if (alpha_size <= 0) {
        return -1.0;
    }
    
    // Convert input matrix to our Matrix class
    Matrix score_matrix(alpha_size, alpha_size);
    score_matrix.copy_from(raw_mat_b);
    
    // Initialize background probabilities with standard amino acid frequencies
    // These are the standard background frequencies used in BLOSUM62
    const double background_freqs[20] = {
        0.078,  // A
        0.051,  // R
        0.041,  // N
        0.052,  // D
        0.024,  // C
        0.034,  // Q
        0.059,  // E
        0.083,  // G
        0.025,  // H
        0.062,  // I
        0.092,  // L
        0.056,  // K
        0.024,  // M
        0.044,  // F
        0.043,  // P
        0.059,  // S
        0.055,  // T
        0.014,  // W
        0.034,  // Y
        0.072   // V
    };
    
    p.resize(alpha_size);
    q.resize(alpha_size);
    
    // Copy the background frequencies
    for (int i = 0; i < std::min(20, alpha_size); ++i) {
        p[i] = background_freqs[i];
        q[i] = background_freqs[i];
    }
    
    // If alpha_size > 20 (e.g., including gap), initialize remaining positions uniformly
    if (alpha_size > 20) {
        double remaining_prob = 1.0;
        for (int i = 0; i < 20; ++i) {
            remaining_prob -= background_freqs[i];
        }
        double uniform_prob = remaining_prob / (alpha_size - 20);
        for (int i = 20; i < alpha_size; ++i) {
            p[i] = uniform_prob;
            q[i] = uniform_prob;
        }
    }
    
    // Normalize to ensure probabilities sum to 1
    double sum_p = 0.0, sum_q = 0.0;
    for (int i = 0; i < alpha_size; ++i) {
        sum_p += p[i];
        sum_q += q[i];
    }
    for (int i = 0; i < alpha_size; ++i) {
        p[i] /= sum_p;
        q[i] /= sum_q;
    }
    
    const double epsilon = 1e-6;
    const int max_iterations = 100;
    int iteration = 0;
    double lambda = 0.0;
    double prev_lambda = 0.0;
    
    // Find maximum absolute score for scaling
    double max_score = 0.0;
    for (int i = 0; i < alpha_size; ++i) {
        for (int j = 0; j < alpha_size; ++j) {
            max_score = std::max(max_score, std::abs(score_matrix.at(i, j)));
        }
    }
    
    // Scale factor to prevent overflow
    const double scale = 1.0 / max_score;
    
    #ifdef LambdaCalculation_DEBUG
    Debug(Debug::ERROR) << "Initial background probabilities:\n";
    for (int i = 0; i < alpha_size; ++i) {
        Debug(Debug::ERROR) << "p[" << i << "]=" << p[i] << ", q[" << i << "]=" << q[i] << "\n";
    }
    #endif
    
    while (iteration < max_iterations) {
        // Find lambda that satisfies sum(p_i * p_j * exp(lambda * s_ij)) = 1
        double lambda_low = 0.0;
        double lambda_high = 1.0;
        double sum;
        
        // Bisection search for lambda
        for (int bisect_iter = 0; bisect_iter < 50; ++bisect_iter) {
            lambda = (lambda_low + lambda_high) / 2.0;
            
            // Calculate sum for current lambda, using scaled values
            sum = 0.0;
            double max_exp = -std::numeric_limits<double>::infinity();
            
            // First pass: find maximum exponent to prevent overflow
            for (int i = 0; i < alpha_size; ++i) {
                for (int j = 0; j < alpha_size; ++j) {
                    double exp_term = lambda * score_matrix.at(i, j) * scale;
                    max_exp = std::max(max_exp, exp_term);
                }
            }
            
            // Second pass: calculate sum with scaling
            for (int i = 0; i < alpha_size; ++i) {
                for (int j = 0; j < alpha_size; ++j) {
                    double exp_term = lambda * score_matrix.at(i, j) * scale;
                    sum += p[i] * q[j] * std::exp(exp_term - max_exp);
                }
            }
            
            // Scale back the sum
            sum *= std::exp(max_exp);
            
            #ifdef LambdaCalculation_DEBUG
            Debug(Debug::ERROR) << "Bisect " << bisect_iter << ": lambda=" << lambda << ", sum=" << sum << "\n";
            #endif
            
            if (std::abs(sum - 1.0) < epsilon) {
                break;
            }
            
            if (sum > 1.0) {
                lambda_high = lambda;
            } else {
                lambda_low = lambda;
            }
        }
        
        // Update background probabilities using current lambda
        std::vector<double> new_p(alpha_size, 0.0);
        std::vector<double> new_q(alpha_size, 0.0);
        
        // Calculate new background probabilities with scaling
        double max_exp = -std::numeric_limits<double>::infinity();
        
        // First pass: find maximum exponent
        for (int i = 0; i < alpha_size; ++i) {
            for (int j = 0; j < alpha_size; ++j) {
                double exp_term = lambda * score_matrix.at(i, j) * scale;
                max_exp = std::max(max_exp, exp_term);
            }
        }
        
        // Second pass: calculate probabilities
        for (int i = 0; i < alpha_size; ++i) {
            for (int j = 0; j < alpha_size; ++j) {
                double exp_term = lambda * score_matrix.at(i, j) * scale;
                double term = p[i] * q[j] * std::exp(exp_term - max_exp);
                new_p[i] += term;
                new_q[j] += term;
            }
        }
        
        // Normalize probabilities
        double sum_p = 0.0, sum_q = 0.0;
        for (int i = 0; i < alpha_size; ++i) {
            sum_p += new_p[i];
            sum_q += new_q[i];
        }
        
        for (int i = 0; i < alpha_size; ++i) {
            p[i] = new_p[i] / sum_p;
            q[i] = new_q[i] / sum_q;
        }
        
        #ifdef LambdaCalculation_DEBUG
        Debug(Debug::ERROR) << "Iteration " << iteration << ": lambda=" << lambda << ", sum=" << sum << "\n";
        #endif
        
        // Check for convergence of lambda
        if (std::abs(lambda - prev_lambda) < epsilon) {
            break;
        }
        
        prev_lambda = lambda;
        iteration++;
    }
    
    return lambda;
}

void Lambda::bisection_search(const Matrix& score_matrix,
                                const Probs p,
                                const Probs q) {
    if (upper_bound < 0.0){
        value = -1.0;
        return;
    }
    // Initialize to the ~middle of the interval
    value = upper_bound / 2.0;
    for (int iters = 0; iters < max_iters; ++iters){
        #ifdef LambdaCalculation_DEBUG
        Debug(Debug::ERROR) << "iters: " << iters + 1 << " LB: " << lower_bound << " Val: " << value << " UB: " << upper_bound << "\n";
        #endif

        // Calculate the value of the restriction condition.
        // TODO: we can short-circuit this by just checking UB, then sign, then and only then checking sum(lb)
        double val_sum = restriction_value(p, q, value, score_matrix);

        if (iters > min_iters && val_sum - 1.0 < epsilon){
            break;
        }
        double ub_sum = restriction_value(p, q, upper_bound, score_matrix);
        double lb_sum = restriction_value(p, q, lower_bound, score_matrix);

        if (sign_change(val_sum, ub_sum)){
            lower_bound = value;
        }
        if (sign_change(val_sum, lb_sum)){
            upper_bound = value;
        }
        value = (lower_bound + upper_bound) / 2.0;
    }
    #ifdef LambdaCalculation_DEBUG
    Debug(Debug::ERROR) << "initial value by bisection: " << value << "\n";
    #endif
}

void Lambda::newton_refine(const Matrix& score_matrix,
                                const Probs p,
                                const Probs q) {
    #ifdef LambdaCalculation_DEBUG
    Debug(Debug::ERROR) << "Running Newton Refinment.\n";
    #endif
    for (int iters = 0; iters < max_iters; ++iters){
        double f_lambda = restriction_value(p, q, value, score_matrix);
        double f_prime_lambda = restriction_value_first_derivative(p, q, value, score_matrix);
        #ifdef LambdaCalculation_DEBUG
        Debug(Debug::ERROR) << "iters: " << iters + 1 << " LB: " << lower_bound << " Val: " << value << " UB: " << upper_bound <<  " F " << f_lambda << " F' " << f_prime_lambda << "\n";
        #endif

        if (std::fabs(f_lambda) < epsilon){
            converged = true;
            #ifdef LambdaCalculation_DEBUG
            Debug(Debug::ERROR) << "converged.\n"; 
            #endif
            break;
        }
        if (std::find(convergence_history.begin(), convergence_history.end(), value) != convergence_history.end()){
            converged = true;
            #ifdef LambdaCalculation_DEBUG
            Debug(Debug::ERROR) << "converged by oscillation.\n";
            #endif
            break;
        }
        convergence_history.insert(value);
        // Protect from big oscillations by capping the convergence rate.
        double step = f_lambda / f_prime_lambda;
        step = signer(step) * std::min(std::fabs(f_lambda / f_prime_lambda), 0.01);
        value = value - step;
    }
    #ifdef LambdaCalculation_DEBUG
    Debug(Debug::ERROR) << "final value: " << value << "\n";
    #endif
}
