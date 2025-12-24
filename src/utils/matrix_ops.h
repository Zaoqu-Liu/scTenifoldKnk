#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

#include <RcppEigen.h>
#include <vector>
#include <algorithm>
#include <random>

// [[Rcpp::depends(RcppEigen)]]

namespace matrix_ops {

// Compute quantile threshold
// Optimized: Use nth_element instead of full sort (O(n) vs O(n log n))
inline double computeQuantile(const Eigen::MatrixXd& X, double q) {
    std::vector<double> values;
    values.reserve(X.size());
    
    // Collect non-zero absolute values
    const double* data = X.data();
    const int total_size = X.size();
    
    for (int i = 0; i < total_size; ++i) {
        double val = std::abs(data[i]);
        if (val > 1e-10) {
            values.push_back(val);
        }
    }
    
    if (values.empty()) return 0.0;
    if (values.size() == 1) return values[0];
    
    // Use nth_element for O(n) complexity instead of O(n log n) sort
    size_t idx = static_cast<size_t>(q * values.size());
    if (idx >= values.size()) idx = values.size() - 1;
    if (idx == 0) idx = 0;
    
    // Partial sort: elements before idx are <= values[idx]
    std::nth_element(values.begin(), values.begin() + idx, values.end());
    
    return values[idx];
}

// Convert dense to sparse matrix with threshold
inline Eigen::SparseMatrix<double> densityToSparse(
    const Eigen::MatrixXd& dense,
    double threshold
) {
    std::vector<Eigen::Triplet<double>> triplets;
    int nRows = dense.rows();
    int nCols = dense.cols();
    
    // Estimate non-zeros
    triplets.reserve(nRows * nCols / 10);
    
    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            double val = dense(i, j);
            if (std::abs(val) >= threshold) {
                triplets.emplace_back(i, j, val);
            }
        }
    }
    
    Eigen::SparseMatrix<double> sparse(nRows, nCols);
    sparse.setFromTriplets(triplets.begin(), triplets.end());
    sparse.makeCompressed();
    
    return sparse;
}

// Random subsample columns (cells) - optimized version
// Uses reservoir sampling for O(nSamples) instead of O(nCols) complexity
template<typename MatrixType>
MatrixType subsampleColumns(
    const MatrixType& X,
    int nSamples,
    unsigned int seed = 42
) {
    int nCols = X.cols();
    if (nSamples >= nCols) return X;
    
    // Optimized: Only generate nSamples random indices instead of shuffling all nCols
    std::vector<int> indices;
    indices.reserve(nSamples);
    
    std::mt19937 gen(seed);
    
    // Use reservoir sampling for efficient random selection
    // This is O(nSamples) instead of O(nCols)
    for (int i = 0; i < nCols; ++i) {
        if (i < nSamples) {
            // Fill reservoir
            indices.push_back(i);
        } else {
            // Randomly replace elements with decreasing probability
            std::uniform_int_distribution<> dis(0, i);
            int j = dis(gen);
            if (j < nSamples) {
                indices[j] = i;
            }
        }
    }
    
    // Sort indices for better cache locality when accessing columns
    std::sort(indices.begin(), indices.end());
    
    MatrixType result(X.rows(), nSamples);
    for (int i = 0; i < nSamples; ++i) {
        result.col(i) = X.col(indices[i]);
    }
    
    return result;
}

// Scale matrix (z-score normalization by rows)
inline void scaleRows(Eigen::MatrixXd& X) {
    #pragma omp parallel for if(X.rows() > 100)
    for (int i = 0; i < X.rows(); ++i) {
        double mean = X.row(i).mean();
        double stdDev = std::sqrt((X.row(i).array() - mean).square().mean());
        if (stdDev > 1e-10) {
            X.row(i).array() = (X.row(i).array() - mean) / stdDev;
        }
    }
}

// Convert dense matrix to sparse (extract non-zero elements)
// Optimized: Single pass with thread-local storage
inline Eigen::SparseMatrix<double> denseToSparseOptimized(
    const Eigen::MatrixXd& dense,
    double threshold = 1e-10
) {
    int nRows = dense.rows();
    int nCols = dense.cols();
    
    std::vector<Eigen::Triplet<double>> triplets;
    
    // Sequential version for small matrices, parallel for large
    if (nRows * nCols < 10000) {
        // Sequential
        triplets.reserve(nRows * nCols / 4);
        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < nCols; ++j) {
                double val = dense(i, j);
                if (std::abs(val) >= threshold) {
                    triplets.emplace_back(i, j, val);
                }
            }
        }
    } else {
        // Parallel with thread-local storage
        #pragma omp parallel
        {
            std::vector<Eigen::Triplet<double>> localTriplets;
            
            #ifdef _OPENMP
            int nThreads = omp_get_num_threads();
            localTriplets.reserve(nRows * nCols / (4 * nThreads));
            #else
            localTriplets.reserve(nRows * nCols / 4);
            #endif
            
            #pragma omp for collapse(2)
            for (int i = 0; i < nRows; ++i) {
                for (int j = 0; j < nCols; ++j) {
                    double val = dense(i, j);
                    if (std::abs(val) >= threshold) {
                        localTriplets.emplace_back(i, j, val);
                    }
                }
            }
            
            #pragma omp critical
            {
                triplets.insert(triplets.end(), localTriplets.begin(), localTriplets.end());
            }
        }
    }
    
    Eigen::SparseMatrix<double> sparse(nRows, nCols);
    sparse.setFromTriplets(triplets.begin(), triplets.end());
    sparse.makeCompressed();
    
    return sparse;
}

} // namespace matrix_ops

#endif // MATRIX_OPS_H

