#include <RcppEigen.h>
#include "utils/matrix_ops.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;

//' Enhance network directionality with C++ acceleration
//' 
//' @param X Network adjacency matrix
//' @param lambda Directionality parameter (0 to 1)
//' @return Enhanced network matrix
//' @export
// [[Rcpp::export]]
Eigen::SparseMatrix<double> strictDirectionCpp(
    SEXP XSEXP,
    double lambda = 1.0
) {
    // Convert SEXP to Eigen matrix
    MatrixXd X = as<MatrixXd>(XSEXP);
    
    if (lambda == 0.0) {
        // Fast path: no directionality enhancement, just convert to sparse
        return matrix_ops::denseToSparseOptimized(X, 1e-10);
    }
    
    int n = X.rows();
    
    // Optimized: Combine all operations in single pass to avoid multiple matrix traversals
    // This reduces 3 loops to 1 loop
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(n * n / 4);  // Estimate: ~25% non-zero after processing
    
    #pragma omp parallel
    {
        // Thread-local storage for triplets
        std::vector<Eigen::Triplet<double>> localTriplets;
        #ifdef _OPENMP
        int nThreads = omp_get_num_threads();
        #else
        int nThreads = 1;
        #endif
        localTriplets.reserve(n * n / (4 * nThreads));
        
        #pragma omp for collapse(2)
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double xij = X(i, j);
                double xji = X(j, i);
                
                // Apply directionality: keep stronger edge, zero weaker edge
                double sij = (std::abs(xij) < std::abs(xji)) ? 0.0 : xij;
                
                // Combine: (1-lambda)*X + lambda*S
                double result = (1.0 - lambda) * xij + lambda * sij;
                
                // Only store if non-zero
                if (std::abs(result) > 1e-10) {
                    localTriplets.emplace_back(i, j, result);
                }
            }
        }
        
        // Merge local triplets into global list
        #pragma omp critical
        {
            triplets.insert(triplets.end(), localTriplets.begin(), localTriplets.end());
        }
    }
    
    Eigen::SparseMatrix<double> sparse(n, n);
    sparse.setFromTriplets(triplets.begin(), triplets.end());
    sparse.makeCompressed();
    
    return sparse;
}

//' Perform virtual knockout on network
//' 
//' @param network Gene regulatory network matrix
//' @param geneIdx Index of gene to knockout (0-based)
//' @return Network with gene knocked out
//' @export
// [[Rcpp::export]]
Eigen::SparseMatrix<double> knockoutGeneCpp(
    SEXP networkSEXP,
    int geneIdx
) {
    // Convert SEXP to Eigen matrix
    MatrixXd network = as<MatrixXd>(networkSEXP);
    int n = network.rows();
    
    if (geneIdx < 0 || geneIdx >= n) {
        stop("Invalid gene index");
    }
    
    MatrixXd ko = network;
    
    // Set outgoing edges to zero (this gene cannot regulate others)
    ko.row(geneIdx).setZero();
    
    // Convert to sparse using optimized utility function
    return matrix_ops::denseToSparseOptimized(ko, 1e-10);
}

//' Convert sparse matrix to dense and transpose
//' 
//' @param X Sparse matrix
//' @return Dense transposed matrix
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd sparseToDenseTranspose(
    SEXP XSEXP
) {
    // Try to convert as sparse matrix first, fallback to dense
    Eigen::MatrixXd dense;
    if (Rf_isS4(XSEXP)) {
        // Sparse matrix
        Eigen::SparseMatrix<double> sparse = as<Eigen::SparseMatrix<double>>(XSEXP);
        dense = MatrixXd(sparse);
    } else {
        // Dense matrix
        dense = as<MatrixXd>(XSEXP);
    }
    return dense.transpose();
}

