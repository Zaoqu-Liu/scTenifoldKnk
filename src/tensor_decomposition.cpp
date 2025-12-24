#include <RcppEigen.h>
#include <vector>
#include <cmath>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::MatrixXd;

// Compute Frobenius norm of difference between tensor and reconstruction
double computeReconstructionError(
    const std::vector<MatrixXd>& networks,
    const MatrixXd& A,  // nGenes x K
    const MatrixXd& B,  // nGenes x K  
    const MatrixXd& C   // nNet x K
) {
    int nNet = networks.size();
    double totalError = 0.0;
    
    #pragma omp parallel for reduction(+:totalError)
    for (int k = 0; k < nNet; ++k) {
        // Reconstruct network k: A * diag(C[k,:]) * B^T
        MatrixXd reconstructed = MatrixXd::Zero(A.rows(), B.rows());
        
        for (int r = 0; r < C.cols(); ++r) {
            reconstructed += C(k, r) * A.col(r) * B.col(r).transpose();
        }
        
        // Compute squared Frobenius norm of difference
        double error = (networks[k] - reconstructed).squaredNorm();
        totalError += error;
    }
    
    // Relative error
    double totalNorm = 0.0;
    for (const auto& net : networks) {
        totalNorm += net.squaredNorm();
    }
    
    return std::sqrt(totalError / totalNorm);
}

// Update one factor matrix in ALS
MatrixXd updateFactor(
    const std::vector<MatrixXd>& networks,
    const MatrixXd& A,
    const MatrixXd& B,
    const MatrixXd& C,
    int mode  // 0 = update A, 1 = update B, 2 = update C
) {
    int K = A.cols();
    
    if (mode == 0) {
        // Update A: solve for A given B, C
        // Correct CP-ALS: A ← X_(1) * (C ⊙ B) * ((C^T*C) ⊙ (B^T*B))^-1
        int nGenes = A.rows();
        MatrixXd newA = MatrixXd::Zero(nGenes, K);
        
        // Precompute Gram matrix: (C^T*C) ⊙ (B^T*B) (element-wise product)
        MatrixXd CtC = C.transpose() * C;  // K x K
        MatrixXd BtB = B.transpose() * B;  // K x K
        MatrixXd gram = CtC.array() * BtB.array();  // Hadamard product
        
        // Add regularization for numerical stability
        gram.diagonal().array() += 1e-10;
        
        // Precompute inverse for thread safety and maximum speed
        // gram is small (K×K, typically 3×3), inverse is very fast
        MatrixXd gram_inv = gram.inverse();
        
        // For each row of A (thread-safe with pre-computed inverse)
        #pragma omp parallel for
        for (int i = 0; i < nGenes; ++i) {
            // Compute right-hand side: sum_k X_k[i,:] * B * diag(C[k,:])
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(K);
            
            for (size_t k = 0; k < networks.size(); ++k) {
                // X_k[i,:] * B weighted by C[k,:]
                for (int r = 0; r < K; ++r) {
                    rhs(r) += C(k, r) * networks[k].row(i).dot(B.col(r));
                }
            }
            
            // Solve: a_i = gram^-1 * rhs (thread-safe matrix multiplication)
            newA.row(i) = gram_inv * rhs;
        }
        
        return newA;
        
    } else if (mode == 1) {
        // Update B: solve for B given A, C
        // Correct CP-ALS: B ← X_(2) * (C ⊙ A) * ((C^T*C) ⊙ (A^T*A))^-1
        int nGenes = B.rows();
        MatrixXd newB = MatrixXd::Zero(nGenes, K);
        
        // Precompute Gram matrix: (C^T*C) ⊙ (A^T*A)
        MatrixXd CtC = C.transpose() * C;  // K x K
        MatrixXd AtA = A.transpose() * A;  // K x K
        MatrixXd gram = CtC.array() * AtA.array();  // Hadamard product
        
        // Add regularization
        gram.diagonal().array() += 1e-10;
        
        // Precompute inverse for thread safety and speed
        MatrixXd gram_inv = gram.inverse();
        
        #pragma omp parallel for
        for (int j = 0; j < nGenes; ++j) {
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(K);
            
            for (size_t k = 0; k < networks.size(); ++k) {
                // X_k[:,j] * A weighted by C[k,:]
                for (int r = 0; r < K; ++r) {
                    rhs(r) += C(k, r) * networks[k].col(j).dot(A.col(r));
                }
            }
            
            newB.row(j) = gram_inv * rhs;
        }
        
        return newB;
        
    } else {
        // Update C: solve for C given A, B
        int nNet = C.rows();
        MatrixXd newC = MatrixXd::Zero(nNet, K);
        
        // Precompute Gramian (diagonal of Khatri-Rao product's Gram matrix)
        Eigen::VectorXd gramian_diag(K);
        for (int r = 0; r < K; ++r) {
            gramian_diag(r) = A.col(r).squaredNorm() * B.col(r).squaredNorm();
        }
        
        // Add regularization for numerical stability
        gramian_diag.array() += 1e-10;
        
        #pragma omp parallel for
        for (int k = 0; k < nNet; ++k) {
            const MatrixXd& Xk = networks[k];
            
            // Compute right-hand side for each factor
            Eigen::VectorXd rhs(K);
            for (int r = 0; r < K; ++r) {
                // Efficient computation: tr(Xk * A_r * B_r^T) = (B_r^T * Xk^T * A_r)
                rhs(r) = (B.col(r).transpose() * Xk.transpose() * A.col(r))(0);
            }
            
            // Solve: newC[k, :] * diag(gramian_diag) = rhs
            // This is equivalent to: newC[k, r] = rhs[r] / gramian_diag[r]
            newC.row(k) = (rhs.array() / gramian_diag.array()).matrix().transpose();
        }
        
        return newC;
    }
}

//' CP Tensor Decomposition with C++ acceleration
//' 
//' @param networkList List of network matrices
//' @param K Rank of decomposition
//' @param maxIter Maximum iterations
//' @param maxError Convergence threshold
//' @param nDecimal Decimal places for output
//' @return List with reconstructed network and factors
//' @export
// [[Rcpp::export]]
List tensorDecompositionCpp(
    const List& networkList,
    int K = 3,
    int maxIter = 1000,
    double maxError = 1e-5,
    int nDecimal = 3
) {
    int nNet = networkList.size();
    
    // Validate inputs
    if (nNet < 1) {
        stop("networkList must contain at least 1 network");
    }
    if (K < 1) {
        stop("K (rank) must be at least 1");
    }
    if (maxIter < 1) {
        stop("maxIter must be at least 1");
    }
    if (maxError <= 0) {
        stop("maxError must be positive");
    }
    
    // Convert R list to C++ vector of matrices
    // Note: We need dense matrices for tensor decomposition operations
    std::vector<MatrixXd> networks(nNet);
    int totalNonZeros = 0;
    int totalElements = 0;
    
    for (int i = 0; i < nNet; ++i) {
        SEXP elem = networkList[i];
        
        // Handle both dense and sparse matrices
        if (Rf_isS4(elem)) {
            // Sparse matrix - convert to dense for ALS operations
            Eigen::SparseMatrix<double> sparse = as<Eigen::SparseMatrix<double>>(elem);
            networks[i] = MatrixXd(sparse);
            totalNonZeros += sparse.nonZeros();
            totalElements += sparse.rows() * sparse.cols();
        } else {
            // Dense matrix
            networks[i] = as<MatrixXd>(elem);
        }
    }
    
    if (totalNonZeros > 0) {
        double sparsity = 100.0 * (1.0 - (double)totalNonZeros / totalElements);
        Rcout << "  Network sparsity: " << sparsity << "% zeros" << std::endl;
    }
    
    int nGenes = networks[0].rows();
    
    // Validate K relative to nGenes
    if (K >= nGenes) {
        Rcout << "  Warning: K (" << K << ") >= nGenes (" << nGenes 
              << "), reducing K to " << (nGenes - 1) << std::endl;
        K = nGenes - 1;
    }
    
    Rcout << "Starting CP decomposition: " << nGenes << " genes, " 
          << nNet << " networks, rank " << K << std::endl;
    
    // Initialize factor matrices randomly
    MatrixXd A = MatrixXd::Random(nGenes, K);
    MatrixXd B = MatrixXd::Random(nGenes, K);
    MatrixXd C = MatrixXd::Random(nNet, K);
    
    // Normalize
    for (int k = 0; k < K; ++k) {
        A.col(k).normalize();
        B.col(k).normalize();
        C.col(k).normalize();
    }
    
    double prevError = 1e10;
    int iter;
    
    // ALS iterations
    for (iter = 0; iter < maxIter; ++iter) {
        // Update each factor matrix
        A = updateFactor(networks, A, B, C, 0);
        B = updateFactor(networks, A, B, C, 1);
        C = updateFactor(networks, A, B, C, 2);
        
        // Check convergence every 10 iterations
        if (iter % 10 == 0) {
            double error = computeReconstructionError(networks, A, B, C);
            
            if (iter % 50 == 0) {
                Rcout << "  Iteration " << iter << ", error: " << error << std::endl;
            }
            
            if (error < maxError) {
                Rcout << "Converged at iteration " << iter << std::endl;
                break;
            }
            
            // Early stopping: if error reduction < 0.1%, stop
            if (iter > 20 && (prevError - error) / prevError < 0.001) {
                Rcout << "Early stopping at iteration " << iter << std::endl;
                break;
            }
            
            prevError = error;
        }
    }
    
    Rcout << "Reconstructing network..." << std::endl;
    
    // Optimized: Precompute all K rank-one matrices A_r * B_r^T
    // This avoids recomputing them nNet times
    std::vector<MatrixXd> rankOneMatrices(K);
    
    #pragma omp parallel for if(K > 1)
    for (int r = 0; r < K; ++r) {
        rankOneMatrices[r].noalias() = A.col(r) * B.col(r).transpose();
    }
    
    // Reconstruct weighted average network
    MatrixXd reconstructed = MatrixXd::Zero(nGenes, nGenes);
    
    for (int k = 0; k < nNet; ++k) {
        double weight = C.row(k).sum() / (K * nNet);  // Normalize by K and nNet together
        
        for (int r = 0; r < K; ++r) {
            reconstructed.noalias() += weight * C(k, r) * rankOneMatrices[r];
        }
    }
    
    // Round to specified decimal places (only if nDecimal >= 4 to avoid losing edges)
    // Note: nDecimal=3 would zero out values < 0.001, losing network structure
    if (nDecimal >= 4) {
        double scale = std::pow(10.0, nDecimal);
        reconstructed = (reconstructed * scale).array().round() / scale;
        Rcout << "  Rounded to " << nDecimal << " decimal places" << std::endl;
    } else if (nDecimal > 0 && nDecimal < 4) {
        Rcout << "  Note: nDecimal=" << nDecimal << " may lose network edges, skipping rounding" << std::endl;
    }
    
    Rcout << "Tensor decomposition completed!" << std::endl;
    
    return List::create(
        Named("X") = reconstructed,
        Named("A") = A,
        Named("B") = B,
        Named("C") = C,
        Named("iterations") = iter,
        Named("error") = prevError
    );
}

