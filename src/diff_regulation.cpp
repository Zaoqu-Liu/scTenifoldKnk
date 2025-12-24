#include <RcppEigen.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Box-Cox transformation with input validation
VectorXd boxCoxTransform(const VectorXd& x) {
    int n = x.size();
    
    // Critical: Box-Cox transformation requires strictly positive values
    double minVal = x.minCoeff();
    VectorXd xPositive = x;
    
    if (minVal <= 0) {
        // Add offset to make all values positive
        double offset = std::abs(minVal) + 1e-6;
        xPositive = x.array() + offset;
        Rcout << "  Note: Applied offset " << offset << " to ensure positive values for Box-Cox" << std::endl;
    }
    
    // Additional check: ensure no values are too close to zero
    if (xPositive.minCoeff() < 1e-10) {
        xPositive = xPositive.array().max(1e-10);
    }
    
    // Find optimal lambda using profile likelihood
    // Search range: -2 to 2
    const int nLambda = 1000;
    VectorXd lambdaValues = VectorXd::LinSpaced(nLambda, -2.0, 2.0);
    VectorXd logLikelihood(nLambda);
    
    double xMean = xPositive.mean();
    
    // Parallel search for optimal lambda
    #pragma omp parallel for
    for (int i = 0; i < nLambda; ++i) {
        double lambda = lambdaValues(i);
        
        if (std::abs(lambda) < 1e-6) {
            lambda = 1e-6;  // Avoid division by zero
        }
        
        VectorXd transformed(n);
        if (lambda < 0) {
            transformed = 1.0 / (xPositive.array().pow(lambda));
        } else {
            transformed = xPositive.array().pow(lambda);
        }
        
        double mean = transformed.mean();
        double variance = (transformed.array() - mean).square().mean();
        
        // Log-likelihood (simplified)
        double ll = -0.5 * n * std::log(variance);
        logLikelihood(i) = ll;
    }
    
    // Find lambda with maximum likelihood
    int maxIdx;
    logLikelihood.maxCoeff(&maxIdx);
    double optimalLambda = lambdaValues(maxIdx);
    
    // Apply optimal transformation
    VectorXd transformed(n);
    if (optimalLambda < 0) {
        transformed = 1.0 / (xPositive.array().pow(optimalLambda));
    } else {
        transformed = xPositive.array().pow(optimalLambda);
    }
    
    // Validate output: check for NaN or Inf
    for (int i = 0; i < n; ++i) {
        if (std::isnan(transformed(i)) || std::isinf(transformed(i))) {
            Rcout << "Warning: Box-Cox produced invalid values, using log transformation" << std::endl;
            return xPositive.array().log();
        }
    }
    
    return transformed;
}

// Compute chi-square p-values
VectorXd chiSquarePvalues(const VectorXd& chiSquareStats, int df = 1) {
    int n = chiSquareStats.size();
    VectorXd pvalues(n);
    
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        // P[X > x] for chi-square distribution
        // Using R's pchisq via Rcpp
        pvalues(i) = R::pchisq(chiSquareStats(i), df, 0, 0);  // lower.tail = FALSE
    }
    
    return pvalues;
}

// FDR correction (Benjamini-Hochberg)
VectorXd fdrCorrection(const VectorXd& pvalues) {
    int n = pvalues.size();
    
    // Create index vector
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    
    // Sort indices by p-values
    std::sort(indices.begin(), indices.end(), [&](int a, int b) {
        return pvalues(a) < pvalues(b);
    });
    
    VectorXd adjusted(n);
    double cumMin = 1.0;
    
    // BH procedure (from largest to smallest)
    for (int i = n - 1; i >= 0; --i) {
        int idx = indices[i];
        double bh = pvalues(idx) * n / (i + 1);
        cumMin = std::min(cumMin, bh);
        adjusted(idx) = std::min(cumMin, 1.0);
    }
    
    return adjusted;
}

//' Differential regulation analysis with C++ acceleration
//' 
//' @param manifoldOutput Matrix with WT and KO gene embeddings (2*nGenes x d)
//' @param geneNames Vector of gene names
//' @param gKO Name of knocked out gene
//' @return DataFrame with differential regulation statistics
//' @export
// [[Rcpp::export]]
DataFrame dRegulationCpp(
    SEXP manifoldOutputSEXP,
    const std::vector<std::string>& geneNames,
    const std::string& gKO
) {
    // Convert SEXP to Eigen matrix
    MatrixXd manifoldOutput = as<MatrixXd>(manifoldOutputSEXP);
    
    int nGenes = geneNames.size();
    int d = manifoldOutput.cols();
    
    // Validate dimensions
    if (manifoldOutput.rows() != 2 * nGenes) {
        stop("manifoldOutput should have 2*nGenes rows");
    }
    
    Rcout << "Computing distances for " << nGenes << " genes..." << std::endl;
    
    // Compute Euclidean distances between WT and KO embeddings
    // This is VECTORIZED - huge speedup over R's sapply!
    VectorXd distances(nGenes);
    
    #pragma omp parallel for
    for (int i = 0; i < nGenes; ++i) {
        // WT embedding (first nGenes rows)
        VectorXd x = manifoldOutput.row(i);
        
        // KO embedding (last nGenes rows)
        VectorXd y = manifoldOutput.row(i + nGenes);
        
        // Euclidean distance
        distances(i) = (x - y).norm();
    }
    
    Rcout << "Applying Box-Cox transformation..." << std::endl;
    
    // Box-Cox transformation
    VectorXd transformed = boxCoxTransform(distances);
    
    // Z-score normalization
    double mean = transformed.mean();
    double stdDev = std::sqrt((transformed.array() - mean).square().mean());
    
    // Prevent division by zero
    if (stdDev < 1e-10) {
        stdDev = 1.0;
        Rcout << "  Warning: Near-zero standard deviation, using stdDev=1" << std::endl;
    }
    
    VectorXd Z = (transformed.array() - mean) / stdDev;
    
    Rcout << "Computing fold changes and p-values..." << std::endl;
    
    // Sort genes by distance (descending)
    std::vector<int> sortedIdx(nGenes);
    std::iota(sortedIdx.begin(), sortedIdx.end(), 0);
    std::sort(sortedIdx.begin(), sortedIdx.end(), [&](int a, int b) {
        return distances(a) > distances(b);
    });
    
    // Find knockout gene index
    int koIdx = -1;
    for (int i = 0; i < nGenes; ++i) {
        if (geneNames[i] == gKO) {
            koIdx = i;
            break;
        }
    }
    
    // Compute FC: distance^2 / mean(other distances^2)
    VectorXd distSquared = distances.array().square();
    
    // Compute mean excluding KO gene
    double sumSquared = distSquared.sum();
    int denominator = nGenes - 1;
    
    if (koIdx >= 0) {
        sumSquared -= distSquared(koIdx);
        denominator = nGenes - 1;
    } else {
        denominator = nGenes;
    }
    
    double meanSquared = sumSquared / denominator;
    
    // Protect against division by zero
    if (meanSquared < 1e-10) {
        Rcout << "  Warning: All genes have near-zero distances, setting meanSquared=1" << std::endl;
        meanSquared = 1.0;
    }
    
    VectorXd FC = distSquared / meanSquared;
    
    // P-values using chi-square distribution
    VectorXd pvalues = chiSquarePvalues(FC, 1);
    
    // FDR correction
    VectorXd padj = fdrCorrection(pvalues);
    
    Rcout << "Creating output DataFrame..." << std::endl;
    
    // Create output DataFrame (sorted by distance)
    std::vector<std::string> sortedGenes(nGenes);
    VectorXd sortedDist(nGenes);
    VectorXd sortedZ(nGenes);
    VectorXd sortedFC(nGenes);
    VectorXd sortedPval(nGenes);
    VectorXd sortedPadj(nGenes);
    
    for (int i = 0; i < nGenes; ++i) {
        int idx = sortedIdx[i];
        sortedGenes[i] = geneNames[idx];
        sortedDist(i) = distances(idx);
        sortedZ(i) = Z(idx);
        sortedFC(i) = FC(idx);
        sortedPval(i) = pvalues(idx);
        sortedPadj(i) = padj(idx);
    }
    
    DataFrame result = DataFrame::create(
        Named("gene") = sortedGenes,
        Named("distance") = sortedDist,
        Named("Z") = sortedZ,
        Named("FC") = sortedFC,
        Named("p.value") = sortedPval,
        Named("p.adj") = sortedPadj
    );
    
    Rcout << "Differential regulation analysis completed!" << std::endl;
    
    return result;
}

