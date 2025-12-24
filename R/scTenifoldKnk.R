#' @export scTenifoldKnk
#' @title Optimized Virtual Knockout Experiment
#' @description Performs in-silico gene knockout experiments using optimized algorithms
#' @param countMatrix Gene expression count matrix (genes x cells)
#' @param gKO Gene(s) to knockout
#' @param qc_mtThreshold Maximum mitochondrial read ratio
#' @param qc_minLSize Minimum library size
#' @param nc_lambda Lambda parameter for strict directionality
#' @param nc_nNet Number of networks to generate
#' @param nc_nCells Number of cells per network
#' @param nc_nComp Number of principal components
#' @param nc_scaleScores Whether to scale network scores
#' @param nc_symmetric Whether to make network symmetric
#' @param nc_q Quantile threshold for edge filtering
#' @param td_K Number of tensor components
#' @param td_maxIter Maximum tensor decomposition iterations
#' @param td_maxError Tensor decomposition error tolerance
#' @param td_nDecimal Number of decimal places
#' @param ma_nDim Number of manifold dimensions
#' @param nCores Number of cores for parallel processing
#' @param verbose Whether to print progress information
#' @return A list containing tensor networks, manifold alignment, and differential regulation results
#' @examples 
#' # Loading single-cell data
#' scRNAseq <- system.file("single-cell/example.csv", package="scTenifoldKnk")
#' scRNAseq <- read.csv(scRNAseq, row.names = 1)
#' 
#' # Running scTenifoldKnk
#' result <- scTenifoldKnk(countMatrix = scRNAseq, gKO = 'G100', qc_minLSize = 0)
scTenifoldKnk <- function(countMatrix, 
                           gKO = NULL, 
                           qc_mtThreshold = 0.1, 
                           qc_minLSize = 1000, 
                           nc_lambda = 0, 
                           nc_nNet = 10, 
                           nc_nCells = 500, 
                           nc_nComp = 3,
                           nc_scaleScores = TRUE, 
                           nc_symmetric = FALSE, 
                           nc_q = 0.9, 
                           td_K = 3, 
                           td_maxIter = 1000,
                           td_maxError = 1e-05, 
                           td_nDecimal = 3, 
                           ma_nDim = 2,
                           nCores = NULL,
                           verbose = TRUE){
  
  if(verbose) cat("=== scTenifoldKnk: Virtual Knockout Analysis ===\n\n")
  
  # Quality control
  if(verbose) cat("Step 1/6: Quality control...\n")
  countMatrix <- scQC(countMatrix, mtThreshold = qc_mtThreshold, minLSize = qc_minLSize)
  
  # Gene filtering
  if(ncol(countMatrix) > 500){
    countMatrix <- countMatrix[Matrix::rowMeans(countMatrix != 0) >= 0.05, , drop = FALSE]
  } else {
    countMatrix <- countMatrix[Matrix::rowSums(countMatrix != 0) >= 25, , drop = FALSE]
  }
  
  if(verbose) cat(sprintf("  Retained %d genes and %d cells after QC\n\n", 
                         nrow(countMatrix), ncol(countMatrix)))
  
  # Network construction
  if(verbose) cat("Step 2/6: Constructing gene regulatory networks (optimized & parallelized)...\n")
  WT <- makeNetworksFast(
    X = countMatrix, 
    q = nc_q, 
    nNet = nc_nNet, 
    nCells = nc_nCells, 
    scaleScores = nc_scaleScores, 
    symmetric = nc_symmetric, 
    nComp = nc_nComp,
    nCores = nCores,
    verbose = verbose
  )
  
  if(verbose) cat("\n")
  
  # Tensor decomposition
  if(verbose) cat("Step 3/6: Tensor decomposition...\n")
  WT <- tensorDecomposition(
    xList = WT, 
    K = td_K, 
    maxError = td_maxError, 
    maxIter = td_maxIter, 
    nDecimal = td_nDecimal
  )
  
  if(verbose) cat("  Tensor decomposition complete\n\n")
  
  # Perform knockout
  if(verbose) cat("Step 4/6: Performing virtual knockout...\n")
  WT <- WT$X
  WT <- strictDirection(WT, lambda = nc_lambda)
  WT <- as.matrix(WT)
  diag(WT) <- 0
  WT <- t(WT)
  
  KO <- WT
  if(is.null(gKO)){
    stop("gKO must be specified (gene to knockout)")
  }
  
  if(!gKO %in% rownames(KO)){
    stop(sprintf("Gene '%s' not found in the network. Available genes: %s", 
                gKO, paste(head(rownames(KO)), collapse = ", ")))
  }
  
  KO[gKO, ] <- 0
  
  if(verbose) cat(sprintf("  Gene '%s' knocked out\n\n", gKO))
  
  # Manifold alignment
  if(verbose) cat("Step 5/6: Manifold alignment...\n")
  MA <- manifoldAlignment(WT, KO, d = ma_nDim)
  if(verbose) cat("  Manifold alignment complete\n\n")
  
  # Differential regulation
  if(verbose) cat("Step 6/6: Differential regulation analysis...\n")
  DR <- dRegulation(MA, gKO)
  if(verbose){
    n_sig <- sum(DR$p.adj < 0.05)
    cat(sprintf("  Found %d significantly affected genes (FDR < 0.05)\n\n", n_sig))
  }
  
  # Prepare output
  outputList <- list()
  outputList$tensorNetworks$WT <- Matrix::Matrix(WT, sparse = TRUE)
  outputList$tensorNetworks$KO <- Matrix::Matrix(KO, sparse = TRUE)
  outputList$manifoldAlignment <- MA
  outputList$diffRegulation <- DR
  
  if(verbose) cat("=== Analysis complete ===\n")
  
  return(outputList)
}

