#' @export manifoldAlignment
#' @importFrom RSpectra eigs
#' @title Performs non-linear manifold alignment of two gene regulatory networks
#' @description Build comparable low-dimensional features for two weight-averaged denoised 
#' single-cell gene regulatory networks using non-linear network embedding
#' @param X A gene regulatory network
#' @param Y A gene regulatory network
#' @param d The dimension of the low-dimensional feature space
#' @param nCores An integer value. Defines the number of cores to be used (currently not used, kept for compatibility)
#' @return A low-dimensional projection for the two gene regulatory networks used as input
manifoldAlignment <- function(X, Y, d = 30, nCores = parallel::detectCores()){
  
  sharedGenes <- intersect(rownames(X), rownames(Y))
  
  if(length(sharedGenes) == 0){
    stop("No shared genes between the two networks")
  }
  
  X <- X[sharedGenes, sharedGenes]
  Y <- Y[sharedGenes, sharedGenes]
  
  nGenes <- length(sharedGenes)
  L <- diag(nGenes)
  
  # Build alignment matrix
  wX <- as.matrix(X) + 1
  wY <- as.matrix(Y) + 1
  wXY <- 0.9 * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  
  # Construct combined Laplacian
  W <- rbind(
    cbind(wX, wXY), 
    cbind(t(wXY), wY)
  )
  
  W <- -W
  diag(W) <- 0
  diag(W) <- -rowSums(W)
  
  # Compute eigenvectors
  E <- suppressWarnings(RSpectra::eigs(W, k = d * 2, which = 'SR'))
  E$values <- suppressWarnings(as.numeric(E$values))
  E$vectors <- suppressWarnings(apply(E$vectors, 2, as.numeric))
  
  # Sort by eigenvalue
  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- E$vectors[, newOrder, drop = FALSE]
  
  # Filter out near-zero eigenvalues
  E$vectors <- E$vectors[, E$values > 1e-8, drop = FALSE]
  
  # Select top d dimensions
  if(ncol(E$vectors) < d){
    warning(sprintf("Only %d valid dimensions found (requested %d)", ncol(E$vectors), d))
    d <- ncol(E$vectors)
  }
  
  alignedNet <- E$vectors[, 1:d, drop = FALSE]
  colnames(alignedNet) <- paste0('NLMA_', 1:d)
  rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
  
  return(alignedNet)
}

