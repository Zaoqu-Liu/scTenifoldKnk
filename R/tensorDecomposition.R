#' @export tensorDecomposition
#' @title Performs CANDECOMP/PARAFAC (CP) Tensor Decomposition
#' @description Generate weight-averaged denoised gene regulatory networks using CANDECOMP/PARAFAC (CP) Tensor Decomposition
#' @param xList A list of gene regulatory networks
#' @param yList Optional. A list of gene regulatory networks
#' @param nDecimal An integer value indicating the number of decimal places to be used
#' @param K The number of rank-one tensors used to approximate the data
#' @param maxError A decimal value between 0 and 1. Defines the relative Frobenius norm error tolerance
#' @param maxIter An integer value. Defines the maximum number of iterations
#' @return A list of weight-averaged denoised gene regulatory networks
tensorDecomposition <- function(xList, yList = NULL, nDecimal = 1, K = 5, maxError = 1e-5, maxIter = 1e3) {
  xNets <- length(xList)
  if (!is.null(yList)) {
    yNets <- length(yList)
    if (xNets != yNets) {
      stop("Same number of networks are required in both cases")
    }
    nNet <- unique(c(xNets, yNets))
    xGenes <- unique(unlist(lapply(xList, rownames)))
    yGenes <- unique(unlist(lapply(yList, rownames)))
    sGenes <- intersect(xGenes, yGenes)
  } else {
    nNet <- xNets
    xGenes <- unique(unlist(lapply(xList, rownames)))
    sGenes <- xGenes
  }
  nGenes <- length(sGenes)

  # Build tensors
  tensorX <- array(data = 0, dim = c(nGenes, nGenes, 1, nNet))
  if (!is.null(yList)) {
    tensorY <- array(data = 0, dim = c(nGenes, nGenes, 1, nNet))
  }

  for (i in seq_len(nNet)) {
    tempX <- matrix(0, nGenes, nGenes)
    rownames(tempX) <- colnames(tempX) <- sGenes
    temp <- as.matrix(xList[[i]])
    tGenes <- sGenes[sGenes %in% rownames(temp)]
    tempX[tGenes, tGenes] <- temp[tGenes, tGenes]
    tensorX[, , , i] <- tempX

    if (!is.null(yList)) {
      tempY <- matrix(0, nGenes, nGenes)
      rownames(tempY) <- colnames(tempY) <- sGenes
      temp <- as.matrix(yList[[i]])
      tGenes <- sGenes[sGenes %in% rownames(temp)]
      tempY[tGenes, tGenes] <- temp[tGenes, tGenes]
      tensorY[, , , i] <- tempY
    }
  }

  set.seed(1)
  tensorX <- as.tensor(tensorX)
  tensorX <- cpDecomposition(tnsr = tensorX, num_components = K, max_iter = maxIter, tol = maxError)
  tX <- tensorX$est$data[, , , 1]
  for (i in seq_len(nNet)[-1]) {
    tX <- tX + tensorX$est$data[, , , i]
  }
  tX <- tX / nNet
  tX <- tX / max(abs(tX))
  tX <- round(tX, nDecimal)
  tX <- Matrix::Matrix(tX, sparse = TRUE)
  rownames(tX) <- colnames(tX) <- sGenes

  if (!is.null(yList)) {
    set.seed(1)
    tensorY <- as.tensor(tensorY)
    tensorY <- cpDecomposition(tnsr = tensorY, num_components = K, max_iter = maxIter, tol = maxError)
    tY <- tensorY$est$data[, , , 1]
    for (i in seq_len(nNet)[-1]) {
      tY <- tY + tensorY$est$data[, , , i]
    }
    tY <- tY / nNet
    tY <- tY / max(abs(tY))
    tY <- round(tY, nDecimal)
    tY <- Matrix::Matrix(tY, sparse = TRUE)
    rownames(tY) <- colnames(tY) <- sGenes
  }

  tensorOutput <- list()
  tensorOutput$X <- tX
  if (!is.null(yList)) {
    tensorOutput$Y <- tY
  }
  return(tensorOutput)
}

# Helper function to create tensor object
as.tensor <- function(arr) {
  list(
    data = arr,
    num_modes = length(dim(arr)),
    modes = dim(arr)
  )
}
