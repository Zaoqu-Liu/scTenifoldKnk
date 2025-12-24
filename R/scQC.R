#' @importFrom Matrix Matrix colSums
#' @importFrom stats lm predict
#' @title Single-cell quality control
#' @description Performs quality control filtering on single-cell RNA-seq data
#' @param X A count matrix or Seurat object
#' @param mtThreshold Maximum mitochondrial read ratio threshold
#' @param minLSize Minimum library size threshold
#' @return Filtered count matrix
scQC <- function(X, mtThreshold = 0.1, minLSize = 1000) {
  if (inherits(X, "Seurat")) {
    countMatrix <- X@assays$RNA@counts
  } else {
    countMatrix <- X
  }

  librarySize <- Matrix::colSums(countMatrix)
  countMatrix <- countMatrix[, librarySize >= minLSize, drop = FALSE]
  librarySize <- Matrix::colSums(countMatrix)

  mtGenes <- grep("^MT-", toupper(rownames(countMatrix)))
  nGenes <- Matrix::colSums(countMatrix != 0)

  genesLM <- lm(nGenes ~ librarySize)
  genesLM <- as.data.frame(predict(genesLM, data.frame(librarySize), interval = "prediction"))

  if (isTRUE(length(mtGenes) > 0)) {
    mtCounts <- Matrix::colSums(countMatrix[grep("^MT-", toupper(rownames(countMatrix))), , drop = FALSE])
    mtProportion <- mtCounts / librarySize
    mtLM <- lm(mtCounts ~ librarySize)
    mtLM <- as.data.frame(predict(mtLM, data.frame(librarySize), interval = "prediction"))
    selectedCells <- mtCounts > mtLM$lwr & mtCounts < mtLM$upr &
      nGenes > genesLM$lwr & nGenes < genesLM$upr &
      mtProportion <= mtThreshold & librarySize < 2 * mean(librarySize)
  } else {
    selectedCells <- nGenes > genesLM$lwr & nGenes < genesLM$upr &
      librarySize < 2 * mean(librarySize)
  }

  selectedCells <- colnames(countMatrix)[selectedCells]

  if (inherits(X, "Seurat")) {
    X <- subset(X, cells = selectedCells)
  } else {
    X <- countMatrix[, selectedCells, drop = FALSE]
  }

  return(X)
}
