test_that("pcNetFast produces valid network", {
  # Generate test data
  set.seed(123)
  nCells <- 100
  nGenes <- 50
  X <- matrix(rnbinom(n = nGenes * nCells, size = 20, prob = 0.9),
    ncol = nCells
  )
  rownames(X) <- paste0("gene", 1:nGenes)
  X <- X[rowSums(X) > 0, ]

  # Test network construction
  network <- pcNetFast(X,
    nComp = 3, scaleScores = TRUE,
    symmetric = FALSE, q = 0.9, verbose = FALSE
  )

  # Check dimensions
  expect_equal(nrow(network), nrow(X))
  expect_equal(ncol(network), nrow(X))

  # Check sparsity
  expect_true(sum(network != 0) < nrow(X) * nrow(X))

  # Check diagonal is zero
  expect_equal(sum(diag(network)), 0)

  # Check it's a sparse matrix
  expect_s4_class(network, "dgCMatrix")

  # Check gene names preserved
  expect_equal(rownames(network), rownames(X))
})

test_that("pcNetFast handles symmetric option", {
  set.seed(123)
  nCells <- 100
  nGenes <- 30
  X <- matrix(rnbinom(n = nGenes * nCells, size = 20, prob = 0.9),
    ncol = nCells
  )
  rownames(X) <- paste0("gene", 1:nGenes)
  X <- X[rowSums(X) > 0, ]

  network <- pcNetFast(X, nComp = 3, symmetric = TRUE, verbose = FALSE)

  # Check symmetry
  expect_equal(as.matrix(network), t(as.matrix(network)), tolerance = 1e-10)
})

test_that("pcNetFast input validation works", {
  X <- matrix(1:100, ncol = 10)
  rownames(X) <- paste0("gene", 1:10)

  # Test nComp too large
  expect_error(pcNetFast(X, nComp = 15), "nComp should be")

  # Test nComp too small
  expect_error(pcNetFast(X, nComp = 1), "nComp should be")
})
