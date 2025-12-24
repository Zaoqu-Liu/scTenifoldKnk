test_that("makeNetworksFast produces correct number of networks", {
  set.seed(123)
  nCells <- 200
  nGenes <- 50
  X <- matrix(rnbinom(n = nGenes * nCells, size = 20, prob = 0.9),
    ncol = nCells
  )
  rownames(X) <- paste0("gene", 1:nGenes)
  X <- X[rowSums(X) > 0, ]

  networks <- makeNetworksFast(X,
    nNet = 5, nCells = 100, nComp = 3,
    nCores = 1, verbose = FALSE
  )

  expect_equal(length(networks), 5)
  expect_true(all(sapply(networks, function(x) inherits(x, "dgCMatrix"))))
})

test_that("makeNetworksFast parallel vs sequential gives same results", {
  set.seed(123)
  nCells <- 200
  nGenes <- 50
  X <- matrix(rnbinom(n = nGenes * nCells, size = 20, prob = 0.9),
    ncol = nCells
  )
  rownames(X) <- paste0("gene", 1:nGenes)
  X <- X[rowSums(X) > 0, ]

  # Sequential
  set.seed(456)
  nets_seq <- makeNetworksFast(X,
    nNet = 3, nCells = 100, nComp = 3,
    nCores = 1, verbose = FALSE, seed = 456
  )

  # Parallel (if cores available)
  if (parallel::detectCores() > 1) {
    set.seed(456)
    nets_par <- makeNetworksFast(X,
      nNet = 3, nCells = 100, nComp = 3,
      nCores = 2, verbose = FALSE, seed = 456
    )

    # Should produce identical results with same seed
    for (i in 1:3) {
      expect_equal(sum(abs(nets_seq[[i]] - nets_par[[i]])), 0)
    }
  }
})
