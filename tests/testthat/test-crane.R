context("test CRANE")

test_that("craneBipartite() returns edge list from adjacency matrix", {
  set.seed(42)
  A <- matrix(runif(12), nrow = 3, ncol = 4,
              dimnames = list(paste0("TF", 1:3), paste0("G", 1:4)))

  result <- craneBipartite(A, alpha = 0.1)

  expect_true(is.data.frame(result))
  expect_equal(ncol(result), 3)
  expect_equal(nrow(result), 12)  # 3 * 4 edges
})

test_that("craneBipartite() with getAdj returns matrix", {
  set.seed(42)
  A <- matrix(runif(12), nrow = 3, ncol = 4,
              dimnames = list(paste0("TF", 1:3), paste0("G", 1:4)))

  result <- craneBipartite(A, alpha = 0.1, getAdj = TRUE)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 4))
})

test_that("craneBipartite() preserves row and column sums approximately", {
  set.seed(42)
  A <- matrix(runif(12), nrow = 3, ncol = 4,
              dimnames = list(paste0("TF", 1:3), paste0("G", 1:4)))

  result <- craneBipartite(A, alpha = 0.0, getAdj = TRUE)

  # With alpha=0, node strengths should be preserved exactly
  expect_equal(rowSums(result), rowSums(A), tolerance = 1e-6)
  expect_equal(colSums(result), colSums(A), tolerance = 1e-6)
})

test_that("craneUnipartite() returns perturbed matrix", {
  set.seed(42)
  n <- 5
  A <- matrix(runif(n * n), n, n,
              dimnames = list(paste0("N", 1:n), paste0("N", 1:n)))
  A <- (A + t(A)) / 2  # symmetric

  result <- craneUnipartite(A, alpha = 0.1)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(n, n))
})

test_that("craneUnipartite() preserves node strengths with alpha=0", {
  set.seed(42)
  n <- 5
  A <- matrix(runif(n * n), n, n,
              dimnames = list(paste0("N", 1:n), paste0("N", 1:n)))
  A <- (A + t(A)) / 2

  result <- craneUnipartite(A, alpha = 0.0)

  expect_equal(rowSums(result), rowSums(A), tolerance = 1e-6)
})

test_that("craneBipartite() accepts edge list input", {
  set.seed(42)
  el <- data.frame(
    TF = rep(paste0("TF", 1:3), each = 4),
    Gene = rep(paste0("G", 1:4), 3),
    Weight = runif(12)
  )

  result <- craneBipartite(el, alpha = 0.1)

  expect_true(is.data.frame(result))
  expect_equal(ncol(result), 3)
})
