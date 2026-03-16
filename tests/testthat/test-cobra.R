context("test COBRA")

test_that("cobra() returns correct structure with pearson method", {
  skip_if_not_installed("rARPACK")
  set.seed(42)
  g <- 20
  n <- 8
  q <- 2
  X <- cbind(rep(1, n), rbinom(n, 1, 0.5))
  expressionData <- matrix(rnorm(g * n), nrow = g, ncol = n)

  result <- cobra(X, expressionData, method = "pearson")

  expect_type(result, "list")
  expect_named(result, c("psi", "Q", "D", "G"))
  expect_equal(nrow(result$psi), q)
  expect_equal(ncol(result$psi), n)
  expect_equal(nrow(result$G), g)
  expect_equal(ncol(result$G), n)
})

test_that("cobra() works with pcorsh method", {
  skip_if_not_installed("rARPACK")
  skip_if_not_installed("corpcor")
  set.seed(42)
  g <- 20
  n <- 8
  q <- 2
  X <- cbind(rep(1, n), rbinom(n, 1, 0.5))
  expressionData <- matrix(rnorm(g * n), nrow = g, ncol = n)

  result <- cobra(X, expressionData, method = "pcorsh")

  expect_type(result, "list")
  expect_named(result, c("psi", "Q", "D", "G"))
})

test_that("cobra() psi dimensions match design matrix", {
  skip_if_not_installed("rARPACK")
  set.seed(123)
  g <- 15
  n <- 10
  q <- 3
  X <- cbind(rep(1, n), rnorm(n), rbinom(n, 1, 0.5))
  expressionData <- matrix(rnorm(g * n), nrow = g, ncol = n)

  result <- cobra(X, expressionData)

  expect_equal(nrow(result$psi), q)
  expect_equal(ncol(result$psi), n)
})

test_that("cobra() with dragon method works", {
  skip_if_not_installed("rARPACK")
  set.seed(42)
  g1 <- 10
  g2 <- 8
  n <- 12
  q <- 2
  X <- cbind(rep(1, n), rbinom(n, 1, 0.5))
  layer1 <- matrix(rnorm(g1 * n), nrow = g1, ncol = n)
  layer2 <- matrix(rnorm(g2 * n), nrow = g2, ncol = n)

  result <- cobra(X, list(layer1, layer2), method = "dragon")

  expect_type(result, "list")
  expect_named(result, c("psi", "Q", "D", "G"))
  expect_equal(nrow(result$psi), q)
  expect_equal(ncol(result$psi), n)
  # G should have g1+g2 rows
  expect_equal(nrow(result$G), g1 + g2)
})

test_that("cobra() errors with invalid method", {
  skip_if_not_installed("rARPACK")
  set.seed(42)
  X <- cbind(rep(1, 8), rbinom(8, 1, 0.5))
  expressionData <- matrix(rnorm(160), nrow = 20, ncol = 8)

  expect_error(
    cobra(X, expressionData, method = "invalid"),
    "Only Pearson and pcor methods are supported"
  )
})

test_that("cobra() dragon method errors with wrong input", {
  skip_if_not_installed("rARPACK")
  set.seed(42)
  X <- cbind(rep(1, 8), rbinom(8, 1, 0.5))
  expressionData <- matrix(rnorm(160), nrow = 20, ncol = 8)

  expect_error(
    cobra(X, expressionData, method = "dragon"),
    "Dragon needs two layers"
  )
})

test_that("cobra() dragon method errors with mismatched samples", {
  skip_if_not_installed("rARPACK")
  set.seed(42)
  X <- cbind(rep(1, 8), rbinom(8, 1, 0.5))
  layer1 <- matrix(rnorm(80), nrow = 10, ncol = 8)
  layer2 <- matrix(rnorm(60), nrow = 10, ncol = 6)

  expect_error(
    cobra(X, list(layer1, layer2), method = "dragon"),
    "same number of samples"
  )
})

test_that("cobra() Q eigenvectors have correct dimensions", {
  skip_if_not_installed("rARPACK")
  set.seed(42)
  g <- 20
  n <- 8
  X <- cbind(rep(1, n), rbinom(n, 1, 0.5))
  expressionData <- matrix(rnorm(g * n), nrow = g, ncol = n)

  result <- cobra(X, expressionData, method = "pearson")

  expect_equal(nrow(result$Q), g)
  expect_equal(length(result$D), n)
})
