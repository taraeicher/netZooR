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
